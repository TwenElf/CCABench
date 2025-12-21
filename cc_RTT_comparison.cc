#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/tcp-l4-protocol.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

using namespace ns3;

// =================== Utilities ===================

static void EnsureDir (const std::string& dir)
{
  std::string cmd = "mkdir -p " + dir;
  std::system (cmd.c_str ());
}

static std::string ToSafeName (const std::string& s)
{
  // e.g. "ns3::TcpVegas" -> "TcpVegas"
  auto pos = s.find ("::");
  if (pos != std::string::npos) return s.substr (pos + 2);
  return s;
}

static TypeId MustLookupTypeId (const std::string& typeName)
{
  TypeId tid;
  bool ok = TypeId::LookupByNameFailSafe (typeName, &tid);
  if (!ok)
    {
      std::cerr << "[FATAL] Cannot find TypeId for: " << typeName << "\n";
      std::cerr << "        Check your ns-3 build has this TCP variant.\n";
      std::exit (1);
    }
  return tid;
}

// =================== Topology (N flows dumbbell) ===================
//
// sender[i] -- r0 --(bottleneck)-- r1 -- receiver[i]
//
struct NFlowTopology
{
  NodeContainer senders;    // N
  NodeContainer routers;    // 2
  NodeContainer receivers;  // N
  std::vector<Ipv4Address> recvIps; // receiver[i] IP
};

static Time ComputeBottleDelay (double targetRttSec, Time accessDelay)
{
  // One-way path delay = access + bottle + access = 2*access + bottle
  // RTT = 2*(2*access + bottle) = 4*access + 2*bottle
  // => bottle = RTT/2 - 2*access
  double bottle = targetRttSec / 2.0 - 2.0 * accessDelay.GetSeconds ();
  if (bottle < 0.0) bottle = 0.0;
  return Seconds (bottle);
}

static uint32_t ComputeQueuePktsBdp (const std::string& bottleRateStr,
                                    double targetRttSec,
                                    uint32_t segmentSizeBytes,
                                    double bdpFactor)
{
  // BDP(bytes) = rate(bps) * RTT(s) / 8
  DataRate rate (bottleRateStr);
  double bdpBytes = (rate.GetBitRate () * targetRttSec) / 8.0;
  double targetBytes = bdpBytes * bdpFactor;

  uint32_t pkts = static_cast<uint32_t> (std::ceil (targetBytes / segmentSizeBytes));
  if (pkts < 1) pkts = 1;
  return pkts;
}

static NFlowTopology BuildNFlowDumbbell (uint32_t nFlows,
                                        const std::string& bottleRate,
                                        double targetRttSec,
                                        const std::string& accessRate,
                                        Time accessDelay,
                                        uint32_t queuePkts)
{
  NFlowTopology topo;
  topo.senders.Create (nFlows);
  topo.receivers.Create (nFlows);
  topo.routers.Create (2);

  InternetStackHelper stack;
  stack.Install (topo.senders);
  stack.Install (topo.receivers);
  stack.Install (topo.routers);

  PointToPointHelper access;
  access.SetDeviceAttribute ("DataRate", StringValue (accessRate));
  access.SetChannelAttribute ("Delay", TimeValue (accessDelay));

  PointToPointHelper bottle;
  bottle.SetDeviceAttribute ("DataRate", StringValue (bottleRate));
  bottle.SetChannelAttribute ("Delay", TimeValue (ComputeBottleDelay (targetRttSec, accessDelay)));
  bottle.SetQueue ("ns3::DropTailQueue<Packet>",
                   "MaxSize", QueueSizeValue (QueueSize (std::to_string (queuePkts) + "p")));

  // bottleneck r0 <-> r1
  NetDeviceContainer devBottle = bottle.Install (topo.routers.Get (0), topo.routers.Get (1));

  Ipv4AddressHelper addr;

  // r0-r1
  addr.SetBase ("10.0.0.0", "255.255.255.0");
  addr.Assign (devBottle);

  topo.recvIps.resize (nFlows);

  // Each flow gets its own access subnets
  for (uint32_t i = 0; i < nFlows; ++i)
    {
      // sender[i] <-> r0
      {
        NetDeviceContainer d = access.Install (topo.senders.Get (i), topo.routers.Get (0));
        std::ostringstream net;
        net << "10.1." << (i + 1) << ".0";
        addr.SetBase (net.str ().c_str (), "255.255.255.0");
        addr.Assign (d);
      }

      // r1 <-> receiver[i]
      {
        NetDeviceContainer d = access.Install (topo.routers.Get (1), topo.receivers.Get (i));
        std::ostringstream net;
        net << "10.2." << (i + 1) << ".0";
        addr.SetBase (net.str ().c_str (), "255.255.255.0");
        Ipv4InterfaceContainer ifc = addr.Assign (d);
        topo.recvIps[i] = ifc.GetAddress (1); // receiver[i] side IP
      }
    }

  Ipv4GlobalRoutingHelper::PopulateRoutingTables ();
  return topo;
}

// =================== Apps ===================

static Ptr<PacketSink> InstallSink (Ptr<Node> node, uint16_t port, double start, double stop)
{
  Address local (InetSocketAddress (Ipv4Address::GetAny (), port));
  PacketSinkHelper sinkHelper ("ns3::TcpSocketFactory", local);
  ApplicationContainer apps = sinkHelper.Install (node);
  apps.Start (Seconds (start));
  apps.Stop (Seconds (stop));
  return apps.Get (0)->GetObject<PacketSink> ();
}

static ApplicationContainer InstallBulk (Ptr<Node> node,
                                        Ipv4Address dstIp,
                                        uint16_t port,
                                        uint32_t sendSize,
                                        double start,
                                        double stop)
{
  Address remote (InetSocketAddress (dstIp, port));
  BulkSendHelper bulk ("ns3::TcpSocketFactory", remote);
  bulk.SetAttribute ("MaxBytes", UintegerValue (0));
  bulk.SetAttribute ("SendSize", UintegerValue (sendSize));

  ApplicationContainer apps = bulk.Install (node);
  apps.Start (Seconds (start));
  apps.Stop (Seconds (stop));
  return apps;
}

// =================== Convergence monitoring (same semantics as your old code) ===================
//
// Metrics:
//  1) T_bar: average throughput in last tailWindow seconds of sim (steady-state proxy)
//  2) relDiff: |T_last - T_prev| / T_last  (two tail windows)
//     approxStable = 1 if relDiff < 0.05
//  3) Tconv_abs: first time where sliding avg (windowLen) >= eps * T_bar
//  4) Tconv_rel = Tconv_abs - appStart
//
class ConvergenceMonitor
{
public:
  ConvergenceMonitor (Ptr<PacketSink> sink,
                      Ptr<OutputStreamWrapper> stream,
                      double sampleInterval,
                      double measureStart,
                      double simTime,
                      double appStart,
                      double tailWindow,
                      double eps,
                      double windowLen)
    : m_sink (sink),
      m_stream (stream),
      m_sampleInterval (sampleInterval),
      m_measureStart (measureStart),
      m_simTime (simTime),
      m_appStart (appStart),
      m_tailWindow (tailWindow),
      m_eps (eps),
      m_windowLen (windowLen),
      m_initialized (false),
      m_lastTotalRx (0)
  {
    *m_stream->GetStream () << "# time_sec, throughput_Mbps" << std::endl;
  }

  void Start ()
  {
    Simulator::Schedule (Seconds (m_sampleInterval),
                         &ConvergenceMonitor::Sample, this);
  }

  void Finish ()
  {
    ComputeStatsAndPrint ();
  }

private:
  void Sample ()
  {
    double now = Simulator::Now ().GetSeconds ();
    uint64_t curRx = m_sink->GetTotalRx ();

    if (!m_initialized)
      {
        m_lastTotalRx = curRx;
        m_initialized = true;
        Simulator::Schedule (Seconds (m_sampleInterval),
                             &ConvergenceMonitor::Sample, this);
        return;
      }

    double mbps =
        (curRx - m_lastTotalRx) * 8.0 / (m_sampleInterval * 1e6);
    m_lastTotalRx = curRx;

    if (now >= m_measureStart)
      {
        *m_stream->GetStream () << std::fixed << std::setprecision (6)
                                << now << "," << mbps << std::endl;
        m_times.push_back (now);
        m_mbps.push_back (mbps);
      }

    Simulator::Schedule (Seconds (m_sampleInterval),
                         &ConvergenceMonitor::Sample, this);
  }

  void ComputeStatsAndPrint ()
  {
    std::ostream &o = *m_stream->GetStream ();

    if (m_times.empty ())
      {
        o << "# No samples collected" << std::endl;
        return;
      }

    // ---- 1) T_bar from last tailWindow ----
    double tailStartLast = m_simTime - m_tailWindow;
    if (tailStartLast < m_measureStart)
      {
        tailStartLast = m_measureStart;
      }

    double tailStartPrev = m_simTime - 2.0 * m_tailWindow;
    if (tailStartPrev < m_measureStart)
      {
        tailStartPrev = m_measureStart;
      }

    double sumPrev = 0.0;
    int    countPrev = 0;
    double sumLast = 0.0;
    int    countLast = 0;

    for (size_t i = 0; i < m_times.size (); ++i)
      {
        double t = m_times[i];
        double mbps = m_mbps[i];

        if (t >= tailStartPrev && t < tailStartLast)
          {
            sumPrev += mbps;
            countPrev++;
          }
        if (t >= tailStartLast)
          {
            sumLast += mbps;
            countLast++;
          }
      }

    double T_prev = (countPrev > 0 ? sumPrev / countPrev : 0.0);
    double T_last = (countLast > 0 ? sumLast / countLast : 0.0);
    double T_bar  = T_last;

    // ---- 2) relDiff + approxStable ----
    double relDiff = 0.0;
    bool   approxStable = false;
    if (T_last > 0.0 && countPrev > 0)
      {
        relDiff = std::fabs (T_last - T_prev) / T_last;
        approxStable = (relDiff < 0.05);
      }

    // ---- 3) epsilon convergence time ----
    double Tconv_abs = -1.0;

    if (T_bar > 0.0)
      {
        double threshold = m_eps * T_bar;
        int windowSize =
            std::max (1, static_cast<int> (m_windowLen / m_sampleInterval));

        for (size_t i = 0; i < m_mbps.size (); ++i)
          {
            if (i + 1 < static_cast<size_t> (windowSize))
              {
                continue;
              }

            double sumWin = 0.0;
            for (int k = 0; k < windowSize; ++k)
              {
                sumWin += m_mbps[i - k];
              }
            double That = sumWin / windowSize;

            if (That >= threshold)
              {
                Tconv_abs = m_times[i];
                break;
              }
          }
      }

    double Tconv_rel = (Tconv_abs >= 0.0 ? Tconv_abs - m_appStart : -1.0);

    // ---- 4) Output to CSV (comment line) ----
    o << "# T_bar_Mbps=" << T_bar
      << ", T_prev_Mbps=" << T_prev
      << ", T_last_Mbps=" << T_last
      << ", relDiff=" << relDiff
      << ", approxStable=" << (approxStable ? 1 : 0)
      << ", eps=" << m_eps
      << ", windowLen_s=" << m_windowLen
      << ", Tconv_abs_s=" << Tconv_abs
      << ", Tconv_rel_s=" << Tconv_rel
      << std::endl;

    // ---- 5) Print console (match your screenshot style) ----
    std::cout << "[Convergence] T_bar=" << T_bar
              << " Mbps, relDiff=" << relDiff
              << " (approxStable=" << (approxStable ? 1 : 0)
              << "), Tconv_abs=" << Tconv_abs
              << " s, Tconv_rel=" << Tconv_rel
              << " s since flow start"
              << std::endl;
  }

private:
  Ptr<PacketSink>          m_sink;
  Ptr<OutputStreamWrapper> m_stream;

  double m_sampleInterval;
  double m_measureStart;
  double m_simTime;
  double m_appStart;

  double m_tailWindow;
  double m_eps;
  double m_windowLen;

  bool      m_initialized;
  uint64_t  m_lastTotalRx;

  std::vector<double> m_times;
  std::vector<double> m_mbps;
};

// =================== Experiment Runners ===================

static void RunSingleFlow (const std::string& outDir,
                           const std::string& ccaTypeName,
                           const std::string& bottleRate,
                           double rttSec,
                           const std::string& accessRate,
                           Time accessDelay,
                           uint32_t segmentSize,
                           uint32_t sendSize,
                           uint32_t queuePkts,
                           double simTime,
                           double sinkStart,
                           double appStart,
                           double sampleInterval,
                           double tailWindow,
                           double eps,
                           double windowLen)
{
  std::cout << "\n============================\n";
  std::cout << " Single-flow: " << ToSafeName (ccaTypeName) << "\n";
  std::cout << "============================\n";
  std::cout << "\n--- RTT = " << rttSec * 1000.0 << " ms ---\n";

  // Global default CCA (single flow)
  Config::SetDefault ("ns3::TcpL4Protocol::SocketType",
                      TypeIdValue (MustLookupTypeId (ccaTypeName)));

  // Common TCP tuning
  Config::SetDefault ("ns3::TcpSocket::SegmentSize", UintegerValue (segmentSize));
  Config::SetDefault ("ns3::TcpSocket::DelAckCount", UintegerValue (1));
  Config::SetDefault ("ns3::TcpSocket::SndBufSize", UintegerValue (1 << 20));
  Config::SetDefault ("ns3::TcpSocket::RcvBufSize", UintegerValue (1 << 20));

  NFlowTopology topo =
    BuildNFlowDumbbell (1, bottleRate, rttSec, accessRate, accessDelay, queuePkts);

  uint16_t port = 5000;

  Ptr<PacketSink> sink =
    InstallSink (topo.receivers.Get (0), port, sinkStart, simTime);

  InstallBulk (topo.senders.Get (0), topo.recvIps[0], port, sendSize, appStart, simTime);

  AsciiTraceHelper ascii;
  std::ostringstream fn;
  int rttMs = static_cast<int> (std::round (rttSec * 1000.0));

  fn << outDir << "/throughput_" << ToSafeName (ccaTypeName)
     << "_single_rtt_" << rttMs << "ms.csv";

  Ptr<OutputStreamWrapper> stream = ascii.CreateFileStream (fn.str ());

  ConvergenceMonitor mon (sink, stream,
                          sampleInterval,
                          appStart,     // measureStart
                          simTime,
                          appStart,
                          tailWindow,
                          eps,
                          windowLen);
  mon.Start ();

  Simulator::Stop (Seconds (simTime));
  Simulator::Run ();

  mon.Finish ();
  Simulator::Destroy ();

  std::cout << "Finished single-flow CCA=" << ToSafeName (ccaTypeName)
            << ", RTT=" << rttMs << " ms\n";
}

static void RunTwoFlowPair (const std::string& outDir,
                            const std::string& ccaA,
                            const std::string& ccaB,
                            const std::string& bottleRate,
                            double rttSec,
                            const std::string& accessRate,
                            Time accessDelay,
                            uint32_t segmentSize,
                            uint32_t sendSize,
                            uint32_t queuePkts,
                            double simTime,
                            double sinkStart,
                            double appStart,
                            double sampleInterval,
                            double tailWindow,
                            double eps,
                            double windowLen)
{
  std::cout << "\n=======================================\n";
  std::cout << " Two-flow competition, RTT = " << rttSec * 1000.0 << " ms\n";
  std::cout << " Flow 0: " << ToSafeName (ccaA) << ", Flow 1: " << ToSafeName (ccaB) << "\n";
  std::cout << "=======================================\n";

  NFlowTopology topo =
    BuildNFlowDumbbell (2, bottleRate, rttSec, accessRate, accessDelay, queuePkts);

  // Per-node TCP variants
  Ptr<TcpL4Protocol> tcp0 = topo.senders.Get (0)->GetObject<TcpL4Protocol> ();
  Ptr<TcpL4Protocol> tcp1 = topo.senders.Get (1)->GetObject<TcpL4Protocol> ();
  tcp0->SetAttribute ("SocketType", TypeIdValue (MustLookupTypeId (ccaA)));
  tcp1->SetAttribute ("SocketType", TypeIdValue (MustLookupTypeId (ccaB)));

  // Common TCP tuning
  Config::SetDefault ("ns3::TcpSocket::SegmentSize", UintegerValue (segmentSize));
  Config::SetDefault ("ns3::TcpSocket::DelAckCount", UintegerValue (1));
  Config::SetDefault ("ns3::TcpSocket::SndBufSize", UintegerValue (1 << 20));
  Config::SetDefault ("ns3::TcpSocket::RcvBufSize", UintegerValue (1 << 20));

  uint16_t port0 = 5000;
  uint16_t port1 = 5001;

  Ptr<PacketSink> sink0 =
    InstallSink (topo.receivers.Get (0), port0, sinkStart, simTime);
  Ptr<PacketSink> sink1 =
    InstallSink (topo.receivers.Get (1), port1, sinkStart, simTime);

  InstallBulk (topo.senders.Get (0), topo.recvIps[0], port0, sendSize, appStart, simTime);
  InstallBulk (topo.senders.Get (1), topo.recvIps[1], port1, sendSize, appStart, simTime);

  AsciiTraceHelper ascii;
  int rttMs = static_cast<int> (std::round (rttSec * 1000.0));

  std::ostringstream f0, f1;
  f0 << outDir << "/throughput_"
     << ToSafeName (ccaA) << "_vs_" << ToSafeName (ccaB)
     << "_flow0_rtt_" << rttMs << "ms.csv";

  f1 << outDir << "/throughput_"
     << ToSafeName (ccaA) << "_vs_" << ToSafeName (ccaB)
     << "_flow1_rtt_" << rttMs << "ms.csv";

  Ptr<OutputStreamWrapper> s0 = ascii.CreateFileStream (f0.str ());
  Ptr<OutputStreamWrapper> s1 = ascii.CreateFileStream (f1.str ());

  ConvergenceMonitor m0 (sink0, s0, sampleInterval, appStart, simTime, appStart,
                         tailWindow, eps, windowLen);
  ConvergenceMonitor m1 (sink1, s1, sampleInterval, appStart, simTime, appStart,
                         tailWindow, eps, windowLen);

  m0.Start ();
  m1.Start ();

  Simulator::Stop (Seconds (simTime));
  Simulator::Run ();

  m0.Finish ();
  m1.Finish ();

  Simulator::Destroy ();

  std::cout << "Finished two-flow " << ToSafeName (ccaA) << " vs " << ToSafeName (ccaB)
            << ", RTT=" << rttMs << " ms\n";
}

// =================== main ===================

int main (int argc, char* argv[])
{
  // ---- defaults ----
  std::string outDir      = "results";
  std::string bottleRate  = "25Mbps";
  std::string accessRate  = "100Mbps";

  // RTT list (seconds)
  std::vector<double> rtts = {0.002, 0.020, 0.040};

  // timing
  double simTime   = 120.0;
  double sinkStart = 0.5;
  double appStart  = 1.0;

  // packet sizing
  uint32_t segmentSize = 1448; // MSS-like
  uint32_t sendSize    = 1500; // BulkSend SendSize

  // measurement (matches your old code)
  double sampleInterval = 0.1;
  double tailWindow     = 20.0;
  double eps            = 0.8;
  double windowLen      = 2.0;

  // access delay (set 0ms to match your old code exactly)
  double accessDelayMs = 0.0;  // <-- if you want to avoid 0ms artifacts, set 1.0

  // queue sizing
  double bdpFactor = 1.0; // 1x BDP

  // CCA names (change bbrType to ns3::TcpBbrV3 if your tree uses that)
  std::string vegasType = "ns3::TcpVegas";
  std::string cubicType = "ns3::TcpCubic";
  std::string bbrType   = "ns3::TcpBbr"; // or "ns3::TcpBbrV3"

  CommandLine cmd;
  cmd.AddValue ("outDir", "Output directory", outDir);
  cmd.AddValue ("bottleRate", "Bottleneck rate", bottleRate);
  cmd.AddValue ("accessRate", "Access link rate", accessRate);

  cmd.AddValue ("simTime", "Simulation time (s)", simTime);
  cmd.AddValue ("sinkStart", "Sink start time (s)", sinkStart);
  cmd.AddValue ("appStart", "App start time (s)", appStart);

  cmd.AddValue ("segmentSize", "TCP SegmentSize (bytes)", segmentSize);
  cmd.AddValue ("sendSize", "BulkSend SendSize (bytes)", sendSize);

  cmd.AddValue ("sampleInterval", "Throughput sampling interval (s)", sampleInterval);
  cmd.AddValue ("tailWindow", "Tail window (s) for T_bar", tailWindow);
  cmd.AddValue ("eps", "Epsilon threshold for convergence", eps);
  cmd.AddValue ("windowLen", "Sliding window length (s)", windowLen);

  cmd.AddValue ("accessDelayMs", "Access link one-way delay (ms)", accessDelayMs);
  cmd.AddValue ("bdpFactor", "Queue size factor * BDP (packets)", bdpFactor);

  cmd.AddValue ("vegasType", "Vegas TypeId name", vegasType);
  cmd.AddValue ("cubicType", "Cubic TypeId name", cubicType);
  cmd.AddValue ("bbrType",   "BBR/BBRv3 TypeId name", bbrType);

  cmd.Parse (argc, argv);

  EnsureDir (outDir);

  Time accessDelay = MilliSeconds (accessDelayMs);

  // Fail fast if CCA types are missing
  MustLookupTypeId (vegasType);
  MustLookupTypeId (cubicType);
  MustLookupTypeId (bbrType);

  // ---- Single flow: Vegas, Cubic, BBR ----
  std::vector<std::string> singleCcas = {vegasType, cubicType, bbrType};

  for (double rtt : rtts)
    {
      uint32_t queuePkts = ComputeQueuePktsBdp (bottleRate, rtt, segmentSize, bdpFactor);

      for (const auto& cca : singleCcas)
        {
          RunSingleFlow (outDir, cca,
                         bottleRate, rtt,
                         accessRate, accessDelay,
                         segmentSize, sendSize,
                         queuePkts,
                         simTime, sinkStart, appStart,
                         sampleInterval,
                         tailWindow, eps, windowLen);
        }
    }

  // ---- Two-flow pairwise competitions ----
  std::vector<std::pair<std::string, std::string>> pairs = {
    {vegasType, cubicType},
    {vegasType, bbrType},
    {cubicType, bbrType}
  };

  for (double rtt : rtts)
    {
      uint32_t queuePkts = ComputeQueuePktsBdp (bottleRate, rtt, segmentSize, bdpFactor);

      for (const auto& p : pairs)
        {
          RunTwoFlowPair (outDir, p.first, p.second,
                          bottleRate, rtt,
                          accessRate, accessDelay,
                          segmentSize, sendSize,
                          queuePkts,
                          simTime, sinkStart, appStart,
                          sampleInterval,
                          tailWindow, eps, windowLen);
        }
    }

  std::cout << "\nAll experiments finished. CSVs in: " << outDir << "\n";
  return 0;
}
