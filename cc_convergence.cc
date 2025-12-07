#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/tcp-l4-protocol.h"
#include "ns3/tcp-cubic.h"
#include "ns3/tcp-vegas.h"

#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <iostream>

using namespace ns3;

// ------------------ Global constants ------------------

static const uint32_t PKT_SIZE   = 1500;   // TCP payload bytes
static const double   SIM_TIME   = 120.0;  // Total simulation time (s)
static const double   START_TIME = 0.5;    // Sink start time (s)
static const double   APP_START  = 1.0;    // BulkSend start time (s)

static const std::string DATA_RATE = "25Mbps";              // Bottleneck bandwidth
static const double RTT_LIST[]     = {0.002, 0.020, 0.040}; // 2, 20, 40 ms

static const std::string RESULTS_FOLDER = "results";        // Output directory

// ------------------ Helper: ensure results directory exists ------------------

void
EnsureResultsDir ()
{
  std::string cmd = "mkdir -p " + RESULTS_FOLDER;
  std::system (cmd.c_str ());   // Warning can be ignored
}

// ------------------ Single-flow topology ------------------

/**
 * Single-flow topology:
 *   sender0 -- r0 --(bottleneck)-- r1 -- receiver0
 */
struct SingleTopology
{
  NodeContainer senders;   // 1 sender
  NodeContainer routers;   // 2 routers (r0, r1)
  NodeContainer receivers; // 1 receiver

  Ipv4Address recvIp;      // IP address of receiver0
};

SingleTopology
BuildSingleDumbbell (const std::string &bottleRate, double targetRtt)
{
  SingleTopology topo;

  topo.senders.Create (1);
  topo.receivers.Create (1);
  topo.routers.Create (2);

  InternetStackHelper stack;
  stack.Install (topo.senders);
  stack.Install (topo.receivers);
  stack.Install (topo.routers);

  // Access links: 100 Mbps, 0 ms delay
  PointToPointHelper access;
  access.SetDeviceAttribute ("DataRate", StringValue ("100Mbps"));
  access.SetChannelAttribute ("Delay", StringValue ("0ms"));

  // Bottleneck link: DATA_RATE, 2 * delay ≈ targetRtt
  PointToPointHelper bottle;
  bottle.SetDeviceAttribute ("DataRate", StringValue (bottleRate));

  // RTT/2 → in milliseconds
  int bottleDelayMs = std::max (1,
      static_cast<int> (targetRtt * 1000.0 / 2.0));
  bottle.SetChannelAttribute ("Delay",
      StringValue (std::to_string (bottleDelayMs) + "ms"));

  // sender0 ↔ r0
  NetDeviceContainer devS0R0 =
      access.Install (topo.senders.Get (0), topo.routers.Get (0));

  // r1 ↔ receiver0
  NetDeviceContainer devR1D0 =
      access.Install (topo.routers.Get (1), topo.receivers.Get (0));

  // r0 ↔ r1 bottleneck
  NetDeviceContainer devBottle =
      bottle.Install (topo.routers.Get (0), topo.routers.Get (1));

  Ipv4AddressHelper address;

  // sender0 - r0
  address.SetBase ("10.1.1.0", "255.255.255.0");
  address.Assign (devS0R0);

  // r0 - r1 bottleneck
  address.SetBase ("10.1.2.0", "255.255.255.0");
  address.Assign (devBottle);

  // r1 - receiver0
  Ipv4InterfaceContainer ifR1D0;
  address.SetBase ("10.1.3.0", "255.255.255.0");
  ifR1D0 = address.Assign (devR1D0);

  topo.recvIp = ifR1D0.GetAddress (1); // receiver0 IP

  Ipv4GlobalRoutingHelper::PopulateRoutingTables ();

  return topo;
}

// ------------------ Two-flow topology ------------------

/**
 * Two-flow competing topology:
 *
 * sender0 --\
 *           r0 ---- bottleneck ---- r1 ---- receiver0  (flow1: Cubic)
 * sender1 --/                            \- receiver1  (flow2: Vegas)
 */
struct MultiTopology
{
  NodeContainer senders;   // 2 senders
  NodeContainer routers;   // 2 routers
  NodeContainer receivers; // 2 receivers

  Ipv4Address recvIp0;     // IP address of receiver0
  Ipv4Address recvIp1;     // IP address of receiver1
};

MultiTopology
BuildMultiDumbbell (const std::string &bottleRate, double targetRtt)
{
  MultiTopology topo;

  topo.senders.Create (2);
  topo.receivers.Create (2);
  topo.routers.Create (2);

  InternetStackHelper stack;
  stack.Install (topo.senders);
  stack.Install (topo.receivers);
  stack.Install (topo.routers);

  PointToPointHelper access;
  access.SetDeviceAttribute ("DataRate", StringValue ("100Mbps"));
  access.SetChannelAttribute ("Delay", StringValue ("0ms"));

  PointToPointHelper bottle;
  bottle.SetDeviceAttribute ("DataRate", StringValue (bottleRate));
  int bottleDelayMs = std::max (1,
      static_cast<int> (targetRtt * 1000.0 / 2.0));
  bottle.SetChannelAttribute ("Delay",
      StringValue (std::to_string (bottleDelayMs) + "ms"));

  // sender0 ↔ r0
  NetDeviceContainer devS0R0 =
      access.Install (topo.senders.Get (0), topo.routers.Get (0));

  // sender1 ↔ r0
  NetDeviceContainer devS1R0 =
      access.Install (topo.senders.Get (1), topo.routers.Get (0));

  // r0 ↔ r1 bottleneck
  NetDeviceContainer devBottle =
      bottle.Install (topo.routers.Get (0), topo.routers.Get (1));

  // r1 ↔ receiver0
  NetDeviceContainer devR1D0 =
      access.Install (topo.routers.Get (1), topo.receivers.Get (0));

  // r1 ↔ receiver1
  NetDeviceContainer devR1D1 =
      access.Install (topo.routers.Get (1), topo.receivers.Get (1));

  Ipv4AddressHelper address;

  // sender0 - r0
  address.SetBase ("10.2.1.0", "255.255.255.0");
  address.Assign (devS0R0);

  // sender1 - r0
  address.SetBase ("10.2.2.0", "255.255.255.0");
  address.Assign (devS1R0);

  // r0 - r1 bottleneck
  address.SetBase ("10.2.3.0", "255.255.255.0");
  address.Assign (devBottle);

  // r1 - receiver0
  Ipv4InterfaceContainer ifR1D0;
  address.SetBase ("10.2.4.0", "255.255.255.0");
  ifR1D0 = address.Assign (devR1D0);

  // r1 - receiver1
  Ipv4InterfaceContainer ifR1D1;
  address.SetBase ("10.2.5.0", "255.255.255.0");
  ifR1D1 = address.Assign (devR1D1);

  topo.recvIp0 = ifR1D0.GetAddress (1);
  topo.recvIp1 = ifR1D1.GetAddress (1);

  Ipv4GlobalRoutingHelper::PopulateRoutingTables ();

  return topo;
}

// ------------------ Application installation ------------------

Ptr<PacketSink>
InstallSink (Ptr<Node> node, uint16_t port)
{
  Address local (InetSocketAddress (Ipv4Address::GetAny (), port));
  PacketSinkHelper sinkHelper ("ns3::TcpSocketFactory", local);
  ApplicationContainer apps = sinkHelper.Install (node);
  apps.Start (Seconds (START_TIME));
  apps.Stop (Seconds (SIM_TIME));
  return apps.Get (0)->GetObject<PacketSink> ();
}

ApplicationContainer
InstallBulk (Ptr<Node> node, Ipv4Address dstIp, uint16_t port)
{
  Address remote (InetSocketAddress (dstIp, port));
  BulkSendHelper bulk ("ns3::TcpSocketFactory", remote);
  bulk.SetAttribute ("MaxBytes", UintegerValue (0));        // Infinite-length flow
  bulk.SetAttribute ("SendSize", UintegerValue (PKT_SIZE)); // Packet size 1500 B

  ApplicationContainer apps = bulk.Install (node);
  apps.Start (Seconds (APP_START));
  apps.Stop (Seconds (SIM_TIME));
  return apps;
}

// ------------------ Convergence monitoring (ε-convergence time) ------------------

class ConvergenceMonitor
{
public:
  /**
   * @param sink           PacketSink to monitor
   * @param stream         Output file (time, throughput)
   * @param sampleInterval Sampling interval Δt (s)
   * @param measureStart   Time to start recording samples (s)
   * @param tailWindow     Use the last tailWindow seconds to estimate long-run average T_bar
   * @param eps            ε, e.g. 0.8
   * @param windowLen      Length of short-term averaging window (s)
   */
  ConvergenceMonitor (Ptr<PacketSink> sink,
                      Ptr<OutputStreamWrapper> stream,
                      double sampleInterval,
                      double measureStart,
                      double tailWindow,
                      double eps,
                      double windowLen)
    : m_sink (sink),
      m_stream (stream),
      m_sampleInterval (sampleInterval),
      m_measureStart (measureStart),
      m_tailWindow (tailWindow),
      m_eps (eps),
      m_windowLen (windowLen),
      m_initialized (false),
      m_lastTotalRx (0)
  {
    *m_stream->GetStream () << "# time_sec, throughput_Mbps" << std::endl;
  }

  void
  Start ()
  {
    Simulator::Schedule (Seconds (m_sampleInterval),
                         &ConvergenceMonitor::Sample, this);
  }

  void
  Finish ()
  {
    ComputeStatsAndPrint ();
  }

private:
  void
  Sample ()
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
        *m_stream->GetStream () << now << "," << mbps << std::endl;
        m_times.push_back (now);
        m_mbps.push_back (mbps);
      }

    Simulator::Schedule (Seconds (m_sampleInterval),
                         &ConvergenceMonitor::Sample, this);
  }

  void
  ComputeStatsAndPrint ()
  {
    std::ostream &o = *m_stream->GetStream ();

    if (m_times.empty ())
      {
      o << "# No samples collected" << std::endl;
      return;
      }

    // ---- 1) Estimate long-run average T_bar ----
    double tailStartLast = SIM_TIME - m_tailWindow;
    if (tailStartLast < m_measureStart)
      {
        tailStartLast = m_measureStart;
      }

    double tailStartPrev = SIM_TIME - 2.0 * m_tailWindow;
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
    double Tbar   = T_last;

    double relDiff = 0.0;
    bool   approxStable = false;
    if (T_last > 0.0 && countPrev > 0)
      {
        // Within 5% difference, consider approximately stable
        relDiff = std::fabs (T_last - T_prev) / T_last;
        approxStable = (relDiff < 0.05);
      }

    // ---- 2) Use sliding window to compute T_conv(ε) ----
    double Tconv = -1.0;

    if (Tbar > 0.0)
      {
        double threshold = m_eps * Tbar;
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
                Tconv = m_times[i];
                break;
              }
          }
      }

    double TconvRel = (Tconv >= 0.0 ? Tconv - APP_START : -1.0);

    // ---- 3) Output ----
    o << "# T_bar_Mbps=" << Tbar
      << ", T_prev_Mbps=" << T_prev
      << ", T_last_Mbps=" << T_last
      << ", relDiff=" << relDiff
      << ", approxStable=" << (approxStable ? 1 : 0)
      << ", eps=" << m_eps
      << ", windowLen_s=" << m_windowLen
      << ", Tconv_abs_s=" << Tconv
      << ", Tconv_rel_s=" << TconvRel
      << std::endl;

    std::cout << "[Convergence] T_bar=" << Tbar
              << " Mbps, relDiff=" << relDiff
              << " (approxStable=" << (approxStable ? 1 : 0)
              << "), Tconv_abs=" << Tconv
              << " s, Tconv_rel=" << TconvRel
              << " s since flow start"
              << std::endl;
  }

private:
  Ptr<PacketSink>          m_sink;
  Ptr<OutputStreamWrapper> m_stream;

  double m_sampleInterval;
  double m_measureStart;
  double m_tailWindow;
  double m_eps;
  double m_windowLen;

  bool      m_initialized;
  uint64_t  m_lastTotalRx;

  std::vector<double> m_times;
  std::vector<double> m_mbps;
};

// ------------------ main: first single-flow, then two-flow competition ------------------

int
main (int argc, char *argv[])
{
  CommandLine cmd;
  cmd.Parse (argc, argv);

  EnsureResultsDir ();

  // ============ Part 1: Single-flow experiments (Vegas / Cubic) ============

  struct CcaCfg
  {
    std::string shortName; // "vegas" / "cubic"
    TypeId      typeId;
  };

  std::vector<CcaCfg> ccas = {
    {"vegas", TcpVegas::GetTypeId ()},
    {"cubic", TcpCubic::GetTypeId ()}
  };

  for (const auto &cfg : ccas)
    {
      std::cout << "\n============================\n";
      std::cout << "  Single-flow: Tcp" << cfg.shortName << "\n";
      std::cout << "============================\n";

      // Set global default CCA
      Config::SetDefault ("ns3::TcpL4Protocol::SocketType",
                          TypeIdValue (cfg.typeId));

      for (double rtt : RTT_LIST)
        {
          std::cout << "\n--- RTT = " << rtt * 1000.0 << " ms ---\n";

          // 1) Build single-flow topology
          SingleTopology topo = BuildSingleDumbbell (DATA_RATE, rtt);

          uint16_t port = 5000;

          // 2) Install sink and BulkSend
          Ptr<PacketSink> sink =
              InstallSink (topo.receivers.Get (0), port);

          ApplicationContainer app =
              InstallBulk (topo.senders.Get (0), topo.recvIp, port);

          // 3) Output file
          AsciiTraceHelper ascii;
          std::ostringstream oss;
          int rttMs = int (rtt * 1000.0);

          oss << RESULTS_FOLDER << "/throughput_"
              << cfg.shortName << "_single_rtt_"
              << rttMs << "ms.csv";

          Ptr<OutputStreamWrapper> stream =
              ascii.CreateFileStream (oss.str ());

          // 4) Configure ε-convergence monitor
          double sampleInterval = 0.1;   // 0.1s
          double measureStart   = 1.0;   // Start recording as soon as flow starts
          double tailWindow     = 20.0;  // Use last 20s to estimate T_bar
          double eps            = 0.8;   // 80% of T_bar
          double windowLen      = 2.0;   // 2s sliding window

          ConvergenceMonitor mon (sink, stream,
                                  sampleInterval,
                                  measureStart,
                                  tailWindow,
                                  eps,
                                  windowLen);
          mon.Start ();

          // 5) Run simulation
          Simulator::Stop (Seconds (SIM_TIME));
          Simulator::Run ();

          mon.Finish ();
          Simulator::Destroy ();

          std::cout << "Finished single-flow CCA="
                    << cfg.shortName
                    << ", RTT=" << rttMs << " ms\n";
        }
    }

  // ============ Part 2: Two-flow competition (Flow1=Cubic, Flow2=Vegas) ============

  for (double rtt : RTT_LIST)
    {
      std::cout << "\n=======================================\n";
      std::cout << "  Two-flow competition, RTT = "
                << rtt * 1000.0 << " ms\n";
      std::cout << "   Flow 1: Cubic, Flow 2: Vegas\n";
      std::cout << "=======================================\n";

      // 1) Build two-flow topology
      MultiTopology topo = BuildMultiDumbbell (DATA_RATE, rtt);

      // 2) Set different CCAs for the two senders
      Ptr<TcpL4Protocol> tcp0 =
          topo.senders.Get (0)->GetObject<TcpL4Protocol> ();
      Ptr<TcpL4Protocol> tcp1 =
          topo.senders.Get (1)->GetObject<TcpL4Protocol> ();

      tcp0->SetAttribute ("SocketType",
          TypeIdValue (TcpCubic::GetTypeId ())); // Flow1 = Cubic
      tcp1->SetAttribute ("SocketType",
          TypeIdValue (TcpVegas::GetTypeId ())); // Flow2 = Vegas

      uint16_t port1 = 5000; // Cubic flow
      uint16_t port2 = 5001; // Vegas flow

      // 3) Sinks
      Ptr<PacketSink> sink1 =
          InstallSink (topo.receivers.Get (0), port1);
      Ptr<PacketSink> sink2 =
          InstallSink (topo.receivers.Get (1), port2);

      // 4) BulkSend apps
      ApplicationContainer app1 =
          InstallBulk (topo.senders.Get (0), topo.recvIp0, port1);
      ApplicationContainer app2 =
          InstallBulk (topo.senders.Get (1), topo.recvIp1, port2);

      // 5) Output files
      AsciiTraceHelper ascii;
      std::ostringstream oss1, oss2;
      int rttMs = int (rtt * 1000.0);

      oss1 << RESULTS_FOLDER << "/throughput_cubic_compete_rtt_"
           << rttMs << "ms.csv";
      oss2 << RESULTS_FOLDER << "/throughput_vegas_compete_rtt_"
           << rttMs << "ms.csv";

      Ptr<OutputStreamWrapper> stream1 =
          ascii.CreateFileStream (oss1.str ());
      Ptr<OutputStreamWrapper> stream2 =
          ascii.CreateFileStream (oss2.str ());

      // 6) ε-convergence monitors
      double sampleInterval = 0.1;
      double measureStart   = 1.0;
      double tailWindow     = 20.0;
      double eps            = 0.8;
      double windowLen      = 2.0;

      ConvergenceMonitor mon1 (sink1, stream1,
                               sampleInterval,
                               measureStart,
                               tailWindow,
                               eps,
                               windowLen);

      ConvergenceMonitor mon2 (sink2, stream2,
                               sampleInterval,
                               measureStart,
                              tailWindow,
                               eps,
                               windowLen);

      mon1.Start ();
      mon2.Start ();

      // 7) Run simulation
      Simulator::Stop (Seconds (SIM_TIME));
      Simulator::Run ();

      mon1.Finish ();
      mon2.Finish ();

      Simulator::Destroy ();

      std::cout << "Finished two-flow competition, RTT="
                << rttMs << " ms\n";
    }

  return 0;
}
