#include "ns3/core-module.h"
#include "ns3/internet-module.h"
#include "ns3/network-module.h"
#include "ns3/applications-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/csma-module.h"
#include "ns3/flow-monitor-module.h"
#include <iostream>
#include <vector>
#include <limits>

using namespace ns3;

static const uint32_t PKT_SIZE = 1500;
static const double SIM_TIME = 120.0;
static const double START_TIME = 1.0;
static const double APP_START = 2.0;
static const std::string DATA_RATE = "25Mbps";
static const double BUF_FACTOR = 2.0;
static const double RTT_LIST[3] = {0.002, 0.020, 0.040}; // in seconds
static const std::string RESULTS_FOLDER = "./results/";

static std::map<uint32_t, bool> firstCwnd;                      //!< First congestion window.
static std::map<uint32_t, bool> firstSshThr;                    //!< First SlowStart threshold.
static std::map<uint32_t, Ptr<OutputStreamWrapper>> cWndStream; //!< Congstion window output stream.
static std::map<uint32_t, uint32_t> cWndValue;                      //!< congestion window value.
static std::map<uint32_t, Ptr<OutputStreamWrapper>>
    ssThreshStream; //!< SlowStart threshold output stream.
static std::map<uint32_t, uint32_t> ssThreshValue;                  //!< SlowStart threshold value.



Ptr<PacketSink> InstallSink(Ptr<Node> node, Ipv4Address addr, uint16_t port)
{
    Address sinkAddr = InetSocketAddress(addr, port);
    PacketSinkHelper helper("ns3::TcpSocketFactory", sinkAddr);
    ApplicationContainer apps = helper.Install(node);
    apps.Start(Seconds(START_TIME));
    apps.Stop(Seconds(SIM_TIME));
    return apps.Get(0)->GetObject<PacketSink>();
}

ApplicationContainer InstallBulk(Ptr<Node> node, Ipv4Address addr, uint16_t port)
{
    Address sinkAddr = InetSocketAddress(addr, port);
    BulkSendHelper bulk("ns3::TcpSocketFactory", sinkAddr);
    bulk.SetAttribute("MaxBytes", UintegerValue(0));
    bulk.SetAttribute("SendSize", UintegerValue(PKT_SIZE));
    ApplicationContainer apps = bulk.Install(node);
    apps.Start(Seconds(APP_START));
    apps.Stop(Seconds(SIM_TIME));
    return apps;
}

void
ThroughputTracer (Ptr<PacketSink> sink, uint64_t &lastTotalRx, Ptr<OutputStreamWrapper> stream)
{
    double now = Simulator::Now().GetSeconds ();

    uint64_t cur = sink->GetTotalRx ();
    double mbps = (cur - lastTotalRx) * 8.0 / 1e6;   // per-second Mbps
    lastTotalRx = cur;

    //std::cout << "t=" << now << "s  throughput=" << mbps << " Mbps" << std::endl;
    *stream->GetStream() << now << "," << mbps << std::endl;

    Simulator::Schedule (Seconds(1.0), &ThroughputTracer, sink, std::ref(lastTotalRx), stream);
}

struct Topology
{
    NodeContainer senders, routers, receivers;
    Ipv4Address ip_n2, ip_n3;
};

Topology BuildDumbbell(std::string rate, double rtt)
{
    Topology topo;

    topo.senders.Create(3);   // n0, n1, n2
    topo.routers.Create(2);   // R0, R1
    topo.receivers.Create(3); // n3, n4, n5

    InternetStackHelper stack;
    stack.Install(topo.senders);
    stack.Install(topo.receivers);
    stack.Install(topo.routers);

    // Divide RTT: to 3 parts on the bottleneck link and 1 part on each side
    double rtt_ptp = rtt/5*3 ;
    double rtt_csma = rtt/5*1;
    std::cout << "rtt_ptp=" << rtt_ptp*500 << "ms rtt_csma=" << rtt_csma*500 << "ms" << std::endl;

    // Non-bottleneck CSMA links
    CsmaHelper csma;
    csma.SetChannelAttribute("DataRate", StringValue(rate));
    csma.SetChannelAttribute("Delay", StringValue(std::to_string(int(rtt_csma * 500)) + "ms"));

    // Bottleneck P2P link
    PointToPointHelper bottleneck;
    bottleneck.SetDeviceAttribute("DataRate", StringValue(rate));
    bottleneck.SetChannelAttribute("Delay",
                                   StringValue(std::to_string(int(rtt_ptp * 500)) + "ms"));

    Ipv4AddressHelper addr;

    // senders → R0
    NodeContainer left;
    left.Add(topo.routers.Get(0));
    left.Add(topo.senders);
    NetDeviceContainer dev0 = csma.Install(left);
    addr.SetBase("10.1.1.0", "255.255.255.0");
    addr.Assign(dev0);

    // R1 → receivers
    NodeContainer right;
    right.Add(topo.routers.Get(1));
    right.Add(topo.receivers);
    NetDeviceContainer dev1 = csma.Install(right);
    addr.SetBase("10.1.3.0", "255.255.255.0");
    Ipv4InterfaceContainer recIf = addr.Assign(dev1);

    topo.ip_n2 = recIf.GetAddress(1);
    topo.ip_n3 = recIf.GetAddress(2);

    // R0 ↔ R1 bottleneck
    addr.SetBase("10.1.2.0", "255.255.255.0");
    NetDeviceContainer devb = bottleneck.Install(topo.routers.Get(0), topo.routers.Get(1));
    addr.Assign(devb);

    Ipv4GlobalRoutingHelper::PopulateRoutingTables();
    return topo;
}

class ConvergenceMonitor
{
public:
  ConvergenceMonitor (Ptr<PacketSink> sink,
                      Ptr<OutputStreamWrapper> stream,
                      std::vector<double> boundsMbps,
                      double sampleIntervalSec,
                      double measureStartSec,
                      double continuousWindowSec)
    : m_sink(sink),
      m_stream(stream),
      m_bounds(boundsMbps),
      m_sampleInterval(sampleIntervalSec),
      m_measureStart(measureStartSec),
      m_continuousWindow(continuousWindowSec),
      m_initialized(false)
  {
    size_t n = m_bounds.size();
    m_reached.assign(n, false);
    m_firstReach.assign(n, std::numeric_limits<double>::infinity());
    m_timeAbove.assign(n, 0.0);
    m_inside.assign(n, false);
    m_contStart.assign(n, -1.0);
    m_settleTime.assign(n, std::numeric_limits<double>::infinity());
  }

  void
  Start ()
  {
    // initialize lastTotalRx on first sample inside Sample()
    Simulator::Schedule (Seconds(m_sampleInterval), &ConvergenceMonitor::Sample, this);
    // schedule a print at the end (adjust if SIM_TIME variable name differs)
    Simulator::Schedule (Seconds(SIM_TIME + 0.001), &ConvergenceMonitor::PrintResults, this);
  }

  void
  Sample ()
  {
    double now = Simulator::Now ().GetSeconds ();
    uint64_t cur = m_sink->GetTotalRx ();

    if (!m_initialized)
      {
        // set baseline and don't compute throughput until next sample
        m_lastTotalRx = cur;
        m_initialized = true;
        Simulator::Schedule (Seconds(m_sampleInterval), &ConvergenceMonitor::Sample, this);
        return;
      }

    // compute Mbps over the sample interval
    double mbps = (double)(cur - m_lastTotalRx) * 8.0 / (m_sampleInterval * 1e6);
    m_lastTotalRx = cur;

    // optionally write per-sample trace
    if (now >= m_measureStart)
      {
        *m_stream->GetStream () << now << "," << mbps << std::endl;
      }

    // update per-bound stats (only count time after measureStart)
    if (now >= m_measureStart)
      {
        for (size_t i = 0; i < m_bounds.size (); ++i)
          {
            double B = m_bounds[i];

            // mark first reach time
            if (!m_reached[i] && mbps >= B)
              {
                m_reached[i] = true;
                m_firstReach[i] = now;
              }

            // accumulate total time above bound (sampleInterval resolution)
            if (mbps >= B)
              {
                m_timeAbove[i] += m_sampleInterval;
              }

            // continuous (settling) detection
            if (mbps >= B)
              {
                if (m_contStart[i] < 0) m_contStart[i] = now; // start new continuous segment
                double contDur = now - m_contStart[i];
                if (contDur >= m_continuousWindow && m_settleTime[i] == std::numeric_limits<double>::infinity ())
                  {
                    // first time we have a continuous window above B
                    m_settleTime[i] = m_contStart[i] + m_continuousWindow;
                  }
              }
            else
              {
                // break continuous segment
                m_contStart[i] = -1.0;
              }
          }
      }

    // schedule next sample
    Simulator::Schedule (Seconds(m_sampleInterval), &ConvergenceMonitor::Sample, this);
  }

  void
  PrintResults ()
  {
    double measDur = SIM_TIME - m_measureStart;
    if (measDur < 0) measDur = 0.0;

    std::ostream &o = *m_stream->GetStream ();

    o << "# Convergence summary (measure start = " << m_measureStart
      << " s, sampleInterval = " << m_sampleInterval << " s, continuousWindow = "
      << m_continuousWindow << " s)" << std::endl;

    o << "bound_mbps,first_reach_s,settle_time_s,total_time_above_s,percent_time_above" << std::endl;

    for (size_t i = 0; i < m_bounds.size (); ++i)
      {
        double firstReach = (m_firstReach[i] == std::numeric_limits<double>::infinity()) ? -1.0 : m_firstReach[i];
        double settle = (m_settleTime[i] == std::numeric_limits<double>::infinity()) ? -1.0 : m_settleTime[i];
        double timeAbove = m_timeAbove[i];
        double pct = (measDur > 0.0 ? (timeAbove / measDur * 100.0) : 0.0);

        o << m_bounds[i] << "," << firstReach << "," << settle << "," << timeAbove << "," << pct << std::endl;
      }
  }

private:
  Ptr<PacketSink> m_sink;
  Ptr<OutputStreamWrapper> m_stream;
  std::vector<double> m_bounds;

  double m_sampleInterval;
  double m_measureStart;
  double m_continuousWindow;

  // state
  bool m_initialized;
  uint64_t m_lastTotalRx;

  std::vector<bool> m_reached;
  std::vector<double> m_firstReach;
  std::vector<double> m_timeAbove;
  std::vector<bool> m_inside;
  std::vector<double> m_contStart;   // -1 if not inside a continuous segment
  std::vector<double> m_settleTime;
};

std::vector<double> bounds = {5.0, 10.0, 15.0, 20.0, 25.0}; // Mbps - pick your 5 bounds
double sampleInterval = 0.1;     // 100 ms
double measureStart = 5.0;       // ignore first 5s warm-up
double continuousWindow = 2.0;   // require 2s continuous above bound to consider "settled"

int main(int argc, char *argv[])
{
    // CCAs to test
    std::vector<std::string> ccas = {
        "ns3::TcpVegas",
        "ns3::TcpCubic"
    };

    // Loop over CCAs
    for (const auto &cca : ccas)
    {
        std::cout << "\n============================" << std::endl;
        std::cout << " Testing CCA: " << cca << std::endl;
        std::cout << "============================\n" << std::endl;

        // Set CCA
        Config::SetDefault("ns3::TcpL4Protocol::SocketType", StringValue(cca));

        // Loop over RTT values
        for (double rtt : RTT_LIST)
        {
            std::cout << "\n=== Running RTT = "
                      << rtt * 1000 << " ms for " << cca << " ==="
                      << std::endl;

            // Build topology
            Topology topo = BuildDumbbell(DATA_RATE, rtt);

            uint16_t port = 5000;

            // Install sink (receiver side)
            Ptr<PacketSink> sink = InstallSink(
                topo.receivers.Get(0),
                topo.ip_n2,
                port
            );

            // Install sender (left side)
            ApplicationContainer senderApp = InstallBulk(
                topo.senders.Get(0),
                topo.ip_n2,
                port
            );

            // Create directory + output filename
            std::string filename =
                RESULTS_FOLDER +
                "conv_" +
                (cca.find("Vegas") != std::string::npos ? "vegas" : "cubic") +
                "_rtt_" + std::to_string(int(rtt * 1000)) + "ms.csv";

            AsciiTraceHelper ascii;
            Ptr<OutputStreamWrapper> stream = ascii.CreateFileStream(filename);

            // Define convergence bounds in Mbps
            std::vector<double> bounds = {5.0, 10.0, 15.0, 20.0, 25.0};

            double sampleInterval   = 0.1;  // seconds
            double measureStart     = 5.0;  // ignore transient at start
            double continuousWindow = 2.0;  // required stability window

            // Create convergence monitor
            ConvergenceMonitor *mon = new ConvergenceMonitor(
                sink,
                stream,
                bounds,
                sampleInterval,
                measureStart,
                continuousWindow
            );

            mon->Start();

            // Run simulation
            Simulator::Stop(Seconds(SIM_TIME + 0.1));
            Simulator::Run();
            Simulator::Destroy();

            std::cout << "Finished RTT " << rtt * 1000 << " ms for " << cca << "\n";
        }
    }

    return 0;
}
