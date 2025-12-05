#include "ns3/core-module.h"
#include "ns3/internet-module.h"
#include "ns3/network-module.h"
#include "ns3/applications-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/csma-module.h"
#include "ns3/flow-monitor-module.h"
#include <iostream>

using namespace ns3;


// #
// # Experiment goals:
// #   (A) Single-flow experiment
// #       - Topology: dumbbell
// #       - One long-lived flow: n0 -> n2
// #       - CCAs: CUBIC / Vegas (run separately)
// #       - Variable: RTT = 10 / 20 / 40 ms
// #
// #   (B) Two-flows competition experiment
// #       - Same dumbbell topology
// #       - Flow 1: CUBIC (n0 -> n2), starts at t = 1s
// #       - Flow 2: Vegas (n1 -> n3), joins at t = 5s
// #       - Both flows share the same RTT (same bottleneck)
// #       - Variable: RTT = 10 / 20 / 40 ms
// #
// # Other assumptions:
// #   - Bottleneck bandwidth fixed at 25Mbps
// #   - buffer = 1 × BDP
// #   - Simulation length = 40s, applications start from t = 1s

static const uint32_t PKT_SIZE = 1500;
static const double SIM_TIME = 120.0;
static const double START_TIME = 1.0;
static const double APP_START = 2.0;
static const std::string DATA_RATE = "20Mbps";
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

// static void CwndTracer(std::string context ,uint32_t oldval, uint32_t newval)
// {
//     std::cout << "Moving cwnd from " << oldval << " to " << newval << std::endl;
//     return;
// }

// /**
//  * Get the Node Id From Context.
//  *
//  * @param context The context.
//  * @return the node ID.
//  */
// static uint32_t
// GetNodeIdFromContext(std::string context)
// {
//     const std::size_t n1 = context.find_first_of('/', 1);
//     const std::size_t n2 = context.find_first_of('/', n1 + 1);
//     return std::stoul(context.substr(n1 + 1, n2 - n1 - 1));
// }

/**
 * Congestion window tracer.
 *
 * @param context The context.
 * @param oldval Old value.
 * @param newval New value.
 */
// static void
// CwndTracer2(std::string context, uint32_t oldval, uint32_t newval)
// {
//     uint32_t nodeId = GetNodeIdFromContext(context);

//     if (firstCwnd[nodeId])
//     {
//         *cWndStream[nodeId]->GetStream() << "0.0 " << oldval << std::endl;
//         firstCwnd[nodeId] = false;
//     }
//     *cWndStream[nodeId]->GetStream() << Simulator::Now().GetSeconds() << " " << newval << std::endl;
//     cWndValue[nodeId] = newval;

//     if (!firstSshThr[nodeId])
//     {
//         *ssThreshStream[nodeId]->GetStream()
//             << Simulator::Now().GetSeconds() << " " << ssThreshValue[nodeId] << std::endl;
//     }
// }

/**
//  * Congestion window trace connection.
//  *
//  * @param cwnd_tr_file_name Congestion window trace file name.
//  * @param nodeId Node ID.
//  */
// static void
// TraceCwnd(std::string cwnd_tr_file_name, uint32_t nodeId)
// {
//     AsciiTraceHelper ascii;
//     cWndStream[nodeId] = ascii.CreateFileStream(cwnd_tr_file_name);
//     Config::Connect("/NodeList/" + std::to_string(nodeId) +
//                         "/$ns3::TcpL4Protocol/SocketList/0/CongestionWindow",
//                     MakeCallback(&CwndTracer));
// }

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



// BDP (in packets)
// (Compute BDP in units of packets.
int CalcBdpPkts(const std::string &rate, double rtt)
{
    DataRate dr(rate);
    double bits = dr.GetBitRate() * rtt;
    double bytes = bits / 8.0;
    return bytes / PKT_SIZE;
}

// -----------------------------------------------------------------------------
//                         Build Dumbbell Topology (your diagram)
//              _      1 rtt           1 rtt      _
//             |    n0------+        +------n3     |
//             |            | 3 rtt  |
//     Senders |    n1------R1------R2-----n4      |    Receivers
//             |            |        |             |
//             |_   n2------+        +------n5    _|
// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
//                             Single-flow experiment
// -----------------------------------------------------------------------------
// """
// (A) Single-flow experiment:
//     - dumbbell topology
//     - one long-lived flow: n0 -> n2
//     - all nodes use same CCA (Cubic/Vegas)
//     - variable: RTT
// """
void RunSingleFlow(double rtt, std::string ccaName, std::string ccaType)
{
    // Set system-wide CCA
    Config::SetDefault("ns3::TcpL4Protocol::SocketType",
                       TypeIdValue(TypeId::LookupByName(ccaType)));

    int qpkts = CalcBdpPkts(DATA_RATE, rtt) * BUF_FACTOR;
    Config::SetDefault("ns3::DropTailQueue<Packet>::MaxSize",
                       StringValue(std::to_string(qpkts) + "p"));

    Topology topo = BuildDumbbell(DATA_RATE, rtt);

    Ptr<PacketSink> sink =
        InstallSink(topo.receivers.Get(0), topo.ip_n2, 8080);

    // Track throughput per second
    AsciiTraceHelper ascii;
    Ptr<OutputStreamWrapper> stream =
                        ascii.CreateFileStream(RESULTS_FOLDER+"single_" + ccaName +
                        "_rtt" + std::to_string(int(rtt * 1000)) +
                        "ms_throughput.csv");

    *stream->GetStream() << "time,throughput_mbps" << std::endl;

    uint64_t lastTotalRx = 0;
    Simulator::Schedule(Seconds(START_TIME), &ThroughputTracer, sink, std::ref(lastTotalRx),stream);

    InstallBulk(topo.senders.Get(0), topo.ip_n2, 8080);

    // Ptr<TcpL4Protocol> tcp = topo.senders.Get(0)->GetObject<TcpL4Protocol>();
    // Ptr<Socket> sock = tcp->FindSocket(InetSocketAddress(topo.ip_n2, 8080));
    // Ptr<TcpSocketBase> tcpSock = DynamicCast<TcpSocketBase>(sock);

    // tcpSock->TraceConnectWithoutContext(
    //     "CongestionWindow",
    //     MakeCallback(&CwndTracer));


    // int index = 0;
    // std::string prefix_file_name = "./results/cc_rtt_single_" + ccaName + "_rtt" + std::to_string(int(rtt * 1000)) + "ms";
    // Simulator::Schedule(Seconds(START_TIME + 0.00001),
    //                         &TraceCwnd,
    //                         prefix_file_name + "-cwnd.data",
    //                         index + 1);

    Simulator::Stop(Seconds(SIM_TIME));
    Simulator::Run();

    double rx = sink->GetTotalRx();
    double thr = (rx * 8.0) / (SIM_TIME - APP_START) / 1e6;
    std::cout << "[Single] CCA=" << ccaName
              << " RTT=" << int(rtt * 1000)
              << "ms Thr=" << thr << " Mbps\n";

    Simulator::Destroy();
}

// -----------------------------------------------------------------------------
//                             Competing-flows experiment
// -----------------------------------------------------------------------------
void RunCompeting(double rtt)
{
    int qpkts = CalcBdpPkts(DATA_RATE, rtt) * BUF_FACTOR;
    Config::SetDefault("ns3::DropTailQueue<Packet>::MaxSize",
                       StringValue(std::to_string(qpkts) + "p"));

    Topology topo = BuildDumbbell(DATA_RATE, rtt);

    // Assign CCAs to nodes
    Config::Set("/NodeList/0/$ns3::TcpL4Protocol/SocketType",
                TypeIdValue(TypeId::LookupByName("ns3::TcpCubic")));
    Config::Set("/NodeList/1/$ns3::TcpL4Protocol/SocketType",
                TypeIdValue(TypeId::LookupByName("ns3::TcpVegas")));

    Ptr<PacketSink> sink1 =
        InstallSink(topo.receivers.Get(0), topo.ip_n2, 8080);
    InstallBulk(topo.senders.Get(0), topo.ip_n2, 8080);


    Ptr<PacketSink> sink2 =
        InstallSink(topo.receivers.Get(1), topo.ip_n3, 8081);
    ApplicationContainer bulk2 =
        InstallBulk(topo.senders.Get(1), topo.ip_n3, 8081);
    bulk2.Get(0)->SetStartTime(Seconds(5.0)); // starts later
    
    // Track throughput per second for both flows
    uint64_t lastRx2 = 0;
    uint64_t lastRx1 = 0;

    AsciiTraceHelper ascii1;
    Ptr<OutputStreamWrapper> stream1 =
                    ascii1.CreateFileStream(RESULTS_FOLDER+"compete_Cubic_rtt" + 
                    std::to_string(int(rtt * 1000)) +
                    "ms_throughput.csv");
    *stream1->GetStream() << "time,throughput_mbps" << std::endl;

    AsciiTraceHelper ascii2;
    Ptr<OutputStreamWrapper> stream2 =
        ascii2.CreateFileStream(RESULTS_FOLDER+"compete_Vegas_rtt" +
                                std::to_string(int(rtt * 1000)) +
                                "ms_throughput.csv");

    *stream2->GetStream() << "time,throughput_mbps" << std::endl;


    Simulator::Schedule(Seconds(START_TIME+0.0001),
                    &ThroughputTracer,
                    sink1, std::ref(lastRx1), stream1);
    Simulator::Schedule(Seconds(START_TIME+0.0001),
                    &ThroughputTracer,
                    sink2, std::ref(lastRx2), stream2);


    Simulator::Stop(Seconds(SIM_TIME));
    Simulator::Run();

    double rx1 = sink1->GetTotalRx();
    double rx2 = sink2->GetTotalRx();
    double thr1 = (rx1 * 8.0) / (SIM_TIME - APP_START) / 1e6;
    double thr2 = (rx2 * 8.0) / (SIM_TIME - APP_START) / 1e6;

    std::cout << "[Compete] RTT=" << int(rtt * 1000)
              << "ms  Cubic=" << thr1 << " Mbps"
              << "  Vegas=" << thr2 << " Mbps\n";

    Simulator::Destroy();
}

// -----------------------------------------------------------------------------
//                                       MAIN
// -----------------------------------------------------------------------------
int main()
{
    //LogComponentEnable("TcpL4Protocol", LOG_LEVEL_FUNCTION);
    LogComponentEnable("TcpL4Protocol", LOG_LEVEL_WARN);

    std::cout << "=== (A) Single-flow experiments ===\n";

    for (double rtt : RTT_LIST)
    {
        RunSingleFlow(rtt, "Cubic", "ns3::TcpCubic");
        RunSingleFlow(rtt, "Vegas", "ns3::TcpVegas");
    }


    std::cout << "\n=== (B) Competition experiments ===\n";

    for (double rtt : RTT_LIST)
        RunCompeting(rtt);

    return 0;
}