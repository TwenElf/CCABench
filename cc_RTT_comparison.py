# cc_comparison.py
#
# Experiment goals:
#   (A) Single-flow experiment
#       - Topology: dumbbell
#       - One long-lived flow: n0 -> n2
#       - CCAs: CUBIC / Vegas (run separately)
#       - Variable: RTT = 10 / 20 / 40 ms
#
#   (B) Two-flows competition experiment
#       - Same dumbbell topology
#       - Flow 1: CUBIC (n0 -> n2), starts at t = 1s
#       - Flow 2: Vegas (n1 -> n3), joins at t = 5s
#       - Both flows share the same RTT (same bottleneck)
#       - Variable: RTT = 10 / 20 / 40 ms
#
# Other assumptions:
#   - Bottleneck bandwidth fixed at 25Mbps
#   - buffer = 1 × BDP
#   - Simulation length = 40s, applications start from t = 1s

try:
    from ns import ns
except ModuleNotFoundError:
    raise SystemExit(
        "Error: ns3 Python module not found;"
        " Python bindings may not be enabled"
        " or your PYTHONPATH might not be properly configured"
    )
from ctypes import c_bool, c_int

verbose = c_bool(True)
if verbose.value:
    ns.LogComponentEnable("UdpEchoClientApplication", ns.LOG_LEVEL_INFO)
    ns.LogComponentEnable("UdpEchoServerApplication", ns.LOG_LEVEL_INFO)

# --------- Global parameters -------------
PKT_SIZE = 1500         # bytes
SIM_TIME = 40.0         # simulation time (s)
APP_START = 1.0         # application start time (s)

DATA_RATE = "25Mbps"    # bottleneck bandwidth

BUF_FACTOR = 1.0        # buffer size = BDP × this factor

RTT_LIST = [0.010, 0.020, 0.040]  # seconds

CCA_TYPES_SINGLE = [
    ("Cubic", "ns3::TcpCubic"),
    ("Vegas", "ns3::TcpVegas"),
]


def set_tcp_cca(cca_type_str: str):
    """Set global default TCP CCA."""
    tid = ns.TypeId.LookupByName(cca_type_str)
    ns.Config.SetDefault("ns3::TcpL4Protocol::SocketType", ns.TypeIdValue(tid))


def set_node_tcp_cca(node_index: int, cca_type_str: str):
    """Set TCP CCA for one specific node (used in competition experiment)."""
    tid = ns.TypeId.LookupByName(cca_type_str)
    path = f"/NodeList/{node_index}/$ns3::TcpL4Protocol/SocketType"
    ns.Config.Set(path, ns.TypeIdValue(tid))


def calc_bdp_pkts(data_rate_str: str, rtt_sec: float, pkt_size: int) -> float:
    """Compute BDP in units of packets."""
    dr = ns.DataRate(data_rate_str)
    bw_bps = dr.GetBitRate()
    bdp_bits = bw_bps * rtt_sec
    bdp_bytes = bdp_bits / 8.0
    bdp_pkts = bdp_bytes / pkt_size
    return bdp_pkts


def build_dumbbell_topology(data_rate_str: str, rtt_sec: float):
def build_dumbbell_topology2(data_rate_str: str, rtt_sec: float):
    """
    Construct a dumbbell topology:

             _                                 _
            |    n1------+        +------n4     |
            |            |        |             |
    Senders |    n2------R1------R2-----n5      |    Receivers
            |            |        |             |
            |_   n3------+        +------n6    _|
  

    Node indices:
        0: n0 (left sender 1)
        1: n1 (left sender 2)
        2: R0 (left router)
        3: R1 (right router)
        4: n2 (right receiver 1)
        5: n3 (right receiver 2)

    RTT setting:
        - Non-bottleneck links use 0ms delay
        - Bottleneck one-way delay = RTT/2
          so end-to-end RTT ≈ RTT

    Returns:
        nodes: NodeContainer with 6 nodes
        (ip_n2, ip_n3): IPv4 addresses of n2 and n3
    """
    networkTop = '''
    Network Topology:
             _                                 _
            |    n1------+        +------n4     |
            |            |        |             |
    Senders |    n2------R1------R2-----n5      |    Receivers
            |            |        |             |
            |_   n3------+        +------n6    _|
    '''
    #print(networkTop)


    # nodes = ns.NodeContainer()
    # nodes.Create(6)

    senders = ns.NodeContainer()
    senders.Create(3)

    routers = ns.NodeContainer()
    routers.Create(2)

    receivers = ns.NodeContainer()
    receivers.Create(3)

    stack = ns.InternetStackHelper()
    stack.Install(senders)
    stack.Install(receivers)
    stack.Install(routers)

    # Non-bottleneck links: fast, 0ms delay
    # print("setting up non-bottleneck links")
    edge = ns.CsmaHelper()
    edge.SetChannelAttribute("DataRate", ns.StringValue("1Mbps"))
    edge.SetChannelAttribute("Delay", ns.StringValue("2ms"))


    # Bottleneck: controlled RTT
    # print("setting up bottleneck links")
    bottleneck = ns.PointToPointHelper()
    bottleneck.SetDeviceAttribute("DataRate", ns.StringValue(data_rate_str))

    one_way_ms = int((rtt_sec * 1000.0) / 2.0)
    bottleneck.SetChannelAttribute("Delay", ns.StringValue(f"{one_way_ms}ms"))

    address = ns.Ipv4AddressHelper()
    # print("setting up c0")

    # n1,n2,n3 -- R0
    c0 = ns.NodeContainer()
    c0.Add(routers.Get(0)) # router R1
    c0.Add(senders) # n1,n2,n3
    dev0 = edge.Install(c0)
    address.SetBase(ns.Ipv4Address("10.1.1.0"), ns.Ipv4Mask("255.255.255.0"))
    address.Assign(dev0)

    # print("setting up c1")
    # R1 -- n4, n5, n6
    c1 = ns.NodeContainer()
    c1.Add(routers.Get(1)) # router R2
    c1.Add(receivers) # n4 n5 n6
    dev1 = edge.Install(c1)
    address.SetBase(ns.Ipv4Address("10.1.3.0"), ns.Ipv4Mask("255.255.255.0"))
    recAddresses = address.Assign(dev1)

    # print("setting up R0-R1")
    # R0 -- R1 (bottleneck)
    devb = bottleneck.Install(routers.Get(0), routers.Get(1)) # R0-R1
    address.SetBase(ns.Ipv4Address("10.1.2.0"), ns.Ipv4Mask("255.255.255.0"))
    address.Assign(devb)

    ns.Ipv4GlobalRoutingHelper.PopulateRoutingTables()

    ip_n2 = recAddresses.GetAddress(2)
    ip_n3 = recAddresses.GetAddress(3)

    return (senders, routers, receivers), (ip_n2, ip_n3), edge


def run_single_flow(data_rate_str, rtt_sec, cca_name, cca_type_str):
    """
    (A) Single-flow experiment:
        - dumbbell topology
        - one long-lived flow: n0 -> n2
        - all nodes use same CCA (Cubic/Vegas)
        - variable: RTT
    """
    set_tcp_cca(cca_type_str)
    bdp_pkts = calc_bdp_pkts(data_rate_str, rtt_sec, PKT_SIZE)
    queue_pkts = int(round(bdp_pkts * BUF_FACTOR))

    ns.Config.SetDefault(
        "ns3::DropTailQueue<Packet>::MaxSize",
        ns.StringValue(f"{queue_pkts}p")
    )

    (senders, routers, receivers), (ip_n2, _), edge = build_dumbbell_topology2(data_rate_str, rtt_sec)

    port = 8080

    # Receiver on n2
    sock_addr = ns.InetSocketAddress(ip_n2, port)
    sink_addr = sock_addr.ConvertTo()
    sink_helper = ns.PacketSinkHelper("ns3::TcpSocketFactory", sink_addr)
    sink_apps = sink_helper.Install(receivers.Get(1))
    sink_apps.Start(ns.Seconds(0.0))
    sink_apps.Stop(ns.Seconds(SIM_TIME))

    # Sender on n0
    bulk = ns.BulkSendHelper("ns3::TcpSocketFactory", sink_addr)
    bulk.SetAttribute("MaxBytes", ns.UintegerValue(0))
    bulk.SetAttribute("SendSize", ns.UintegerValue(PKT_SIZE))

    client = bulk.Install(senders.Get(0))
    client.Start(ns.Seconds(APP_START))
    client.Stop(ns.Seconds(SIM_TIME))


    trace_file = ascii.CreateFileStream("cc_rtt_single_" + cca_name + f"_rtt{int(rtt_sec*1000)}ms.tr")
    edge.EnableAsciiAll(trace_file)



    ns.Simulator.Stop(ns.Seconds(SIM_TIME))
    ns.Simulator.Run()


    total_rx = sink_apps.Get(0).GetTotalRx()

    duration = SIM_TIME - APP_START
    throughput_mbps = (total_rx * 8.0) / duration / 1e6

    rtt_ms = int(rtt_sec * 1000)

    print(f"[Single ] CCA={cca_name:6s}, RTT={rtt_ms:2d}ms, "
          f"BW={data_rate_str:7s}, Buffer={BUF_FACTOR:.1f}×BDP "
          f"(≈ {queue_pkts:3d} pkts), TotalRx={total_rx} bytes, "
          f"Throughput≈{throughput_mbps:.3f} Mbps")

    ns.Simulator.Destroy()


def run_competing_flows(data_rate_str, rtt_sec, start2=5.0):
    """
    (B) Two-flows competition experiment:

        - Flow 1: CUBIC (n0 -> n2), starts at t=1s
        - Flow 2: Vegas (n1 -> n3), starts at t=5s
        - Same RTT for both flows
    """
    bdp_pkts = calc_bdp_pkts(data_rate_str, rtt_sec, PKT_SIZE)
    queue_pkts = int(round(bdp_pkts * BUF_FACTOR))

    ns.Config.SetDefault(
        "ns3::DropTailQueue<Packet>::MaxSize",
        ns.StringValue(f"{queue_pkts}p")
    )

    (senders, routers, receivers), (ip_n2, ip_n3) = build_dumbbell_topology2(data_rate_str, rtt_sec)

    set_node_tcp_cca(0, "ns3::TcpCubic")
    set_node_tcp_cca(1, "ns3::TcpVegas")

    port1 = 8080
    port2 = 8081

    # Flow 1 receiver
    addr1 = ns.InetSocketAddress(ip_n2, port1)
    sink_addr1 = addr1.ConvertTo()
    sink_h1 = ns.PacketSinkHelper("ns3::TcpSocketFactory", sink_addr1)
    sink_apps1 = sink_h1.Install(receivers.Get(0))
    sink_apps1.Start(ns.Seconds(0.0))
    sink_apps1.Stop(ns.Seconds(SIM_TIME))

    # Flow 1 sender
    bulk1 = ns.BulkSendHelper("ns3::TcpSocketFactory", sink_addr1)
    bulk1.SetAttribute("MaxBytes", ns.UintegerValue(0))
    bulk1.SetAttribute("SendSize", ns.UintegerValue(PKT_SIZE))
    client1 = bulk1.Install(senders.Get(0))
    client1.Start(ns.Seconds(APP_START))
    client1.Stop(ns.Seconds(SIM_TIME))

    # Flow 2 receiver
    addr2 = ns.InetSocketAddress(ip_n3, port2)
    sink_addr2 = addr2.ConvertTo()
    sink_h2 = ns.PacketSinkHelper("ns3::TcpSocketFactory", sink_addr2)
    sink_apps2 = sink_h2.Install(receivers.Get(2))
    sink_apps2.Start(ns.Seconds(0.0))
    sink_apps2.Stop(ns.Seconds(SIM_TIME))

    # Flow 2 sender
    bulk2 = ns.BulkSendHelper("ns3::TcpSocketFactory", sink_addr2)
    bulk2.SetAttribute("MaxBytes", ns.UintegerValue(0))
    bulk2.SetAttribute("SendSize", ns.UintegerValue(PKT_SIZE))
    client2 = bulk2.Install(senders.Get(1))
    client2.Start(ns.Seconds(start2))
    client2.Stop(ns.Seconds(SIM_TIME))

    ns.Simulator.Stop(ns.Seconds(SIM_TIME))
    ns.Simulator.Run()

    rx1 = sink_apps1.Get(0).GetTotalRx()
    rx2 = sink_apps2.Get(0).GetTotalRx()

    duration = SIM_TIME - APP_START
    thr1 = (rx1 * 8.0) / duration / 1e6
    thr2 = (rx2 * 8.0) / duration / 1e6

    rtt_ms = int(rtt_sec * 1000)

    print(f"[Compete] RTT={rtt_ms:2d}ms, BW={data_rate_str:7s}, "
          f"Buffer={BUF_FACTOR:.1f}×BDP (≈{queue_pkts:3d} pkts), "
          f"Cubic: Rx={rx1}, Thr≈{thr1:.3f} Mbps; "
          f"Vegas: Rx={rx2}, Thr≈{thr2:.3f} Mbps")

    ns.Simulator.Destroy()


def main():
    print("=== (A) Single-flow experiment: CUBIC / Vegas, RTT in {10,20,40} ms ===")
    for rtt_sec in RTT_LIST:
        for cca_name, cca_type in CCA_TYPES_SINGLE:
            run_single_flow(DATA_RATE, rtt_sec, cca_name, cca_type)

    print("\n=== (B) Two-flows competition: Cubic vs Vegas, RTT in {10,20,40} ms ===")
    for rtt_sec in RTT_LIST:
        run_competing_flows(DATA_RATE, rtt_sec, start2=5.0)


if __name__ == "__main__":
    main()
