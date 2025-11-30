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

from ns import ns

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
    """
    Construct a dumbbell topology:

        n0 ----\
                \
                 R0 ====== R1
                /          \
        n1 ----/            \
                              n2
                              n3

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
    nodes = ns.NodeContainer()
    nodes.Create(6)

    stack = ns.InternetStackHelper()
    stack.Install(nodes)

    # Non-bottleneck links: fast, 0ms delay
    edge = ns.PointToPointHelper()
    edge.SetDeviceAttribute("DataRate", ns.StringValue("1Gbps"))
    edge.SetChannelAttribute("Delay", ns.StringValue("0ms"))

    # Bottleneck: controlled RTT
    bottleneck = ns.PointToPointHelper()
    bottleneck.SetDeviceAttribute("DataRate", ns.StringValue(data_rate_str))

    one_way_ms = int((rtt_sec * 1000.0) / 2.0)
    bottleneck.SetChannelAttribute("Delay", ns.StringValue(f"{one_way_ms}ms"))

    address = ns.Ipv4AddressHelper()

    # n0 -- R0
    c0 = ns.NodeContainer(nodes.Get(0), nodes.Get(2))
    dev0 = edge.Install(c0)
    address.SetBase(ns.Ipv4Address("10.1.1.0"), ns.Ipv4Mask("255.255.255.0"))
    _ = address.Assign(dev0)

    # n1 -- R0
    c1 = ns.NodeContainer(nodes.Get(1), nodes.Get(2))
    dev1 = edge.Install(c1)
    address.SetBase(ns.Ipv4Address("10.1.2.0"), ns.Ipv4Mask("255.255.255.0"))
    _ = address.Assign(dev1)

    # R0 -- R1 (bottleneck)
    cb = ns.NodeContainer(nodes.Get(2), nodes.Get(3))
    devb = bottleneck.Install(cb)
    address.SetBase(ns.Ipv4Address("10.1.3.0"), ns.Ipv4Mask("255.255.255.0"))
    _ = address.Assign(devb)

    # R1 -- n2
    c2 = ns.NodeContainer(nodes.Get(3), nodes.Get(4))
    dev2 = edge.Install(c2)
    address.SetBase(ns.Ipv4Address("10.1.4.0"), ns.Ipv4Mask("255.255.255.0"))
    iface2 = address.Assign(dev2)

    # R1 -- n3
    c3 = ns.NodeContainer(nodes.Get(3), nodes.Get(5))
    dev3 = edge.Install(c3)
    address.SetBase(ns.Ipv4Address("10.1.5.0"), ns.Ipv4Mask("255.255.255.0"))
    iface3 = address.Assign(dev3)

    ns.Ipv4GlobalRoutingHelper.PopulateRoutingTables()

    ip_n2 = iface2.GetAddress(1)
    ip_n3 = iface3.GetAddress(1)

    return nodes, (ip_n2, ip_n3)


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

    nodes, (ip_n2, _) = build_dumbbell_topology(data_rate_str, rtt_sec)

    port = 8080

    # Receiver on n2
    sock_addr = ns.InetSocketAddress(ip_n2, port)
    sink_addr = sock_addr.ConvertTo()
    sink_helper = ns.PacketSinkHelper("ns3::TcpSocketFactory", sink_addr)
    sink_apps = sink_helper.Install(nodes.Get(4))
    sink_apps.Start(ns.Seconds(0.0))
    sink_apps.Stop(ns.Seconds(SIM_TIME))

    # Sender on n0
    bulk = ns.BulkSendHelper("ns3::TcpSocketFactory", sink_addr)
    bulk.SetAttribute("MaxBytes", ns.UintegerValue(0))
    bulk.SetAttribute("SendSize", ns.UintegerValue(PKT_SIZE))

    client = bulk.Install(nodes.Get(0))
    client.Start(ns.Seconds(APP_START))
    client.Stop(ns.Seconds(SIM_TIME))

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

    nodes, (ip_n2, ip_n3) = build_dumbbell_topology(data_rate_str, rtt_sec)

    set_node_tcp_cca(0, "ns3::TcpCubic")
    set_node_tcp_cca(1, "ns3::TcpVegas")

    port1 = 8080
    port2 = 8081

    # Flow 1 receiver
    addr1 = ns.InetSocketAddress(ip_n2, port1)
    sink_addr1 = addr1.ConvertTo()
    sink_h1 = ns.PacketSinkHelper("ns3::TcpSocketFactory", sink_addr1)
    sink_apps1 = sink_h1.Install(nodes.Get(4))
    sink_apps1.Start(ns.Seconds(0.0))
    sink_apps1.Stop(ns.Seconds(SIM_TIME))

    bulk1 = ns.BulkSendHelper("ns3::TcpSocketFactory", sink_addr1)
    bulk1.SetAttribute("MaxBytes", ns.UintegerValue(0))
    bulk1.SetAttribute("SendSize", ns.UintegerValue(PKT_SIZE))
    client1 = bulk1.Install(nodes.Get(0))
    client1.Start(ns.Seconds(APP_START))
    client1.Stop(ns.Seconds(SIM_TIME))

    # Flow 2 receiver
    addr2 = ns.InetSocketAddress(ip_n3, port2)
    sink_addr2 = addr2.ConvertTo()
    sink_h2 = ns.PacketSinkHelper("ns3::TcpSocketFactory", sink_addr2)
    sink_apps2 = sink_h2.Install(nodes.Get(5))
    sink_apps2.Start(ns.Seconds(0.0))
    sink_apps2.Stop(ns.Seconds(SIM_TIME))

    bulk2 = ns.BulkSendHelper("ns3::TcpSocketFactory", sink_addr2)
    bulk2.SetAttribute("MaxBytes", ns.UintegerValue(0))
    bulk2.SetAttribute("SendSize", ns.UintegerValue(PKT_SIZE))
    client2 = bulk2.Install(nodes.Get(1))
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
