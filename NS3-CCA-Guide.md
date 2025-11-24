# NS-3 CCA Benchmarking Guide  
A step-by-step guide for running Congestion Control Algorithm (CCA) experiments in ns-3.

---

## 1. Download and Build ns-3

### Recommended version: ns-3.40  
Download: https://www.nsnam.org/releases/ns-3.40/

Extract the package:
tar xf ns-allinone-3.40.tar.bz2

Enter the ns-3 directory:
cd ns-allinone-3.40/ns-3.40

Configure and build:
./ns3 configure --enable-examples --enable-tests --enable-modules=flow-monitor
./ns3 build

Using:  nano examples/tcp/CMakeLists.txt
Check whether the CMakeLists.txt include :
build_example(
  NAME tcp-variants-comparison
  SOURCE_FILES tcp-variants-comparison.cc
  LIBRARIES_TO_LINK
    ${libpoint-to-point}
    ${libinternet}
    ${libapplications}
    ${libtraffic-control}
    ${libnetwork}
    ${libinternet-apps}
    ${libflow-monitor}
)

iF YOU DO CHANGE REBUILD IT
Rebuild:
./ns3 configure --enable-examples --enable-tests --enable-modules=flow-monitor
./ns3 build


After building, the executable will be at:
build/examples/tcp/ns3.40-tcp-variants-comparison-default


Check it works:
./build/examples/tcp/ns3.40-tcp-variants-comparison-default --help


2. Running an Experiment
Example experiment using TCP NewReno, 10 Mbps bottleneck, 5 ms delay, 2 flows:
./build/examples/tcp/ns3.40-tcp-variants-comparison-default \
  --transport_prot=TcpNewReno \
  --bandwidth=10Mbps \
  --delay=5ms \
  --access_bandwidth=100Mbps \
  --access_delay=0ms \
  --num_flows=2 \
  --duration=60 \
  --flow_monitor=true


3. Retrieving the FlowMonitor Output
The simulation writes:
TcpVariantsComparison.flowmonitor


Find it:
find . -name "*.flowmonitor"
