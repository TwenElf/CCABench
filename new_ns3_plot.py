import pandas as pd
import matplotlib.pyplot as plt
import os

# ====== 结果文件所在路径 ======
base = r"A:\ns3_result"

# ====== 输出图片目录 ======
output_dir = r"A:\ns3-diagram"
os.makedirs(output_dir, exist_ok=True)   # 自动创建目录

# RTT 列表
RTTs = ["2ms", "20ms", "40ms"]

# ================================
# 文件名模板（对齐你现在目录里的名字）
# ================================

# 单流：Vegas / Cubic / BBR
single_vegas = "throughput_TcpVegas_single_rtt_{}"
single_cubic = "throughput_TcpCubic_single_rtt_{}"
single_bbr   = "throughput_TcpBbr_single_rtt_{}"

# Vegas vs BBR 竞争（flow0=Vegas, flow1=BBR）
comp_vb_flow0 = "throughput_TcpVegas_vs_TcpBbr_flow0_rtt_{}"
comp_vb_flow1 = "throughput_TcpVegas_vs_TcpBbr_flow1_rtt_{}"

# Cubic vs BBR 竞争（flow0=Cubic, flow1=BBR）
comp_cb_flow0 = "throughput_TcpCubic_vs_TcpBbr_flow0_rtt_{}"
comp_cb_flow1 = "throughput_TcpCubic_vs_TcpBbr_flow1_rtt_{}"

# Vegas vs Cubic 竞争（flow0=Vegas, flow1=Cubic）
comp_vc_flow0 = "throughput_TcpVegas_vs_TcpCubic_flow0_rtt_{}"
comp_vc_flow1 = "throughput_TcpVegas_vs_TcpCubic_flow1_rtt_{}"


# matplotlib 设置（保留你原本风格）
plt.style.use("seaborn-v0_8")
plt.rcParams["figure.figsize"] = (10, 5)
plt.rcParams["font.size"] = 12


# ======================================================
# 读取 CSV（忽略 # 注释行）
# ======================================================
def load_throughput_csv(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"文件不存在: {path}")
    df = pd.read_csv(
        path,
        comment="#",
        header=None,
        names=["time_sec", "throughput_Mbps"]
    )
    # 防御：去掉坏行/空行
    df["time_sec"] = pd.to_numeric(df["time_sec"], errors="coerce")
    df["throughput_Mbps"] = pd.to_numeric(df["throughput_Mbps"], errors="coerce")
    df = df.dropna().reset_index(drop=True)
    return df


# ======================================================
# 画单流图：Cubic + Vegas + BBR
# ======================================================
def plot_single(rtt: str):
    file_cubic = os.path.join(base, single_cubic.format(rtt) + ".csv")
    file_vegas = os.path.join(base, single_vegas.format(rtt) + ".csv")
    file_bbr   = os.path.join(base, single_bbr.format(rtt) + ".csv")

    df_cubic = load_throughput_csv(file_cubic)
    df_vegas = load_throughput_csv(file_vegas)
    df_bbr   = load_throughput_csv(file_bbr)

    plt.figure()
    plt.plot(df_cubic["time_sec"], df_cubic["throughput_Mbps"],
             label="TcpCubic (single)", linewidth=2)
    plt.plot(df_vegas["time_sec"], df_vegas["throughput_Mbps"],
             label="TcpVegas (single)", linewidth=2)
    plt.plot(df_bbr["time_sec"], df_bbr["throughput_Mbps"],
             label="TcpBbr (single)", linewidth=2)

    plt.xlabel("Time (s)")
    plt.ylabel("Throughput (Mbps)")
    plt.title(f"Single-flow Throughput (RTT = {rtt})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"single_rtt_{rtt}.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved: {save_path}")


# ======================================================
# 画竞争图：Vegas vs BBR
# ======================================================
def plot_compete_vegas_bbr(rtt: str):
    file_flow0 = os.path.join(base, comp_vb_flow0.format(rtt) + ".csv")
    file_flow1 = os.path.join(base, comp_vb_flow1.format(rtt) + ".csv")

    df0 = load_throughput_csv(file_flow0)
    df1 = load_throughput_csv(file_flow1)

    plt.figure()
    plt.plot(df0["time_sec"], df0["throughput_Mbps"],
             label="flow0 (TcpVegas)", linewidth=2)
    plt.plot(df1["time_sec"], df1["throughput_Mbps"],
             label="flow1 (TcpBbr)", linewidth=2)

    plt.xlabel("Time (s)")
    plt.ylabel("Throughput (Mbps)")
    plt.title(f"Two-flow Competition: TcpVegas vs TcpBbr (RTT = {rtt})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"competition_vegas_bbr_rtt_{rtt}.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved: {save_path}")


# ======================================================
# 画竞争图：Cubic vs BBR
# ======================================================
def plot_compete_cubic_bbr(rtt: str):
    file_flow0 = os.path.join(base, comp_cb_flow0.format(rtt) + ".csv")
    file_flow1 = os.path.join(base, comp_cb_flow1.format(rtt) + ".csv")

    df0 = load_throughput_csv(file_flow0)
    df1 = load_throughput_csv(file_flow1)

    plt.figure()
    plt.plot(df0["time_sec"], df0["throughput_Mbps"],
             label="flow0 (TcpCubic)", linewidth=2)
    plt.plot(df1["time_sec"], df1["throughput_Mbps"],
             label="flow1 (TcpBbr)", linewidth=2)

    plt.xlabel("Time (s)")
    plt.ylabel("Throughput (Mbps)")
    plt.title(f"Two-flow Competition: TcpCubic vs TcpBbr (RTT = {rtt})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"competition_cubic_bbr_rtt_{rtt}.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved: {save_path}")


# ======================================================
# 画竞争图：Vegas vs Cubic
# ======================================================
def plot_compete_vegas_cubic(rtt: str):
    file_flow0 = os.path.join(base, comp_vc_flow0.format(rtt) + ".csv")
    file_flow1 = os.path.join(base, comp_vc_flow1.format(rtt) + ".csv")

    df0 = load_throughput_csv(file_flow0)
    df1 = load_throughput_csv(file_flow1)

    plt.figure()
    plt.plot(df0["time_sec"], df0["throughput_Mbps"],
             label="flow0 (TcpVegas)", linewidth=2)
    plt.plot(df1["time_sec"], df1["throughput_Mbps"],
             label="flow1 (TcpCubic)", linewidth=2)

    plt.xlabel("Time (s)")
    plt.ylabel("Throughput (Mbps)")
    plt.title(f"Two-flow Competition: TcpVegas vs TcpCubic (RTT = {rtt})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"competition_vegas_cubic_rtt_{rtt}.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved: {save_path}")


# ======================================================
# 主程序
# ======================================================
if __name__ == "__main__":
    for rtt in RTTs:
        print(f"\n==== 绘制单流 RTT = {rtt} ====")
        plot_single(rtt)

        print(f"\n==== 绘制 Vegas vs BBR 竞争 RTT = {rtt} ====")
        plot_compete_vegas_bbr(rtt)

        print(f"\n==== 绘制 Cubic vs BBR 竞争 RTT = {rtt} ====")
        plot_compete_cubic_bbr(rtt)

        print(f"\n==== 绘制 Vegas vs Cubic 竞争 RTT = {rtt} ====")
        plot_compete_vegas_cubic(rtt)

    print(f"\n✔ 所有图已生成，已保存到：{output_dir}")
