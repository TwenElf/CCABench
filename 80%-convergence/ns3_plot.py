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

# 文件名模板
single_cubic = "throughput_cubic_single_rtt_{}"
single_vegas = "throughput_vegas_single_rtt_{}"
comp_cubic   = "throughput_cubic_compete_rtt_{}"
comp_vegas   = "throughput_vegas_compete_rtt_{}"

# matplotlib 设置
plt.style.use("seaborn-v0_8")
plt.rcParams["figure.figsize"] = (10, 5)
plt.rcParams["font.size"] = 12


# ======================================================
# 读取 CSV（忽略 # 注释行）
# ======================================================
def load_throughput_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(
        path,
        comment="#",
        header=None,
        names=["time_sec", "throughput_Mbps"]
    )


# ======================================================
# 画单流图
# ======================================================
def plot_single(rtt: str):
    file_cubic = os.path.join(base, single_cubic.format(rtt) + ".csv")
    file_vegas = os.path.join(base, single_vegas.format(rtt) + ".csv")

    df_cubic = load_throughput_csv(file_cubic)
    df_vegas = load_throughput_csv(file_vegas)

    plt.figure()
    plt.plot(df_cubic["time_sec"], df_cubic["throughput_Mbps"], label="Cubic (single)", linewidth=2)
    plt.plot(df_vegas["time_sec"], df_vegas["throughput_Mbps"], label="Vegas (single)", linewidth=2)

    plt.xlabel("Time (s)")
    plt.ylabel("Throughput (Mbps)")
    plt.title(f"Single-flow Throughput (RTT = {rtt})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # ====== 保存到 A:\ns3-diagram ======
    save_path = os.path.join(output_dir, f"single_rtt_{rtt}.png")
    plt.savefig(save_path, dpi=300)
    plt.close()


# ======================================================
# 画竞争图
# ======================================================
def plot_compete(rtt: str):
    file_cubic = os.path.join(base, comp_cubic.format(rtt) + ".csv")
    file_vegas = os.path.join(base, comp_vegas.format(rtt) + ".csv")

    df_cubic = load_throughput_csv(file_cubic)
    df_vegas = load_throughput_csv(file_vegas)

    plt.figure()
    plt.plot(df_cubic["time_sec"], df_cubic["throughput_Mbps"], label="Cubic (compete)", linewidth=2)
    plt.plot(df_vegas["time_sec"], df_vegas["throughput_Mbps"], label="Vegas (compete)", linewidth=2)

    plt.xlabel("Time (s)")
    plt.ylabel("Throughput (Mbps)")
    plt.title(f"Two-flow Competition (RTT = {rtt})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # ====== 保存到 A:\ns3-diagram ======
    save_path = os.path.join(output_dir, f"competition_rtt_{rtt}.png")
    plt.savefig(save_path, dpi=300)
    plt.close()


# ======================================================
# 主程序
# ======================================================
if __name__ == "__main__":
    for rtt in RTTs:
        print(f"\n==== 绘制单流 RTT = {rtt} ====")
        plot_single(rtt)

        print(f"\n==== 绘制竞争 RTT = {rtt} ====")
        plot_compete(rtt)

    print(f"\n✔ 所有图已生成（共 6 张 PNG），已保存到：{output_dir}")
