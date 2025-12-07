import pandas as pd
import matplotlib.pyplot as plt
import os

# ====== Path to result CSV files ======
base = r"A:\ns3_result"

# ====== Output directory for figures ======
output_dir = r"A:\ns3-diagram"
os.makedirs(output_dir, exist_ok=True)   # Automatically create directory if not exists

# RTT list
RTTs = ["2ms", "20ms", "40ms"]

# Filename templates
single_cubic = "throughput_cubic_single_rtt_{}"
single_vegas = "throughput_vegas_single_rtt_{}"
comp_cubic   = "throughput_cubic_compete_rtt_{}"
comp_vegas   = "throughput_vegas_compete_rtt_{}"

# matplotlib settings
plt.style.use("seaborn-v0_8")
plt.rcParams["figure.figsize"] = (10, 5)
plt.rcParams["font.size"] = 12


# ======================================================
# Load throughput CSV (ignoring comment lines starting with #)
# ======================================================
def load_throughput_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(
        path,
        comment="#",
        header=None,
        names=["time_sec", "throughput_Mbps"]
    )


# ======================================================
# Plot single-flow throughput
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

    # ====== Save figure to A:\ns3-diagram ======
    save_path = os.path.join(output_dir, f"single_rtt_{rtt}.png")
    plt.savefig(save_path, dpi=300)
    plt.close()


# ======================================================
# Plot two-flow competition throughput
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

    # ====== Save figure to A:\ns3-diagram ======
    save_path = os.path.join(output_dir, f"competition_rtt_{rtt}.png")
    plt.savefig(save_path, dpi=300)
    plt.close()


# ======================================================
# Main script
# ======================================================
if __name__ == "__main__":
    for rtt in RTTs:
        print(f"\n==== Plotting single-flow RTT = {rtt} ====")
        plot_single(rtt)

        print(f"\n==== Plotting competition RTT = {rtt} ====")
        plot_compete(rtt)

    print(f"\n✔ All figures generated (6 PNG files) and saved to: {output_dir}")
