import pandas as pd
import matplotlib.pyplot as plt
import os 
import subprocess
import re

files = [f for f in os.listdir('.') if re.match(r'.*_throughput\.csv', f)]


for f in files:
    if f.startswith("compete_Cubic"):
        continue
    elif f.startswith("compete_Vegas"):
        df = pd.read_csv(f)
        df2 = pd.read_csv(f.replace("compete_Vegas", "compete_Cubic"))
        name = f[:-len('_throughput.csv')].replace("compete_Vegas", "compare_Cubic_vs_Vegas")
        plt.figure()
        plt.plot(df["time"], df["throughput_mbps"],label='Vegas')
        plt.plot(df2["time"], df2["throughput_mbps"],label='Cubic')
        plt.legend()
        plt.xlabel("Time (s)")
        plt.ylabel("Throughput (Mbps)")
        plt.title(name)
        plt.grid()
        plt.savefig(name+".png", bbox_inches='tight')
        plt.close()
    else:
        df = pd.read_csv(f)
        name = f[:-len('_throughput.csv')].replace("compete_Vegas", "compare_Cubic_vs_Vegas")
        plt.figure()
        plt.plot(df["time"], df["throughput_mbps"])
        plt.xlabel("Time (s)")
        plt.ylabel("Throughput (Mbps)")
        plt.title(name)
        plt.grid()
        plt.savefig(name+".png", bbox_inches='tight')
        plt.close()