#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO   # <-- FIX: use built-in StringIO

TIME_COLUMN = "time_s"
THROUGHPUT_COLUMN = "throughput_mbps"

# "Comment block" prefix in your CSVs
COMMENT_PREFIX = "# Convergence summary"

def plot_file(filename):
    print(f"=== Processing {filename} ===")

    # Read the raw file
    with open(filename, "r") as f:
        lines = f.readlines()

    # Find the start of the convergence summary section
    summary_index = None
    for i, line in enumerate(lines):
        if line.startswith(COMMENT_PREFIX):
            summary_index = i
            break

    if summary_index is None:
        raise RuntimeError("No convergence summary found in file: " + filename)

    # --- DATA PART (top of file) ---
    data_lines = lines[:summary_index]

    # Convert CSV text → DataFrame using built-in StringIO
    df = pd.read_csv(
        StringIO("".join(data_lines)),
        names=[TIME_COLUMN, THROUGHPUT_COLUMN]
    )

    # --- SUMMARY PART (bottom of file) ---
    summary_lines = lines[summary_index + 1:]
    summary_df = pd.read_csv(
        StringIO("".join(summary_lines))
    )

    # Now plot
    plt.figure(figsize=(10, 5))
    plt.plot(df[TIME_COLUMN], df[THROUGHPUT_COLUMN], label="Throughput")

    # Draw threshold lines (bound_mbps)
    for _, row in summary_df.iterrows():
        bound = row["bound_mbps"]
        percent = row["percent_time_above"]
        plt.axhline(bound, linestyle="--", alpha=0.5, label=f"{bound} Mbps ({percent:.1f}%)")

    plt.title(os.path.basename(filename))
    plt.xlabel("Time (s)")
    plt.ylabel("Throughput (Mbps)")
    plt.legend()
    plt.grid(True)

    out_png = filename.replace(".csv", ".png")
    plt.savefig(out_png, dpi=150)
    plt.close()

    print(f"Saved {out_png}")


def main():
    cwd = os.getcwd()
    print("Scanning directory:", cwd)

    for f in sorted(os.listdir(".")):
        if f.endswith(".csv"):
            plot_file(f)


if __name__ == "__main__":
    main()
