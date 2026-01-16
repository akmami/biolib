import re
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from collections import defaultdict
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, FuncFormatter

out_dir = Path("../out/fa-hmin/plots")
out_dir.mkdir(exist_ok=True)

# -----------------------------
# Configuration
# -----------------------------
ks = [6, 9, 12, 15]
ws = [4, 6, 8, 10]
data_dir = Path("../out/fa-hmin")  # directory containing the files

k_markers = {
    6: "o",    # triangle
    9: "^",    # square
    12: "s",   # pentagon
    15: "p",   # hexagon
}

w_colors = {
    4: "tab:blue",
    6: "tab:orange",
    8: "tab:green",
    10: "tab:red",
}

length_by_w = defaultdict(list)

# -----------------------------
# Helpers
# -----------------------------
def log_comma_formatter(x, pos):
    if x >= 1:
        return f"{int(x):,}"
    return ""

def clean_line(line):
    line = line.replace("\\\\", "")
    line = re.sub(r"\\(midrule|bottomrule|toprule)", "", line)
    return line.strip()

def parse_row(line):
    line = clean_line(line)

    values = []
    for token in line.split("&")[1:]:
        token = token.strip().replace(",", "")
        try:
            values.append(float(token))
        except ValueError:
            pass
    return values

def parse_row_int(line):
    line = clean_line(line)

    values = []
    for token in line.split("&")[1:]:
        token = token.strip().replace(",", "")
        try:
            values.append(int(token))
        except ValueError:
            pass
    return values

def parse_row_str(line):
    line = clean_line(line)

    values = []
    for token in line.split("&")[1:]:
        token = token.strip().replace(",", "")
        try:
            values.append(token)
        except ValueError:
            pass
    return values

# -----------------------------
# Containers for plots
# -----------------------------
core_data = []
length_data = []
distance_data = []

def parse_file(filepath):
    stats = {}

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if line.startswith("Hierarchy level"):
                stats["level"] = parse_row_str(line)
            elif line.startswith("Total \\# Cores"):
                stats["total_cores"] = parse_row_int(line)
            elif line.startswith("Unique Cores"):
                stats["unique_cores"] = parse_row_int(line)
            elif line.startswith("Avg Distance"):
                stats["avg_dist"] = parse_row(line)
            elif line.startswith("StdDev Distance"):
                stats["std_dist"] = parse_row(line)
            elif line.startswith("Avg Length"):
                stats["avg_len"] = parse_row(line)
            elif line.startswith("StdDev Length"):
                stats["std_len"] = parse_row(line)
    return stats

for k in ks:
    for w in ws:
        file = data_dir / f"minimizers-{k}-{w}-output.txt"
        if not file.exists():
            print(f"Skipping missing file: {file}")
            continue

        stats = parse_file(file)
        level = np.array(stats["level"])

        label = f"k={k}, w={w}"

        core_data.append((level, stats["total_cores"], stats["unique_cores"], k, w, label))
        length_data.append((level, stats["avg_len"], stats["std_len"], k, w, label))
        distance_data.append((level, stats["avg_dist"], stats["std_dist"], k, w, label))

# -----------------------------
# Plot 1: Core counts
# -----------------------------
plt.figure(figsize=(8, 6))
for level, total, unique, k, w, label in core_data:
    plt.plot(level, total, color=w_colors[w], marker=k_markers[k], linestyle="--", label=f"Total ({label})")
    # plt.plot(level, unique, color=w_colors[w], marker=k_markers[k], linestyle="--", label=f"Unique ({label})")

plt.xlabel("Hierarchy level")
plt.ylabel("Core count")
plt.title("Core counts vs Hierarchy level")
plt.legend(fontsize=8, ncol=2)
plt.grid(True)

plt.gca().set_yscale("log")
plt.gca().yaxis.set_major_locator(LogLocator(base=10))
plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x):,}"))

k_legend = [
    Line2D([0], [0], marker=marker, color="black", linestyle="None", markersize=8, label=f": {k}")
    for k, marker in k_markers.items()
]

w_legend = [
    Line2D([0], [0], color=color, lw=2, label=f": {w}")
    for w, color in w_colors.items()
]

legend1 = plt.legend(
    handles=k_legend,
    title="k-mer size",
    loc="upper right",
    bbox_to_anchor=(0.98, 1),
)

plt.gca().add_artist(legend1)

plt.legend(
    handles=w_legend,
    title="window size",
    loc="upper right",
    bbox_to_anchor=(0.85, 1),
)

plt.tight_layout()
plt.savefig(out_dir / "core_counts_vs_level.png", dpi=300)
plt.close()

# -----------------------------
# Plot 2: Core length
# -----------------------------
plt.figure(figsize=(8, 6))
for level, avg, std, k, w, label in length_data:
    plt.errorbar(level, avg, yerr=std, color=w_colors[w], marker=k_markers[k], linestyle="-", capsize=3, label=label)

plt.xlabel("Hierarchy level")
plt.ylabel("Core length")
plt.title("Core length vs Hierarchy level")
plt.legend(fontsize=8, ncol=2)
plt.grid(True)

plt.gca().set_yscale("log")
plt.gca().yaxis.set_major_locator(LogLocator(base=10))
plt.gca().yaxis.set_major_formatter(
    FuncFormatter(lambda x, pos: f"{int(x):,}")
)

k_legend = [
    Line2D([0], [0], marker=marker, color="black", linestyle="None",
           markersize=8, label=f": {k}")
    for k, marker in k_markers.items()
]

w_legend = [
    Line2D([0], [0], color=color, lw=2, label=f": {w}")
    for w, color in w_colors.items()
]

legend1 = plt.legend(
    handles=k_legend,
    title="k-mer size",
    loc="upper left",
    bbox_to_anchor=(0.02, 1),
)

plt.gca().add_artist(legend1)

plt.legend(
    handles=w_legend,
    title="window size",
    loc="upper left",
    bbox_to_anchor=(0.15, 1),
)

plt.tight_layout()
plt.savefig(out_dir / "core_length_vs_level.png", dpi=300)
plt.close()

# -----------------------------
# Plot 3: Distance
# -----------------------------
plt.figure(figsize=(8, 6))
for level, avg, std, k, w, label in distance_data:
    plt.errorbar(level, avg, yerr=std, color=w_colors[w], marker=k_markers[k], linestyle="-", capsize=3, label=label)

plt.xlabel("Hierarchy level")
plt.ylabel("Distance")
plt.title("Distances vs Hierarchy level")
plt.legend(fontsize=8, ncol=2)
plt.grid(True)

plt.gca().set_yscale("log")
plt.gca().yaxis.set_major_locator(LogLocator(base=10))
plt.gca().yaxis.set_major_formatter(
    FuncFormatter(lambda x, pos: f"{int(x):,}")
)

k_legend = [
    Line2D([0], [0], marker=marker, color="black", linestyle="None",
           markersize=8, label=f": {k}")
    for k, marker in k_markers.items()
]

w_legend = [
    Line2D([0], [0], color=color, lw=2, label=f": {w}")
    for w, color in w_colors.items()
]

legend1 = plt.legend(
    handles=k_legend,
    title="k-mer size",
    loc="upper left",
    bbox_to_anchor=(0.02, 1),
)

plt.gca().add_artist(legend1)

plt.legend(
    handles=w_legend,
    title="window size",
    loc="upper left",
    bbox_to_anchor=(0.17, 1),
)

plt.tight_layout()
plt.savefig(out_dir / "distance_vs_level.png", dpi=300)
plt.close()

print("Plots saved in ../out/fa-hmin/plots/")
