import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.scale import FuncScale
from matplotlib.ticker import FixedLocator
from adjustText import adjust_text
import matplotlib as mpl

from helpers.Constants import *


new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none',
"font.family": "sans-serif",
"font.sans-serif": "Arial"
}
mpl.rcParams.update(new_rc_params)
# ------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------
df = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/DongAhnSD_Var/DongAnSD_PaperStats.tsv', sep="\t")

euchro = df.loc[(df['region'] == 'euchromatin') & (~df['cyto'].isin(['4q28.3', '4q32.2']))]
euchro.sort_values(by='max_min_diff', ascending=False, inplace=True)
euchro.set_index('region_id', inplace=True)

# ------------------------------------------------------------------
# 2. Style mappings
# ------------------------------------------------------------------
color_map = {"euchromatin": "tab:blue", "SubTel": "tab:orange", "CenSat": "tab:green"}
marker_map = {"HasINV": "s", "NoINV": "o"}

# Draw HasINV (squares) on top of NoINV (circles) since it's the rarer group
inv_order = ["NoINV", "HasINV"]

# ------------------------------------------------------------------
# 3. Build a short label from cyto + genes
# ------------------------------------------------------------------
MAX_GENES_SHOWN = 3


def make_label(row):
    cyto = row["cyto"] if isinstance(row["cyto"], str) and row["cyto"] else "?"
    genes_raw = row["genes"] if isinstance(row["genes"], str) else ""
    gene_list = [g for g in genes_raw.split(",") if g and g != "."]
    if not gene_list:
        gene_str = "no genes"
    elif len(gene_list) <= MAX_GENES_SHOWN:
        gene_str = ", ".join(gene_list)
    else:
        gene_str = ", ".join(gene_list[:MAX_GENES_SHOWN]) + f" +{len(gene_list) - MAX_GENES_SHOWN} more"
    return f"{cyto}\n{gene_str}"


# ------------------------------------------------------------------
# 4. Custom "broken but continuous" y-scale
# ------------------------------------------------------------------
y_break = 1.0  # values below this are linear/full-size
compress = 0.12  # values above this are compressed to `compress` x their span


def forward(y):
    y = np.asarray(y, dtype=float)
    return np.where(y <= y_break, y, y_break + (y - y_break) * compress)


def inverse(y):
    y = np.asarray(y, dtype=float)
    return np.where(y <= y_break, y, y_break + (y - y_break) / compress)


# ------------------------------------------------------------------
# 5. Plot (single axis, continuous data, no clipping)
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(23, 9))

for inv_tag in inv_order:
    marker = marker_map[inv_tag]
    for region, color in color_map.items():
        sub = df[(df["region"] == region) & (df["inv_tag"] == inv_tag)]
        ax.scatter(
            sub["median"], sub["max_min_diff"],
            c=color, marker=marker, s=140, rasterized=True,
            alpha=0.55, edgecolors="none",
        )

ax.set_xscale("log")
ax.set_yscale(FuncScale(ax, (forward, inverse)))

y_max = df["max_min_diff"].max() * 1.08
ax.set_ylim(-0.02, y_max)

# Explicit ticks: dense spacing below the break, sparser above it
low_ticks = np.arange(0, y_break + 1e-9, 0.2)
high_ticks = np.arange(2, np.ceil(y_max) + 1, 1)
ax.yaxis.set_major_locator(FixedLocator(np.concatenate([low_ticks, high_ticks])))

ax.set_xlabel("median (log scale)")
ax.set_ylabel(f"max_min_diff (Mbp) -- compressed above {y_break} Mbp")
ax.set_title("max_min_diff vs median, by region and inversion status")

# Mark where the compression starts with a subtle reference line
ax.axhline(y_break, color="gray", lw=0.8, ls="--", alpha=0.6)

# ------------------------------------------------------------------
# 6. Label points with max_min_diff > y_break
# ------------------------------------------------------------------
outliers = df[df["max_min_diff"] > y_break]

texts = []
for _, row in outliers.iterrows():
    if row['region_id'] in euchro.index and row['max_min_diff']>=1:
        texts.append(
            ax.text(
                row["median"], row["max_min_diff"], make_label(row),
                fontsize=16, color="black",
            )
        )

adjust_text(
    texts, ax=ax,
    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
    expand=(1.3, 1.8),
    force_text=(0.5, 0.8),
    force_points=(0.3, 0.5),
    max_move=200,
)

# ------------------------------------------------------------------
# 7. Combined legend (region = color, inv_tag = marker shape)
# ------------------------------------------------------------------
legend_elements = (
        [Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=9, label=r)
         for r, c in color_map.items()]
        + [Line2D([0], [0], marker=m, color="k", markerfacecolor="k", markersize=9, label=t)
           for t, m in marker_map.items()]
)
ax.legend(handles=legend_elements, loc="upper left", frameon=True, fontsize=18)
ax.spines[['right', 'top']].set_visible(False)
plt.tight_layout()
fig.savefig(f'{VOL28}/SVREF/hq_asm/CHM13/DongAhnSD_Var/sites_MaxMedian3.svg')
plt.show()