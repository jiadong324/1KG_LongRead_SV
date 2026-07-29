#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 7/21/26
'''

"""
Visualize `max_min_diff` (Mbp) vs `ctg_brk_num` NORMALIZED BY REGION LENGTH,
from DongAnSD_PaperStats.tsv.

Normalization: region length = end - start (bp). `ctg_brk_num` is divided
by region length in Mbp, giving "contig breaks per Mbp" -- this puts
regions of very different sizes on a comparable scale, instead of large
regions trivially accumulating more breaks just because they're bigger.

- `max_min_diff` uses a single, continuous y-axis with a custom nonlinear
  scale: values below `y_break` are linear/full-size, values above it are
  compressed. This avoids squashing the dense low-value cluster while
  still showing the sparse high outliers, with no marker cropping since
  it's one continuous Axes (not two split panels).
- `ctg_brk_num_norm` (breaks per Mbp) is on the x-axis, linear by default
  -- flip USE_LOG_X to True if it turns out to be heavily right-skewed.
- Points with max_min_diff > y_break are labeled with their cytoband
  (`cyto`) and a short gene summary.

Color = region, marker shape = inv_tag.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

# ------------------------------------------------------------------
# 2. Normalize ctg_brk_num by region length (breaks per Mbp)
# ------------------------------------------------------------------
region_length_mbp = (df["end"] - df["start"]) / 1e6
df["ctg_brk_num_norm"] = df["ctg_brk_num"] / region_length_mbp

# ------------------------------------------------------------------
# 2b. Recalculate the allele-length difference directly from
#     min_allele/max_allele (bp -> Mbp), instead of relying on the
#     precomputed max_min_diff column.
# ------------------------------------------------------------------
df["allele_diff_recalc"] = (df["max_allele"] - df["min_allele"]) / 1e6

# ------------------------------------------------------------------
# 3. Style mappings
# ------------------------------------------------------------------
color_map = {"euchromatin": "tab:blue", "SubTel": "tab:orange", "CenSat": "tab:green"}
region_order = ["euchromatin", "SubTel", "CenSat"]
inv_order = ["NoINV", "HasINV"]
inv_hatch = {"NoINV": None, "HasINV": "//"}


def log_hist_panel(ax, column, xlabel, title):
    """Plot overlaid per-region histograms of `column` on a log x-axis.
    Non-positive values (<=0) can't be log-scaled and are dropped from
    this view -- report how many were dropped so nothing goes missing
    silently.
    """
    vals_all = df[column]
    n_dropped = (vals_all <= 0).sum()

    positive = df[df[column] > 0]
    bins = np.logspace(
        np.log10(positive[column].min()),
        np.log10(positive[column].max()),
        30,
    )

    for region in region_order:
        vals = positive.loc[positive["region"] == region, column]
        if len(vals) == 0:
            continue
        ax.hist(
            vals, bins=bins, color=color_map[region], alpha=0.5,
            label=region, edgecolor="none",
        )

    ax.set_xscale("log")
    suffix = f"\n({n_dropped} zero/negative values excluded)" if n_dropped else ""
    ax.set_xlabel(xlabel + suffix)
    ax.set_ylabel("count")
    ax.set_title(title)
    ax.legend(frameon=True, fontsize=8)


def linear_hist_panel(ax, column, xlabel, title, n_bins=30):
    """Plot overlaid per-region histograms of `column` on a LINEAR x-axis
    (no exclusions needed -- zeros and negatives are fine on a linear scale).
    """
    bins = np.linspace(df[column].min(), df[column].max(), n_bins)

    for region in region_order:
        vals = df.loc[df["region"] == region, column]
        if len(vals) == 0:
            continue
        ax.hist(
            vals, bins=bins, color=color_map[region], alpha=0.5,
            label=region, edgecolor="none",
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")
    ax.set_title(title)
    ax.legend(frameon=True, fontsize=8)


# ------------------------------------------------------------------
# 4. Figure with three panels
# ------------------------------------------------------------------
fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(19, 6))

# --- Panel A: grouped bar counts of region x inv_tag ---
counts = df.groupby(["region", "inv_tag"]).size().unstack(fill_value=0)
counts = counts.reindex(index=region_order, columns=inv_order, fill_value=0)

x = np.arange(len(region_order))
bar_width = 0.35

for i, inv_tag in enumerate(inv_order):
    offset = (i - 0.5) * bar_width
    bars = ax_a.bar(
        x + offset, counts[inv_tag], width=bar_width,
        color=[color_map[r] for r in region_order],
        hatch=inv_hatch[inv_tag], edgecolor="black", linewidth=0.6,
        label=inv_tag,
    )
    ax_a.bar_label(bars, padding=2, fontsize=8)

ax_a.set_xticks(x)
ax_a.set_xticklabels(region_order)
ax_a.set_ylabel("count")
ax_a.set_title("Region counts by inversion status")

hatch_legend = [
    Patch(facecolor="white", edgecolor="black", hatch=inv_hatch[t], label=t)
    for t in inv_order
]
ax_a.legend(handles=hatch_legend, loc="upper right", frameon=True)

# --- Panel B: distribution of ctg_brk_num_norm by region ---
log_hist_panel(
    ax_b, "ctg_brk_num_norm",
    "ctg_brk_num normalized by region length (breaks / Mbp, log scale)",
    "Distribution of contig breaks per Mbp, by region",
)

# --- Panel C: distribution of recalculated allele-length diff by region
#     (log scale) ---
log_hist_panel(
    ax_c, "allele_diff_recalc",
    "max_allele - min_allele (Mbp, log scale)",
    "Distribution of recalculated allele-length difference, by region",
)

plt.tight_layout()
plt.savefig(f'{VOL28}/SVREF/hq_asm/CHM13/DongAhnSD_Var/sites_MaxBrks.svg')
plt.show()
