#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 5/26/26
'''

import gzip
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from helpers.Constants import *

EXPR_FILE = f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/Iris_eQTL/TMM.expression.bed.gz"
VCF_FILE = f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/Iris_eQTL/strict_set.maf05.vcf.gz"

# GENE_ID = "ENSG00000184227.8"
# GENE_NAME = "ACOT1"

# GENE_ID = "ENSG00000119673.15"
# GENE_NAME = "ACOT2"

# GENE_ID = "ENSG00000177465.5"
# GENE_NAME = "ACOT4"

# GENE_ID = "ENSG00000268500.6"
# GENE_NAME = "SIGLEC5"

# GENE_ID = "ENSG00000254415.3"
# GENE_NAME = "SIGLEC14"

# GENE_ID = "ENSG00000128383.13"
# GENE_NAME = "APOBEC3A"

GENE_ID = "ENSG00000179750.16"
GENE_NAME = "APOBEC3B"

# SV_ID = "chr14-73529736-DEL-27727"
# SV_ID = "chr19-51630516-DEL-16424"
SV_ID = "chr22-38962276-DEL-29936"

OUT_FILE = f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/Iris_eQTL/Figures/{GENE_NAME}.{SV_ID}.boxplot"

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none',
"font.family": "sans-serif",
"font.sans-serif": "Arial"
}
mpl.rcParams.update(new_rc_params)

def get_expression_for_gene(expr_file, gene_id):
    with gzip.open(expr_file, "rt") as f:
        header = f.readline().strip().split("\t")
        sample_names = header[4:]

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if fields[3] == gene_id:
                expr_values = fields[4:]
                expr_dict = {}
                for s, v in zip(sample_names, expr_values):
                    try:
                        expr_dict[s] = float(v)
                    except ValueError:
                        expr_dict[s] = np.nan
                return expr_dict

    raise ValueError(f"Gene ID not found: {gene_id}")


def classify_genotype(gt_str):
    if gt_str is None:
        return None

    gt = gt_str.split(":")[0].replace("|", "/")

    if gt in {"./.", ".", ".|."}:
        return None
    if gt == "0/0":
        return "AA"
    if gt in {"0/1", "1/0"}:
        return "AB"
    if gt == "1/1":
        return "BB"
    return None


def get_sv_genotypes(vcf_file, sv_id):
    sample_names = None

    with gzip.open(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                continue

            fields = line.rstrip("\n").split("\t")
            if fields[2] != sv_id:
                continue

            ref = fields[3]
            alt = fields[4]
            sample_fields = fields[9:]

            gt_dict = {}
            for sample, sample_info in zip(sample_names, sample_fields):
                gt_class = classify_genotype(sample_info)
                if gt_class is not None:
                    gt_dict[sample] = gt_class

            return ref, alt, gt_dict

    raise ValueError(f"SV ID not found: {sv_id}")


def make_plot(expr_dict, gt_dict, gene_name, out_file):
    groups = {"AA": [], "AB": [], "BB": []}

    for sample, gt_class in gt_dict.items():
        if sample in expr_dict and not np.isnan(expr_dict[sample]):
            groups[gt_class].append(expr_dict[sample])

    order = ["AA", "AB", "BB"]
    data = [groups[g] for g in order]

    fig, ax = plt.subplots(figsize=(3.5, 4))

    # match the style in your example more closely
    box_width = 0.42
    jitter_width = 0.075  # keep points clearly inside the box

    ax.boxplot(
        data,
        positions=[1, 2, 3],
        widths=box_width,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="#f4a259", linewidth=1.4),
        boxprops=dict(facecolor="white", edgecolor="black", linewidth=1.1),
        whiskerprops=dict(color="black", linewidth=1.0),
        capprops=dict(color="black", linewidth=1.0),
    )

    rng = np.random.default_rng(123)
    colors = {"AA": "green", "AB": "orange", "BB": "blue"}

    for i, g in enumerate(order, start=1):
        y = np.array(groups[g], dtype=float)
        if len(y) == 0:
            continue

        x = i + rng.uniform(-jitter_width, jitter_width, size=len(y))

        ax.scatter(
            x,
            y,
            s=5,
            color=colors[g],
            alpha=0.8,
            edgecolors="none",
            zorder=3,
            rasterized = True,
        )
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["AA", "AB", "BB"], fontsize=12)
    ax.set_xlabel("Genotype", fontsize=13)
    ax.set_ylabel("Expression", fontsize=13)
    ax.set_title(gene_name, fontsize=14, pad=10)

    ax.tick_params(axis="y", labelsize=12)
    ax.tick_params(axis="x", labelsize=12)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    plt.tight_layout()
    plt.savefig(f'{out_file}.png', dpi=600, bbox_inches="tight")
    plt.savefig(f'{out_file}.svg')
    plt.close()


def main():
    expr_dict = get_expression_for_gene(EXPR_FILE, GENE_ID)
    ref, alt, gt_dict = get_sv_genotypes(VCF_FILE, SV_ID)

    make_plot(
        expr_dict=expr_dict,
        gt_dict=gt_dict,
        gene_name=GENE_NAME,
        out_file=OUT_FILE
    )


if __name__ == "__main__":
    main()

if __name__ == '__main__':
    main()
