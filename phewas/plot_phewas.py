#!/usr/bin/env python3
"""
Generate a PheWAS Manhattan plot from a PheTK-ready results TSV.
Labels all points passing Bonferroni correction (0.05 / n_tests).
"""
import argparse
import matplotlib as mpl
import numpy as np
import pandas as pd
from phetk.plot import Plot


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--phewas-results", required=True,
                   help="PheTK-ready TSV produced by format_results.py")
    p.add_argument("--output", default=None,
                   help="Output image path; defaults to {gene}_{SVID}.png or {SVID}.png")
    p.add_argument("--sv-info", default=None,
                   help="Optional TSV with columns: ID, AF, gene, loc")
    p.add_argument("--title", default=None,
                   help="Plot title; overrides the title derived from --sv-info/SVID")
    p.add_argument("--label-only-bonferroni", action="store_true",
                   help="Label only Bonferroni-significant points; takes precedence over --manhattan-label-count")
    p.add_argument("--manhattan-label-count", type=int, default=15,
                   help="Number of top associations to label when --label-only-bonferroni is not set")
    p.add_argument("--min-case-carrier-ct", type=int, default=5,
                   help="Minimum CASE_HET_A1_CT + CASE_HOM_A1_CT to include in plot")
    p.add_argument("--min-carrier-ct-above-bonferroni-only", action="store_true",
                   help="Apply --min-case-carrier-ct only to points passing the Bonferroni "
                        "threshold; points below the threshold are not filtered by carrier count")
    return p.parse_args()


def main():
    args = parse_args()

    import os

    mpl.rcParams.update(
        {
            'text.usetex': False, 
            'svg.fonttype': 'none', 
            'font.family': 'sans-serif', 
            'font.sans-serif': 'Arial'
        }
    )
    sv_id = os.path.basename(args.phewas_results)
    sv_id = sv_id[:sv_id.rindex('.')] # strip ext

    gene = None
    title = sv_id
    if args.sv_info:
        sv_info_df = pd.read_csv(args.sv_info, sep="\t", dtype=str)
        row = sv_info_df[sv_info_df["ID"] == sv_id]
        if len(row) == 1:
            gene = row.iloc[0]["gene"]
            loc  = row.iloc[0]["loc"]
            hom_ref = float(row.iloc[0]["0/0"])
            het     = float(row.iloc[0]["1/0"])
            hom_alt = float(row.iloc[0]["1/1"])
            total   = hom_ref + het + hom_alt
            af = f"{(het + 2 * hom_alt) / (2 * total):.4f}" if total > 0 else "NA"
            title = f"{gene} {sv_id} ({loc}, AF={af})"

    if args.title:
        title = args.title

    output = args.output or (f"{gene}_{sv_id}.png" if gene else f"{sv_id}.png")

    df = pd.read_csv(args.phewas_results, sep="\t")

    # Bonferroni uses all phenotypes with a valid category (before carrier filter)
    df_all = df.dropna(subset=["phecode_category"])
    n_tests = len(df_all)
    bonferroni_threshold = 0.05 / n_tests if n_tests > 0 else 0.05

    # Points shown in plot also require sufficient alt-allele cases
    carrier_ct = df_all["CASE_HET_A1_CT"] + df_all["CASE_HOM_A1_CT"]
    if args.min_carrier_ct_above_bonferroni_only:
        # No carrier-count filter below the Bonferroni threshold (min 0); only
        # significant points must meet --min-case-carrier-ct
        passes_bonferroni = df_all["p_value"] < bonferroni_threshold
        df_plot = df_all[~passes_bonferroni | (carrier_ct >= args.min_case_carrier_ct)].copy()
    else:
        df_plot = df_all[carrier_ct >= args.min_case_carrier_ct].copy()

    # Annotate labels with OR and case/control carrier counts
    def fmt_label(row):
        try:
            or_val = f"OR={np.exp(row['beta']):.2f}" if pd.notna(row.get("beta")) else "OR=NA"
        except Exception:
            or_val = "OR=NA"
        def ct(col):
            v = row.get(col)
            return int(v) if pd.notna(v) else 0
        cases = ct("CASE_HET_A1_CT") + ct("CASE_HOM_A1_CT") + ct("CASE_NON_A1_CT")
        controls = ct("CTRL_HET_A1_CT") + ct("CTRL_HOM_A1_CT") + ct("CTRL_NON_A1_CT")
        return f"{row['phecode_string']} ({or_val}, cases={cases}, controls={controls})"

    df_plot["phecode_string"] = df_plot.apply(fmt_label, axis=1)

    df_plot.to_csv("phewas_results_for_plot.tsv", sep="\t", index=False)

    neg_log_bonferroni = -np.log10(bonferroni_threshold)
    # Ensure y-axis extends at least to the Bonferroni line (+ buffer)
    max_neg_log_p = df_plot["neg_log_p_value"].dropna().max()
    max_neg_log_p = max_neg_log_p if pd.notna(max_neg_log_p) else 0
    y_limit = int(max(max_neg_log_p, neg_log_bonferroni) * 1.1) + 1

    if args.label_only_bonferroni:
        sig = df_plot[df_plot["p_value"] < bonferroni_threshold]
        label_values = sig["phecode"].dropna().astype(str).tolist()
        label_count = len(label_values)
        print(f"n_tests={n_tests}, Bonferroni threshold={bonferroni_threshold:.2e}, "
              f"significant in plot={label_count}", flush=True)
    else:
        label_values = "p_value"
        label_count = args.manhattan_label_count
        print(f"n_tests={n_tests}, Bonferroni threshold={bonferroni_threshold:.2e}, "
              f"labelling top {label_count}", flush=True)

    phewas_plot = Plot("phewas_results_for_plot.tsv", bonferroni=neg_log_bonferroni)
    phewas_plot.manhattan(
        title=title,
	title_text_size=14,
        label_values=label_values,
        label_count=label_count,
        y_limit=y_limit,
        output_file_path=output,
        save_plot=True,
        label_size=12,
	axis_text_size=13,
	show_legend=False,
	#label_split_threshold=10
    )
    print(f"Manhattan plot written to {output}", flush=True)


if __name__ == "__main__":
    main()
