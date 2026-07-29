#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 7/20/26
'''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import matplotlib as mpl

from helpers.Constants import *

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none',
"font.family": "sans-serif",
"font.sans-serif": "Arial"
}
mpl.rcParams.update(new_rc_params)

# plt.rcParams['svg.fonttype'] = 'none'
df = pd.read_csv('ASCC1_chr10-72100639-DEL-2740_all_firth_results_SV_snps_minAC100.tsv', sep='\t')

# Keep text as real, editable text elements (not outlined paths) in exported SVGs
plt.rcParams['svg.fonttype'] = 'none'

sv_id = 'chr10-72100639-DEL-2740'
gwas_snp_id = 'chr10_72100529_T_A'

sv_row = df[df['ID'] == sv_id].iloc[0]
sv_pos = sv_row['POS']

df['logP'] = -np.log10(df['P'])
df['dist_kb'] = (df['POS'] - sv_pos) / 1000.0

snps = df[df['ID'] != sv_id].copy()
gwas_row = df[df['ID'] == gwas_snp_id].iloc[0]
gwas_dist = (gwas_row['POS'] - sv_pos) / 1000.0

fig, ax = plt.subplots(figsize=(7, 5))

cmap = plt.cm.viridis
norm = Normalize(vmin=snps['logP'].min(), vmax=snps['logP'].max())

sc = ax.scatter(snps['dist_kb'], snps['UNPHASED_R2'], c=snps['logP'], cmap=cmap, norm=norm, rasterized=True,
                s=40, edgecolor='black', linewidth=0.25, alpha=0.85, zorder=3)

# SV itself at distance 0, R2=1
ax.scatter(0, 1.0, marker='D', s=260, c='r', edgecolor='r', linewidth=1.2, rasterized=True,zorder=6)
ax.annotate(sv_id, (0, 1.0), xytext=(12, -18), textcoords='offset points',
            fontsize=9, fontweight='bold', color='black',
            arrowprops=dict(arrowstyle='-', color='black', lw=0.8))

# GWAS SNP highlighted
ax.scatter(gwas_dist, gwas_row['UNPHASED_R2'], marker='*', s=550, facecolor='none',rasterized=True,
           edgecolor='limegreen', linewidth=2.2, zorder=7)
ax.annotate(gwas_snp_id, (gwas_dist, gwas_row['UNPHASED_R2']), xytext=(12, 10),
            textcoords='offset points', fontsize=9, fontweight='bold', color='darkgreen',
            arrowprops=dict(arrowstyle='-', color='darkgreen', lw=0.8))

# Other high-LD SNPs (R2 >= 0.8), excluding GWAS SNP
high_ld = snps[(snps['UNPHASED_R2'] >= 0.8) & (snps['ID'] != gwas_snp_id)].sort_values('UNPHASED_R2', ascending=False)
offsets = [(15, 6), (15, -16), (-90, 8), (-90, -16), (15, 24), (-90, 24)]
for i, (_, row) in enumerate(high_ld.iterrows()):
    ax.scatter(row['dist_kb'], row['UNPHASED_R2'], s=90, facecolor='none', edgecolor='purple',
               linewidth=1.6, zorder=5)
    off = offsets[i % len(offsets)]
    ax.annotate(row['ID'], (row['dist_kb'], row['UNPHASED_R2']), xytext=off, textcoords='offset points',
                fontsize=7.5, color='purple', arrowprops=dict(arrowstyle='-', color='purple', lw=0.6))

# R2 = 0.8 threshold line
ax.axhline(0.8, color='grey', linestyle='--', linewidth=1, zorder=1)
ax.text(ax.get_xlim()[1] if ax.get_xlim()[1] else 100, 0.8, '  $R^2$=0.8', va='bottom',
        ha='right', fontsize=8, color='grey')

cbar = fig.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label(r'$-\log_{10}(P)$', fontsize=10)

ax.set_xlabel('Distance to SV', fontsize=11)
ax.set_ylabel(r'LD with SV ($R^2$)', fontsize=11)
# ax.set_title('LD decay with distance from SV chr16-90023596-DEL-364\n(phecode_172, chr16 locus, colored by significance)', fontsize=12)
ax.set_ylim(-0.03, 1.05)

legend_elements = [
    Line2D([0], [0], marker='D', color='w', markerfacecolor='r', markeredgecolor='white',
           markersize=12, label='SV (2740 bp DEL)'),
    Line2D([0], [0], marker='*', color='w', markerfacecolor='none', markeredgecolor='limegreen',
           markersize=18, markeredgewidth=2, label='Tagging SNP (rs186180298)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='none', markeredgecolor='purple',
           markersize=10, markeredgewidth=1.6, label=r'Other SNPs $R^2\geq0.8$'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=8.5, framealpha=0.9)

ax.grid(alpha=0.2, zorder=0)
plt.tight_layout()
plt.savefig(f'{VOL28}/SVREF/locityper_inAoU/aou_sr_gts/phewas_results/ASCC1_snp_plink_chr10-72100639-DEL-2740.png', dpi=300)
plt.savefig(f'{VOL28}/SVREF/locityper_inAoU/aou_sr_gts/phewas_results/ASCC1_snp_plink_chr10-72100639-DEL-2740.svg')
print("saved")
plt.show()
