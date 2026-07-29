#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 7/24/26
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

df = pd.read_csv(f'{VOL28}/SVREF/locityper_inAoU/aou_sr_gts/phewas_results/INSR_chr19-7258352-DEL-930_all_firth_results_SV_snps_minAC100.tsv', sep='\t')

sv_id = 'chr19-7258352-DEL-930'
gwas_snp_id = 'chr19_7257979_C_G'

df['logP'] = -np.log10(df['P'])
df['POS_Mb'] = df['POS'] / 1e6

sv_row = df[df['ID'] == sv_id].iloc[0]
gwas_row = df[df['ID'] == gwas_snp_id].iloc[0]

# Split SNPs vs SV
is_sv = df['ID'] == sv_id
snps = df[~is_sv].copy()

fig, ax = plt.subplots(figsize=(11, 7))

# Color by R2 (LD with the SV)
cmap = plt.cm.RdYlBu_r
norm = Normalize(vmin=0, vmax=1)

sc = ax.scatter(snps['POS_Mb'], snps['logP'], c=snps['UNPHASED_R2'], cmap=cmap, norm=norm, rasterized=True,
                s=45, edgecolor='black', linewidth=0.3, zorder=3, alpha=0.9)

# Plot the SV itself as a diamond, larger, distinctly marked
ax.scatter(sv_row['POS_Mb'], sv_row['logP'], marker='D', s=260, c='r', rasterized=True,
           edgecolor='white', linewidth=1.2, zorder=6, label='Structural variant (SV)')

# Highlight the GWAS SNP
ax.scatter(gwas_row['POS_Mb'], gwas_row['logP'], marker='*', s=550, facecolor='none', rasterized=True,
           edgecolor='limegreen', linewidth=2.2, zorder=7)

# Label the GWAS SNP
ax.annotate(gwas_snp_id, (gwas_row['POS_Mb'], gwas_row['logP']),
            xytext=(15, 18), textcoords='offset points', fontsize=9, fontweight='bold',
            color='darkgreen', arrowprops=dict(arrowstyle='-', color='darkgreen', lw=0.8))

# Label the SV
ax.annotate(sv_id, (sv_row['POS_Mb'], sv_row['logP']),
            xytext=(15, -22), textcoords='offset points', fontsize=9, fontweight='bold',
            color='black', arrowprops=dict(arrowstyle='-', color='black', lw=0.8))

# Label other SNPs in high LD (R2 >= 0.8) with the SV, excluding the GWAS SNP (already labeled) and SV itself
high_ld = snps[(snps['UNPHASED_R2'] >= 0.8) & (snps['ID'] != gwas_snp_id)].copy()
# high_ld = snps[(snps['UNPHASED_R2'] == 1) & (snps['ID'] != gwas_snp_id)].copy()
high_ld = high_ld.sort_values('logP', ascending=False)

# To avoid overplotting labels, offset them progressively
offsets = [(20, 6), (20, -14), (-70, 10), (-70, -14), (20, 26), (-70, 26), (20, -30), (-70, -30)]
for i, (_, row) in enumerate(high_ld.iterrows()):
    ax.scatter(row['POS_Mb'], row['logP'], s=90, facecolor='none', edgecolor='purple',rasterized=True,
               linewidth=1.6, zorder=5)
    off = offsets[i % len(offsets)]
    # ax.annotate(row['ID'], (row['POS_Mb'], row['logP']),
    #             xytext=off, textcoords='offset points', fontsize=7.5, color='purple',
    #             arrowprops=dict(arrowstyle='-', color='purple', lw=0.6))

print(f"Number of SNPs with R2>=0.8 with SV (excluding GWAS SNP labeled separately): {len(high_ld)}")

# Colorbar for R2
cbar = fig.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label(r'LD ($R^2$) with SV', fontsize=10)

ax.set_xlabel('Position (Mb)', fontsize=11)
ax.set_ylabel(r'$-\log_{10}(P)$', fontsize=11)
# ax.set_title('Regional association plot: phecode_172, chr16 locus\ncolored by LD ($R^2$) with SV chr16-90023596-DEL-364', fontsize=12)

# Genome-wide significance line (5e-8)
gw_line = -np.log10(5e-8)
ax.axhline(gw_line, color='grey', linestyle='--', linewidth=1, zorder=1)
ax.text(ax.get_xlim()[1], gw_line, '  P=5e-8', va='center', ha='left', fontsize=8, color='grey')

legend_elements = [
    Line2D([0], [0], marker='D', color='r', markerfacecolor='r', markeredgecolor='r',
           markersize=12, label=f'SV ({sv_id})'),
    Line2D([0], [0], marker='*', color='w', markerfacecolor='none', markeredgecolor='limegreen',
           markersize=18, markeredgewidth=2, label=f'Tagging SNP ({gwas_snp_id})'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='none', markeredgecolor='purple',
           markersize=10, markeredgewidth=1.6, label=r'Other SNPs in strong LD with SV'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=8.5, framealpha=0.9)

ax.grid(alpha=0.2, zorder=0)
plt.tight_layout()
plt.savefig(f'{VOL28}/SVREF/locityper_inAoU/aou_sr_gts/phewas_results/INSR_snp_plink_chr19-7258352-DEL-930_pval.png', dpi=300)
plt.savefig(f'{VOL28}/SVREF/locityper_inAoU/aou_sr_gts/phewas_results/INSR_snp_plink_chr19-7258352-DEL-930_pval.svg')
print("saved")
plt.show()