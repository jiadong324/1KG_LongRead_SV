#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 3/2/26
'''
import json
import gzip
import math
import joblib as jl
import pysam
import numpy as np
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from functools import reduce
from adjustText import adjust_text
import matplotlib as mpl
from scipy import stats
from intervaltree import IntervalTree

from helpers.Constants import *
from helpers.VcfFunc import read_svsize_bed
from helpers.Functions import nonoverlap_intervals

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none',
"font.family": "sans-serif",
"font.sans-serif": "Arial"
}
mpl.rcParams.update(new_rc_params)


def final_hete_pcrt_fig3a(workdir, input_pcrt):
    bi_hete = pd.read_csv(f'{workdir}/SV_HETE/{input_pcrt}', sep='\t', header=0)
    bi_hete['group'] = ['Biallelic' for _ in range(len(bi_hete))]

    print('#Biallelic', len(bi_hete))
    print('#Genic biallelic', len(bi_hete.loc[bi_hete['gene'] != '.']))
    print('#CDS biallelic', len(bi_hete.loc[bi_hete['location'] == 'CDS']))
    print('#UTR biallelic', len(bi_hete.loc[bi_hete['location'].str.contains('UTR')]))

    legends_1 = [Line2D([0], [0], label='CDS', color='#984ea3', ls='None', marker='o'),
                 Line2D([0], [0], label='3UTR', color='#377eb8', ls='None', marker='o'),
                 Line2D([0], [0], label='5UTR', color='#f781bf', ls='None', marker='o'),
                 # Line2D([0], [0], label='REG', color='#bf812d', ls='None', marker='o'),
                 Line2D([0], [0], label='Intergenic SV', markeredgecolor='grey', markerfacecolor='None', ls='None',
                        marker='o'), ]

    # pathogenic, cds_tr_hete, utr = load_tr_hete()
    # pathogenic, cds_tr_hete, utr = load_hg38_tr_hete()

    pathogenic = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/final_patho_stats.tsv',
                             sep='\t', header=0)
    cds_tr_hete = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/final_cds_stats.tsv',
                              sep='\t', header=0)
    utr = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/final_utr_stats.tsv', sep='\t',
                      header=0)

    trs = pd.concat([cds_tr_hete, utr])

    sd_hete = load_sd_gene_hete()

    print('#SD genes', len(sd_hete))
    print('\t#SD genes (%HET>0.5)', len(sd_hete.loc[sd_hete['all_het'] > 0.5]),
          len(sd_hete.loc[sd_hete['all_het'] > 0.5]) / len(sd_hete))
    print('#Pathogenic sites', len(pathogenic))
    print('\t#Pathogenic (%HET>0.5)', len(pathogenic.loc[pathogenic['all_het'] > 0.5]),
          len(pathogenic.loc[pathogenic['all_het'] > 0.5]) / len(pathogenic))
    print('#CDS sites', len(cds_tr_hete.loc[cds_tr_hete['context'] == 'CDS']))
    print('\t#CDS (%HET>0.5)', len(cds_tr_hete.loc[(cds_tr_hete['all_het'] > 0.5) & (cds_tr_hete['context'] == 'CDS')]),
          len(cds_tr_hete.loc[(cds_tr_hete['all_het'] > 0.5) & (cds_tr_hete['context'] == 'CDS')]) / len(cds_tr_hete))
    print('#UTR sites', len(utr))
    print('\t#UTR (%HET>0.5)', len(utr.loc[utr['all_het'] > 0.5]), len(utr.loc[utr['all_het'] > 0.5]) / len(utr))

    fig, axes = plt.subplots(1, 6, sharey=True, figsize=(14, 7))

    bi_hete = bi_hete.loc[~bi_hete['chrom'].isin(['chrX', 'chrY'])]

    for ele, color in zip(['3UTR', '5UTR', 'CDS'], ['#377eb8', '#f781bf', '#984ea3']):
        # print(ele, len(bi_hete.loc[bi_hete['location'] == ele]))
        if ele in ['.', 'REG']:
            sns.stripplot(data=bi_hete.loc[bi_hete['location'] == ele], y='all_het', x='group',
                          facecolors='none', jitter=True, rasterized=True, edgecolor='grey', linewidth=1, ax=axes[0])
        else:
            sns.stripplot(data=bi_hete.loc[bi_hete['location'] == ele], y='all_het', x='group',
                          color=color, jitter=True, rasterized=True, ax=axes[0])

    # sns.stripplot(data=bi_hete.loc[(bi_hete['gene'] != '.')&(~bi_hete['location'].isin(['CDS', '5UTR', '3UTR']))], y='all_het', x='group',
    #               facecolors='none', jitter=True, rasterized=True, edgecolor='#ff7f00', linewidth=1, ax=axes[0])

    sns.stripplot(data=pathogenic, x='context', y='all_het', rasterized=True, color='#e41a1c', jitter=True, ax=axes[1])

    text_path = []
    for collection in axes[1].collections:
        offsets = collection.get_offsets()
        for (x, y), label, het_pcrt in zip(offsets, pathogenic["locus"], pathogenic['all_het']):
            if het_pcrt > 0.5:
                text_path.append(axes[1].text(x, y, label.split('_')[0], fontsize=9))
    adjust_text(text_path, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='k', alpha=0.5), style='italic', ax=axes[1])

    for i, ele in enumerate(['CDS', '3UTR', '5UTR']):
        i += 2
        this_tr = trs.loc[trs['context'] == ele]
        print(ele, len(this_tr))

        sns.stripplot(data=this_tr, x='context', y='all_het', rasterized=True, color='#4daf4a',
                      # hue='gene_type', hue_order=cds_order, palette=cds_color,
                      ax=axes[i], jitter=True)

        text_tr = []
        for collection in axes[i].collections:
            offsets = collection.get_offsets()
            for (x, y), label, het_pcrt in zip(offsets, this_tr["gene"], this_tr['all_het']):
                if ele == 'CDS' and het_pcrt >= 0.5:
                    text_tr.append(axes[i].text(x, y, label, fontsize=9))
                elif het_pcrt > 0.9:
                    text_tr.append(axes[i].text(x, y, label, fontsize=9))
        adjust_text(text_tr, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='k', alpha=0.5), style='italic',
                    ax=axes[i])

    sns.stripplot(data=sd_hete, x='context', y='all_het', color='#4daf4a', rasterized=True, ax=axes[5], jitter=True)

    text_sd = []
    for collection in axes[5].collections:
        offsets = collection.get_offsets()
        for (x, y), label, het_pcrt in zip(offsets, sd_hete["locus"], sd_hete['all_het']):
            if het_pcrt > 0.65:
                text_sd.append(axes[5].text(x, y, label.split('_')[0], fontsize=9))
    adjust_text(text_sd, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='k', alpha=0.5), style='italic', ax=axes[5])

    axes[0].set_xlabel('Biallelic SVs', fontsize=13)
    axes[0].legend(handles=legends_1)

    for i, ax in enumerate(axes):
        ax.spines[['right', 'top']].set_visible(False)
        ax.set_xlabel('')
        if i > 1:
            ax.legend('', frameon=False)
            ax.spines[['left', 'right', 'top']].set_visible(False)

        ax.axhline(y=.5, ls='--', lw=1.5)
        ax.set_ylabel('Observed heterozygosity (%)', fontsize=14)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1], labels=[0, 0.25, 0.5, 0.75, 1])

    fig.tight_layout()
    fig.savefig(f'{VOL28}/SVREF/all_cohorts/Figures/sv_BpClust_het_pcrt_hg38.svg')
    plt.show()

def genetic_diversity(counts):
    counts = np.array(counts)
    freqs = counts / counts.sum()
    # freqs = np.array(freqs)
    freqs = freqs[freqs > 0]

    # Expected heterozygosity
    He = 1 - np.sum(freqs ** 2)

    # Effective allele number
    # Ne = 1 / np.sum(freqs ** 2)

    return He

def normalized_shannon(counts):
    counts = np.array(counts)
    p = counts / counts.sum()
    p = p[p > 0]  # remove zeros

    H = -np.sum(p * np.log(p))
    H_max = np.log(len(p))

    return H, H / H_max


def get_biallelic_sv_hete_pcrt(workdir, input_vcf, output_pcrt):
    # vcf = pysam.VariantFile(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/yang_annot/unique_sv_fill_tags.vcf.gz', 'rb')
    vcf = pysam.VariantFile(f'{workdir}/Annot/{input_vcf}')
    # sv_annot = pd.read_csv(f'{workdir}/Annot/gt_af_anno.tsv.gz',sep='\t',index_col=['ID'])
    # sv_ele_annot = pd.read_csv(f'{workdir}/yang_annot/allreg_sv_T2T.txt',sep='\t',index_col=['ID'])
    sv_ele_annot = pd.read_csv(f'{workdir}/yang_annot/allreg_sv.txt', sep='\t', index_col=['ID'])
    # sv_annot.fillna('.', inplace=True)


    # sample_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/sample_group.txt', sep='\t', names=['sample', 'group'])
    sample_tbl = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/sample_group.txt', sep='\t', names=['sample', 'group'])
    afr_samples = len(sample_tbl.loc[sample_tbl['group']=='AFR'])
    nonafr_samples = len(sample_tbl.loc[sample_tbl['group'] == 'NonAFR'])

    het_pcrt_list = []
    for rec in vcf.fetch():
        af_ac_het = rec.info['AC_Het_AFR'][0]
        nonaf_ac_het = rec.info['AC_Het_NonAFR'][0]

        # gene = sv_annot.at[rec.id, 'Gene']
        context = sv_ele_annot.at[rec.id, 'Location_GENCODE']
        gene = sv_ele_annot.at[rec.id, 'gene']
        missing_al = 0

        for sample, val in rec.samples.items():
            h1, h2 = val.get('GT')
            if h1 is None:
                missing_al += 1

            if h2 is None:
                missing_al += 1

        if missing_al / len(rec.samples) * 2 > 0.2:
            continue

        het_pcrt_list.append([rec.chrom, rec.start, rec.id, rec.info['SVTYPE'], rec.info['SVLEN'], gene, context, rec.info['AF'][0], rec.info['MAF'],
                              af_ac_het / afr_samples, nonaf_ac_het / nonafr_samples, (af_ac_het + nonaf_ac_het)/(afr_samples + nonafr_samples)])

    df_het = pd.DataFrame(het_pcrt_list, columns=['chrom', 'start', 'svid', 'svtype', 'svlen', 'gene', 'location', 'af', 'maf', 'afr_het_pcrt', 'nonafr_het_pcrt', 'all_het'])
    # df_het.to_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/SV_HETE/unique_sv_hete_pcrt.tsv', sep='\t', header=True, index=False)
    # df_het.to_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/SV_HETE/biallelic_sv_hete_pcrt.tsv', sep='\t', header=True, index=False)
    df_het.to_csv(f'{workdir}/SV_HETE/{output_pcrt}',sep='\t', header=True, index=False)

def load_sd_gene_hete():
    sd_hete = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/SV_InSD_Blck/results/ly_gene_het.tsv',sep='\t', header=0)
    sd_hete = sd_hete.loc[~sd_hete['locus'].str.contains('LOC')]
    sd_hete['context'] = ['SD Genes' for _ in range(len(sd_hete))]
    return sd_hete

def update_patho_tr_type(row):
    tag = 'STR'
    if row['unit_size'] == 3:
        tag = 'STR-Triplet'
    elif row['unit_size'] > 6:
        tag = 'VNTR'
    return tag

def main():
    final_hete_pcrt_fig3a(f'{VOL28}/SVREF/hq_asm/GRCh38/pct09_dist500/complete',
                   'biallelic_BpClust_hete_pcrt.tsv')

    ## Calculate the het percent for biallelic SVs
    get_biallelic_sv_hete_pcrt(f'{VOL28}/SVREF/hq_asm/GRCh38/pct09_dist500/complete',
                               'biallelic_BpClust.vcf', 'biallelic_BpClust_hete_pcrt.tsv')

if __name__ == '__main__':
    main()
