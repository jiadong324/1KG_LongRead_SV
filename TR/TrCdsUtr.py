#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 3/4/26
'''


import math
import matplotlib as mpl
from adjustText import adjust_text
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

VOL28 = './'

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none',
"font.family": "sans-serif",
"font.sans-serif": "Arial"
}
mpl.rcParams.update(new_rc_params)

def str_length_cv_fig3b():
    patho, cds, utr, cv_hue, cv_color = load_variation_db('length')
    print('#CDS sites:', len(cds))
    print('#UTR sites:', len(utr))
    print('#Pathogenic sites:', len(patho))

    utr_str = utr.loc[(utr['TR_TAG'].str.contains('STR')) & (utr['weighted_cv'] > 0) & (utr['context']=='5UTR')]
    utr_str['plot_class'] = utr_str.apply(lambda row: 'Outlier' if row['weighted_cv']>=30 else 'Normal', axis=1)
    patho_str = patho.loc[(patho['TR_TAG'].str.contains('STR')) & (patho['weighted_cv'] > 0)]
    patho_str['plot_class'] = ['Patho' for _ in range(len(patho_str))]

    combined = pd.concat([utr_str, patho_str])
    combined['log_score'] = combined.apply(lambda row: math.log10(row['weighted_cv'] + 1), axis=1)

    fig, ax = plt.subplots(figsize=(7, 5))
    sns.scatterplot(data=combined, x='all_het', y='log_score', hue='plot_class', hue_order=['Patho', 'Normal'],
                    palette=['#fdb462', '#66c2a5'], s=70, rasterized=True, ax=ax)

    sns.scatterplot(data=combined.loc[combined['plot_class']=='Outlier'], x='all_het', y='log_score',
                    color='#e41a1c', s=120, rasterized=True, ax=ax)
    text = []
    for idx, row in combined.iterrows():
        if row['weighted_cv'] >= 30:
            text.append(ax.text(row['all_het'], row['log_score'], row['gene'],
                                                 horizontalalignment='center', verticalalignment='bottom', size='small',
                                                 color='k'))
    adjust_text(text, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=ax)

    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/Figures/unstable_str_patho_variation_main.svg')

    plt.show()

def all_length_cv():
    patho, cds, utr, cv_hue, cv_color = load_variation_db('length')

    cds[['locus', 'length_std', 'length_avg', 'pct99', 'pct50', 'max_len', 'min_len', 'unique_alleles', 'allele_range', 'context', 'gene',
         'SegDup','TR_TAG', 'all_het']].to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/final_cds_stats.tsv',
                    sep='\t', header=True, index=False)

    utr[['trid', 'length_std', 'length_avg', 'pct99', 'pct50', 'max_len', 'min_len', 'unique_alleles', 'allele_range', 'context', 'gene',
         'SegDup', 'TR_TAG', 'all_het']].to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/final_utr_stats.tsv',
        sep='\t', header=True, index=False)

    patho['context'] = ['Patho' for _ in range(len(patho))]
    patho[['locus', 'length_std', 'length_avg', 'pct99', 'pct50', 'max_len', 'min_len', 'unique_alleles', 'allele_range', 'context', 'gene', 'TR_TAG', 'all_het']].to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/final_patho_stats.tsv',
        sep='\t', header=True, index=False)

    print('#CDS sites:', len(cds))
    print('#UTR sites:', len(utr))
    print('#Pathogenic sites:', len(patho))


    fig, axes = plt.subplots(2, 3, sharex=True, figsize=(15, 9))
    utr5 = utr.loc[utr['context']=='5UTR']
    fig_i = 0
    for ele, tmp_df in zip(['CDS', "5'UTR", 'Pathogenic'], [cds, utr5, patho]):
        tmp_df['log_score'] = tmp_df.apply(lambda row: math.log10(row['weighted_cv'] + 1), axis=1)

        print(ele)
        legends = []
        for color, tr_type in zip(['#bf812d', '#7570b3'], ['STR', 'VNTR']):
            # this_tr = tmp_df.loc[(tmp_df['TR_TAG'] == tr_type)&(tmp_df['weighted_cv'] > 0)]
            this_tr = tmp_df.loc[(tmp_df['TR_TAG'].str.contains(tr_type)) & (tmp_df['weighted_cv'] > 0)]
            legends.append(Line2D([0], [0], marker='o', color=color, ls='none', label=tr_type))

            this_tr_corr, _ = pearsonr(this_tr['all_het'], this_tr['log_score'])
            slope, intercept, r_value, p_value, std_err = stats.linregress(this_tr['all_het'], this_tr['log_score'])

            if tr_type != 'VNTR':
                sns.regplot(data=this_tr, x="all_het", y="log_score", color=color, ci=None, ax=axes[0][fig_i])
            else:
                sns.regplot(data=this_tr, x="all_het", y="log_score", color=color, ci=None, ax=axes[1][fig_i])

            # axes[fig_i].text(0, 4.5, f'Pearson = {corr:.2f}', fontsize=12)
            print(f'\t{tr_type}: {this_tr_corr}')
            print(f'\t{tr_type} R-squared: {r_value**2}')

        axes[0][fig_i].set_ylabel('')
        axes[1][fig_i].set_ylabel('')

        if fig_i == 0:
            text_vntr, text_str = [], []
            for idx, row in tmp_df.iterrows():
                if row['log_score'] > 1.25:
                    if row['TR_TAG'] != 'VNTR':
                        text_str.append(axes[0][fig_i].text(row['all_het'], row['log_score'], row['gene'], horizontalalignment='center',
                                        verticalalignment='bottom', size='small', color='k'))
                    else:
                        text_vntr.append(axes[1][fig_i].text(row['all_het'], row['log_score'], row['gene'],
                                                         horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))
            adjust_text(text_str, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=axes[0][fig_i])
            adjust_text(text_vntr, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=axes[1][fig_i])
            axes[0][fig_i].set_ylabel('Allele length variability score', fontsize=14)
            axes[1][fig_i].set_ylabel('Allele length variability score', fontsize=14)

        elif fig_i == 1:
            text_vntr, text_str = [], []
            for idx, row in tmp_df.iterrows():
                if row['log_score'] > 1.5:
                    if row['TR_TAG'] != 'VNTR':
                        text_str.append(axes[0][fig_i].text(row['all_het'], row['log_score'], row['gene'], horizontalalignment='center',
                                         verticalalignment='bottom', size='small', color='k'))
                    else:
                        text_vntr.append(axes[1][fig_i].text(row['all_het'], row['log_score'], row['gene'],
                                                         horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))
            adjust_text(text_str, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5),ax=axes[0][fig_i])
            adjust_text(text_vntr, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=axes[1][fig_i])

        elif fig_i == 2:
            text_vntr, text_str = [], []
            for idx, row in tmp_df.iterrows():
                if row['all_het'] > 0.05:
                    if row['TR_TAG'] != 'VNTR':
                        text_str.append(axes[0][fig_i].text(row['all_het'], row['log_score'], row['gene'], horizontalalignment='center',
                                         verticalalignment='bottom', size='small', color='k'))
                    else:
                        text_vntr.append(axes[1][fig_i].text(row['all_het'], row['log_score'], row['gene'],
                                                         horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))
            adjust_text(text_str, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5),ax=axes[0][fig_i])
            adjust_text(text_vntr, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=axes[1][fig_i])

        if fig_i == 0:
            axes[0][fig_i].legend(handles=legends, loc='lower right')

        axes[0][fig_i].spines[['top', 'right']].set_visible(False)
        axes[1][fig_i].spines[['top', 'right']].set_visible(False)
        axes[0][fig_i].set_xlabel('')
        axes[1][fig_i].set_xlabel('Structural heterozygosity', fontsize=14)
        fig_i += 1

    axes[0][0].set_yticks([0, 0.5, 1, 1.5, 2], labels=[0, 0.5, 1, 1.5, 2])
    axes[0][1].set_yticks([0, 0.5, 1, 1.5, 2, 2.5], labels=[0, 0.5, 1, 1.5, 2, 2.5])
    axes[1][0].set_yticks([0, 0.5, 1, 1.5, 2], labels=[0, 0.5, 1, 1.5, 2])
    axes[1][1].set_yticks([0, 0.5, 1, 1.5, 2, 2.5], labels=[0, 0.5, 1, 1.5, 2, 2.5])

    axes[0][2].set_yticks([0, 1, 2, 3, 4], labels=[0, 1, 2, 3, 4])
    axes[1][2].set_yticks([0, 1, 2, 3], labels=[0, 1, 2, 3])

    fig.tight_layout()
    # fig.savefig(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/Figures/het_var_corr.png', dpi=300)
    # fig.savefig(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/Figures/het_var_corr.svg')
    # fig.savefig(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/Figures/het_var_corr2.svg')
    plt.show()

    # cds_sort_extreme = cds.sort_values(by='extreme', ascending=False)
    # utr_sort_extreme = utr.sort_values(by='extreme', ascending=False)
    #
    cds.sort_values(by='weighted_cv', inplace=True, ascending=False)
    utr.sort_values(by='weighted_cv', inplace=True, ascending=False)


    unstable_sites = pd.concat([cds.loc[(cds['weighted_cv']>0)&(cds['unique_alleles']>2)], utr.loc[(utr['weighted_cv']>0)&(utr['unique_alleles']>2)]])

    unstable_sites.to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/unstable_cds_utr.tsv',
                    sep='\t', header=True, index=False)

    unstable_5utr = utr.loc[(utr['weighted_cv'] > 0)&(utr['unique_alleles']>2)&(utr['context']=='5UTR')]

    print('CDS unstable', len(cds.loc[(cds['weighted_cv'] > 0)&(cds['unique_alleles']>2)]))
    print('UTR unstable', len(utr.loc[(utr['weighted_cv'] > 0)&(utr['unique_alleles']>2)]))
    print('5UTR unstable', len(unstable_5utr))

    unstable_sites['single_gene'] = unstable_sites.apply(lambda row: row['gene'].split(',')[0], axis=1)
    unstable_sites[['single_gene']].to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/unstable_cds_utr_genes.tsv',
                    sep='\t', header=False, index=False)

    hyper_var = pd.concat([cds.loc[cds['weighted_cv']>=6], utr.loc[utr['weighted_cv']>=6]])
    # hyper_var.to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/hypervariable_cds_utr.tsv',
    #                 sep='\t', header=True, index=False)

    hyper_cds = hyper_var.loc[(hyper_var['context']=='CDS')&(hyper_var['unique_alleles']>2)]
    top_cds = hyper_cds[0: int(len(hyper_cds)*0.4)]
    top_cds[['trid', 'motifs','unique_alleles', 'extreme', 'allele_range', 'gene', 'context', 'TR_TAG']].to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/hypervariable_top_cds.tsv',
                     sep='\t', header=True, index=False)

    hyper_utrs = hyper_var.loc[(hyper_var['context'].str.contains('UTR'))&(hyper_var['unique_alleles']>2)]
    top10_utr = hyper_utrs[0: int(len(hyper_utrs)*0.2)]
    top10_utr[['trid', 'motifs','unique_alleles', 'extreme', 'allele_range', 'gene', 'context', 'TR_TAG']].to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/hypervariable_top_utr.tsv',
        sep='\t', header=True, index=False)

    top1_unstable_5utr = unstable_5utr[0: int(len(unstable_5utr) * 0.05)]
    top1_unstable_5utr[['trid', 'motifs', 'unique_alleles', 'extreme', 'allele_range', 'gene', 'context', 'TR_TAG']].to_csv(
        f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/unstable_top5pcrt_5utr.tsv',
        sep='\t', header=True, index=False)


    # hyper_var['single_gene'] = hyper_var.apply(lambda row: row['gene'].split(',')[0], axis=1)
    # hyper_var[['single_gene']].to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/hypervariable_cds_utr_genes.tsv',
    #                 sep='\t', header=False, index=False)

    print('CDS hypervariable', len(cds.loc[cds['weighted_cv'] >= 6]))
    print('UTR hypervariable', len(utr.loc[utr['weighted_cv'] >= 6]))

    print('CDS triplets hypervariable', len(cds.loc[(cds['weighted_cv'] >= 6)&(cds['TR_TAG']=='STR-Triplet')]))
    print('UTR triplets hypervariable', len(utr.loc[(utr['weighted_cv'] >= 6)&(utr['TR_TAG']=='STR-Triplet')]))

def update_context(row, bdry_cds):
    if 'CDS' in row['context']:
        if row['trid'] not in bdry_cds:
            return 'CDS'
        else:
            return 'BdryCDS'
    else:
        if '5UTR_PROM' in row['context'].split(',')[0]:
            return '5UTR'
        elif '3UTR' in row['context'].split(',')[0]:
            return '3UTR'
        else:
            return 'NA'

def update_tr_type(row, triplet_sites):
    tag = 'STR-Triplet' if row['trid'] in triplet_sites else 'STR'

    if tag != 'STR-Triplet':
        for ele in row['motif'].split(','):
            if len(ele) > 6:
                tag = 'VNTR'
    return tag

def update_patho_tr_type(row):
    tag = 'STR'
    if row['unit_size'] == 3:
        tag = 'STR-Triplet'
    elif row['unit_size'] > 6:
        tag = 'VNTR'
    return tag

def update_weight_cv(row):
    cv = row['length_std'] / row['length_avg']
    # weight = 1 + (row['max_len'] - row['pct50'])/row['pct50'] if row['pct50'] != 0 else 1
    # div = row['cn_diversity']
    # ex_index = (row['max_len'] - row['pct50'])/row['pct50'] if row['pct50'] != 0 else 0
    weight = 1 + abs(row['max_len'] - row['pct99']) / row['pct99'] if row['pct99'] != 0 else 1

    return weight * cv * 100
    # return weight * div * 100


def load_variation_db(stats_type):
    cds_tr_hete = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/svtr_trgt/final_genic_sites/locus_het_pcrt2.tsv', sep='\t', header=0)
    cds_tr_hete.rename(columns={"pcrt_99th": "pct99", "pcrt_50th": "pct50", "max_allele": 'max_allele_cn', "min_allele": 'min_allele_cn'},inplace=True)
    cds_tr_hete['gene'] = cds_tr_hete.apply(lambda row: row['trid'], axis=1)
    # cds_tr_hete['cv'] = cds_tr_hete.apply(lambda row: row['length_std'] * 100 / row['length_avg'], axis=1)
    # cds_tr_hete['tr'] = cds_tr_hete.apply(lambda row: row[''] / row['pct50'], axis=1)

    tmp_patho = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/svtr_trgt/final_genic_sites/pathogenic_het_pcrt.tsv',
        sep='\t', header=0)
    tmp_patho.rename(columns={"pcrt_99th": "pct99", "pcrt_50th": "pct50", "max_allele": 'max_len', "min_allele": 'min_len', "all_unique_alleles": "unique_alleles"}, inplace=True)

    tmp_patho = tmp_patho.loc[tmp_patho['context'] == 'Patho']
    tmp_patho['gene'] = tmp_patho.apply(lambda row: row['locus'].split('_')[0], axis=1)

    patho_vntr = ['DUX4', 'LPA', 'PLIN4', 'DRD4', 'CEL', 'ACAN', 'DSPP', 'TENT5A', 'PER3', 'IRF5', 'PHGR1', 'IVL', 'TCHH', 'CLEC4M', 'GP1BA', 'MUC1']
    unit_size_list = [3297, 5606, 99, 48, 33, 57, 9, 15, 54, 15, 33, 30, 18, 69, 39,60]

    mucus_genes = ['MUC2', 'MUC4', 'MUC6', 'MUC12', 'MUC17', 'MUC21', 'MUC22']
    mucus_genes_unit = [69, 48, 507, 84, 177, 45, 30]

    for gene, unit_size in zip(mucus_genes, mucus_genes_unit):
        cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'allele_range'] = cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'allele_range'] * unit_size
        cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'length_std'] = cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'length_std'] * unit_size
        cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'length_avg'] = cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'length_avg'] * unit_size
        cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'pct99'] = cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'pct99'] * unit_size
        cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'pct50'] = cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'pct50'] * unit_size
        cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'max_len'] = cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'max_allele_cn'] * unit_size
        cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'min_len'] = cds_tr_hete.loc[cds_tr_hete['gene'] == gene, 'min_allele_cn'] * unit_size

    cds_tr_hete.rename(columns={"all_unique_alleles": "unique_alleles"}, inplace=True)
    cds_tr_hete['weighted_cv'] = cds_tr_hete.apply(lambda row: update_weight_cv(row), axis=1)
    cds_tr_hete['extreme'] = cds_tr_hete.apply(lambda row: row['max_len'] - row['pct99'], axis=1)

    mucus = cds_tr_hete.loc[cds_tr_hete['gene'].isin(mucus_genes)]
    mucus['context'] = ['CDS' for _ in range(len(mucus))]
    mucus['TR_TAG'] = ['VNTR' for _ in range(len(mucus))]

    cds_sites = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_cds_stats/cds_{stats_type}_annot.tsv',sep='\t', header=0)

    cds_sites.rename(columns={"std": "length_std", "avg": "length_avg", "trid": "locus"}, inplace=True)
    # cds_sites['tr'] = cds_sites.apply(lambda row: row['pct99'] * 100 / row['pct50'], axis=1)
    cds_sites['weighted_cv'] = cds_sites.apply(lambda row: update_weight_cv(row), axis=1)
    cds_sites['cv'] = cds_sites.apply(lambda row: row['length_std'] / row['length_avg'], axis=1)
    cds_sites['extreme'] = cds_sites.apply(lambda row: row['max_len'] - row['pct99'], axis=1)

    ## add other VNTR CDS
    tmp_patho2 = cds_tr_hete.loc[cds_tr_hete['trid'].isin(patho_vntr)]

    for gene, unit_size in zip(patho_vntr, unit_size_list):
        tmp_patho2.loc[tmp_patho2['gene'] == gene, 'allele_range'] = tmp_patho2.loc[tmp_patho2['gene'] == gene, 'allele_range'] * unit_size
        tmp_patho2.loc[tmp_patho2['gene'] == gene, 'length_std'] = tmp_patho2.loc[tmp_patho2['gene'] == gene, 'length_std'] * unit_size
        tmp_patho2.loc[tmp_patho2['gene'] == gene, 'length_avg'] = tmp_patho2.loc[tmp_patho2['gene'] == gene, 'length_avg'] * unit_size
        tmp_patho2.loc[tmp_patho2['gene'] == gene, 'pct99'] = tmp_patho2.loc[tmp_patho2['gene'] == gene, 'pct99'] * unit_size
        tmp_patho2.loc[tmp_patho2['gene'] == gene, 'pct50'] = tmp_patho2.loc[tmp_patho2['gene'] == gene, 'pct50'] * unit_size
        tmp_patho2.loc[tmp_patho2['gene'] == gene, 'max_len'] = tmp_patho2.loc[tmp_patho2['gene'] == gene, 'max_allele_cn'] * unit_size
        tmp_patho2.loc[tmp_patho2['gene'] == gene, 'min_len'] = tmp_patho2.loc[tmp_patho2['gene'] == gene, 'min_allele_cn'] * unit_size

    ## add EP400
    tmp_patho3 = cds_sites.loc[cds_sites['locus'] == 'chr12_132062523_132062666_trsolve']
    tmp_patho3['unit_size'] = 3
    tmp_patho3['locus'] = 'EP400_chr12_132062523_132062666_trsolve'

    patho_columns = ['locus', 'unit_size', 'unique_alleles', 'allele_range', 'gene', 'length_std', 'length_avg', 'pct99', 'pct50', 'max_len', 'min_len', 'all_het']
    pathogenic = pd.concat([tmp_patho[patho_columns], tmp_patho2[patho_columns], tmp_patho3[patho_columns]])

    pathogenic['TR_TAG'] = pathogenic.apply(lambda row: update_patho_tr_type(row), axis=1)
    pathogenic['cv'] = pathogenic.apply(lambda row: row['length_std'] / row['length_avg'], axis=1)
    # pathogenic['tr'] = pathogenic.apply(lambda row: , axis=1)
    pathogenic['weighted_cv'] = pathogenic.apply(lambda row: update_weight_cv(row), axis=1)
    pathogenic['extreme'] = pathogenic.apply(lambda row: row['max_len'] - row['pct99'], axis=1)

    cds_sites = cds_sites.loc[(cds_sites['SegDup']=='OutSD')
                              &(~cds_sites['gene'].str.contains('ENSG'))
                              &(~cds_sites['gene'].str.contains('MUC'))
                              & (~cds_sites['context'].str.contains('UTR'))
                              & (cds_sites['context']=='CDS')
                              &(~cds_sites['gene'].isin(patho_vntr + ['EP400', 'TPGS2', 'LHX4', 'UBAP2', 'RYR3']))]

    # cds_sites.sort_values(by='allele_range')
    utr_sites = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_utr_stats/utr_{stats_type}_annot.tsv',sep='\t', header=0)
    utr_sites['tr'] = utr_sites.apply(lambda row: row['pct99'] * 100 / row['pct50'], axis=1)
    utr_sites.rename(columns={"std": "length_std", "avg": "length_avg"}, inplace=True)


    ## Exclude sites
    utr_sites = utr_sites.loc[(utr_sites['SegDup']=='OutSD')
                              &(~utr_sites['gene'].str.contains('ENSG'))]

    utr_sites['cv'] = utr_sites.apply(lambda row: row['length_std'] * 100 / row['length_avg'], axis=1)
    utr_sites['weighted_cv'] = utr_sites.apply(lambda row: update_weight_cv(row), axis=1)
    utr_sites['extreme'] = utr_sites.apply(lambda row: row['max_len'] - row['pct99'], axis=1)

    final_cds_sites = pd.concat([cds_sites, mucus])

    utr_sites = utr_sites.loc[(utr_sites['SegDup'] == 'OutSD')
                              & (~utr_sites['gene'].str.contains('ENSG'))
                              & (~utr_sites['gene'].str.contains('MUC'))
                              & (~utr_sites['gene'].isin(patho_vntr))]

    max_score = max(cds_sites['weighted_cv'].max(), utr_sites['weighted_cv'].max(), pathogenic['weighted_cv'].max())

    bins = [0, 10, 50, 500, max_score]
    names = ['#fee391', '#fe9929', '#cc4c02', '#990000']

    ## assign colors
    pathogenic['color'] = pd.cut(pathogenic['weighted_cv'], bins, labels=names)
    pathogenic['color'] = pathogenic.apply(lambda row: '#fff7bc' if row['cv']==0 else row['color'], axis=1)

    final_cds_sites['color'] = pd.cut(final_cds_sites['weighted_cv'], bins, labels=names)
    final_cds_sites['color'] = final_cds_sites.apply(lambda row: '#fff7bc' if row['cv'] == 0 else row['color'], axis=1)

    utr_sites['color'] = pd.cut(utr_sites['weighted_cv'], bins, labels=names)
    utr_sites['color'] = utr_sites.apply(lambda row: '#fff7bc' if row['cv'] == 0 else row['color'], axis=1)

    return pathogenic, final_cds_sites, utr_sites, bins, ['#fff7bc'] + names

def load_utr_cds_genes():
    sd_utr_cds = pd.read_csv(f'{VOL28}/ref_diff/hg38/gencode_v45/platinumTRs-v1.0.trgt.cds_utr_InSD.bed',
                             sep='\t', usecols=[3, 4, 5, 6], names=['context', 'gene', 'trid', 'motif'])

    sd_utr_cds.set_index('trid', inplace=True)

    # utr_cds_genes = pd.concat([utr_sites, cds_sites])
    # utr_cds_genes = pd.read_csv(f'{VOL28}/ref_diff/hg38/gencode_v45/refined_cds_utr/platinumTRs-v1.0.trgt.cds_utr.bed',
    #                             sep='\t', names=['chrom', 'start', 'end', 'context', 'gene', 'trid', 'motif'])

    utr_cds_genes = pd.read_csv(f'{VOL28}/ref_diff/hg38/gencode_v45/refined_cds_utr/platinumTRs-v1.0.trgt.cds_utr_TRE.bed',
                                sep='\t', names=['chrom', 'start', 'end', 'context', 'gene', 'trid', 'motif'])

    utr_triple = pd.read_csv(f'{VOL28}/ref_diff/hg38/gencode_v45/refined_cds_utr/platinumTRs-v1.0.trgt.utr_triple.tsv',
                             sep='\t', header=0)
    cds_triple = pd.read_csv(f'{VOL28}/ref_diff/hg38/gencode_v45/refined_cds_utr/platinumTRs-v1.0.trgt.cds_triple.tsv',
                             sep='\t', header=0)

    hg38_bdry_cds = pd.read_csv(f'{VOL28}/ref_diff/hg38/gencode_v45/refined_cds_utr/platinumTRs-v1.0.trgt.bdry_cds.tsv', sep='\t', index_col=['trid'])

    all_triplet = pd.concat([utr_triple, cds_triple])

    utr_cds_grp = utr_cds_genes.groupby('trid').agg({
                    'context': lambda x: ','.join(x.unique()),
                    'gene': lambda x: ','.join(x.unique()),
                    'motif': lambda x: ','.join(x.unique()),
                }).reset_index()

    utr_cds_grp = utr_cds_grp.loc[~utr_cds_grp['gene'].str.contains('ENSG')]

    utr_cds_grp['final_context'] = utr_cds_grp.apply(lambda row: update_context(row, hg38_bdry_cds.index), axis=1)
    utr_cds_grp['SegDup'] = utr_cds_grp.apply(lambda row: 'InSD' if row['trid'] in sd_utr_cds.index else 'OutSD', axis=1)
    utr_cds_grp['TR_TAG'] = utr_cds_grp.apply(lambda row: update_tr_type(row, all_triplet['trid'].tolist()), axis=1)
    utr_cds_grp.loc[utr_cds_grp['gene'] == 'UBXN11', 'TR_TAG'] = 'VNTR'
    # utr_cds_grp.set_index('trid', inplace=True)

    utr_cds_genes['tr_tag'] = utr_cds_genes.apply(lambda row: update_tr_type(row, all_triplet['trid'].tolist()), axis=1)

    return utr_cds_genes, utr_cds_grp, all_triplet


def length_variation_exfig2a(stats_type):
    unit = 'bp' if stats_type == 'length' else 'unit'

    pathogenic, cds_sites, utr_sites, cv_hue, cv_color = load_variation_db(stats_type)
    pathogenic = pathogenic.loc[pathogenic['unique_alleles'] >= 2]
    utr_sites = utr_sites.loc[(utr_sites['context'] == '5UTR') & (utr_sites['unique_alleles'] > 2)]
    cds_sites = cds_sites.loc[cds_sites['unique_alleles'] > 2]

    print('====== raw sites count ======')
    print('\tCDS: ', len(cds_sites))
    print('\t5UTR: ', len(utr_sites))

    print('====== after filtering sites count ======')
    print('\tCDS: ', len(cds_sites))
    print('\t5UTR: ', len(utr_sites))

    fig, axes = plt.subplots(3, 3, sharey='row', sharex=True, figsize=(14, 10))

    for row_i, ele in enumerate(['STR-Triplet', 'STR', 'VNTR']):
        text1, text2, text3 = [], [], []

        tmp_cds = cds_sites.loc[cds_sites['TR_TAG']==ele]
        tmp_utr = utr_sites.loc[utr_sites['TR_TAG']==ele]
        tmp_patho = pathogenic.loc[pathogenic['TR_TAG']==ele]

        print(f'======== {ele} ===========')
        print(f'\t#CDS: ', len(tmp_cds))
        print(f'\t#UTR: ', len(tmp_utr))
        print(f'\t#Pathogenic: ', len(tmp_patho))

        sns.scatterplot(data=tmp_cds, x='unique_alleles',
                        # color='#7fbc41',
                        rasterized=True, hue='color',
                        hue_order=cv_color, palette=cv_color, edgecolor="grey",
                        y='allele_range', ax=axes[row_i][0])

        sns.scatterplot(data=tmp_utr, x='unique_alleles',rasterized=True,
                        # hue='context', hue_order=['five_prime_UTR', 'three_prime_UTR'], palette=['#35978f', '#bf812d'],
                        hue='color', hue_order=cv_color, palette=cv_color, edgecolor="grey",
                        y='allele_range', ax=axes[row_i][1])

        sns.scatterplot(data=tmp_patho, x='unique_alleles',rasterized=True,
                        hue='color', hue_order=cv_color, palette=cv_color, edgecolor="grey",
                        y='allele_range', ax=axes[row_i][2])

        for i in range(3):
            ax = axes[row_i][i]
            if row_i == 0:
                ax.set_ylabel(f'Triplet length difference ({unit})', fontsize=13)
            elif row_i == 1:
                ax.set_ylabel(f'STR length difference ({unit})', fontsize=13)
            else:
                ax.set_ylabel(f'VNTR length difference ({unit})', fontsize=13)

            ax.set_yticks([2, 10, 100], labels=[2, 10, 100])
            ax.set_yscale('log')
            ax.set_xscale('log')
            # ax.set_xlabel('Observed heterozygousity', fontsize=13)
            ax.set_xlabel('#Unique alleles', fontsize=13)

            # if i != 2:
            #     ax.axhline(y=8, ls='--', lw=1, color='grey')
            #     ax.axvline(x=8, ls='--', lw=1, color='grey')

            ax.spines[['top', 'right']].set_visible(False)


        ## Add text labels
        for idx, row in cds_sites.loc[cds_sites['TR_TAG']==ele].iterrows():
            # if ele == 'VNTR':
            #     if row['allele_range'] >= 8 and row['unique_alleles'] >= 10:
            #         text1.append(axes[row_i][0].text(row['unique_alleles'], row['allele_range'], row['gene'],
            #                               horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))
            # else:
            #     if row['allele_range'] >= 50 and row['unique_alleles'] >= 5:
            #         text1.append(axes[row_i][0].text(row['unique_alleles'], row['allele_range'], row['gene'],
            #                               horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))
            if row['color'] in ['#fe9929', '#cc4c02', '#990000']:
                text1.append(axes[row_i][0].text(row['unique_alleles'], row['allele_range'], row['gene'],
                                       horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))


        for idx, row in utr_sites.loc[utr_sites['TR_TAG']==ele].iterrows():
            # if ele == 'VNTR':
            #     if row['allele_range'] >= 8 and row['unique_alleles'] >= 10:
            #         text2.append(axes[row_i][1].text(row['unique_alleles'], row['allele_range'], row['gene'],
            #                               horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))
            # else:
            #     if stats_type == 'repnum' and row['allele_range'] >= 50 and row['unique_alleles'] >= 5:
            #         text2.append(axes[row_i][1].text(row['unique_alleles'], row['allele_range'], row['gene'],
            #                               horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))
            #     if stats_type == 'length' and row['allele_range'] >= 100 and row['unique_alleles'] >= 5:
            #         text2.append(axes[row_i][1].text(row['unique_alleles'], row['allele_range'], row['gene'],
            #                               horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))
            if row['color'] in ['#fe9929', '#cc4c02', '#990000']:
                text2.append(axes[row_i][1].text(row['unique_alleles'], row['allele_range'], row['gene'],
                                       horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))


        for idx, row in pathogenic.loc[pathogenic['TR_TAG']==ele].iterrows():

            if row['unique_alleles'] >= 2:
                text3.append(axes[row_i][2].text(row['unique_alleles'], row['allele_range'], row['locus'].split('_')[0],
                                          horizontalalignment='center', verticalalignment='bottom', size='small', color='k'))

        adjust_text(text1, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=axes[row_i][0])
        adjust_text(text2, expand_points=(1.2, 1.2), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=axes[row_i][1])
        adjust_text(text3, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=axes[row_i][2])

    fig.tight_layout()
    fig.savefig(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/Figures/all_CdsUtr_{stats_type}_variation.svg')
    plt.show()



def main():
    str_length_cv_fig3b()
    all_length_cv()
    length_variation_exfig2a('length')

if __name__ == '__main__':
    main()
