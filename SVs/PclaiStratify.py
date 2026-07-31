#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 6/22/26
'''

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pickle
import gzip, json
import matplotlib as mpl
from adjustText import adjust_text


new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none',
"font.family": "sans-serif",
"font.sans-serif": "Arial"
}
mpl.rcParams.update(new_rc_params)

VOL28 = './'
POPORDER = ['C1', 'C2', 'C3', 'C4']

def parse_vcf_info_column(info_str):
    info_tokens = info_str.split(";")
    info_dict = {}

    for token in info_tokens:
        if "=" not in token:
            continue
        info_dict[token.split('=')[0]] = token.split('=')[1]

    return info_dict

def group_sv_fst(workdir, pop_list):

    fst_list = []
    ac_value_bins = list(range(0, 100, 1))
    pop_fst_bin_dict = {}

    total_count = 0
    for pop in pop_list:
        bin_count = {i: 0 for i in ac_value_bins}
        fst_tbl = pd.read_csv(f'{workdir}/FST_PCLAI/{pop}/{pop}_list.weir.fst', sep='\t', header=[0])
        fst_tbl.dropna(inplace=True)
        # fst_tbl = fst_tbl.loc[fst_tbl['WEIR_AND_COCKERHAM_FST'] > 0]
        total_count += len(fst_tbl.loc[fst_tbl['WEIR_AND_COCKERHAM_FST'] >= 0.15])

        for idx, row in fst_tbl.iterrows():
            fst_val = abs(row['WEIR_AND_COCKERHAM_FST'])
            bin_idx = fst_val // 0.01
            bin_count[bin_idx] += 1
        pop_fst_bin_dict[pop] = bin_count

        fst_tbl['Label'] = [pop for _ in range(len(fst_tbl))]
        fst_tbl['ID'] = fst_tbl.apply(lambda row: "{0}-{1}".format(row['CHROM'], row['POS']), axis=1)
        fst_list.append(fst_tbl)


    with open(f'{workdir}/FST_PCLAI/fst_freq.pickle', 'wb') as f:
        pickle.dump(pop_fst_bin_dict, f)
    print(f'Total stratified: {total_count}')

    fst_all = pd.concat(fst_list)
    fst_all.set_index('ID', inplace=True)

    grouped_svs = []
    counter = 0
    for line in gzip.open(f'{workdir}/plcai_pop_fill_tags.vcf.gz', 'rt'):
        if line.startswith('#'):
            continue
        entries = line.strip().split('\t')
        svid = f'{entries[0]}-{entries[1]}'
        info_dict = parse_vcf_info_column(entries[7])
        svtype, svlen = info_dict['SVTYPE'], int(info_dict['SVLEN'])
        end = int(entries[1]) + 1 if svtype == 'INS' else int(entries[1]) + abs(svlen)
        if svid in fst_all.index:
            fst_vals = fst_all.loc[svid]
            pop_fst = ['.', '.', '.', '.']
            pop_ac = []
            for i, pop in enumerate(POPORDER[::-1]):
                pop_ac.append(info_dict[f'AC_{pop}'])
                if type(fst_vals) != pd.Series:
                    val = fst_vals.loc[fst_vals['Label']==pop]['WEIR_AND_COCKERHAM_FST'].values
                    if len(val) > 0:
                        pop_fst[i] = max(val)
                else:
                    this_pop = POPORDER[::-1].index(fst_vals['Label'])
                    pop_fst[this_pop] = fst_vals['WEIR_AND_COCKERHAM_FST']

            grouped_svs.append([entries[0], entries[1], end, entries[2], svtype, svlen] + pop_fst + pop_ac)

        if counter % 10000 == 0:
            print('Processed SVs: ', counter)

        counter += 1

    headers = ['CHROM', 'POS', 'END', 'SVID', 'SVTYPE', 'SVLEN']
    for pop in POPORDER[::-1]:
        headers.append(f'{pop}_fst')
    for pop in POPORDER[::-1]:
        headers.append(f'{pop}_ac')
    df_grouped = pd.DataFrame(grouped_svs, columns=headers)
    df_grouped.to_csv(f'{workdir}/FST_PCLAI/sv_fst_ac.tsv', sep='\t', header=True, index=False)


def annot_chm13_stratified_svs():
    biallelics = [line.strip() for line in open(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/Annot/allele_class/biallelic_insdel.txt')]

    sd_svs = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/Annot/gt_af_sd.tsv.gz', sep='\t', usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'svid'])
    sd_svs.set_index('svid', inplace=True)
    tr_svs = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/InTR/svs_popaf.tsv.gz',sep='\t',
                         usecols=[0, 1, 2, 3, 25], names=['chrom', 'start', 'end', 'svid', 'motif'])
    tr_svs.set_index('svid', inplace=True)

    gene_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/Annot/gt_af_gene.tsv', sep='\t', usecols=[0,1,2,3,9,10], names=['chrom', 'start', 'end', 'svid', 'gene', 'context'])
    gene_tbl.set_index('svid', inplace=True)

    reg_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/Annot/gt_af_reg.tsv',sep='\t', usecols=[0, 1, 2, 3, 8],names=['chrom', 'start', 'end', 'svid', 'element'])
    reg_tbl.set_index('svid', inplace=True)

    fst_svs = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/FST_PCLAI/sv_fst_ac.tsv', sep='\t', header=0)
    sv_af_db = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/fill_tags.tsv.gz', sep='\t', index_col=['ID'])

    ## Chenlu introgression results
    # cl_ai = pd.read_csv('/Volumes/eichler-vol28/home/cdi6/nobackups/SV/introgression/sv_fst_ac.AI_annotated.tsv', sep='\t', index_col=['SVID'])
    # cl_ai.fillna('.', inplace=True)

    annot_out = []

    for idx, row in fst_svs.iterrows():
        svid = row['SVID']
        sv_af = sv_af_db.at[svid, 'AF_Human']

        sd_tag = 'InSD' if svid in sd_svs.index else 'OutSD'
        tr_tag = tr_svs.at[svid, 'motif'] if svid in tr_svs.index else '.'

        bi_tag = 'Yes' if svid in biallelics else 'No'

        # deni = cl_ai.at[svid, 'Denisovan_AI_MaLAdapt'] if svid in cl_ai.index else '.'
        # nean = cl_ai.at[svid, 'Neanderthal_AI_MaLAdapt'] if svid in cl_ai.index else '.'

        deni = '.'
        nean = '.'

        fsts = []
        for pop in POPORDER:
            if row[f'{pop}_fst'] != '.' and float(row[f'{pop}_fst']) >= 0.15:
                fsts.append([pop, float(row[f'{pop}_fst'])])

        if len(fsts) > 0:
            sorted_fst = sorted(fsts, key=lambda x:x[0], reverse=True)
            stratified_pop, fst = sorted_fst[0]

            gene, key_struct, struct = '.', '.', '.'
            if svid in gene_tbl.index:
                gene_info = gene_tbl.loc[svid]
                if isinstance(gene_info, pd.DataFrame):
                    struct_list = gene_info['context'].unique().tolist()
                    gene = ';'.join(gene_info['gene'].unique().tolist())
                    struct = ';'.join(gene_info['context'].unique().tolist())
                    if 'UTR' in struct_list:
                        key_struct = 'UTR'
                    if 'CDS' in struct_list:
                        key_struct = 'CDS'
                else:
                    gene = gene_info['gene']
                    struct = gene_info['context']
                    key_struct = struct

            reg = '.'
            if svid in reg_tbl.index:
                reg_info = reg_tbl.loc[svid]
                if isinstance(reg_info, pd.DataFrame):
                    reg = ';'.join(reg_info['element'].unique().tolist())
                else:
                    reg = reg_info['element']

            if reg != '.' and key_struct == '.':
                key_struct = 'REG'
            location = 'Genic' if gene != '.' else 'Intergenic'
            annot_out.append(row.tolist() + [bi_tag, sv_af, gene, key_struct, location, struct, reg, tr_tag, sd_tag, stratified_pop, fst, deni, nean])

    annot_df = pd.DataFrame(annot_out, columns=list(fst_svs.columns) + ['Biallelic', 'AF', 'GENE', 'CONTEXT', 'LOCATION',  'INFO', 'REG', 'TR', 'SD', 'Stratified_POP', 'FST', 'AI_Denisovan', 'AI_Neanderthal'])
    annot_df.to_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/FST_PCLAI/stratified_svs/sv_sites_func.tsv', sep='\t', header=True, index=False)

def plot_binned_fst_val_fig4c(workdir, pop_list):
    reported = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/peterE_reported_stratify_gene.txt', sep='\t', index_col=['GENE'])

    fst_annot = pd.read_csv(f'{workdir}/stratified_svs/sv_sites_func.tsv', sep='\t', header=0)
    stratified = fst_annot.loc[fst_annot['Stratified_POP']!='.']

    fst_grps = stratified.groupby('Stratified_POP').size().reset_index(name='Count')
    print('All stratified SVs')
    print(fst_grps.head())

    biall = stratified.loc[stratified['Biallelic']=='Yes']
    biall_grps = biall.groupby('Stratified_POP').size().reset_index(name='Count')
    print('Biallelic stratified SVs')
    print(biall_grps.head())

    biall_genic = stratified.loc[(stratified['Biallelic']=='Yes')&(stratified['GENE']!='.')]
    biall_genic['Novel'] = biall_genic.apply(lambda row: 'No' if row['GENE'] in reported.index else 'Yes', axis=1)
    biall_genic_grps = biall_genic.groupby('Stratified_POP').size().reset_index(name='Count')
    print(biall_genic_grps.head())
    print('Stratified biallelic genic SVs: ', biall_genic_grps['Count'].sum())

    biall_genic.to_csv(f'{workdir}/stratified_svs/biallelic_genic_svs.tsv', sep='\t', header=True, index=False)

    biall_cds_utr = biall.loc[(biall['CONTEXT'].str.contains('UTR'))|(biall['CONTEXT']=='CDS')]
    biall_cds_utr_grps = biall_cds_utr.groupby('Stratified_POP').size().reset_index(name='Count')
    print(biall_cds_utr_grps.head())
    print('Stratified biallelic CDS/UTR SVs: ', biall_cds_utr_grps['Count'].sum())
    biall_cds_utr.to_csv(f'{workdir}/stratified_svs/biallelic_cds_utr_svs.tsv', sep='\t', header=True, index=False)

    biall_cds_utr_reg = biall.loc[(biall['CONTEXT'].str.contains('UTR')) | (biall['CONTEXT'] == 'CDS')| (biall['CONTEXT'] == 'REG')]
    biall_cds_utr_reg_grps = biall_cds_utr_reg.groupby('Stratified_POP').size().reset_index(name='Count')
    print(biall_cds_utr_reg_grps.head())
    print('Stratified biallelic CDS/UT/REG SVs: ', biall_cds_utr_reg_grps['Count'].sum())
    biall_cds_utr.to_csv(f'{workdir}/stratified_svs/biallelic_cds_utr_reg_svs.tsv', sep='\t', header=True, index=False)


    with open(f'{workdir}/fst_freq.pickle', 'rb') as f:
        pop_fst_bin_dict = pickle.load(f)

    fig, ax = plt.subplots(figsize=(7, 4))
    # xticks = np.arange(len(ac_value_bins))
    # ax.set_title(f'{func.upper()} SVs (N={len(func_fst)})', fontsize=14)
    val_list = []
    for pop, val_dict in pop_fst_bin_dict.items():
        for fst, val in val_dict.items():
            val_list.append([fst, val, pop])
        # ax.plot( xticks, list([ele + 1 for ele in val_dict.values()]),color=POPCOLOR[pop], label=pop)

    df_val = pd.DataFrame(val_list, columns=['fst', 'val', 'pop'])
    sns.lineplot(data=df_val, x='fst', y='val', hue='pop', hue_order=pop_list,
                 palette=[POPCOLOR[ele] for ele in pop_list], ax=ax, errorbar=None)

    ax.set_xticks([0, 20, 40, 60, 80, 100], labels=[0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=12)
    ax.axvline(15, ls='--', color='grey')
    # ax.text(17, 500, f'N={total_num} (Fst>=0.15)')
    ax.set_yscale('log')
    ax.set_ylabel('NO. of SVs', fontsize=13)
    ax.set_xlabel('Fst', fontsize=13)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(f'{workdir}/sv_fst.png', dpi=300)
    fig.savefig(f'{workdir}/sv_fst.svg')
    # plt.show()


def fst_svlen_af_scatter_fig4d(workdir):
    reported = pd.read_csv(f'{workdir}/peterE_reported_stratify_gene.txt', sep='\t', header=0)

    sv_af = pd.read_csv(f'{workdir}/pclai_grp_fill_tags.tsv', sep='\t', header=0)
    sv_tbl = pd.read_csv(f'{workdir}/FST_PCLAI/stratified_svs/sv_sites_func.tsv', sep='\t', header=0)

    # biall = sv_tbl.loc[(sv_tbl['Biallelic']=='Yes')&(sv_tbl['CONTEXT'].isin(['REG', '5UTR', '3UTR', 'CDS']))]
    # biall = sv_tbl.loc[(sv_tbl['Biallelic']=='Yes')&(sv_tbl['CONTEXT'].isin(['REG', '5UTR', '3UTR', 'CDS']))]
    biall = sv_tbl.loc[(sv_tbl['Biallelic'] == 'Yes') & (sv_tbl['LOCATION']=='Genic')]
    biall['ABS_SVLEN'] = biall.apply(lambda row: abs(row['SVLEN']), axis=1)
    print(biall['CONTEXT'].unique().tolist())

    fig, ax = plt.subplots(figsize=(8, 7))

    # sns.scatterplot(data=biall, x='ABS_SVLEN', y='FST', style='SVTYPE',
    #                 style_order=['INS', 'DEL'], hue='Stratified_POP', hue_order=['C1', 'C2', 'C3', 'C4'],
    #                 palette=[POPCOLOR['C1'], POPCOLOR['C2'], POPCOLOR['C3'], POPCOLOR['C4']], markers=True, rasterized=True, ax=ax)

    sns.scatterplot(data=biall, x='ABS_SVLEN', y='FST', style='SVTYPE',
                    style_order=['INS', 'DEL'], hue='CONTEXT', hue_order=['.', 'intron'], markers=['o', 'D'],
                    s=50, palette=['#bababa', '#bababa'], ec='face',
                    rasterized=True, ax=ax)
    sns.scatterplot(data=biall, x='ABS_SVLEN', y='FST', style='SVTYPE', ec='face',
                    style_order=['INS', 'DEL'], hue='CONTEXT', hue_order=['REG', 'CDS', '3UTR', '5UTR'], markers=['o', 'D'], s=150,
                    rasterized=True, ax=ax)


    ax.set_xscale('log')
    text = []
    novel = 0
    for idx, row in biall.iterrows():
        if row['GENE']!='.' and row['GENE'] not in reported['GENE'].tolist():
            novel += 1
        if row['CONTEXT'] in ['CDS', '3UTR', '5UTR', 'REG']:
            text.append(ax.text(float(row['ABS_SVLEN']), float(row['FST']), row['GENE'].split(';')[0], horizontalalignment='center',
                                         verticalalignment='bottom', fontsize=12, color='k'))
    adjust_text(text, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=ax)
    fig.tight_layout()
    fig.savefig(f'{workdir}/FST_PCLAI/stratified_svs/fst_svlen_plot.svg')

    print('#Novel SVs compared to Ebert: ', novel)

    biall_af = pd.merge(biall, sv_af[['SVID', 'AF_C1', 'AF_NonC1']], on='SVID', how='left')
    biall_af['MARKER_SIZE'] = biall_af.apply(lambda row: row['ABS_SVLEN'] / 100, axis=1)


    svs_num = biall_af['SVID'].nunique()
    fig, ax = plt.subplots(figsize=(8, 7))
    sns.scatterplot(data=biall_af, x='AF_C1', y='AF_NonC1', hue='Stratified_POP', style='SVTYPE',
                    hue_order=['C1', 'C2', 'C3', 'C4'], palette=[POPCOLOR['C1'], POPCOLOR['C2'], POPCOLOR['C3'], POPCOLOR['C4']],
                    style_order=['INS', 'DEL'], markers=True, s=100, ax=ax)


    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1], labels=[0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=12)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1], labels=[0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=12)
    ax.set_title(f'#Stratified biallelic CDS/UTR/REG SVs: {svs_num}', fontsize=14)
    ax.set_xlabel('C1 allele frequency', fontsize=13)
    ax.set_ylabel('NonC1 allele frequency', fontsize=13)
    fig.tight_layout()
    # fig.savefig(f'{workdir}/FST_PCLAI/stratified_svs/fst_freq_plot.svg')
    plt.show()

def main():
    group_sv_fst(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete', ['C1', 'C2', 'C3', 'C4'])

    ## Full set annotation of stratified SVs
    annot_chm13_stratified_svs()

    plot_binned_fst_val_fig4c(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete/FST_PCLAI', ['C1', 'C2', 'C3', 'C4'])

    fst_svlen_af_scatter_fig4d(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete')

if __name__ == '__main__':
    main()
