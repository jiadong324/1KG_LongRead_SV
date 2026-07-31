#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 7/16/26
'''
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib as mpl
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none',
"font.family": "sans-serif",
"font.sans-serif": "Arial"
}
mpl.rcParams.update(new_rc_params)

VOL28 = './'

def summarize_positive(workdir):
    sv_af = pd.read_csv(f'{workdir}/plcai_pop_fill_tags.tsv', sep='\t', index_col=['ID'])
    sv_tbl = pd.read_csv(f'{workdir}/FST_PCLAI/stratified_svs/sv_sites_func.tsv', sep='\t', index_col=['SVID'])
    slt_res = pd.read_csv(f'{workdir}/FST_PCLAI/stratified_svs/DongAhn_Selection/positive_region.tsv', sep='\t', header=0)

    merged_annt = []
    for idx, row in slt_res.iterrows():
        td_bot = row['TD_bottom1_of_interest']
        ihs_top = row['iHS_top1_of_interest']
        for svid in row['variant_id'].split(','):
            strati_pop = sv_tbl.at[svid, 'Stratified_POP']
            svaf_info = sv_af.loc[svid][['AF_C1', 'AF_C2', 'AF_C3', 'AF_C4']].tolist()
            # if svid == 'chr18-40874829-DEL-6857':
            #     print('sss')
            af_diff = 0 if strati_pop == 'C1' else float(sv_af.at[svid, f'AF_{strati_pop}']) - float(sv_af.at[svid, 'AF_C1'])
            svinfo = [svid, td_bot, ihs_top] + sv_tbl.loc[svid][['AF', 'GENE', 'Stratified_POP', 'CONTEXT', 'Biallelic']].tolist() + svaf_info + [af_diff]
            merged_annt.append(svinfo)

    annt_df = pd.DataFrame(merged_annt, columns=['SVID', 'TD_bottom1_of_interest', 'iHS_top1_of_interest', 'AF', 'GENE', 'Stratified_POP', 'CONTEXT', 'Biallelic', 'AF_C1', 'AF_C2', 'AF_C3', 'AF_C4', 'AF_DIFF'])

    annt_df.to_csv(f'{workdir}/FST_PCLAI/stratified_svs/DongAhn_Selection/positive_svs_info.tsv', sep='\t', header=True, index=False)

    # print(annt_df.head())


def mahattan_selection_af_diff_fig4e(workdir):
    tmp_tbl = pd.read_csv(f'{workdir}/FST_PCLAI/stratified_svs/DongAhn_Selection/positive_svs_info.tsv', sep='\t', header=0)
    # selt_tbl = selt_tbl.loc[selt_tbl['Biallelic']=='Yes']

    lct_info = pd.DataFrame([{'SVID': 'chr2-136274028-DEL-3951', 'GENE': 'LCT', 'CONTEXT': 'intron', 'Stratified_POP': 'C3',
                              'Biallelic': 'Yes', 'AF_C1': 0.190217, 'AF_C3': 0.548433, 'AF_DIFF': 0.358216}])

    selt_tbl = pd.concat([tmp_tbl, lct_info])

    print(selt_tbl['SVID'].nunique())

    ## biallelic genic SVs
    bia_genic_tbl = selt_tbl.loc[(selt_tbl['Biallelic']=='Yes')&(selt_tbl['GENE']!='.')]
    bia_genic_tbl.set_index('SVID', inplace=True)

    print(len(bia_genic_tbl.index))

    selt_tbl['CHROM'] = selt_tbl.apply(lambda row: row['SVID'].split('-')[0], axis=1)
    order_map = {val: i for i, val in enumerate(AUTOSOMES)}

    selt_tbl.sort_values(by='CHROM', key=lambda x: x.map(order_map), inplace=True)

    selt_tbl['ind'] = range(len(selt_tbl))
    selt_tbl['ABS_AF_DIFF'] = selt_tbl.apply(lambda row: abs(row['AF_DIFF']), axis=1)
    df_grouped = selt_tbl.groupby('CHROM')
    x_labels = []
    x_labels_pos = []
    x_vline_pos = set()

    legends = [Line2D([0], [0], ls='none', marker='v', label='Decrease', markeredgecolor="k", markerfacecolor='none'),
               Line2D([0], [0], ls='none', marker='^', label='Increase', markeredgecolor="k", markerfacecolor='none'),]

    for pop in ['C2', 'C3', 'C4']:
        legends.append(Patch(label=pop, color=POPCOLOR[pop]))

    # fig, axes = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(24, 10))
    fig, ax = plt.subplots(figsize=(24, 8))

    texts = []
    for num, (name, group) in enumerate(df_grouped):

        for i, pop in enumerate(['C2', 'C3', 'C4']):
            # group = selt_tbl.loc[selt_tbl['Stratified_POP']==pop]
            minus = group.loc[(group['AF_DIFF'] < 0)&(group['Stratified_POP']==pop)]
            plus = group.loc[(group['AF_DIFF'] > 0)&(group['Stratified_POP']==pop)]
            # sns.scatterplot(data=minus, x='ind', y='ABS_AF_DIFF', marker='v', color=POPCOLOR[pop], rasterized=True, s=80, ax=axes[i])
            # sns.scatterplot(data=plus, x='ind', y='ABS_AF_DIFF', marker='^', color=POPCOLOR[pop], rasterized=True, s=80, ax=axes[i])
            sns.scatterplot(data=minus.loc[minus['SVID'].isin(bia_genic_tbl.index)], x='ind', y='ABS_AF_DIFF', marker='v', color=POPCOLOR[pop], rasterized=True, s=300, ax=ax)
            sns.scatterplot(data=minus.loc[~minus['SVID'].isin(bia_genic_tbl.index)], x='ind', y='ABS_AF_DIFF',
                            marker='v', color='#bababa', rasterized=True, s=100, ax=ax)

            sns.scatterplot(data=plus.loc[plus['SVID'].isin(bia_genic_tbl.index)], x='ind', y='ABS_AF_DIFF', marker='^', color=POPCOLOR[pop], rasterized=True, s=300, ax=ax)
            sns.scatterplot(data=plus.loc[~plus['SVID'].isin(bia_genic_tbl.index)], x='ind', y='ABS_AF_DIFF',
                            marker='^', color='#bababa', rasterized=True, s=100, ax=ax)

        x_labels.append(name)
        x_vline_pos.add(group['ind'].iloc[0])
        # x_vline_pos.add(group['ind'].iloc[-1])
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))


    for idx, row in selt_tbl.iterrows():
        if row['SVID'] in bia_genic_tbl.index and row['ABS_AF_DIFF'] > 0.01:
            texts.append(ax.text(row['ind'], row['ABS_AF_DIFF'], row['GENE'], horizontalalignment='center',
                                      verticalalignment='bottom', size='small', color='black'))

    adjust_text(texts, expand_points=(0.1, 0.1), arrowprops=dict(arrowstyle='-', color='grey', alpha=0.5), ax=ax)

    for val in list(x_vline_pos):
        ax.axvline(x=val, ls='--', color='grey')
    ax.spines[['right', 'top']].set_visible(False)
    plt.margins(x=0.01)
    plt.margins(y=0.05)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1], labels=[0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=12)
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, fontsize=14, rotation=90)
    ax.set_xlabel('')
    ax.legend(handles=legends)
    fig.tight_layout()
    fig.savefig(f'{workdir}/FST_PCLAI/stratified_svs/DongAhn_Selection/positive_svs_info_addLCT.svg')
    plt.show()

def main():
    summarize_positive(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete')
    mahattan_selection_af_diff_fig4e(f'{VOL28}/SVREF/all_cohorts/CHM13_0711/pct09_dist500/complete')

if __name__ == '__main__':
    main()
