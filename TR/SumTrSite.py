#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 7/30/26
'''
import json
import os.path
from scipy.stats import kurtosis
import numpy as np
import pandas as pd
from scipy.stats import kstest

VOL28 = './'


def group_trgt_outputs(ref, trgt_type, trtype):

    ## workdir for all CDS/UTR TR sites
    workdir = f'/net/eichler/vol28/projects/medical_reference/nobackups/SVREF/all_cohorts/{ref}/pct09_dist500/matt_trgt'

    ## workdir for TR sites interrupt by SVs
    if trtype == 'svtr':
        workdir = f'/net/eichler/vol28/projects/medical_reference/nobackups/SVREF/all_cohorts/{ref}/pct09_dist500/svtr_trgt'

    if not os.path.exists(f'{workdir}/tr_{trgt_type}_stats'):
        os.mkdir(f'{workdir}/tr_{trgt_type}_stats')

    trgt_tbl = pd.read_csv(f'{workdir}/hgsvc_hprc_trgt.tab', sep='\t', index_col=['SAMPLE'])
    pop_dict = json.load(open('/net/eichler/vol28/projects/medical_reference/nobackups/Samples/igsr_samples_pop.json'))

    all_trs = []
    for sample in trgt_tbl.index:
        tr_info_file = f'{workdir}/tr_{trgt_type}_info/{sample}_trid_info.tsv'

        if not os.path.exists(tr_info_file) or os.path.getsize(tr_info_file) == 0:
            print(sample)
            continue
        allele_tbl = pd.read_csv(tr_info_file, sep='\t', header=[0])

        allele_tbl['POP'] = [pop_dict[sample] for _ in range(len(allele_tbl))]
        all_trs.append(allele_tbl)

    all_df = pd.concat(all_trs)
    grouped_trs = all_df.groupby('TRID')
    allele_info = []
    allele_gt = []
    missed_allele = {}
    allele_len_perct = []
    allele_num_perct = []
    sample_set = set()
    for trid, group in grouped_trs:
        tr_lens = []
        tr_nums = []
        tr_len_bypop = {ele: [] for ele in ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']}
        tr_num_bypop = {ele: [] for ele in ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']}
        motifs = group['MOTIFS'].values[0]

        this_tr_gt = [trid]
        for idx, row in group.iterrows():
            if row['AL_LENGTH'] == '.':
                if trid not in missed_allele:
                    missed_allele[trid] = 0
                missed_allele[trid] += 1
                continue

            al_sd1, al_sd2 = row['AL_SD'].split(',')
            al_len1, al_len2 = row['AL_LENGTH'].split(',')
            al_mc1, al_mc2 = row['AL_MC'].split(',')
            al_num1, al_num2 = sum([int(ele) for ele in al_mc1.split('_')]), sum([int(ele) for ele in al_mc2.split('_')])

            if int(al_sd1) >= 5:
                tr_lens.append(int(al_len1))
                tr_nums.append(al_num1)
                allele_info.append([trid, row['MOTIFS'], row['SAMPLE'], row['POP'], int(al_len1), al_mc1])
                tr_len_bypop[row['POP']].append(int(al_len1))
                tr_num_bypop[row['POP']].append(al_num1)

            if int(al_sd2) >= 5:
                tr_lens.append(int(al_len2))
                tr_nums.append(al_num2)
                allele_info.append([trid, row['MOTIFS'], row['SAMPLE'], row['POP'], int(al_len2), al_mc2])
                tr_len_bypop[row['POP']].append(int(al_len2))
                tr_num_bypop[row['POP']].append(al_num2)

            # tr_lens.extend([int(al_len1), int(al_len2)])
            # tr_nums.extend([sum([int(ele) for ele in al_mc1.split('_')]), sum([int(ele) for ele in al_mc2.split('_')])])

            # allele_info.append([trid, row['MOTIFS'], row['SAMPLE'], row['POP'], int(al_len1), al_mc1])
            # allele_info.append([trid, row['MOTIFS'], row['SAMPLE'], row['POP'], int(al_len2), al_mc2])
            # tr_len_bypop[row['POP']].extend([int(al_len1), int(al_len2)])
            # tr_num_bypop[row['POP']].extend([sum([int(ele) for ele in al_mc1.split('_')]), sum([int(ele) for ele in al_mc2.split('_')])])

            this_tr_gt.append(row['GT'])
            sample_set.add(row['SAMPLE'])

        allele_gt.append(this_tr_gt)

        if len(tr_lens) > 0:
            al_len_std = np.std(tr_lens)
            al_len_avg = np.mean(tr_lens)
            pct99 = np.percentile(tr_lens, 99)
            pct95 = np.percentile(tr_lens, 95)
            pct90 = np.percentile(tr_lens, 90)
            pct50 = np.percentile(tr_lens, 50)
            pct10 = np.percentile(tr_lens, 10)
            pct5 = np.percentile(tr_lens, 5)
            max_len = max(tr_lens)
            min_len = min(tr_lens)

            unique_len = len(set(tr_lens))
            kval = kurtosis(tr_lens)
            div = diversity(tr_lens)
            div_num = diversity(tr_nums)

            pop_al_pcrt = [trid, motifs, al_len_std, al_len_avg, int(pct99), int(pct50), int(pct10), int(pct5), max_len, min_len, unique_len, kval, div, div_num]

            afr, non_afr = tr_len_bypop['AFR'], tr_len_bypop['EAS'] + tr_len_bypop['SAS'] + tr_len_bypop['EUR'] + tr_len_bypop['AMR']
            eas, non_eas = tr_len_bypop['EAS'], tr_len_bypop['AFR'] + tr_len_bypop['SAS'] + tr_len_bypop['EUR'] + tr_len_bypop['AMR']
            amr, non_amr = tr_len_bypop['AMR'], tr_len_bypop['EAS'] + tr_len_bypop['SAS'] + tr_len_bypop['EUR'] + tr_len_bypop['AFR']
            eur, non_eur = tr_len_bypop['EUR'], tr_len_bypop['EAS'] + tr_len_bypop['SAS'] + tr_len_bypop['AFR'] + tr_len_bypop['AMR']
            sas, non_sas = tr_len_bypop['SAS'], tr_len_bypop['EAS'] + tr_len_bypop['AFR'] + tr_len_bypop['EUR'] + tr_len_bypop['AMR']

            for pop, tr_list in tr_len_bypop.items():
                this_pop_pct99 = np.percentile(tr_list, 99) if len(tr_list) > 0 else 0
                this_pop_pct50 = np.percentile(tr_list, 50) if len(tr_list) > 0 else 0
                pop_al_pcrt.extend([int(this_pop_pct99), int(this_pop_pct50)])

            for pop1, pop_rest in [[afr, non_afr], [eas, non_eas], [amr, non_amr], [eur, non_eur], [sas, non_sas]]:
                if len(pop1) > 0 and len(pop_rest) > 0:
                    statistics, pvalue = kstest(np.array(pop1), np.array(pop_rest))
                    pop_al_pcrt.append(pvalue)
                else:
                    pop_al_pcrt.append('.')

            allele_len_perct.append(pop_al_pcrt)

        if len(tr_nums) > 0:
            al_num_std = np.std(tr_nums)
            al_num_avg = np.mean(tr_nums)
            pct99 = np.percentile(tr_nums, 99)
            # pct95 = np.percentile(tr_nums, 95)
            # pct90 = np.percentile(tr_nums, 90)
            pct50 = np.percentile(tr_nums, 50)
            pct10 = np.percentile(tr_nums, 10)
            pct5 = np.percentile(tr_nums, 5)
            max_len = max(tr_nums)
            min_len = min(tr_nums)


            unique_len = len(set(tr_nums))
            pop_num_pcrt = [trid, motifs, al_num_std, al_num_avg, int(pct99), int(pct50), int(pct10), int(pct5), max_len, min_len, unique_len]
            for pop, tr_list in tr_num_bypop.items():
                this_pop_pct99 = np.percentile(tr_list, 99) if len(tr_list) > 0 else 0
                # this_pop_pct75 = np.percentile(tr_list, 75) if len(tr_list) > 0 else 0
                this_pop_pct50 = np.percentile(tr_list, 50) if len(tr_list) > 0 else 0
                pop_num_pcrt.extend([int(this_pop_pct99), int(this_pop_pct50)])

            afr, non_afr = tr_num_bypop['AFR'], tr_num_bypop['EAS'] + tr_num_bypop['SAS'] + tr_num_bypop['EUR'] + tr_num_bypop['AMR']
            eas, non_eas = tr_num_bypop['EAS'], tr_num_bypop['AFR'] + tr_num_bypop['SAS'] + tr_num_bypop['EUR'] + tr_num_bypop['AMR']
            amr, non_amr = tr_num_bypop['AMR'], tr_num_bypop['EAS'] + tr_num_bypop['SAS'] + tr_num_bypop['EUR'] + tr_num_bypop['AFR']
            eur, non_eur = tr_num_bypop['EUR'], tr_num_bypop['EAS'] + tr_num_bypop['SAS'] + tr_num_bypop['AFR'] + tr_num_bypop['AMR']
            sas, non_sas = tr_num_bypop['SAS'], tr_num_bypop['EAS'] + tr_num_bypop['AFR'] + tr_num_bypop['EUR'] + tr_num_bypop['AMR']

            for pop1, pop_rest in [[afr, non_afr], [eas, non_eas], [amr, non_amr], [eur, non_eur], [sas, non_sas]]:
                if len(pop1) > 0 and len(pop_rest) > 0:
                    statistics, pvalue = kstest(np.array(pop1), np.array(pop_rest))
                    pop_num_pcrt.append(pvalue)
                else:
                    pop_num_pcrt.append('.')

            allele_num_perct.append(pop_num_pcrt)


    fout = open(f'{workdir}/tr_{trgt_type}_stats/missed_trid_count.tsv', 'w')
    for trid, count in missed_allele.items():
        print(f'{trid}\t{count}', file=fout)

    # allele_gt_df = pd.DataFrame(allele_gt, columns=['trid'] + list(sample_set))
    # allele_gt_df.to_csv(f'{workdir}/allele_gt.tsv', sep='\t', header=True, index=False)

    allele_df = pd.DataFrame(allele_info, columns=['trid', 'motifs', 'sample', 'pop', 'length', 'mc'])
    allele_df.to_csv(f'{workdir}/tr_{trgt_type}_stats/allele_length.tsv', sep='\t', header=True, index=False)

    headers = ['trid', 'motifs', 'std', 'avg', 'pct99', 'pct50', 'pct10', 'pct5', 'max_len', 'min_len', 'unique_alleles']
    for pop in ['AFR', 'AMR', 'EUR', 'EAS', 'SAS']:
        for pect in ['pct99', 'pct50']:
            headers.append(f'{pop}_{pect}')
    headers.extend(['AFR_nonAFR', 'EAS_nonEAS', 'AMR_nonAMR', 'EUR_nonEUR', 'SAS_nonSAS'])

    allele_pcrt = pd.DataFrame(allele_len_perct, columns=headers + ['kurtosis', 'diversity', 'cn_diversity'])
    allele_pcrt.to_csv(f'{workdir}/tr_{trgt_type}_stats/allele_length_percentiles.tsv', sep='\t', header=True, index=False)

    allele_num = pd.DataFrame(allele_num_perct, columns=headers)
    allele_num.to_csv(f'{workdir}/tr_{trgt_type}_stats/allele_repnum_percentiles.tsv', sep='\t', header=True, index=False)

def diversity(alleles):
    alleles = np.array(alleles)
    n = len(alleles)

    unique, counts = np.unique(alleles, return_counts=True)
    freqs = counts / n

    # Diversity
    H = 1 - np.sum(freqs ** 2)
    return H

def classify_gt(s):
    return 'HETE' if s.max() - s.min() != 0 else 'HOMO'

def tr_heterozygosity(workdir):

    allele_tbl = pd.read_csv(f'{workdir}/allele_length.tsv', sep='\t', header=0)
    allele_tbl['cn'] = allele_tbl.apply(lambda row: sum([int(ele) for ele in row['mc'].split('_')]), axis=1)

    allele_sample_groups = allele_tbl.groupby(['trid', 'sample']).agg(CN_gt=('cn',classify_gt),
                                                                      LEN_gt=('length',classify_gt),
                                                                      CN=('cn',lambda x: '|'.join(x.astype(str))),
                                                                      LEN=('length', lambda x: '|'.join(x.astype(str))),
                                                                      ).reset_index()

    allele_gts = allele_sample_groups.groupby(['trid', 'CN_gt']).size().reset_index(name='Count')
    allele_gt_all = allele_sample_groups.groupby('trid').size().reset_index(name='All_count')

    merged_allele = pd.merge(allele_gts, allele_gt_all, on=['trid'], how='left')
    merged_allele['all_het'] = merged_allele.apply(lambda row: row['Count']/row['All_count'], axis=1)

    merged_allele.loc[merged_allele['CN_gt']=='HETE'].to_csv(f'{workdir}/allele_cn_hete.tsv', sep='\t', header=True, index=False)

    allele_gts_len = allele_sample_groups.groupby(['trid', 'LEN_gt']).size().reset_index(name='Count')
    allele_gt_all_len = allele_sample_groups.groupby('trid').size().reset_index(name='All_count')

    merged_allele_len = pd.merge(allele_gts_len, allele_gt_all_len, on=['trid'], how='left')
    merged_allele_len['all_het'] = merged_allele_len.apply(lambda row: row['Count'] / row['All_count'], axis=1)

    merged_allele_len.loc[merged_allele_len['LEN_gt'] == 'HETE'].to_csv(f'{workdir}/allele_len_hete.tsv', sep='\t', header=True, index=False)

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

def complete_sites_annot(stats_type):
    sd_tr = pd.read_csv(f'{VOL28}/ref_diff/hg38/hg38_sd_tr.bed', sep='\t',names=['chrom', 'start', 'end', 'motif_id', 'motif'])
    sd_tr.set_index('motif_id', inplace=True)

    svtr_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/InTR/tr_svs.stats.tsv', sep='\t', index_col=['motif_id'])

    ## UTR annotation mixed with CDS, then prioritize CDS

    utr_cds_genes, utr_cds_grp, all_triplet = load_utr_cds_genes()

    cds_final = utr_cds_grp.loc[utr_cds_grp['final_context'].str.contains('CDS')]
    utr_final = utr_cds_grp.loc[utr_cds_grp['final_context'].str.contains('UTR')]

    utr_cds_grp.set_index('trid', inplace=True)
    cds_allele_het = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_cds_stats/allele_hete.tsv', sep='\t', index_col=['trid'])

    utr_allele_het = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_utr_stats/allele_hete.tsv', sep='\t', index_col=['trid'])

    cds_sites = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_cds_stats/allele_{stats_type}_percentiles.tsv', sep='\t', header=0)
    cds_sites = cds_sites.loc[(~cds_sites['trid'].str.contains('pathogenic'))&(cds_sites['trid'].isin(cds_final['trid'].tolist()))]

    print('#CDS sites', cds_sites['trid'].nunique())

    cds_sites['allele_range'] = cds_sites.apply(lambda row: row['max_len'] - row['min_len'], axis=1)
    cds_sites['gene'] = cds_sites.apply(lambda row: utr_cds_grp.at[row['trid'], 'gene'], axis=1)
    cds_sites['context'] = cds_sites.apply(lambda row: utr_cds_grp.at[row['trid'], 'final_context'], axis=1)
    cds_sites['SegDup'] = cds_sites.apply(lambda row: utr_cds_grp.at[row['trid'], 'SegDup'], axis=1)
    cds_sites['TR_TAG'] = cds_sites.apply(lambda row: utr_cds_grp.at[row['trid'], 'TR_TAG'], axis=1)
    cds_sites['all_het'] = cds_sites.apply(lambda row: cds_allele_het.at[row['trid'], 'all_het'] if row['trid'] in cds_allele_het.index else 0, axis=1)
    cds_sites['SV_TAG'] = cds_sites.apply(lambda row: 'TR_SV' if row['trid'] in svtr_tbl.index else 'TR', axis=1)

    cds_sites.to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_cds_stats/cds_{stats_type}_annot.tsv',sep='\t', header=True, index=False)


    utr_sites = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_utr_stats/allele_{stats_type}_percentiles.tsv', sep='\t', header=0)
    utr_sites = utr_sites.loc[(~utr_sites['trid'].str.contains('pathogenic'))&(utr_sites['trid'].isin(utr_final['trid'].tolist()))]

    print('#UTR sites', utr_sites['trid'].nunique())
    utr_sites['allele_range'] = utr_sites.apply(lambda row: row['max_len'] - row['min_len'], axis=1)
    utr_sites['gene'] = utr_sites.apply(lambda row: utr_cds_grp.at[row['trid'], 'gene'], axis=1)
    utr_sites['context'] = utr_sites.apply(lambda row: utr_cds_grp.at[row['trid'], 'final_context'], axis=1)
    utr_sites['SegDup'] = utr_sites.apply(lambda row: utr_cds_grp.at[row['trid'], 'SegDup'], axis=1)
    utr_sites['TR_TAG'] = utr_sites.apply(lambda row: utr_cds_grp.at[row['trid'], 'TR_TAG'], axis=1)
    utr_sites['all_het'] = utr_sites.apply(lambda row: utr_allele_het.at[row['trid'], 'all_het'] if row['trid'] in utr_allele_het.index else 0, axis=1)
    utr_sites['SV_TAG'] = utr_sites.apply(lambda row: 'TR_SV' if row['trid'] in svtr_tbl.index else 'TR', axis=1)
    utr_sites.to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_utr_stats/utr_{stats_type}_annot.tsv', sep='\t',header=True, index=False)

    cds_sites = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_cds_stats/cds_{stats_type}_annot.tsv',
        sep='\t', header=0)
    utr_sites = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_utr_stats/utr_{stats_type}_annot.tsv',
        sep='\t', header=0)

    all_sites = pd.concat([cds_sites, utr_sites])


    print(all_sites.groupby('context').size().reset_index(name='count'))

def main():
    group_trgt_outputs('GRCh38_0711', 'utr', 'All')
    group_trgt_outputs('GRCh38_0711', 'cds', 'All')

    tr_heterozygosity(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_cds_stats')
    tr_heterozygosity(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/matt_trgt/tr_utr_stats')

    complete_sites_annot('length')


if __name__ == '__main__':
    main()
