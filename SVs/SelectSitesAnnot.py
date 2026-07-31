#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 2/24/26
'''
import matplotlib.pyplot as plt
import pandas as pd
from venny4py.venny4py import *
from venn import venn
import numpy as np

from helpers.Constants import *

def update_reg(row):
    ele = row['ENCODE_cCRE']
    # print(ele)
    if 'PLS' in ele:
        return 'Promoter'
    elif 'CTCF-only' in ele:
        return 'CTCF'
    elif 'pELS' in ele:
        return 'pELS'
    elif 'dELS' in ele:
        return 'dELS'
    else:
        return '.'
def update_stratify(row):
    amr_fst = float(row['AMR_fst']) if row['AMR_fst'] != '.' else 0
    afr_fst = float(row['AFR_fst']) if row['AFR_fst'] != '.' else 0
    eas_fst = float(row['EAS_fst']) if row['EAS_fst'] !='.' else 0
    eur_fst = float(row['EUR_fst']) if row['EUR_fst'] != '.' else 0
    sas_fst = float(row['SAS_fst']) if row['SAS_fst'] != '.' else 0

    fst_dict = {'EAS': eas_fst,'EUR' :eur_fst, 'SAS': sas_fst, 'AFR': afr_fst, 'AMR': amr_fst}
    sorted_fst = sorted(fst_dict.items(), key=lambda x:x[1], reverse=True)

    stratified_pop = '.'
    if sorted_fst[0][1] >= 0.15:
        stratified_pop = sorted_fst[0][0]
    return stratified_pop


## locityper 1KG test
sv_tested = pd.read_csv('/Volumes/eichler-vol28/home/iwong1/nobackups/aou/svs_sr_locityper/svs_sr_locityper.vcf',sep='\t',comment='#',usecols=[0,1,2], names=['CHROM', 'POS', 'ID'])
sv_tested.set_index('ID', inplace=True)

## Biallelic SVs
bi_insdel = [line.strip() for line in open(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/allele_class/biallelic_insdel.txt')]

## load TR catalog SVs >= 50% covered by repeats
tr_svs = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/gt_af_tr.tsv.gz', sep='\t', index_col=['ID'])

## load trf SVs >= 50% covered by repeats
trf_svs = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/trf_simpleRepeat_SVs.tsv', sep='\t', names=['chrom', 'start', 'end', 'ID'])
trf_svs.set_index('ID', inplace=True)

## load SV_SNP_LD of disease phenotypes
sv_snp_disease_traits = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/gwas_traits/SV_SNP_pairs_disease_trait.tsv', sep='\t', header=0,)
strong_ld = sv_snp_disease_traits.loc[sv_snp_disease_traits['R2']>=0.8]
sv_strong_ld_groups = strong_ld.groupby('SVID').agg({'TRAITS': ';'.join, 'SNP': ';'.join}).reset_index()
sv_strong_ld_groups.rename(columns={'SVID': 'ID'}, inplace=True)

## load ALL SV_SNP_LD
sv_snp_traits = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/gwas_traits/SV_SNP_pairs.tsv', sep='\t', header=0,)
sv_snp_traits = sv_snp_traits.astype({'SNP': str, 'TRAITS': str})
sv_traits_strong_ld_groups = sv_snp_traits.loc[sv_snp_traits['R2']>=0.8].groupby('SVID').agg({'TRAITS': ';'.join, 'SNP': ';'.join}).reset_index()
sv_traits_strong_ld_groups.rename(columns={'SVID': 'ID'}, inplace=True)

## Load SV-eQTL sites
eqtl_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/Iris_eQTL/eQTL_hh_refpanel_annot.tsv', sep='\t', header=0)
eqtl_tbl['eGenes'] = eqtl_tbl.groupby('SV')['geneSymbol'].transform(','.join)
eqtl_tbl.drop_duplicates(subset=['SV'], keep='first', inplace=True)

eqtl_ids = eqtl_tbl['SV'].unique().tolist()
eqtl_tbl.set_index('SV', inplace=True)

## Load SV Fst result
# fst_db = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/FST/sv_fst_ac.tsv', sep='\t', header=0)
fst_db = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/FST_PCLAI/stratified_svs/sv_sites_func.tsv', sep='\t', header=0)
# fst_db['Stratified_Pop'] = fst_db.apply(lambda row: update_stratify(row), axis=1)
fst_db.set_index('SVID', inplace=True)

## Locityper batch2 AoU run
lt_batch2 = [line.strip() for line in open(f"{VOL28}/SVREF/locityper_inAoU/batch2_SVs/aou_1kg_tier1.txt")]


## Load Yang's compelete annotations
cadd_sv_tbl = pd.read_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/combine_anno_all_gnomad_V3.txt.gz", sep='\t', index_col=['ID'])
annot_tbl = pd.read_csv("/Volumes/eichler-vol28/projects/autism_genome_assembly/nobackups/yangsui/JD/JD_final/maxCDSallUTR/allreg_sv.txt", sep='\t', header=0)

## exclude sex chromosomes
annot_tbl = annot_tbl.loc[(~annot_tbl['ID'].str.contains('chrX'))&(~annot_tbl['ID'].str.contains('chrY'))]
annot_tbl.fillna({'ENCODE_cCRE': '.'}, inplace=True)
annot_tbl['ENCODE_REG'] = annot_tbl.apply(lambda row: update_reg(row), axis=1)

annot_tbl['locityper_batch'] = annot_tbl.apply(lambda row: 'Batch2' if row['ID'] in lt_batch2 else '.', axis=1)
# annot_tbl['eQTL'] = annot_tbl.apply(lambda row: 'Yes' if row['ID'] in eqtl_ids else '.', axis=1)
annot_tbl['eGene'] = annot_tbl.apply(lambda row: eqtl_tbl.at[row['ID'], 'eGenes'] if row['ID'] in eqtl_ids else '.', axis=1)
annot_tbl['Stratified_Pop'] = annot_tbl.apply(lambda row: fst_db.at[row['ID'], 'Stratified_POP'] if row['ID'] in fst_db.index else '.', axis=1)
annot_tbl['CADD-SV_PHRED-scorev1.1.2'] = annot_tbl.apply(lambda row: cadd_sv_tbl.at[row['ID'], 'CADD-SV_PHRED-scorev1.1.2'], axis=1)
annot_tbl['LocationSV_AnnotSV'] = annot_tbl.apply(lambda row: cadd_sv_tbl.at[row['ID'], 'LocationSV_AnnotSV'], axis=1)
annot_tbl['Gene_AnnotSV'] = annot_tbl.apply(lambda row: cadd_sv_tbl.at[row['ID'], 'Gene_AnnotSV'], axis=1)



print(annot_tbl['CMRG'].nunique())

header_to_remove = ['#CHROM', 'POS', 'END']
# annot_names = list(annot_tbl.columns)
annot_names = [item for item in list(annot_tbl.columns) if item not in header_to_remove]

## Add AF/AC information to Yang's annotation table
kg_info = pd.read_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/svs_ac_af.tsv", sep='\t', header=0)
kg_info_new_names = {}
for i, ele in enumerate(kg_info.columns):
    if i > 5:
        kg_info_new_names[ele] = f'1KG_LR_{ele}'
kg_info.rename(columns=kg_info_new_names, inplace=True)

annot_tbl = pd.merge(annot_tbl, kg_info, on='ID', how='left')
annot_tbl.fillna('.', inplace=True)


## Output SV-eQTL sites
eqtl_sites = annot_tbl.loc[annot_tbl['ID'].isin(eqtl_ids)]
eqtl_sites['isNovel'] = eqtl_sites.apply(lambda row: eqtl_tbl.at[row['ID'], 'New_eGenes'], axis=1)
print('SV-eQTL sites: ', len(eqtl_sites))

out_columns = annot_names + ['1KG_LR_AC', '1KG_LR_AF', '1KG_LR_AF_AFR', '1KG_LR_AF_AMR', '1KG_LR_AF_EUR', '1KG_LR_AF_SAS', '1KG_LR_AF_EAS']
eqtl_sites[out_columns + ['isNovel']].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/eqtl_sv_annot.tsv", sep='\t', header=True, index=False)

## CADD>=20
subset_cadd = annot_tbl.loc[annot_tbl['CADD-SV_PHRED-scorev1.1.2']!='.']
subset_cadd['CADD-SV_PHRED-scorev1.1.2'] = subset_cadd['CADD-SV_PHRED-scorev1.1.2'].astype(float)
cadd_out = subset_cadd.loc[subset_cadd['CADD-SV_PHRED-scorev1.1.2'] >= 20]
cadd_out[out_columns].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/high_caddsv.tsv", sep='\t', header=True, index=False)

## CDS sites
cds_sites = annot_tbl.loc[annot_tbl['Location_GENCODE']=='CDS']
print('CDS sites: ', len(cds_sites))
cds_sites[out_columns].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/cds_sv_annot.tsv", sep='\t', header=True, index=False)

## SV in LD of GWAS disease phenotype traits
annot_tbl_disease_ld_merge = pd.merge(annot_tbl, sv_strong_ld_groups, on='ID', how='left')
annot_tbl_disease_ld_merge.fillna('.', inplace=True)

tagged_svs = annot_tbl_disease_ld_merge.loc[annot_tbl_disease_ld_merge['SNP']!='.']
out_columns = out_columns + ['TRAITS', 'SNP']
tagged_svs[out_columns].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/disease_ld_sv_annot.tsv", sep='\t', header=True, index=False)

## SV in LD with GWAS traits
annot_tbl_ld_merge = pd.merge(annot_tbl, sv_traits_strong_ld_groups, on='ID', how='left')
annot_tbl_ld_merge.fillna('.', inplace=True)
tagged_traits_svs = annot_tbl_ld_merge.loc[annot_tbl_ld_merge['SNP']!='.']
tagged_traits_svs['Stratified_Pop'] = tagged_traits_svs.apply(lambda row: fst_db.at[row['ID'], 'Stratified_POP'] if row['ID'] in fst_db.index else '.', axis=1)
tagged_traits_svs[out_columns].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/traits_ld_sv_annot.tsv", sep='\t', header=True, index=False)

## Exonic (CDS+UTR)
# exonic = annot_tbl_ld_merge.loc[(annot_tbl_ld_merge['LocationSV_AnnotSV'] == 'Exonic')&(annot_tbl_ld_merge['gene']!='.')]
exonic = annot_tbl_ld_merge.loc[(annot_tbl_ld_merge['Location_GENCODE'] == 'CDS')|(annot_tbl_ld_merge['Location_GENCODE'].str.contains('UTR'))]
print('Exonic sites: ', len(exonic))
exonic[out_columns].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/exonic_sv_annot.tsv", sep='\t', header=True, index=False)

## REG SVs, excluded CDS/UTR
reg_svs = annot_tbl_ld_merge.loc[(annot_tbl_ld_merge['ENCODE_cCRE'] != '.')&(annot_tbl_ld_merge['LocationSV_AnnotSV'] != 'Exonic')]
reg_svs[out_columns].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/reg_sv_annot.tsv", sep='\t', header=True, index=False)

## integrated table for genotyping candidates
eqtl_sites = eqtl_sites.loc[eqtl_sites['Location_GENCODE']!='.']
cadd_out = cadd_out.loc[cadd_out['Location_GENCODE']!='.']
tagged_svs = tagged_svs.loc[tagged_svs['Gene_AnnotSV']!='.']
stratified_genic = annot_tbl.loc[(annot_tbl['Stratified_Pop'].isin(['AFR','AMR']))&(annot_tbl['Gene_AnnotSV']!='.')&(annot_tbl['1KG_LR_AC']>1)]

print('#Stratified AFR/AMR ', len(stratified_genic))

## Save final annotation table
annot_tbl.to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/Final_annot_tbl.tsv', sep='\t', header=True, index=False)

eqtl_sites['source'] = ['SV-eQTL' for _ in range(len(eqtl_sites))]
cadd_out['source'] = ['CADD-SV' for _ in range(len(cadd_out))]
tagged_svs['source'] = ['LD-SNP' for _ in range(len(tagged_svs))]
exonic['source'] = ['CDS-UTR' for _ in range(len(exonic))]
# stratified_genic['source'] = ['Stratified' for _ in range(len(stratified_genic))]

# candidates = pd.concat([eqtl_sites, cadd_out, tagged_svs, exonic, stratified_genic])
candidates = pd.concat([eqtl_sites, cadd_out, tagged_svs, exonic])

## Exclude singletons and inside TR/SD.
# candidates = candidates.loc[(candidates['1KG_LR_AC'] > 1) & (candidates['platinumTRs']=='.')&(candidates['SegDup']=='.')&(candidates['SimpleRepeat']=='.')&(~candidates['RepeatMasker'].str.contains('Simple_repeat'))]
# candidates = candidates.loc[(candidates['1KG_LR_AC'] > 1) & (~candidates['ID'].isin(tr_svs.index)) & (~candidates['ID'].isin(trf_svs.index)) & (candidates['SegDup']=='.')]

## Get biallelic SVs
# candidates = candidates.loc[(candidates['1KG_LR_AC'] > 1) & (candidates['ID'].isin(bi_insdel))]
candidates = candidates.loc[(candidates['ID'].isin(bi_insdel))& (candidates['Location_GENCODE']!='.')]

# candidates[out_columns + ['source']].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/svs_sr_locityper.tsv", sep='\t', header=True, index=False)

## Load SR test results
sr_res = pd.read_csv(f'{VOL28}/SVREF/locityper_inAoU/isaac_loci_eval.tsv',sep='\t', index_col=['ID'])
sr_res.replace({np.nan: 'NA'}, inplace=True)

## Load batch3 SVs
# sr_res2 = pd.read_csv('/Volumes/eichler-vol28/home/iwong1/nobackups/aou/svs_sr_locityper/svs_sr_locityper.vcf', sep='\t', comment='#', usecols=[2], names=['SVID'])
# sr_res2.set_index('SVID', inplace=True)

batch2_3_svs = pd.read_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/svs_sr_locityper.tsv", sep='\t', index_col=['ID'])
batch1_svs = pd.read_csv(f"{VOL28}/SVREF/locityper_inAoU/batch1_1KG_SVID.txt", sep='\t', index_col=['locus'])

all_batches = list(batch1_svs.index) + list(batch2_3_svs.index)
# candidates['locityper_batch'] = candidates.apply(lambda row: 'batch1' if row['ID'] in batch1_svs.index else '.', axis=1)
candidates['locityper_batch'] = candidates.apply(lambda row: 'batch1_2_3' if row['ID'] in all_batches else '.', axis=1)
candidates['F1'] = candidates.apply(lambda row: sr_res.at[row['ID'], 'F1'] if row['ID'] in sr_res.index else '.', axis=1)

# candidates['SR_TestBatch'] = candidates.apply(lambda row: 'Tested' if row['ID'] in sr_res2.index else '.', axis=1)

not_tested = candidates.loc[(~candidates['ID'].isin(sv_tested.index))&(candidates['1KG_LR_AC'] > 1)]
# not_tested[['ID']].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/allele_class/SVtoTest_0608.tsv", sep='\t', header=True, index=False)
print('#SVs to test', len(not_tested))
# candidates[out_columns + ['F1','source']].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/allele_class/svs_sr_locityper.tsv", sep='\t', header=True, index=False)
candidates[out_columns + ['F1','source']].to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/allele_class/svs_sr_locityper_pclai.tsv", sep='\t', header=True, index=False)

print('#SVs: ', candidates['ID'].nunique())

## check stratified CDS/UTR to genotype
# stratified_cds_utr = pd.read_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/FST_PCLAI/stratified_svs/biallelic_cds_utr_svs.tsv",sep='\t', header=0)
# cds_utr_toGT = stratified_cds_utr.loc[~stratified_cds_utr['SVID'].isin(candidates['ID'].tolist())]
# cds_utr_toGT.to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/allele_class/cds_utr_stratified_toTest.tsv", sep='\t',header=True, index=False)

# batch3_svs = candidates.loc[candidates['locityper_batch']=='.']
# pd.Series(batch3_svs['ID'].unique()).to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/locityper_inAoU/batch3_SVs/aou_1kg_batch3.txt', index=False, header=False)

print(len(set(lt_batch2) - set(candidates['ID'].tolist())))

ac2_candidates = candidates.loc[(candidates['1KG_LR_AC'] > 1)&(candidates['Location_GENCODE']!='.')]
ac2_candidates.to_csv(f"{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/allele_class/svs_sr_locityper_pclai_finalTbl.tsv", sep='\t', header=True,index=False)

duplicate_pairs = ac2_candidates[ac2_candidates.duplicated(subset=['ID'])]
duplicate_pairs_grp = duplicate_pairs.groupby('ID').agg(source_list=('source', ','.join),
    source_num=('source', 'count')).reset_index()

print(duplicate_pairs_grp)
print(len(duplicate_pairs_grp))

# sets = {
#     'SV-eQTL': set(ac2_candidates.loc[ac2_candidates['source']=='SV-eQTL']['ID'].tolist()),
#     'CADD-SV': set(ac2_candidates.loc[ac2_candidates['source']=='CADD-SV']['ID'].tolist()),
#     'LD-SNP': set(ac2_candidates.loc[ac2_candidates['source']=='LD-SNP']['ID'].tolist()),
#     'CDS/UTR': set(ac2_candidates.loc[ac2_candidates['source']=='CDS-UTR']['ID'].tolist()),
#     # 'Stratified': set(candidates.loc[candidates['source']=='Stratified']['ID'].tolist())
# }
#
# # venny4py(sets=sets,out=f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/svs_sr_venn')
# venny4py(sets=sets,ext='svg',out=f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/allele_class/svs_ac2_venn')

# fig, ax = plt.subplots(figsize=(7, 6))
# venn(sets, ax=ax)
# fig.tight_layout()
# fig.savefig(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/yang_annot/svs_sr_venn5.png', dpi=300)