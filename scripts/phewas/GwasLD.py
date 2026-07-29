#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 12/16/25
'''
import json
import gzip
import math
import joblib as jl
import pysam
import unicodedata
import re
from collections import Counter
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pickle
import matplotlib as mpl

from helpers.Constants import *
from intervaltree import IntervalTree

def normalize(text):
    text = unicodedata.normalize("NFKD", text.lower().strip())
    return re.sub(r'[^\w\s]', '', text)

def get_sig_snp():
    gwas_tbl = pd.read_csv(f'{VOL28}/SVREF/GWAS/gwas-catalog-download-associations-alt-full.tsv', sep='\t', header=0)

    snp_out = open('/Volumes/eichler-vol28/projects/medical_reference/nobackups/SVREF/GWAS/gwas_catalog_v1.0.sig.bed',
                   'w')
    for idx, row in gwas_tbl.iterrows():
        disease = row['MAPPED_TRAIT']
        chrid = 'chr' + str(row['CHR_ID'])
        # if row['CHR_POS'] == np.NaN:
        #     continue
        try:
            pos = int(row['CHR_POS'])
            gene = row['MAPPED_GENE']
            risk_allele = row['STRONGEST SNP-RISK ALLELE']
            snp_id = row['SNPS']
            snp_af = row['RISK ALLELE FREQUENCY']
            context = row['CONTEXT']
            pval = row['P-VALUE']
            if pval < 5e-8:
                print(f'{chrid}\t{pos}\t{pos + 1}\t{gene}\t{risk_allele}\t{snp_id}\t{snp_af}\t{pval}\t{context}\t{disease}',file=snp_out)
        except ValueError:
            print('Position unknown')

def sv_gwas_snp_ld():
    omim_genes = pd.read_csv('/Volumes/eichler-vol28/projects/medical_reference/nobackups/SVREF/OMIM/genemap.txt',
                             sep='\t', header=[0])
    # omim_genes = pd.read_csv(f'/Volumes/eichler-vol28/projects/medical_reference/nobackups/SVREF/OMIM/genemap.txt',sep='\t', header=[0])
    omim_genes.dropna(subset=['Approved Gene Symbol'], inplace=True)
    omim_genes.drop_duplicates(subset=['Approved Gene Symbol'], inplace=True)
    omim_genes.set_index('Approved Gene Symbol', inplace=True)

    # sd_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/gt_af_sd.tsv.gz',sep='\t', usecols=[0,1,2,3], names=['chrom', 'start', 'end', 'svid'])
    # sd_tbl.set_index('svid', inplace=True)

    sv_sample_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/tables/disco_truvari_collapsed.tsv.gz',
                            sep='\t', index_col=['ID'])

    # cds_utr_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/gencode_v45/sv_cds_utr.tsv', sep='\t',
    #     names=['chrom', 'start', 'end', 'svid', 'svtype', 'af', 'element', 'gene'])
    # cds_utr_tbl.set_index('svid', inplace=True)

    sample_pop = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/sample_pop.txt',
                            sep='\t', names=['sample', 'pop', 'cohort'])
    sample_pop.set_index('sample', inplace=True)

    combi_sv_annot_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/final_annot/fill_tags_cadd_pli_gene.tsv', sep='\t', index_col=['ID'])

    # sv_cds_utr_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/Annot/gt_af_gene.tsv', sep='\t', names=['chrom', 'start', 'end', 'svid', 'svtype', 'af', 'gene', 'context'])
    # sv_cds_utr_tbl.set_index('svid', inplace=True)
    gwas_tbl = pd.read_csv(f'{VOL28}/SVREF/GWAS/gwas_catalog_v1.0.sig.srt.bed',
                            sep='\t', names=['chrom', 'pos', 'end', 'gene', 'risk_allele', 'snp_id', 'snp_af', 'pval', 'context', 'trait'])


    # tr_sv_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/complete/InTR/svs_popaf.tsv.gz', sep='\t', index_col=['ID'])

    gwas_tbl_group = gwas_tbl.groupby('snp_id')

    aou_377_intervals = {}
    for line in open(f'{VOL28}/SVREF/AoU/Iris_377_SV.bed'):
        chrom, start, end, svid, svtype, svlen = line.strip().split('\t')
        if chrom not in aou_377_intervals:
            aou_377_intervals[chrom] = IntervalTree()
        aou_377_intervals[chrom][int(start) - 1000: int(end) + 1000] = (int(start) - 1000, int(end) + 1000, svid)


    snp_info_list = {}

    for snp, group in gwas_tbl_group:
        gene_set = set()
        traits = set()
        context = set()
        pos = group['pos'].values[0]
        for idx, row in group.iterrows():
            gene_set.add(row['gene'])
            traits.add(row['trait'])
            context.add(row['context'])

        gene_out = ';'.join([str(ele) for ele in list(gene_set)])
        trait_out = ','.join([str(ele) for ele in list(traits)])
        context_out = ';'.join([str(ele) for ele in list(context)])
        snp_info_list[pos] = [snp, len(traits), gene_out, trait_out, context_out]
        snp_info_list[pos+1] = [snp, len(traits), gene_out, trait_out, context_out]


    sv_snp_pair = []
    sv_counter = []
    traits_list = []
    strong_sv_counter = []
    strong_trait_list = []
    for num in range(22):
        chrom = f'chr{num+1}'
        ld_tbl = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/gwas_traits/{chrom}/r2_output.tsv', sep='\t', header=0)
        for idx, row in ld_tbl.iterrows():
            var_a, var_b = row['SNP_A'], row['SNP_B']
            var_a_sv = var_a.find('INS') !=-1 or var_a.find('DEL') != -1
            var_b_sv = var_b.find('INS') !=-1 or var_b.find('DEL') != -1

            if var_a_sv != var_b_sv:
                sv_af = combi_sv_annot_tbl.at[var_a, 'AF'] if var_a_sv else combi_sv_annot_tbl.at[var_b, 'AF']
                sv_samples = sv_sample_tbl.at[var_a, 'MERGE_SAMPLES'] if var_a_sv else sv_sample_tbl.at[var_b, 'MERGE_SAMPLES']

                sv_pop = set()
                sv_cohort = set()

                for sample in sv_samples.split(','):
                    sv_pop.add(sample_pop.at[sample, 'pop'])
                    sv_cohort.add(sample_pop.at[sample, 'cohort'])


                sv_maf = combi_sv_annot_tbl.at[var_a, 'MAF'] if var_a_sv else combi_sv_annot_tbl.at[var_b, 'MAF']
                sv_cadd = combi_sv_annot_tbl.at[var_a, 'CADD'] if var_a_sv else combi_sv_annot_tbl.at[var_b, 'CADD']
                sv_func = combi_sv_annot_tbl.at[var_a, 'Func'] if var_a_sv else combi_sv_annot_tbl.at[var_b, 'Func']
                sv_reg = combi_sv_annot_tbl.at[var_a, 'REG'] if var_a_sv else combi_sv_annot_tbl.at[var_b, 'REG']
                sv_brain_reg = combi_sv_annot_tbl.at[var_a, 'brainREG'] if var_a_sv else combi_sv_annot_tbl.at[var_b, 'brainREG']
                sv_context = combi_sv_annot_tbl.at[var_a, 'CONTEXT'] if var_a_sv else combi_sv_annot_tbl.at[var_b, 'CONTEXT']

                # if var_a_sv and var_a_sv in cds_utr_tbl.index:
                #     sv_context = cds_utr_tbl.at[var_a, 'context']
                #
                # if var_b_sv and var_b_sv in cds_utr_tbl.index:
                #     sv_context = cds_utr_tbl.at[var_b, 'context']

                sv_gene = combi_sv_annot_tbl.at[var_a, 'GENE'] if var_a_sv else combi_sv_annot_tbl.at[var_b, 'GENE']
                gene_pli = combi_sv_annot_tbl.at[var_a, 'PLI'] if var_a_sv else combi_sv_annot_tbl.at[var_b, 'PLI']

                omim = []
                if sv_gene!='.':
                    for gene in sv_gene.split(';'):
                        if gene in omim_genes.index:
                            omim.append(gene)
                omim_out = ';'.join(omim)
                if not omim:
                    omim_out = '.'


                snp_pos = var_a.split(':')[1] if not var_a_sv else var_b.split(':')[1]
                if int(snp_pos) not in snp_info_list:
                    continue
                snp_info = snp_info_list[int(snp_pos)]
                traits_list.extend(snp_info[3].split(';'))
                if row['R2'] >= 0.8:
                    strong_trait_list.extend(snp_info[3].split(';'))

                snp_id = var_a if not var_a_sv else var_b
                snp_gwas_id = snp_info[0]
                # snp_traits = snp_info[3]

                if var_a_sv:
                    tr_motif = combi_sv_annot_tbl.at[row['SNP_A'], 'TR_MOTIF'] if row['SNP_A'] in combi_sv_annot_tbl.index else '.'
                    tr_id = combi_sv_annot_tbl.at[row['SNP_A'], 'MOTIF_ID'] if row['SNP_A'] in combi_sv_annot_tbl.index else '.'
                    sv_counter.append(row['SNP_A'])
                    is_sd = combi_sv_annot_tbl.at[row['SNP_A'], 'SegDup']
                    chrom, start, svtype, svlen = row['SNP_A'].split('-')

                    ovlps = aou_377_intervals[chrom].overlap(int(start), int(start) + 1)
                    aou_sv = '.'
                    if ovlps:
                        for ovlp in ovlps:
                            aou_sv = ovlp.data[2]

                    if row['R2']>=0.8:
                        strong_sv_counter.append(row['SNP_A'])

                    this_pair = ([row['CHR_A'], row['BP_A'], row['SNP_A'], row['CHR_B'], row['BP_B'], snp_gwas_id, row['R2']] +
                                 [snp_id, aou_sv, sv_af, sv_maf, sv_cadd, sv_gene, sv_context, omim_out, gene_pli, sv_func, sv_reg, sv_brain_reg, is_sd, tr_motif, tr_id, ','.join(list(sv_pop)), ','.join(list(sv_cohort))] + snp_info[2:] + [row['SNP_A'], snp_gwas_id])
                else:
                    tr_motif = combi_sv_annot_tbl.at[row['SNP_B'], 'TR_MOTIF'] if row['SNP_B'] in combi_sv_annot_tbl.index else '.'
                    tr_id = combi_sv_annot_tbl.at[row['SNP_B'], 'MOTIF_ID'] if row['SNP_B'] in combi_sv_annot_tbl.index else '.'
                    sv_counter.append(row['SNP_B'])
                    is_sd = combi_sv_annot_tbl.at[row['SNP_B'], 'SegDup']
                    chrom, start, svtype, svlen = row['SNP_B'].split('-')
                    ovlps = aou_377_intervals[chrom].overlap(int(start), int(start) + 1)
                    aou_sv = '.'
                    if ovlps:
                        for ovlp in ovlps:
                            aou_sv = ovlp.data[2]
                    if row['R2']>=0.8:
                        strong_sv_counter.append(row['SNP_B'])
                    this_pair = [row['CHR_A'], row['BP_A'], snp_gwas_id, row['CHR_B'], row['BP_B'], row['SNP_B'], row['R2']] + [snp_id, aou_sv, sv_af, sv_maf, sv_cadd, sv_gene, sv_context, omim_out, gene_pli, sv_func, sv_reg, sv_brain_reg, is_sd,tr_motif, tr_id, ','.join(list(sv_pop)), ','.join(list(sv_cohort))] + snp_info[2:] + [row['SNP_B'], snp_gwas_id]
                sv_snp_pair.append(this_pair)

    print('Total pairs', len(sv_snp_pair))

    print('\t#Unique SVs:', len(Counter(sv_counter)))
    print('\t#Unique traits:', len(Counter(traits_list)))

    print('\t#Unique SVs (r2>=0.8):', len(Counter(strong_sv_counter)))
    print('\t#Unique traits (r2>=0.8):', len(Counter(strong_trait_list)))

    df_sv_snps = pd.DataFrame(sv_snp_pair, columns=['CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B', 'R2', 'GATK_SNP_ID', 'AoU_Phase1_SVID', 'AF', 'MAF', 'CADD', 'SV_GENE', 'SV_GENE_CONTEXT', 'OMIM', 'SV_GENE_PLI', 'FuncElement', 'REG', 'brainREG', 'SegDup', 'TR_MOTIF', 'MOTIF_ID', 'SV_POP', 'SV_COHORTS', 'GENE', 'TRAITS', 'LOCATION', 'SVID', 'SNP'])
    df_sv_snps.to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/gwas_traits/SV_SNP_pairs.tsv', sep='\t', header=True, index=False)

def find_disease_traits():
    ref_traits = f'{VOL28}/SVREF/AoU/Disease_disorder_summary.txt'
    gwas_traits = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/gwas_traits/SV_SNP_pairs.tsv', sep='\t', header=0)

    with open(ref_traits, 'r') as f:
        norm_ref_traits = [normalize(line.strip()) for line in f if line.strip()]


    ## create regex patterns for reference traits
    regex_patterns = []
    for trait in norm_ref_traits:
        words = trait.split()
        if words:
            words[-1] = rf"{words[-1]}s?"
            pattern = re.compile(rf"\b{' '.join(words)}\b", re.IGNORECASE)
            regex_patterns.append(pattern)

    ## search gwas trait in the reference patterns
    combined_pattern = "|".join(f"(?:{p.pattern})" for p in regex_patterns)
    disease_trait = []
    disease_trait_set = set()
    for idx, row in gwas_traits.iterrows():
        is_disease_trait = False
        for ele in str(row['TRAITS']).split(','):
            if re.findall(combined_pattern, normalize(ele)):
                is_disease_trait = True
                disease_trait_set.add(ele)
        if is_disease_trait:
            sv_id = row['SNP_A'] if 'INS' in row['SNP_A'] or 'DEL' in row['SNP_A'] else row['SNP_B']
            snp_id = row['SNP_A'] if 'rs' in row['SNP_A'] else row['SNP_B']
            disease_trait.append(row.tolist() + [sv_id, snp_id])

    print("#Medically relevant trait:", len(disease_trait_set))

    df_disease_traits = pd.DataFrame(disease_trait, columns=list(gwas_traits.columns) + ['SVID', 'SNP'])
    df_disease_traits.to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/gwas_traits/SV_SNP_pairs_disease_trait.tsv', sep='\t', header=True, index=False)

def check_ld_results():
    all_svs = pd.read_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/gwas_traits/SV_SNP_pairs.tsv', sep='\t', header=0)
    all_svs = all_svs.loc[all_svs['R2']>=0.8]
    cds_utr_svs = all_svs.loc[(all_svs['SV_GENE_CONTEXT'].str.contains('Exonic')) | (all_svs['SV_GENE_CONTEXT'].str.contains('UTR'))]
    noncoding = all_svs.loc[(~all_svs['SV_GENE_CONTEXT'].str.contains('Exonic')) & (~all_svs['SV_GENE_CONTEXT'].str.contains('UTR'))]
    noncoding_reg = noncoding.loc[(noncoding['REG']=='YES')|(noncoding['brainREG']=='YES')]

    final_tmp = all_svs.loc[(all_svs['R2'] >= 0.8) & (all_svs['SV_GENE'] != '.') & (all_svs['MOTIF_ID'] == '.') & (all_svs['SegDup'] == 'OutSD')]
    final_svset = final_tmp.loc[~final_tmp['SV_GENE'].str.contains('LOC')]
    final_svset.to_csv(f'{VOL28}/SVREF/all_cohorts/GRCh38_0711/pct09_dist500/gwas_traits/final_ld_svset.tsv', sep='\t',
                       header=True, index=False)

    print('# CDS_UTR: ', len(cds_utr_svs))
    print('# REG: ', len(noncoding_reg))


def main():
    # sv_gwas_snp_ld()
    find_disease_traits()

    # check_ld_results()


if __name__ == '__main__':
    main()
