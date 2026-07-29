
import os
import pandas as pd
import lib
import json

global SAMPLES
global REF_DICT
global SD_DICT
global TR_DICT
global GENE_DICT

def find_tr(wildcards):
    return TR_DICT[wildcards.refv]

def find_gene(wildcards):
    return GENE_DICT[wildcards.refv]

def find_segdup(wildcards):
    return SD_DICT[wildcards.refv]

rule summarize_svs_info:
    input:
        filt_bed=expand(rules.parse_truvari_collapse.output.filt_bed,sample=SAMPLES.index,refv=REF_DICT.keys()),
        pav_bed=expand(rules.parse_truvari_collapse.output.pav_bed,sample=SAMPLES.index,refv=REF_DICT.keys()),
    output:
        sv_num="sv_info/persample_sv_num.tsv.gz",
        sv_size100="sv_info/persample_sv_size100.tsv",
        sv_size1k="sv_info/persample_sv_size1k.tsv"
    params:
        cohort = config.get('COHORT')
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    run:
        # fout = open(output.alt_yaml,'w')
        # print('samples:',file=fout)
        # for alt_fa in input.filt_alt_fas:
        #     sample_name = alt_fa.split('/')[1]
        #     print(f" {sample_name}_filt: '{os.getcwd()}/{alt_fa}'",file=fout)
        #     os.system(f'samtools faidx {alt_fa}')
        #
        # for alt_fa in input.pav_alt_fas:
        #     sample_name = alt_fa.split('/')[1]
        #     print(f" {sample_name}_pav-only: '{os.getcwd()}/{alt_fa}'",file=fout)
        #     os.system(f'samtools faidx {alt_fa}')

        pop_dict = json.load(open('../config/igsr_samples_pop.json'))
        sv_num_list = []
        size100_list, size1k_list = [], []
        for filt in input.filt_bed:
            sample_name = filt.split('/')[1]
            ref = filt.split('/')[0]
            pop = pop_dict[sample_name] if sample_name in pop_dict else 'Unk'
            df_filt = pd.read_csv(filt,sep='\t',header=[0])
            sv_num_list.append([sample_name, pop, params.cohort, len(df_filt[df_filt['SVTYPE'] == 'INS']),
                                len(df_filt[df_filt['SVTYPE'] == 'DEL']), len(df_filt), 'Filt', ref])
            size100_df, size1k_df = lib.parsers.read_svsize_bed(df_filt,ref,sample_name)
            size100_df['Set'] = ['Filt' for _ in range(len(size100_df))]
            size100_df['Cohort'] = [params.cohort for _ in range(len(size100_df))]
            size1k_df['Set'] = ['Filt' for _ in range(len(size100_df))]
            size1k_df['Cohort'] = [params.cohort for _ in range(len(size100_df))]
            size100_list.append(size100_df)
            size1k_list.append(size1k_df)

        for filt in input.pav_bed:
            sample_name = filt.split('/')[1]
            ref = filt.split('/')[0]
            pop = pop_dict[sample_name] if sample_name in pop_dict else 'Unk'
            df_filt = pd.read_csv(filt,sep='\t',header=[0])
            sv_num_list.append([sample_name, pop, params.cohort, len(df_filt[df_filt['SVTYPE'] == 'INS']),
                                len(df_filt[df_filt['SVTYPE'] == 'DEL']), len(df_filt), 'PAV-only', ref])
            size100_df, size1k_df = lib.parsers.read_svsize_bed(df_filt,ref,sample_name)
            size100_df['Set'] = ['PAV-only' for _ in range(len(size100_df))]
            size100_df['Cohort'] = [params.cohort for _ in range(len(size100_df))]
            size1k_df['Set'] = ['PAV-only' for _ in range(len(size100_df))]
            size1k_df['Cohort'] = [params.cohort for _ in range(len(size100_df))]
            size100_list.append(size100_df)
            size1k_list.append(size1k_df)
        pd.DataFrame(sv_num_list,columns=['Sample', 'Pop', 'Cohort', 'INS_num', 'DEL_num', 'Total_num', 'Set',
                                          'REF']).to_csv(output.sv_num,sep='\t',header=True,index=False,compression='gzip')

        df_data_len100 = pd.concat(size100_list)
        df_data_len100.to_csv(output.sv_size100)
        df_data_len1k = pd.concat(size1k_list)
        df_data_len1k.to_csv(output.sv_size1k)

rule gene_svs:
    input:
        svs = '{refv}/{sample}/caller_merge/truvari_collapsed.insdel.{subset}.bed.gz',
        gene_region = find_gene,
    output:
        annot_bed = '{refv}/{sample}/caller_merge/Annot/{subset}.gene_svs.bed.gz',
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """        
        bedtools intersect -wo -a {input.svs} -b {input.gene_region} | awk -v OFS="\\t" '!a[$13]++{{print $1,$2,$3,$4,$5,$6,$13}}' | bgzip -c > {output.annot_bed}        
        """

rule tr_svs:
    input:
        svs = '{refv}/{sample}/caller_merge/truvari_collapsed.insdel.{subset}.bed.gz',
        tr_region = find_tr,
    output:
        annot_bed = '{refv}/{sample}/caller_merge/Annot/{subset}.tr_svs.bed.gz',
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """
        bedtools intersect -wo -a {input.svs} -b {input.tr_region} -f 0.5 | cut -f 1-8,11,12,16-24 | awk -v OFS='\\t' 'BEGIN {{print "#CHROM\\tPOS\\tEND\\tSVLEN\\tSVTYPE\\tID\\tGT\\tSUPP\\tTR_START\\tTR_END\\tN_MOTIF\\tMIN_MOTIFLEN\\tMAX_MOTIFLEN\\tMOTIFS\\tSTRUCT\\tSegDup\\tTELO\\tCEN\\tCDS"}};{{print $0}}'| bgzip -c > {output.annot_bed}
        """

rule sd_svs:
    input:
        svs = '{refv}/{sample}/caller_merge/truvari_collapsed.insdel.{subset}.bed.gz',
        sd_region = find_segdup,
    output:
        annot_bed = '{refv}/{sample}/caller_merge/Annot/{subset}.sd_svs.bed.gz',
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """
        bedtools intersect -wo -a {input.svs} -b {input.sd_region} -f 0.5 | awk -v OFS='\\t' '!a[$6]++{{print $1,$2,$3,$4,$5,$6,"SegDup"}}' | bgzip -c > {output.annot_bed}
        """

rule annot:
    input:
        svs='{refv}/{sample}/caller_merge/truvari_collapsed.insdel.{subset}.bed.gz',
        sd_annot=rules.sd_svs.output.annot_bed,
        tr_annot=rules.tr_svs.output.annot_bed,
        gene_annot=rules.gene_svs.output.annot_bed,
    output:
        annot_out = '{refv}/{sample}/caller_merge/Annot/{subset}.annot.bed.gz'
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    run:

        def get_gene(sv_id, gene_db):
            this_sv_tr = gene_db.loc[gene_db['ID'] == sv_id]
            if len(this_sv_tr) == 0:
                return '.'
            annot = this_sv_tr['GENE'].tolist()
            return ','.join(annot)

        def get_tr(sv_id, tr_db):
            this_sv_tr = tr_db.loc[tr_db['ID'] == sv_id]
            if len(this_sv_tr) == 0:
                return '.', '.', '.'
            cpx_tag = []
            cpx_tag.append('CDS' if this_sv_tr['CDS'].tolist()[0] != 0 else '')
            cpx_tag.append('SegDup' if this_sv_tr['SegDup'].tolist()[0] != 0 else '')
            cpx_tag.append('TELO' if this_sv_tr['TELO'].tolist()[0] != 0 else '')
            cpx_tag.append('CEN' if this_sv_tr['CEN'].tolist()[0] != 0 else '')
            cpx_out = ','.join(cpx_tag)
            tr_annot = this_sv_tr['MOTIFS'].tolist()
            for ele in tr_annot:
                motif_list = ele.split(',')
                for motif in motif_list:
                    if len(motif) > 6:
                        return 'VNTR', ','.join(motif_list), cpx_out
            motif_list = tr_annot[0]
            return 'STR', ','.join(motif_list), cpx_out

        def get_sd(sv_id, sd_db):
            this_sv_tr = sd_db.loc[sd_db['ID'] == sv_id]
            if len(this_sv_tr) == 0:
                return '.'
            return 'SegDup'

        tr_svs = pd.read_csv(input.tr_annot,sep='\t',header=[0])
        sd_svs = pd.read_csv(input.sd_annot,sep='\t',names=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'REPTYPE'])
        gene_svs = pd.read_csv(input.gene_annot,sep='\t',names=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GENE'])

        all_svs = pd.read_csv(input.svs, sep='\t',header=[0])
        sv_annots = []
        for idx, row in all_svs.iterrows():
            this_tr, tr_motif, cpx_tag = get_tr(row['ID'],tr_svs)
            this_sd = get_sd(row['ID'],sd_svs)
            this_gene = get_gene(row['ID'],gene_svs)
            sv_annots.append(row.tolist()[:-1] + [this_gene, this_sd, this_tr, tr_motif, cpx_tag])
        pd.DataFrame(sv_annots,columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GT', 'SUPP', 'GENE',
                                        'SD', 'TR', 'TR_MOTIF', 'CPX_TAG']).to_csv(output.annot_out,sep='\t',header=True,index=False,compression='gzip')

rule intra_sample_stats:
    input:
        expand(rules.annot.output.annot_out, refv=REF_DICT.keys(), sample=SAMPLES.index, subset=['filt', 'pav-supp', 'pav-only']),
        "sv_info/persample_sv_num.tsv.gz",
        "sv_info/persample_sv_size100.tsv",
        "sv_info/persample_sv_size1k.tsv"
