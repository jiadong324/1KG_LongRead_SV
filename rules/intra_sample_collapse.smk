import pandas as pd
import pysam
import numpy as np
import json
import os
from intervaltree import IntervalTree
import lib

global SAMPLES
global HEADERS_DICT
global REF_DICT
global CPR_DICT
global MERGE_CALLERS

def find_header(wildcards):
    return HEADERS_DICT[wildcards.refv]

def find_ref(wildcards):
    return REF_DICT[wildcards.refv]
def find_cprmask(wildcards):
    return CPR_DICT[wildcards.refv]


rule bcftool_all:
    input:
        vcf="{refv}/{sample}/vcf_list.txt"
    output:
        outvcf="{refv}/{sample}/caller_merge/insdel.tmp.vcf.gz"
    log:
        "log/{refv}/{sample}.bcftool.log",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0"
    threads: 5
    shell:
        """
        bcftools merge --thread {threads} --merge none --force-samples -O z -o {output.outvcf} --file-list {input.vcf}
        tabix -p vcf {output.outvcf}
        """

rule truvari:
    input:
        bcfvcf=rules.bcftool_all.output.outvcf
    output:
        removed="{refv}/{sample}/caller_merge/removed.vcf.gz",
        collapse="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.vcf.gz"
    log:
        "log/{refv}/{sample}.truvari.log",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "truvari/4.3.1"
    shell:
        """
        truvari collapse -i {input.bcfvcf} -c {output.removed} --sizemin 50 --sizemax 100000 --gt het -k first --intra --pctseq 0.90 --pctsize 0.90 --refdist 500 | bcftools sort --max-mem 8G -O z -o {output.collapse}
        tabix -p vcf {output.collapse}
        """


rule parse_truvari_collapse:
    input:
        vcf=rules.truvari.output.collapse,
        headers=find_header,
        excl=find_cprmask
    params:
        callers=config.get("CALLERS"),
        read_callers=config.get("READ_CALLERS"),
        sample=config.get('SAMPLE'),
    output:
        svs_bed="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.bed.gz",
        filt_vcf="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.filt.vcf",
        filt_bed="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.filt.bed.gz",
        pav_supp_vcf="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.vcf",
        pav_supp_bed="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.bed.gz",
        pav_bed="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.pav-only.bed.gz",
        pav_supp_fa="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.fa",
        pav_fa="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.pav-only.fa",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    run:
        import gzip

        AUTOSOMESXY = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                       "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", 'chrY']

        read_based_callers = params.read_callers.split(',')
        sample = params.sample
        caller_list = params.callers.split(',')

        ## Read excluded regions to intervaltree
        excl_dict = {}
        for line in open(input.excl):
            entries = line.strip().split('\t')
            chrom, start, end = entries[0], entries[1], entries[2]
            if chrom not in excl_dict:
                excl_dict[chrom] = IntervalTree()
            excl_dict[chrom][int(start): int(end)] = (int(start), int(end))

        insdel_dict = {'Total': {'INS': 0, 'DEL': 0}}
        caller_dict = {}
        sv_records = []
        svs_bed, pav_supp_bed, pav_only_svs = [], [], []
        filt_svs_vcf, pav_supp_vcf = [], []

        alt_fa_out = open(output.pav_supp_fa,'w')
        alt_asm_unique_fa = open(output.pav_fa,'w')

        for line in gzip.open(input.vcf,'rt'):
            if line.startswith('#'):
                continue
            entries = line.strip().split('\t')
            supp_vec = lib.parsers.decode_truvari_supp(int(entries[9].split(':')[-1]), len(caller_list))
            supp = sum([int(ele) for ele in supp_vec])

            if entries[0] not in AUTOSOMESXY:
                continue
            if supp == 0:
                continue

            info_dict = lib.parsers.parse_vcf_info_column(entries[7])
            svlen = abs(int(info_dict['SVLEN']))
            if svlen < 50 or svlen > 100000:
                continue

            svtype = info_dict['SVTYPE']
            end = int(entries[1]) + 1 if svtype == 'INS' else int(entries[1]) + svlen

            if excl_dict[entries[0]].overlaps(int(entries[1]),end):
                continue

            insdel_dict['Total'][svtype] += 1
            if f'SUPP={supp}' not in insdel_dict:
                insdel_dict[f'SUPP={supp}'] = {'INS': 0, 'DEL': 0}
            insdel_dict[f'SUPP={supp}'][svtype] += 1

            info_dict['SVLEN'] = -svlen if svtype == 'DEL' else svlen
            info_dict['END'] = end

            supp_caller = ''
            is_read_sv = False
            is_asm_sv = False
            for i, val in enumerate(supp_vec):
                if val == '1':
                    supp_caller += f'{caller_list[i]},'
                    if caller_list[i] in read_based_callers:
                        is_read_sv = True
                    if caller_list[i] == 'pav' or caller_list[i] == 'hapdiff':
                        is_asm_sv = True

            if supp_caller[:-1] not in caller_dict:
                caller_dict[supp_caller[:-1]] = {'INS': 0, 'DEL': 0}
            caller_dict[supp_caller[:-1]][svtype] += 1

            info_dict['SUPP_CALLER'] = supp_caller
            new_info_str = ';'.join(['{0}={1}'.format(key,info_dict[key]) for key in
                                     ['SVLEN', 'END', 'SVTYPE', 'SUPP_CALLER']])
            new_sv_id = f'{entries[0]}-{int(entries[1]) + 1}-{svtype}-{abs(svlen)}'
            svs_bed.append([entries[0], entries[1], end, svlen, svtype, entries[2], new_sv_id, entries[9].split(':')[0],
                            supp_caller[:-1]])

            ## Assembly-only SVs
            if supp_caller[:-1] == 'pav':
                phased_gt = entries[9].split(':')[0]
                sv_record = f'{entries[0]}\t{entries[1]}\t{entries[2]}\t{entries[3]}\t{entries[4]}\t{entries[5]}\tPASS\t{new_info_str[:-1]}\tGT\t{phased_gt}'
                filt_svs_vcf.append(sv_record)

                pav_only_svs.append([entries[0], entries[1], end, svlen, svtype, new_sv_id, entries[9].split(':')[0],
                                     supp_caller[:-1], sample])
                if svtype == 'INS':
                    print(f'>{entries[2]}',file=alt_asm_unique_fa)
                    print(f'{entries[4]}',file=alt_asm_unique_fa)
                if svtype == 'DEL':
                    print(f'>{entries[2]}',file=alt_asm_unique_fa)
                    print(f'{entries[3]}',file=alt_asm_unique_fa)

            ## Keep SVs supported by at least two callers
            # if len(supp_caller[:-1].split(',')) >= 2 and 'pav' in supp_caller:
            if is_asm_sv and is_read_sv:

                phased_gt = entries[9].split(':')[0]

                sv_record = f'{entries[0]}\t{entries[1]}\t{new_sv_id}\t{entries[3]}\t{entries[4]}\t{entries[5]}\tPASS\t{new_info_str[:-1]}\tGT\t{phased_gt}'
                pav_supp_vcf.append(sv_record)
                filt_svs_vcf.append(sv_record)
                pav_supp_bed.append([entries[0], entries[1], end, svlen, svtype, new_sv_id, phased_gt, supp_caller[:-1],
                                     sample])

                if svtype == 'INS':
                    print(f'>{new_sv_id}',file=alt_fa_out)
                    print(f'{entries[4]}',file=alt_fa_out)
                if svtype == 'DEL':
                    print(f'>{new_sv_id}',file=alt_fa_out)
                    print(f'{entries[3]}',file=alt_fa_out)

        alt_asm_unique_fa.close()
        alt_fa_out.close()

        filt_sv_writer = open(output.filt_vcf,'w')
        pav_supp_vcf_writer = open(output.pav_supp_vcf,'w')

        for line in open(input.headers):
            print(line.strip(),file=filt_sv_writer)
            print(line.strip(),file=pav_supp_vcf_writer)

        print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}',file=filt_sv_writer)
        print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}',file=pav_supp_vcf_writer)

        for sv_line in filt_svs_vcf:
            print(sv_line,file=filt_sv_writer)

        for sv_line in pav_supp_vcf:
            print(sv_line,file=pav_supp_vcf_writer)

        ## Save all collapsed SVs to BED file
        sorted_bed = sorted(svs_bed,key=lambda x: (x[0], int(x[1])))
        pd.DataFrame(sorted_bed,columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'VCF_ID', 'GT',
                                         'SUPP']).to_csv(output.svs_bed, sep='\t',header=False,index=False,compression='gzip')

        filt_svs_bed = sorted(pav_supp_bed + pav_only_svs,key=lambda x: (x[0], int(x[1])))
        pd.DataFrame(filt_svs_bed,columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GT', 'SUPP',
                                           'MERGE_SAMPLES']).to_csv(output.filt_bed, sep='\t',header=True,index=False,compression='gzip')

        ## Save PAV SVs supported by callers to BED file
        sorted_supp_bed = sorted(pav_supp_bed,key=lambda x: (x[0], int(x[1])))
        pd.DataFrame(sorted_supp_bed,columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GT', 'SUPP',
                                              'MERGE_SAMPLES']).to_csv(output.pav_supp_bed,sep='\t',header=True,index=False,compression='gzip')

        ## Save PAV only SVs to BED file
        sorted_pav_bed = sorted(pav_only_svs,key=lambda x: (x[0], int(x[1])))
        pd.DataFrame(pav_only_svs,columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GT', 'SUPP',
                                           'MERGE_SAMPLES']).to_csv(output.pav_bed,sep='\t',header=True,index=False,compression='gzip')

        ## Save SV counts by number of supporting callers
        with open(f'{wildcards.refv}/{wildcards.sample}/caller_merge/truvari_collapsed.insdel.stats.json','w') as f:
            json.dump(insdel_dict,f)

        ## Save SVs detected by different combination of callers
        with open(f'{wildcards.refv}/{wildcards.sample}/caller_merge/truvari_collapsed.insdel.stats.json','w') as f:
            json.dump(caller_dict,f)


rule index_filt_vcf:
    input:
        filt_vcf="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.filt.vcf",
        pav_supp_vcf="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.vcf"
    output:
        gzip="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.filt.vcf.gz",
        gzip_pav="{refv}/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.vcf.gz",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0"
    shell:
        """
        bcftools sort -o /dev/stdout -O v {input.filt_vcf} | bgzip -c > {output.gzip}
        tabix -p vcf {output.gzip}
        rm {input.filt_vcf}

        bcftools sort -o /dev/stdout -O v {input.pav_supp_vcf} | bgzip -c > {output.gzip_pav}
        tabix -p vcf {output.gzip_pav}
        rm {input.pav_supp_vcf}    
        """

rule intra_sample_collapse:
    input:
        expand(rules.index_filt_vcf.output.gzip, sample=SAMPLES.index, ref=REF_DICT.keys())