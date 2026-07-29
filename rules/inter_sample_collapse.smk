
import pandas as pd


global SAMPLES
global MERGE_CALLERS
global REF_DICT

def assign_disc_class(samples: int, total: int):
    if samples == 1:
        return "SINGLE", 0
    elif 1 <= samples < (total / 2):
        return "POLY", samples/total
    elif (total / 2) <= samples < total:
        return "MAJOR", samples/total
    elif samples == total:
        return "SHARED", samples/total



rule sample_vcf_list:
    input:
        sample_vcf = expand("{refv}/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.vcf.gz",sample=SAMPLES.index,calllers=MERGE_CALLERS.split(','),refv=REF_DICT.key())
    output:
        vcf_list = '{refv}/vcf_list.txt'

rule bcftool:
    input:
        vcf = vcf_list
    output:
        outvcf = "disco_bcftools_merge.normed.vcf.gz"
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
        bcftools merge --threads 5 --merge none --force-samples --file-list {input.vcf} -O z | bcftools norm --threads {threads} --do-not-normalize --multiallelics -any --output-type z -o {output.outvcf}
        tabix -p vcf {output.outvcf}
        """

rule truvari_collapse:
    input:
        vcf = rules.bcftool.output.outvcf,
    output:
        removed_vcf = "disco_truvari_removed.vcf.gz",
        collapsed_vcf = "disco_truvari_collapsed.vcf.gz",
        samples = "samples.out"
    resources:
        mem=64,
        hrs=24,
    threads: 4
    params:
        opts=config.get("truvari_opts")
    shell:
        """
        truvari collapse --input {input.vcf} --collapsed-output {output.removed_vcf} --sizemin 50 --sizemax 100000 --keep common --gt all {params.opts} | bcftools sort --max-mem $( expr {threads} \\* {resources.mem} )G --output-type z > {output.collapsed_vcf}
        tabix -p vcf {output.collapsed_vcf}
        bcftools query -l {output.collapsed_vcf} > {output.samples}
        """

rule collapse_data_table:
    input:
        vcf="disco_truvari_collapsed.vcf.gz",
        samples = "samples.out"
    output:
        bed="tables/disco_truvari_collapsed.bed.gz",
        tsv="tables/disco_truvari_collapsed.tsv.gz",
        haps="tables/disco_truvari_collapsed_haps.bed.gz"
    threads: 1
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    run:
        sample_list = [line.strip() for line in open(input.samples)]
        sample_haps = []
        for sample in sample_list:
            sample_haps.append(f'{sample}_1')
            sample_haps.append(f'{sample}_2')

        var_info = []
        headers = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'DISC_CLASS', 'DISC_FREQ', 'MERGE_SAMPLES', 'MERGE_GT']
        haps_header = ['ID'] + sample_haps
        vcf = pysam.VariantFile(input.vcf,'r')
        all_haps = []
        for rec in vcf.fetch():
            svlen = rec.info['SVLEN'] if type(rec.info['SVLEN']) is int else rec.info['SVLEN'][0]
            end = int(rec.pos) + 1 if rec.info['SVTYPE'] == 'INS' else int(rec.pos) + abs(int(svlen))
            this_rec = [rec.chrom, rec.pos, end, rec.id, rec.info['SVTYPE'], svlen]
            this_rec_haps = [rec.id] + ['.' for i in range(len(sample_haps))]

            gts = []
            merged_samples = []
            for (sample, gt) in rec.samples.items():
                h1 = '.' if gt.get('GT')[0] is None else gt.get('GT')[0]
                h2 = '.' if gt.get('GT')[1] is None else gt.get('GT')[1]
                if f'{h1}|{h2}' != '.|.':
                    gts.append(f'{h1}|{h2}')
                    merged_samples.append(sample)

                h1_idx = sample_haps.index(f'{sample}_1')
                h2_idx = sample_haps.index(f'{sample}_2')
                this_rec_haps[h1_idx + 1] = h1
                this_rec_haps[h2_idx + 1] = h2

            all_haps.append(this_rec_haps)
            disc_class, freq = assign_disc_class(len(merged_samples), len(sample_list))
            this_rec.append(disc_class)
            this_rec.append(round(freq, 4))
            this_rec.append(','.join(merged_samples))
            this_rec.append(','.join(gts))
            var_info.append(this_rec)
        var_out = pd.DataFrame(var_info, columns=headers)
        var_out[['#CHROM', 'POS', 'END', 'SVTYPE', 'ID', 'SVLEN']].to_csv(output.bed, sep='\t', index=False, compression='gzip', header=False)
        var_out.to_csv(output.tsv, sep='\t', index=False, compression='gzip', header=True)

        pd.DataFrame(all_haps, columns=haps_header).to_csv(output.haps, compression='gzip', sep='\t', index=False, header=True)


rule inter_sample:
    input:
        "tables/disco_truvari_collapsed.bed.gz",
