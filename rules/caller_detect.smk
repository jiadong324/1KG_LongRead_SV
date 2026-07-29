

global SAMPLES
global DETECT_CALLERS
global REF_DICT

def find_ref(wildcards):
    return REF_DICT[wildcards.refv]

def get_h1(wildcards):
    # manifest_df = pd.read_csv(MANIFEST[wildcards.ref],sep='\t',index_col=['NAME'])
    return SAMPLES.at[wildcards.sample, 'HAP1']

def get_h2(wildcards):
    # manifest_df = pd.read_csv(MANIFEST[wildcards.ref],sep='\t',index_col=['NAME'])
    return SAMPLES.at[wildcards.sample, 'HAP2']
def find_bam(wildcards):
    return SAMPLES.at[wildcards.sample, 'BAM']

rule sawfish_discover:
    input:
        bam = find_bam,
        ref = find_ref
    output:
        bcf = '{refv}/{sample}/sawfish_disc/candidate.sv.bcf'
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "sawfish/0.12.4",
    resources:
        mem = 50,
        hrs = 24,
        disk_free = 1,
    threads: 4
    shell:
        """
        sawfish discover --threads {threads} --ref {input.ref} --bam {input.bam} --output-dir $( dirname {output.bcf} ) --clobber
        """

rule sawfish_call:
    input:
        bcf = '{refv}/{sample}/sawfish_disc/candidate.sv.bcf'
    output:
        vcf = '{refv}/{sample}/sawfish_call/genotyped.sv.vcf.gz'
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "sawfish/0.12.4",
    resources:
        mem = 50,
        hrs = 24,
        disk_free = 1,
    threads: 4
    shell:
        """
        sawfish joint-call --threads {threads} --sample $( dirname {input.bcf} ) --clobber --output-dir $( dirname {output.vcf} )
        """


rule sniffles:
     input:
        bam = find_bam,
        ref = find_ref
     output:
         vcf= "{refv}/{sample}/{sample}.sniffles.vcf"
     log:
         "log/{refv}/{sample}.sniffles.log",
     envmodules:
         "modules",
         "modules-init",
         "modules-gs/prod",
         "modules-eichler/prod",
         "sniffles/2.2",
     resources:
         mem=10,
         hrs=24,
         disk_free=1,
     threads: 4
     shell:
         """
         sniffles -i {input.bam} --reference {input.ref} --output-rnames -v {output.vcf} -t {threads}
         """

rule delly:
    input:
        bam=find_bam,
        ref=find_ref
    output:
        vcf = "{refv}/{sample}/{sample}.delly.bcf"
    log:
        "log/{refv}/{sample}.cutesv.log",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "delly/1.2.6",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """
        delly lr -y pb -o {output.vcf} -g {input.ref} {input.bam}
        """

rule delly_bcf:
    input:
        bcf = rules.delly.output.vcf
    output:
        vcf = "{refv}/{sample}/{sample}.delly.vcf"
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
        bcftools view {input.bcf} > {output.vcf}
        rm {input.bcf}
        """

rule cuteSV:
    input:
        bam = find_bam,
        ref = find_ref
    output:
        vcf = "{refv}/{sample}/{sample}.cutesv.vcf"
    log:
        "log/{refv}/{sample}.cutesv.log",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "cuteSV/2.1.0",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 8
    shell:
        """
        cuteSV -t {threads} --write_old_sigs --genotype -l 50 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 {input.bam} {input.ref} {output.vcf} {wildcards.refv}/{wildcards.sample}
        """


rule pbsv_discover:
    input:
        bam = find_bam,
        trf = find_ref
    output:
        sig = "{refv}/{sample}/{sample}.svsig.gz",
    log:
        "log/{refv}/{sample}.pbsv.log",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "pbconda/202403",
    resources:
        mem=50,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        pbsv discover --hifi --tandem-repeats {input.trf} {input.bam} {output.sig}
        """


rule pbsv_call:
    input:
        sig = "{refv}/{sample}/{sample}.svsig.gz",
        ref = find_ref
    output:
        vcf= "{refv}/{sample}/{sample}.pbsv.vcf"
    log:
        "log/{refv}/{sample}.pbsv.log",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "pbconda/202403",
    resources:
        mem=50,
        hrs=24
    threads: 4
    shell:
        """
        pbsv call --hifi -m 50 -j {threads} {input.ref} {input.sig} {output.vcf}
        rm {input.sig}
        """

rule svision:
    input:
        bam=find_bam,
        ref=find_ref,
    output:
        vcf= "{refv}/{sample}/{sample}.svision.vcf"
    log:
        "log/{refv}/{sample}.svision.log",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "SVision/1.4",
    resources:
        mem=20,
        hrs=24,
        disk_free=1,
    threads: 8
    shell:
        """
        SVision -o $( dirname {output.vcf} ) -s 10 -b {input.bam} -t {threads} -n {wildcards.sample} -g {input.ref} -m /net/eichler/vol28/projects/medical_reference/nobackups/svision_model/svision-cnn-model.ckpt
        mv {wildcards.refv}/{wildcards.sample}/{wildcards.sample}.svision.*.vcf {output.vcf}
        """

rule svisionpro:
    input:
        bam = find_bam,
        ref = find_ref,
    output:
        vcf="{refv}/{sample}/{sample}.svisionpro.vcf"
    log:
        "log/{refv}/{sample}.svision.log",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "svision-pro/2.3",
    resources:
        mem=30,
        hrs=24,
        disk_free=1,
    threads: 8
    shell:
        """
        SVision-pro --out_path $( dirname {output.vcf} ) --min_supp 10 --preset hifi --target_path {input.bam} --process_num {threads} --sample_name {wildcards.sample} --genome_path {input.ref} --model_path /net/eichler/vol28/projects/medical_reference/nobackups/SVision-pro/src/pre_process/model_liteunet_256_8_16_32_32_32.pth
        mv {wildcards.refv}/{wildcards.sample}/{wildcards.sample}.svision_pro_*.vcf {output.vcf}
        """

rule hapdiff:
    input:
        h1 = get_h1,
        h2 = get_h2,
        ref = find_ref
    output:
        vcf="{ref}/{sample}/hapdiff_phased.vcf.gz",

    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0",
        "hapdiff/0.9"
    resources:
        mem=20,
        hrs=24,
        disk_free=1,
    threads: 6
    shell:
        """
        hapdiff.py --reference {input.ref} --pat {input.h1} --mat {input.h2} --out-dir $( dirname {output.vcf} ) -t {threads} --sample {wildcards.sample} --sv-size 50
        """

rule caller_detect:
    input:
        expand("{refv}/{sample}/{sample}.{var_caller}.vcf", sample=SAMPLES.index, var_caller=DETECT_CALLERS.split(',')),
        expand("{ref}/{sample}/hapdiff_phased.vcf.gz", sample=SAMPLES.index)
    message:
        "Caller detection complete"