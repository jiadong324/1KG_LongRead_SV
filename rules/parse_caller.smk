## Normalize each VCF file for truvari usage

global SAMPLES
global REF_DICT
global MERGE_CALLERS

def find_ref(wildcards):
    return REF_DICT[wildcards.refv]

rule norm_caller_vcf:
    input:
        vcf="{refv}/{sample}/{sample}.{caller}.vcf",
        ref=find_ref
    output:
        insdel="{refv}/{sample}/{sample}.{caller}.insdel.vcf.gz",
        # inv="{refv}/{sample}/{sample}.{caller}.inv.vcf.gz"
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "truvari/4.2.1",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        bcftools norm --multiallelics - --output-type v {input.vcf} | python scripts/resolve.py /dev/stdin {wildcards.caller} {input.ref} |  bcftools norm --check-ref s --fasta-ref {input.ref} -N -m-any | bcftools annotate -x 'INFO/AF,INFO/STRAND' > {wildcards.refv}/{wildcards.sample}/{wildcards.caller}.tmp.vcf
        bcftools view -i "(SVTYPE=='INS'||SVTYPE=='DEL')&FILTER=='PASS'" -O v {wildcards.refv}/{wildcards.sample}/{wildcards.caller}.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}           
        tabix -p vcf {output.insdel}
        rm {wildcards.refv}/{wildcards.sample}/{wildcards.caller}.tmp.vcf
        """

rule parse_hapdiff:
    input:
        vcf = "{ref}/{sample}/hapdiff_phased.vcf.gz",
        ref = find_ref
    output:
        insdel="{ref}/{sample}/{sample}.hapdiff.insdel.vcf.gz",
        # inv="{ref}/{sample}/{sample}.hapdiff.inv.vcf.gz",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "truvari/4.2.1",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        bcftools norm --multiallelics - --output-type v {input.vcf} | python /net/eichler/vol28/projects/medical_reference/nobackups/Scripts/MedRef/parsers/resolve.py /dev/stdin svimasm | bcftools norm --check-ref s --fasta-ref {input.ref} -N -m-any > {wildcards.ref}/{wildcards.sample}/hapdiff.tmp.vcf
        bcftools view -i "SVTYPE=='INS'||SVTYPE=='DEL'" -O v {wildcards.ref}/{wildcards.sample}/hapdiff.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}
        tabix -p vcf {output.insdel}
        rm {wildcards.ref}/{wildcards.sample}/hapdiff.tmp.vcf
        """

rule parse_pav:
    input:
        vcf='{ref}/PAV/pav_{sample}.vcf.gz',
        ref=find_ref
    output:
        insdel="{ref}/{sample}/{sample}.pav.insdel.vcf.gz",
        inv="{ref}/{sample}/{sample}.pav.inv.vcf.gz",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "truvari/4.2.1",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        python ../scripts/Pav2SV.py {input.vcf} {wildcards.sample} | bcftools norm --multiallelics - --output-type v /dev/stdin | python ../scripts/resolve.py /dev/stdin pav {input.ref} |  bcftools norm --check-ref s --fasta-ref {input.ref} -N -m-any > {wildcards.ref}/{wildcards.sample}/pav.tmp.vcf
        bcftools view -i "SVTYPE=='INS'||SVTYPE=='DEL'" -O v {wildcards.ref}/{wildcards.sample}/pav.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}
        bcftools view -i "SVTYPE=='INV'" -O v {wildcards.ref}/{wildcards.sample}/pav.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}
        tabix -p vcf {output.insdel}
        tabix -p vcf {output.inv}
        rm {wildcards.ref}/{wildcards.sample}/pav.tmp.vcf
        """

rule parse_dipcall:
    input:
        vcf='{ref}/PAV/pav_{sample}.vcf.gz',
        ref=find_ref
    output:
        insdel="{ref}/{sample}/{sample}.pav.insdel.vcf.gz",
        inv="{ref}/{sample}/{sample}.pav.inv.vcf.gz",
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
    threads: 1
    shell:
        """
        bcftools norm --multiallelics - --output-type v {input.vcf} | python ../scripts/Dipcall2SV.py /dev/stdin {input.ref} | bcftools sort -o /dev/stdout -O z - > {output.insdel}
        tabix -p vcf {output.insdel}
        """


rule caller_vcf_list:
    input:
        caller_vcf=expand("{refv}/{sample}/{sample}.{caller}.insdel.vcf.gz", sample=SAMPLES.index, calllers=MERGE_CALLERS.split(','), refv=REF_DICT.key())
    output:
        vcf_list = '{refv}/{sample}/vcf_list.txt'
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    run:
        fout = open(output.vcf_list, 'w')
        for line in input.caller_vcf:
            print(line, file=fout)


rule parse_caller:
    input:
        expand("{ref}/{sample}/vcf_list.txt", caller=MERGE_CALLERS.split(','), samples=SAMPLES.index, ref=REF_DICT.keys())
    message:
        "Caller normalization complete"