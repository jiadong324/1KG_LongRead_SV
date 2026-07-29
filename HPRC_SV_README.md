# HPRC R2 SV

This is a readme file of haplotype phased SVs for HPRC R2 [Lucas et al.](https://www.biorxiv.org/content/10.64898/2026.07.21.739710v1). SVs are seperated by INS/DEL and INV. 
SVs are for each genome are available at https://human-pangenomics.s3.amazonaws.com/index.html?prefix=submissions/759B21AD-0ED8-4640-A433-7C92A57EA3D3--UW_EEE_SV_Calls/
If you have any issues please contact jdlin@uw.edu and ee3@uw.edu.


## 08-26-2025 updates

Uploaded truvari collapsed calls ```{sample}.truvari_collapsed.insdel.vcf.gz``` of different callers. The collapse parameters ```--pctseq 0.8 --pctsize 0.8 --refdist 1000 --sizemin 50 --sizemax 100000 --gt het -k maxqual --intra```. 

## 05-21-2025 updates

The ```TIG_REGION``` is added to the INFO column for each SV. This is identified by PAV.

## 04-15-2025 updates

1. Complete SVs for NA20503, NA20762, NA20806, NA20827 and HG06807
2. Uploaded BAM/PAF for assembly alignment to reference genomes. The alignment pipeline is https://github.com/mrvollger/asm-to-reference-alignment.


## 04-08-2025 updates

SVs are detected based on 232 genomes released in https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assemblies_pre_release_v0.6.1.index.csv.
A total of 233 samples are called, including one genome (HG03492) not in the index file but found at NCBI.

##### File descriptions
For each genome, we created the following files:
1. ```{sample}.{caller}.vcf.gz``` is all variants detected by each caller without any filtering.
2. ```pav_{sample}.vcf.gz``` is all variant detected by PAV without any filtering.
3. ```insdel/{sample}.{caller}.insdel.vcf.gz``` is the normalized INS/DEL calls of each caller. They are further used to create the final callset ```insdel/{sample}.insdel.vcf``` for each genome.
4. ```insdel/{sample}.insdel.vcf``` qced INS/DEL for each genome. SVs at complex regions such as centromere, gaps, etc are excluded. These regions are obtained from https://www.biorxiv.org/content/10.1101/2024.09.24.614721v1.

##### Tandem repeat annotation
We used tandem repeat annotation from https://zenodo.org/records/13178746.

## 18-12-2024 updates

A total of 202 samples in ```assembly_manifest.tab``` are detected with the following callers on both GRCh38 and T2T reference genomes:
1. PAV (v2.3.4)
2. dipcall (v0.3)
3. hapdiff (v0.9): this pipeline used SVIM-asm (v1.0.2) for SV calling.
4. pbsv (v2.9.0)
5. delly (v1.2.6)
6. cutesv (v2.1.0)
7. sniffles (v2.2)
8. sawfish (v0.12.4)

Note that this list included 196 genomes released in https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assemblies_pre_release_v0.2.index.csv. Six genomes (HG005, HG00733, HG01243, HG02145, HG03492 and NA19240) are HPRC year1 samples assembled with hifiasm v0.14. HG00733 and NA19240 are also assembled by HGSVC3 using verkko v1.4 but we only uploaded SVs detected from hifiasm assemblies. 

These are all raw SV calls without any filtering. We are still trying to fix the PAV running issue for HG03130, NA187974 and HG03516 on GRCh38 and will upload afterward.

