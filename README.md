# 1KG_LongRead_SV

This repository contains workflows and analysis scripts for creating a SV control set from 1,725 1KG long-read genomes. 
These genomes are sequenced by HPRC ([Lucas et al. Biorxiv 2026](https://www.biorxiv.org/content/10.64898/2026.07.21.739710v1)), 
HGSVC ([Logsdon et al. Nature 2025](https://www.nature.com/articles/s41586-025-09140-6)), UW-ONT and IB-ONT. 

The UW-ONT sequenced a total of 500 genomes, 400 are novel and 100 genomes were published by [Gustafson et al. 2025](https://genome.cshlp.org/content/34/11/2061).
The IB-ONT genomes are published by Siegfried Schloissnig et.al. ([Nature 2025](https://www.nature.com/articles/s41586-025-09290-7)).

<div align=left><img width=20% height=20% src="https://github.com/jiadong324/1KG_LongRead_SV/blob/main/images/HGSVC_logo.png"/> <img width=30% height=30% src="https://github.com/jiadong324/1KG_LongRead_SV/blob/main/images/hprc_log.png"/><img width=40% height=40% src="https://github.com/jiadong324/1KG_LongRead_SV/blob/main/images/1000G-ONT.png"/></div> 


<div align=left><img width=100% height=100% src="https://github.com/jiadong324/1KG_LongRead_SV/blob/main/images/Fig1_Samples.png"/> </div>

## Genomes

The table below lists the source of 1,218 genomes used to create the callset in Lin et al. 
Note that the number below counts for the unique genomes.
The SVs detected from 1,218 genomes have been used to identify potential pathogenic variants and new SV-disease associations in biobanks.


HG002 and HG005 are included as part of the HPRC dataset. We also included NA12877 and NA12878. NA12878 is also one of the genomes sequenced by GIAB and HGSVC. 
NA12877 is sequenced by IB-ONT but we used the assembly and data published in [David Porubsky et al. Nature 2025](https://www.nature.com/articles/s41586-025-08922-2)

| Dataset | Genomes | Platform      |
|---------|---------|---------------|
| HPRC    | 232     | HiFi, UL-ONT  |
| HGSVC   | 61      | HiFi, UL-ONT  |
| 1KG-ONT | 480     | ONT (R9, R10) |
| IB-ONT  | 445     | ONT (R9)      |


## SV callset

### Individual genomes

Multiple caller merged callset for each genome. Please check the index file (TBA) for download.

### Integrated callset

We provide both GRCh38 and T2T-CHM13 callsets (zenodo link, TBA) for 293 HPRC+HGSVC genomes and 1,218 genomes.

**CHM13_INSDEL_HGSVC_HPRC_wAF.vcf.gz:** Integrated SVs from HGSVC/HPRC genomes with estimated allele frequency.

**CHM13_INSDEL_1218_wAF.vcf.gz:** Integrated SVs from all dataset containing 1,218 genomes with estimated allele frequency.

**GRCh38_INSDEL_HGSVC_HPRC_wAF.vcf.gz:** Integrated SVs from HGSVC genomes with estimated allele frequency.

**GRCh38_INSDEL_1218_wAF.vcf.gz:** Integrated SVs from all dataset containing 1,218 genomes with estimated allele frequency.

## Genome assembly

### HPRC

The HPRC genomes are assembled with hifiasm and verkko. Please refer to HPRC [release](https://github.com/human-pangenomics/hprc_intermediate_assembly). 

### HGSVC

Assemblies used in Glennis Logsdon, et.al. were created by Verkko v1.1 and further phased with StrandSeq data using Graphasing. 
With latest Verkko (v2.2.1), we resolved previous genomes with large unassigned contigs using HiC data, including HG00096, HG00514, HG02587, HG03009, HG03683, NA18939, NA19434 and NA19650. 
The rest of the HGSVC genomes were not reassembled. 

### UW-ONT
Consistent with our previous publication ([Gustafson et al. 2025](https://genome.cshlp.org/content/34/11/2061)) on the first 100 genomes, the NAPU pipeline of the same version is used to build the assembly for R9 and R10 data. 
Briefly, this pipeline used shasta to make the haploid assembly and create the pseudo-haplotypes with hapdup. 
Yak (v0.1) and QUAST (v5.2.1) were used to assess the QV and contig N50 of the pseudo-phased assembly. The QV was calculated with Illumina short-read data.


## SV discovery

### Alignment
HiFi reads are aligned with [pbmm2](https://github.com/PacificBiosciences/pbmm2) v1.13.1 ‘--preset HiFi’. 
1KG-ONT reads are aligned with minimap2 in the NAPU pipeline used in previous publication by Jonas Gustafson et.al.
We used IB-ONT alignment directly from the publication for SV discovery and phasing.

Minimap2 v2.28 is used to align assembly to both references. The alignment pipeline is [here](https://github.com/mrvollger/asm-to-reference-alignment).

### SV callers

| Tool        | Input type | Version |
|-------------|------------|---------|
| PAV         | Assembly   | v2.3.4  | 
| Dipcall     | Assembly   | v0.3    | 
| hapdiff     | Assembly   | v0.9    |
| pbsv        | HiFi       | v2.9.0  |
| sawfish     | HiFi       | v0.12.4 |
| sniffles    | HiFi, ONT  | v2.2    |
| delly       | HiFi, ONT  | v1.2.6  |
| cutesv      | HiFi, ONT  | v2.1.0  |
| Nanovar     | ONT        | v1.8.0  |
| debreak     | ONT        | v1.2.0  |
| SVision     | HiFi, ONT  | v1.4    |
| SVision-pro | HiFi, ONT  | v2.3    |

### Per genome SV

For each genome, we prioritize the PAV calling results and identify SVs supported by at least one another caller with [Truvari (v5.2.0)](https://github.com/acenglish/truvari). 
The pipeline is ```rules/intra_sample_collapse.smk```. The output of this pipeline is a multi-caller integrated VCF used the PAV reported breakpoint, sv length, phased genotype etc.
We then run [BoostSV](https://github.com/jiadong324/BoostSV) on the multi-caller integrated VCF for each genome.

We also used the same annotation as [Logsdon et al. Nature 2025](https://www.nature.com/articles/s41586-025-09140-6) to exclude SVs inside complex regions, gaps, etc. 
Briefly, these regions include UCSC gaps and centromere on GRCh38. 
For T2T-CHM13, complex regions include centromere, acrocentric p-arms, satellite regions except for monomeric satellite.

### TR genotyping

Tandem repeat catalogs for GRCh38 and T2T-CHM13 can be found [here](https://zenodo.org/records/13178746).


| Tool   | Input type | Version | Dataset             | Purpose        |
|--------|------------|---------|---------------------|----------------|
| TRGT   | HiFi       | v1.4.1  | HPRC, HGSVC, UW-ONT | Tandem repeats |
| vamos  | Assembly   | v2.1.5  | HPRC, HGSVC, UW-ONT | Tandem repeats |


### Cohort-level integration

#### Callable regions
We first defined the callable regions for each genome. HGSVC/HPRC and 1KG-ONT callable regions were created by PAV based on the assembly to reference alignment. 
For IB-ONT genomes, we splitted the PMDV phased BAM into two read sets. The callable regions for each haplotype were created by merging each read set (>= 3 reads) into non–overlapping intervals with BEDtools merge ‘-d 500’. 

#### Create callset

The non-redundant set integrated SVs from HPRC, HGSVC, UW-ONT and IB-ONT genomes with coverage >= 15x and read N50 >=15 kbp. 
Truvari (v5.2.0) was used to create the non-redundant SV set for both references with ‘--pctseq 0.90 –pctsize 0.90 –refdist 500 –keep common’. 
Briefly, different alleles at the same SV site were collapsed if they share minimum 90% sequence similarity and 90% allele size similarity. 
The calls with the highest quality predicted by BoostSV were used to represent each collapsed SV site. 
Moreover, this integration only considered INS/DEL ranging from 50bp to 100,000bp. 

We then filled the missing genotypes ‘./.’ with reference genotypes ‘0|0’, ‘0|.’ and ‘.|0’ for each genome based on callable regions. 
Note that the missing genotype in the final VCF only suggests there is no confident read or assembly alignments. 
For each integrated SV, we also kept the allele breakpoint position (FORMAT/APOS) and length (FORMAT/AL) from each sample. 
BCFtools (v1.16) plugin function fill-tags is used to calculate the statistics for each SV site, including allele frequency, minor allele frequency, etc (https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html). 

## Phenome-wide associations

We selected 2,752 high-impact SVs and genotyped them in 230 thousands AoU short-read samples to identify SV-disease associations. 
The WDL pipeline for genotyping and phenotype testing is available at [AoU_WDL](https://github.com/EichlerLab/AoU_WDL).

### Annotation

We selected potential high-impact SVs from callset **GRCh38_INSDEL_1218_wAF.vcf.gz:**. 
SVs were annotated using the comREG pipeline adapted from a previously published framework, with minor modifications (https://github.com/EichlerLab/asap, v1.1.0). 
Briefly, SVs aligned to GRCh38 were intersected with CDS and UTR from the same annotation file used for TR sites annotation. Gene annotations were further validated using the RefSeq database through [AnnotSV (v3.4)](https://github.com/lgmgeo/AnnotSV). 
Regulatory regions (denoted as REG in this study) were annotated using datasets from UCSC genome browser, including ENCODE cCREs (candidate cis-regulatory elements), ORegAnno, GeneHancer, H3K27ac, H3K4me1, and H3K4me3, together with brain-derived epigenomic profiles [Sui et al. 2026](https://www.nature.com/articles/s41467-026-68378-4). 
Additional functional annotations included ENCODE TF clusters, UCSC noncoding RNA tracks (tRNA, snRNA, lincRNA, and sno/miRNA), and repetitive regions (UCSC Segmental Duplications, RepeatMasker, Simple Repeats, and Tandem Repeats–Platinum).

To further evaluate potential functional impact, we incorporated the noncoding constraint Gnocchi score from gnomAD (Chen et al. 2024) and the CADD-SV score (v1.1.2) (Kleinert and Kircher 2022). 
Population SV frequencies from 63,046 unrelated short-read genomes were assessed using gnomAD SVs (gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz, https://gnomad.broadinstitute.org/data#v4-structural-variants). 
Comparison with the gnomAD SV callset was performed using Truvari (v5.2.0) with the command ```truvari bench -c {gnomad.v4.1.INSDEL50.non_neuro_controls.sites.vcf.gz} -b {collapsedSV.vcf.gz} --pctsize 0.5 --pctseq 0.0 --pctovl 0.5 --sizefilt 50 --sizemax 100000 -o {output}```

### Genotyping

SVs are first evaluated with matched 1KG SR data. The SV with F1-score>=0.8 are further genotyped in AoU short-read samples with the pipeline in [AoU_WDL](https://github.com/EichlerLab/AoU_WDL)

### Associations

We applied firth regression model in PLINK2 to conduct PheWAS analysis with covariant, phenotype and genotype files as input. 

#### Prepare inputs
For biallelic SVs and flanking SNPs, the variant genotypes are converted to PLINK format via ```plink2 --vcf variant.vcf --make-pgen --out variant```.
The linkage R2 values between every flanking SNP and SV is calculated via ```plink2 --pfile variant --r2 --out svid_ld --ld-window-r2 0 --ld-snp svid```. 
For phenotype testing, variants (SNPs and SVs) are tested separately using samples with valid genotypes (0/0, 1/0, 0/1 and 1/1) because SVs are not always successfully genotyped in the same set of samples. 
The covariant file contains 12 variables, including age, sex, and first ten principal components (PC1-PC10) derived from ancestry inference. 

Phecodes are derived from ICD-10 codes in participant electronic health records. We excluded phenotypes that have less than 20 cases. 
With this filtering, the number of tested phenotypes reduced from 1,846 to 1,621, covering 17 different major groups, such as infectious diseases, neoplasms, neurological, etc. 
The samples’ conditions are coded as ‘1’, ‘2’ or ‘NA’ for control, case and failure in ICD-10 codes conversion, respectively. 
Note that a case should have at least two occurrences of phecode while others are considered as control. 

#### Testing
To test the associations, we excluded related samples and used command ```plink2 --threads 16 --geno-counts --pfile variant --out variant --pheno samples_pheno.tsv --covar samples_covar.tsv --covar-variance-standardize --covar-name age,sex,pc1,pc2,pc3,..,pc10 --glm firth hide-covar cols=+nobs,+a1countcc,+gcountcc’```. 
The parameters ```cols=+nobs,+a1countcc,+gcountcc``` summarizes the SV genotypes in cases (CASE) and controls (CTRL). The genotype counts are in the output columns CASE_NON_A1_CT, CASE_HET_A1_CT, CASE_HOM_A1_CT, CTRL_NON_A1_CT, CTRL_HET_A1_CT and CTRL_HOM_A1_CT. 
To reduce false positives, we only reported associations in which at least five of the case samples carry the variant (i.e. #CASE_HET_A1_CT+#CASE_HOM_A1_CT>=5) and observed in at least 200,000 samples (OBS_CT>=200000).

The Manhattan plot is created by ```scripts/phewas/plot_phewas.py```, which takes the plink test results as input.

To view LD of SVs and nearby SNPs, we plot the LD decay (```scripts/phewas/SvLdDecay.py```) and variant associations (```scripts/phewas/SvLdAssoc.py```).

### Association fine-mapping

To conduct variant fine-mapping, we specifically assessed the SV and all SNPs within a 100kbp range with an allele count of at least 100 and a min case count of at least 5 for the associated phenotype. 
A genotype file and covariate file are created for further fine-mapping. The genotype file contains the SV genotypes from Locityper and DRAGEN genotypes of nearby filtered SNPs among 230 thousand samples. 
The covariate file includes the phenotype as one column and other columns use the same covariates mentioned above (i.e., age, sex and PC1-PC10). 

The R package [susieR](https://stephenslab.github.io/susieR/) (Zou et al. 2022) is used to identify leading associations. 
It takes the genotypes of the variants (i.e. SV and SNP) and covariates as input and tests against a specific phenotype. 
The WDL pipeline for variant fine-mapping is available at https://github.com/EichlerLab/AoU_WDL/tree/main/susie. 
The two major outputs contain the percentage of credible sets and posterior inclusion probability (PIP) for each variant. 
The value for credible set (CS) suggests the probability that at least one variant in this set is causal. 
The highest-PIP variant is likely to be the leading variant compared to other variants in the CS. 

Based on CS and PIP, we classified one SV-disease association into the following different situations:

- SV is the leading variant. 1) Only one CS; 2) PIP of SV is the highest among variants in the CS. 
- SV is an independent signal. 1) More than one CS; 2) SV is leading in one of the CS; 3) SV is not in LD with leading variants in other CS (R2<0.1).
