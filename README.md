# 1KG_LongRead_SV

This repository contains workflows and analysis scripts for creating a SV control set from 1725 1KG long-read genomes. 
These genomes are sequenced by HPRC ([Lucas et al. Biorxiv 2026](https://www.biorxiv.org/content/10.64898/2026.07.21.739710v1)), 
HGSVC ([Logsdon et al. Nature 2025](https://www.nature.com/articles/s41586-025-09140-6)), UW-ONT and IB-ONT. 
The UW-ONT sequenced a total of 500 genomes, 400 are novel and 100 genomes were published by [Gustafson et al. 2024](https://genome.cshlp.org/content/34/11/2061).
The IB-ONT genomes are published by Siegfried Schloissnig et.al. ([Nature 2025](https://www.nature.com/articles/s41586-025-09290-7)).
The SV control set has been used to identify potential pathogenic variants and new SV-disease associations in biobanks.

<div align=left><img width=20% height=20% src="https://github.com/jiadong324/1KG_LongRead_SV/blob/main/images/HGSVC_logo.png"/> <img width=30% height=30% src="https://github.com/jiadong324/1KG_LongRead_SV/blob/main/images/hprc_log.png"/><img width=40% height=40% src="https://github.com/jiadong324/1KG_LongRead_SV/blob/main/images/1000G-ONT.png"/></div> 


## Genomes

<div align=left><img width=50% height=50% src="https://github.com/jiadong324/1KG_LongRead_SV/blob/main/images/Fig1_Samples.png"/> </div>

The table below lists the source of 1,218 genomes used to create the callset in Lin et al. Note that the number below counts for the unique genomes.
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

Multiple caller merged callset for each genome. Please check the index file for download.

### Integrated callset

We provide both GRCh38 and T2T-CHM13 callsets (zenodo link) for 293 HPRC+HGSVC genomes and 1,218 genomes.

**CHM13_INSDEL_HGSVC_HPRC.vcf.gz:** Integrated SVs from HGSVC/HPRC genomes

**CHM13_INSDEL_1218.vcf.gz:** Integrated SVs from all dataset containing 1,218 genomes

**GRCh38_INSDEL_HGSVC_HPRC.vcf.gz:** Integrated SVs from HGSVC genomes

**GRCh38_INSDEL_1218.vcf.gz:** Integrated SVs from all dataset containing 1,218 genomes

## Genome assembly

### HPRC

The HPRC genomes are assembled with hifiasm and verkko. Please refer to HPRC [release](https://github.com/human-pangenomics/hprc_intermediate_assembly). 

### HGSVC

Assemblies used in Glennis Logsdon, et.al. were created by Verkko v1.1 and further phased with StrandSeq data using Graphasing. 
With latest Verkko (v2.2.1), we resolved previous genomes with large unassigned contigs using HiC data, including HG00096, HG00514, HG02587, HG03009, HG03683, NA18939, NA19434 and NA19650. 
The rest of the HGSVC genomes were not reassembled. 

### UW-ONT
Consistent with our previous publication (Gustafson et al.) on the first 100 genomes, the NAPU pipeline of the same version is used to build the assembly for R9 and R10 data. 
Briefly, this pipeline used shasta to make the haploid assembly and create the pseudo-haplotypes with hapdup. 
Yak (v0.1) and QUAST (v5.2.1) were used to assess the QV and contig N50 of the pseudo-phased assembly. The QV was calculated with Illumina short-read data.


## Per genome discovery

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

### Cross-caller SV integration

For each genome, we prioritize the PAV calling results and identify SVs supported by at least one another caller with [Truvari (v5.2.0)](https://github.com/acenglish/truvari). 
The pipeline is ```rules/intra_sample_collapse.smk```. The output of this pipeline is a multi-caller integrated VCF used the PAV reported breakpoint, sv length, phased genotype etc.

We used the same annotation as [Logsdon et al. Nature 2025](https://www.nature.com/articles/s41586-025-09140-6) to exclude SVs inside complex regions, gaps, etc. 
Briefly, these regions include UCSC gaps and centromere on GRCh38. 
For T2T-CHM13, complex regions include centromere, acrocentric p-arms, satellite regions except for monomeric satellite.

### TR genotyping

Tandem repeat catalogs for GRCh38 and T2T-CHM13 can be found [here](https://zenodo.org/records/13178746).


| Tool   | Input type | Version | Dataset             | Purpose        |
|--------|------------|---------|---------------------|----------------|
| TRGT   | HiFi       | v1.4.1  | HPRC, HGSVC, UW-ONT | Tandem repeats |
| vamos  | Assembly   | v2.1.5  | HPRC, HGSVC, UW-ONT | Tandem repeats |



