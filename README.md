# Timema_sex_chr_evol_code


This is the repository for the collected scripts used in the study:

>Parker, D. J., Jaron, K. S., Dumas, Z., Robinson-Rechavi, M., Schwander, T. 2021. X chromosomes show relaxed selection and complete somatic dosage compensation across _Timema_ stick insect species. bioRxiv. https://doi.org/10.1101/2021.11.28.470265

Raw sequence reads generated in this project have been deposited in NCBI’s sequence read archive under the bioproject number PRJNA725673. 


## ID X linked scaffolds

* **1_read_prep.sh** | Female genome reads bioproject accession number: PRJNA670663;  Male genome reads bioproject accession number: PRJNA725673
* **2_genome_read_mapping.sh** | mapping to version 8 genomes (bioproject accession number: PRJEB31411, fasta and gff also here: https://doi.org/10.5281/zenodo.5636226).
* **Male_and_female_coverage.R** | Coverage analyses, heterozygosity.
* **Male_and_female_coverage_nucla.R** | Nucleotide diversity analyses, Ne.
* **Orthologs_on_the_X.R** | looking if X linked genes are the same in each species

## Gene expression analysis

* **3_RNAseq_mapping.sh** | RNAseq read trimming and mapping | RNA-seq reads Bioproject accession number: PRJNA392384
* **Expression_analyses.R** | Dosage compensation and sex-biased genes
* **Expression_analyses_keep_sex_specific_genes.R** | Dosage compensation and sex-biased genes when sex-limited genes are retained

## Selection analyses

* **4_selection_analyses.sh** | Calc GC, look at selection on the X 
* **GC_CUB.R** | Calc GC
* **GC_genomic.R** | plot GC for contigs
* **selection_on_the_X.R** | looking at selection on the X


## Update linkage groups in the Nosil et al genome
* **Nosil_genome_read_mapping.sh** | map reads 
* **Nosil_Male_and_female_coverage.R**  | coverage analyses

## Blast to Bacillus
* **Bacillus_X_blast.sh** | Blast Timema X orthologs to Bacillus genome.

## data/

* **Nosil_mapping** | coverage for each scaffold from Nosil et al (2018).
* **linkage_groups.tar.gz** | linkage groups for each species inferred from Nosil et al (2018), see 2a_Nosil_genome_read_mapping.sh.
* **Sample_coverage** | per base coverage for each sample.
* **coverage_and_heterozygosity** | Contig coverage and heterozygosity estimates for each sample.
* **TBITCETCMTPATPS_HOG_matrix.txt** | 1-to-1 orthologs for the five species. Orthologs were obtained using OrthoDB standalone pipeline (v. 2.4.4) (Kriventseva et al, 2015) as described in Jaron et al. (2021).
* **selection** | Sequence evolution data from Jaron et al. (2021).
* **counts** | RNAseq read counts with sex chromosome information
* **AllIllumina-PEadapters.fa** | Illumina adapter sequence
* **Timema_X_soft_HOG_longest_aa.fa** | amino acid sequence for _Timema_ X orthologs 
* **Selectome_v07_Timema-nt_masked.zip** | Masked alignments from selectome (See [Selectome](https://selectome.org/timema) and Jaron et al. (2021) for more details)


## References

Jaron, K. S*., Parker, D. J*., Anselmetti, Y., Tran Van, P. T., Bast, J., Dumas,  Z., Figuet, E., François, C. M., Hayward, K., Rossier, V., Simion, P., Robinson-Rechavi,  M., Galtier, N., Schwander, T. 2021. Convergent consequences of parthenogenesis on stick insect genomes. bioRxiv. doi: https:// doi.org/10.1101/2020.11.20.391540

Kriventseva, E. V., Tegenfeldt, F., Petty, T. J., Waterhouse, R. M., Simão, F. A., Pozdnyakov, I. A., Ioannidis, P., & Zdobnov, E. M. 2015. OrthoDB v8: update of the hierarchical catalog of orthologs and the underlying free software. Nucleic Acids Research, 43, D250–D256.

Nosil, P., Villoutreix, R., de Carvalho, C. F., Farkas, T. E., Soria-Carrasco, V., Feder, J. L., Crespi, B. J., Gompert, Z. Natural selection and the predictability of evolution in _Timema_ stick insects. Science. 359, 765–770 (2018).
