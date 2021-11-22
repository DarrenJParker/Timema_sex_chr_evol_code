# Timema_sex_chr_evol

* Need to explain all accessory scripts

## how we get to genes
select largest alignment block. most have one large block and many little ones.
If there was > 1 linkage group I make a note
work out rel position of each gene.
NOTE I match the midpoint of the largest alignment block to midpoint of the part of UNIL scaf it aligned to.


## ID X linked scaffolds

* **1_read_prep.sh**
* **2_genome_read_mapping.sh**
* **Male_and_female_coverage.R**
* **Orthologs_on_the_X.R** | looking if X linked genes are the same in each species

## Gene expression analysis

* **3_RNAseq_mapping.sh** | need to add trimming + download reads
* **Expression_analyses.R** | Dosage compensation and sex-biased genes

## Update linkage groups in the Nosil et al genome
* **2a_Nosil_genome_read_mapping.sh** | map reads 
* **Nosil_Male_and_female_coverage.R**  | coverage analyses

## blast to Bacillus
* **Bacillus_X_blast.sh** | Blast Timema X orthologs to Bacillus genome.

## data

* **Nosil_mapping/Tce_Nosil_30aDR_cov.csv** | coverage for each scaffold from Nosil et al (2018).
* **linkage_groups.tar.gz** | linkage groups for each species inferred from Nosil et al (2018), see 2a_Nosil_genome_read_mapping.sh.
* **Sample_coverage** | per base coverage for each sample.
* **coverage_and_heterozygosity** | Contig coverage and heterozygosity estimates for each sample.
* **TBITCETCMTPATPS_HOG_matrix.txt** | 1-to-1 orthologs for the five species. Orthologs were obtained using OrthoDB standalone pipeline (v. 2.4.4) (Kriventseva et al, 2015) as described in Jaron et al. (2020).
* **counts** | RNAseq read counts with sex chromosome information


## References

Jaron, K. S*., Parker, D. J*., Anselmetti, Y., Tran Van, P. T., Bast, J., Dumas,  Z., Figuet, E., François, C. M., Hayward, K., Rossier, V., Simion, P., Robinson-Rechavi,  M., Galtier, N., Schwander, T. 2020. Convergent consequences of parthenogenesis on stick insect genomes. bioRxiv. doi: https:// doi.org/10.1101/2020.11.20.391540

Kriventseva, E. V., Tegenfeldt, F., Petty, T. J., Waterhouse, R. M., Simão, F. A., Pozdnyakov, I. A., Ioannidis, P., & Zdobnov, E. M. 2015. OrthoDB v8: update of the hierarchical catalog of orthologs and the underlying free software. Nucleic Acids Research, 43, D250–D256.

P. Nosil, R. Villoutreix, C. F. de Carvalho, T. E. Farkas, V. Soria-Carrasco, J. L. Feder, B. J. Crespi, Z. Gompert, Natural selection and the predictability of evolution in Timema stick insects. Science. 359, 765–770 (2018).
