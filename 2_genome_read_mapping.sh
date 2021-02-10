## 2_read_mapping_v8_genomes.sh

#####################################################################################################

### get references

cd /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/Genomes
mkdir Genomes
cd Genomes/
wget ftp://ftpmrr.unil.ch/AsexGenomeEvol/timema/final_references/*
gzip -d *

mkdir REFS
mkdir gffs
mv *gff gffs/
mv *fasta REFS/

cd REFS
mv 1_Tdi_b3v08.fasta Tdi_b3v08.fasta
mv 1_Tps_b3v08.fasta Tps_b3v08.fasta
mv 2_Tcm_b3v08.fasta Tcm_b3v08.fasta
mv 2_Tsi_b3v08.fasta Tsi_b3v08.fasta
mv 3_Tce_b3v08.fasta Tce_b3v08.fasta
mv 3_Tms_b3v08.fasta Tms_b3v08.fasta
mv 4_Tbi_b3v08.fasta Tbi_b3v08.fasta
mv 4_Tte_b3v08.fasta Tte_b3v08.fasta
mv 5_Tge_b3v08.fasta Tge_b3v08.fasta
mv 5_Tpa_b3v08.fasta Tpa_b3v08.fasta

cd ../gffs/
mv 1_Tdi_b3v08.max_arth_b2g_droso_b2g.gff Tdi_b3v08.max_arth_b2g_droso_b2g.gff
mv 1_Tps_b3v08.max_arth_b2g_droso_b2g.gff Tps_b3v08.max_arth_b2g_droso_b2g.gff
mv 2_Tcm_b3v08.max_arth_b2g_droso_b2g.gff Tcm_b3v08.max_arth_b2g_droso_b2g.gff
mv 2_Tsi_b3v08.max_arth_b2g_droso_b2g.gff Tsi_b3v08.max_arth_b2g_droso_b2g.gff
mv 3_Tce_b3v08.max_arth_b2g_droso_b2g.gff Tce_b3v08.max_arth_b2g_droso_b2g.gff
mv 3_Tms_b3v08.max_arth_b2g_droso_b2g.gff Tms_b3v08.max_arth_b2g_droso_b2g.gff
mv 4_Tbi_b3v08.max_arth_b2g_droso_b2g.gff Tbi_b3v08.max_arth_b2g_droso_b2g.gff
mv 4_Tte_b3v08.max_arth_b2g_droso_b2g.gff Tte_b3v08.max_arth_b2g_droso_b2g.gff
mv 5_Tge_b3v08.max_arth_b2g_droso_b2g.gff Tge_b3v08.max_arth_b2g_droso_b2g.gff
mv 5_Tpa_b3v08.max_arth_b2g_droso_b2g.gff Tpa_b3v08.max_arth_b2g_droso_b2g.gff



#####################################################################################################################
### LINKAGE GROUPS
# from Jaron, K. S*., Parker, D. J*., Anselmetti, Y., Tran Van, P. T., Bast, J., Dumas,  Z., Figuet, E., François, C. M., Hayward, K., Rossier, V., Simion, P., Robinson-Rechavi,  M., Galtier, N., Schwander, T. 2020. Convergent consequences of parthenogenesis on stick insect genomes. bioRxiv. doi: https:// doi.org/10.1101/2020.11.20.391540
# provided in data/linkage_groups.tar.gz for convenience

# unzip before use

tar -zxf data/linkage_groups.tar.gz


#####################################################################################################
## mapping with BWA

cd /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker
mkdir mapping_v8

### prep refs

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.3

for i in ./Genomes/REFS/*fasta; do
	bwa index $i
done



#########################################################################################################################################################
#### Mapping strategy

## I want bams for:
#			coverage analysis
#			est heterozygosity

### Mapping as paired-end
#### remove duplicate reads, filter by mapq, and make sorted bam
### BWA MEM does not give an XT:A:U tag in MEM, but XA tags are still there.
### HOWEVER XA tags are only there if there are 1-3 alternative alignments.
### However, this should not be a problem since multiple alignments would result in a low mapping_v8 quality (and can thus be removed on this basis).
# A simple grep would work, but that could leave me with broken flags (filter out one read in a pair and not the other without adjusting the flag).
# SO need to ID multi reads the remove them (I do not need to do this for single end)
# note ONLY WORKS with this samtools 1.3 or higher

# I also filter out Supplemental alignments. These are chimeric alignments
# A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
# So it is a bit of a multimapper. Some of these supplemental alignments have very high mapQ values (60)
# They have an SA:Z tag. If these are paired reads I also remove the other pair, as I do with multi-mappers


# I also add read groups
# to ADD read groups - easiest at BWA stage - but can add / change with picard if I need to
#ID = Read group identifier This tag identifies which read group each read belongs to, so each read group's ID must be unique.
#PU = Platform Unit The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}. 
#SM = Sample The name of the sample sequenced in this read group. 
#PL = Platform/technology used to produce the read. Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO.
#LB = DNA preparation library identifier 

# for BWA
# -R "@RG\tID:S1L6\tSM:S1\tPL:ILLUMINA\tLB:FC-140-1086"
# in my case lib = sample

# I then merge bams together by sample and remove PCR duplicates
 # With Picard
 
### I then do coverage analysis on these BAMs

### calc heterozygosity
 ## Additional step of indel realignment (Gatk)
 ## then use angsD on these bams to get heterozyg: http://popgen.dk/angsd/index.php/Heterozygosity




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### map reads as paired reads with BWA

mkdir ./mapping_v8/BWA_out/
mkdir ./mapping_v8/BWA_out/mapped_as_paired
mkdir ./mapping_v8/BWA_out/flagstat_out_paired

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.3

read_dir="/scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/READS/trimmed_reads_by_RG"
ref_dir="/scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/Genomes/REFS"
map_out_dir="./mapping_v8/BWA_out/mapped_as_paired"
flag_out_dir="./mapping_v8/BWA_out/flagstat_out_paired"
mapqfilt="30"

for i in $read_dir/*_R1_qtrimmed.fq.gz; do
		read_file_name_R1=`echo $i`
		read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
        base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
        base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
		badnames=`echo $base_read_name"_badnames.txt"`
        infile=`echo $i`
		sp=`echo $base_read_name | sed 's/_.*//'`
        ref_fa=`echo $ref_dir"/"$sp"_b3v08.fasta"`
        outsam=`echo $map_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA.sam"`
		outbam=`echo $map_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA_mapqfilt_"$mapqfilt".bam"`
		outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA_mapqfilt_"$mapqfilt"_sorted.bam"`
		flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA_flagstat_out.txt"`
		flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$sp"_v8_pe_BWA_mapqfilt_"$mapqfilt"_flagstat_out.txt"`
		IFS='_' read -r -a sp_want_list <<< "$base_read_name"
		readgroup=`echo ${sp_want_list[-1]}`
		readgroup_txt=`echo "@RG\tID:"$readgroup"\tSM:"$base_read_name3"\tPL:ILLUMINA\tLB:"$base_read_name3`
		
		
		echo $read_file_name_R1
		echo $read_file_name_R2
        echo $base_read_name
		echo $base_read_name3
        echo $ref_fa
        echo $outsam
		echo $outbam
		echo $outbam_sorted
		echo $flagstat_out_sam
		echo $flagstat_out_bam
		echo $readgroup
		echo $readgroup_txt
		echo $badnames
		echo ""

		## map
        bwa mem -t 16 -R $readgroup_txt $ref_fa $read_file_name_R1 $read_file_name_R2 > $outsam
		
		#flagstat
		samtools flagstat $outsam > $flagstat_out_sam
		
		# filter ## filter both reads out to avoid broken flags
		samtools view -S  $outsam | fgrep XA | cut -f 1 > $badnames
		samtools view -Sh $outsam | fgrep -vf $badnames | samtools view -Shbq $mapqfilt - > $outbam
		
		# sort bam
		samtools sort $outbam > $outbam_sorted
		
		#flagstat
		samtools flagstat $outbam_sorted > $flagstat_out_bam
		
		#tidyup
		rm $outsam
		rm $outbam
		
done

rm *_badnames.txt


#########################################################################################################################
##### remove supp reads

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.3


for i in ./mapping_v8/BWA_out/mapped_as_paired/*_mapqfilt_30_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30_sorted.bam/_mapqfilt_30a_sorted.bam/'`
    badnames=`echo $i"_badnames.txt"`
    flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`

	echo $i
	echo $outbam
	echo $flagstat_out_bam
	
    samtools view $i | fgrep SA:Z: | cut -f 1 > $badnames
    samtools view -h $i | fgrep -vf $badnames | samtools view -b > $outbam
    samtools flagstat $outbam > $flagstat_out_bam
done

## remove original bams and _badnames.txt

mv ./mapping_v8/BWA_out/mapped_as_single/*_flagstat_out.txt  mapping_v8/BWA_out/flagstat_out_single/
rm ./mapping_v8/BWA_out/mapped_as_single/*_mapqfilt_30_sorted.bam

mv ./mapping_v8/BWA_out/mapped_as_paired/*_flagstat_out.txt  mapping_v8/BWA_out/flagstat_out_paired/
rm ./mapping_v8/BWA_out/mapped_as_paired/*badnames.txt
rm ./mapping_v8/BWA_out/mapped_as_paired/*_mapqfilt_30_sorted.bam


#########################################################################################################################
##### merge bams

mkdir /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/mapping_v8/BWA_out/mapped_as_paired_merged

samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tbi_F_CC86B_to_Tbi_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tbi_F_CC86B*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tbi_F_CC86C_to_Tbi_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tbi_F_CC86C*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tbi_F_CC87B_to_Tbi_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tbi_F_CC87B*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tbi_F_CC87C_to_Tbi_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tbi_F_CC87C*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tbi_F_CC88B_to_Tbi_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tbi_F_CC88B*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tbi_M_13_Tbi_to_Tbi_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tbi_M_13_Tbi*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tbi_M_14_Tbi_to_Tbi_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tbi_M_14_Tbi*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tbi_M_15_Tbi_to_Tbi_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tbi_M_15_Tbi*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tbi_M_16_Tbi_to_Tbi_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tbi_M_16_Tbi*_mapqfilt_30a_sorted.bam

samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tce_F_CC22B_to_Tce_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tce_F_CC22B*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tce_F_CC22C_to_Tce_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tce_F_CC22C*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tce_F_CC24B_to_Tce_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tce_F_CC24B*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tce_F_CC24C_to_Tce_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tce_F_CC24C*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tce_F_CC25B_to_Tce_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tce_F_CC25B*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tce_M_05_HM15_to_Tce_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tce_M_05_HM15*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tce_M_06_HM16_to_Tce_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tce_M_06_HM16*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tce_M_07_HM33_to_Tce_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tce_M_07_HM33*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tce_M_08_HM61_to_Tce_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tce_M_08_HM61*_mapqfilt_30a_sorted.bam

samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tcm_F_HM217_to_Tcm_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tcm_F_HM217*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tcm_F_HM218_to_Tcm_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tcm_F_HM218*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tcm_F_HM219_to_Tcm_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tcm_F_HM219*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tcm_F_HM220_to_Tcm_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tcm_F_HM220*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tcm_F_HM221_to_Tcm_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tcm_F_HM221*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tcm_M_01_HM148_to_Tcm_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tcm_M_01_HM148*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tcm_M_02_HM149_to_Tcm_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tcm_M_02_HM149*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tcm_M_03_HM150_to_Tcm_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tcm_M_03_HM150*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tcm_M_04_HM151_to_Tcm_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tcm_M_04_HM151*_mapqfilt_30a_sorted.bam

samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tpa_F_H54_to_Tpa_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tpa_F_H54*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tpa_F_H56_to_Tpa_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tpa_F_H56*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tpa_F_Pa_AB_to_Tpa_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tpa_F_Pa_AB*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tpa_F_PA_CD_to_Tpa_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tpa_F_PA_CD*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tpa_F_PA_E_to_Tpa_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tpa_F_PA_E*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tpa_M_09_Tpa_to_Tpa_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tpa_M_09_Tpa*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tpa_M_10_Tpa_to_Tpa_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tpa_M_10_Tpa*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tpa_M_11_Tpa_to_Tpa_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tpa_M_11_Tpa*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tpa_M_12_Tpa_to_Tpa_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tpa_M_12_Tpa*_mapqfilt_30a_sorted.bam

samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tps_F_ReSeq_Ps08_to_Tps_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tps_F_ReSeq_Ps08*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tps_F_ReSeq_Ps12_to_Tps_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tps_F_ReSeq_Ps12*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tps_F_ReSeq_Ps14_to_Tps_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tps_F_ReSeq_Ps14*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tps_F_ReSeq_Ps16_to_Tps_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tps_F_ReSeq_Ps16*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tps_F_ReSeq_Ps18_to_Tps_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tps_F_ReSeq_Ps18*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tps_M_17_HM99_to_Tps_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tps_M_17_HM99*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tps_M_18_HM100_to_Tps_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tps_M_18_HM100*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tps_M_19_HM101_to_Tps_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tps_M_19_HM101*_mapqfilt_30a_sorted.bam
samtools merge mapping_v8/BWA_out/mapped_as_paired_merged/Tps_M_20_15255_to_Tps_v8_pe_BWA_mapqfilt_30a_sorted.bam mapping_v8/BWA_out/mapped_as_paired/Tps_M_20_15255*_mapqfilt_30a_sorted.bam



#########################################################################################################################
##### remove PCR duplicates  (Not just mark, as I don't think angsD pays attention to this) (actually it should, but I also need to calc cov.)
### most samples do this in 20GB RAM - 2 Tpa samples failed - stuck up to 60GB - worked

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/picard-tools/2.9.0 
module load UHTS/Analysis/samtools/1.3

for i in  ./mapping_v8/BWA_out/mapped_as_paired_merged/*_mapqfilt_30a_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30a_sorted.bam/_mapqfilt_30aDR_sorted.bam/'`
	flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
	metric_file=`echo $outbam | sed 's/.bam/_metric.txt/'`

	echo $i
	echo $outbam
	echo $metric_file
	echo $flagstat_out_bam

	picard-tools MarkDuplicates REMOVE_DUPLICATES=true \
	INPUT=$i \
    OUTPUT=$outbam \
    METRICS_FILE=$metric_file
	
	samtools flagstat $outbam > $flagstat_out_bam

done

## tidy - remove original bams

mv ./mapping_v8/BWA_out/mapped_as_paired_merged/*_flagstat_out.txt  mapping_v8/BWA_out/flagstat_out_paired_merged/
rm ./mapping_v8/BWA_out/mapped_as_paired_merged/*_mapqfilt_30a_sorted.bam 
rm ./mapping_v8/BWA_out/mapped_as_paired_merged/*_metric.txt


############################################################################################
###### calc cov per scaf

module add  Bioinformatics/Software/vital-it
module load UHTS/Analysis/BEDTools/2.26.0

cd    /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker

module add  Bioinformatics/Software/vital-it
module load UHTS/Analysis/BEDTools/2.26.0
full_ref_dir="/scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/Genomes/REFS"

for i in ./mapping_v8/BWA_out/mapped_as_paired_merged/*_sorted.bam; do
    sp=`echo $i | sed 's/.*\///' | sed 's/_.*//'`
	inname=`echo $i`
    basename=`echo $i | sed 's/_sorted.bam//'`
    out_file=`echo $basename"_coverage.out"`
    ref_fa=`echo $full_ref_dir"/"$sp"_b3v08.fasta"`
    echo $sp
	echo $inname
    echo $basename
    echo $out_file
    echo $ref_fa
	echo ""
    genomeCoverageBed -ibam $inname -g $ref_fa > $out_file
	
done

mkdir ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR
mv mapping_v8/BWA_out/*_merged*/*_coverage.out ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR

### plot coverage - all sites

for i in ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/*_coverage.out; do
	echo $i
	python3 ~/Gen_BioInf/genomeCoverageBed_tidier_wholegenomecov.py -i $i
done

mkdir ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/for_plotting
mv ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/*_genomecov.txt ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/for_plotting

### plot coverage - for all sites on scafs >= 1000bp

for i in mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/*_coverage.out; do
	sp=`echo $i | sed 's/.*\///' | sed 's/_.*//'`
	want_file=`echo "Genomes/REFS/"$sp"_b3v08_1000.names"`
	echo $i
	echo $want_file
	python3 ~/Gen_BioInf/genomeCoverageBed_tidier_select_scafs.py -i $i -w $want_file -e _1000_
done

mv ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/*_coverage_1000_cov.txt ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/for_plotting

### the output from the code above is provided in data/Sample_coverage for convenience

cd /Users/dparker/Documents/University/Lausanne/Sex_chromosomes/mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/for_plotting 

for i in ./*_genomecov.txt; do
	echo $i
	Rscript ~/Documents/Gen_BioInf/plot_genome_cov.R $i
done

for i in ./*_coverage_1000_cov.txt; do
	echo $i
	Rscript ~/Documents/Gen_BioInf/plot_genome_cov.R $i
done

### this gives cov ests for angsD cutoffs.


#######################################################################################################################
#######################################################################################################################
## ALSO Calc cov PER SCAF for coverage analyses 

mkdir ./mapping_v8/mappingcoverage_ests_BWA_out_aDR
cp ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/*_coverage.out ./mapping_v8/mappingcoverage_ests_BWA_out_aDR ## moving to sep folder for ease

for i in ./mapping_v8/mappingcoverage_ests_BWA_out_aDR/*_coverage.out; do
	in_name=`echo $i`
    basename=`echo $i | sed 's/_coverage.out//'`	
	echo $in_name
    echo $basename
	python ~/Gen_BioInf/genomeCoverageBed_tidier.py $in_name max $basename
done

    
#######################################################################################################################
## join together by sp with contig lengths

cd ./mapping_v8/mappingcoverage_ests_BWA_out_aDR

mkdir Tbi_cont_cov
mkdir Tce_cont_cov
mkdir Tcm_cont_cov
mkdir Tpa_cont_cov
mkdir Tps_cont_cov
mkdir Tms_cont_cov
mkdir Tge_cont_cov
mkdir Tdi_cont_cov
mkdir Tsi_cont_cov
mkdir Tte_cont_cov

mv Tbi_*_contig_cov.txt Tbi_cont_cov
mv Tce_*_contig_cov.txt Tce_cont_cov
mv Tcm_*_contig_cov.txt Tcm_cont_cov
mv Tpa_*_contig_cov.txt Tpa_cont_cov
mv Tps_*_contig_cov.txt Tps_cont_cov
mv Tms_*_contig_cov.txt Tms_cont_cov
mv Tge_*_contig_cov.txt Tge_cont_cov
mv Tdi_*_contig_cov.txt Tdi_cont_cov
mv Tsi_*_contig_cov.txt Tsi_cont_cov
mv Tte_*_contig_cov.txt Tte_cont_cov

cd /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/

## Filter down to smallest contig being 1000 bp 
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tbi_cont_cov -f Genomes/REFS/Tbi_b3v08.fasta -m 1000 -o Tbi -e 30aDR_contig_cov.txt -P
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tce_cont_cov -f Genomes/REFS/Tce_b3v08.fasta -m 1000 -o Tce -e 30aDR_contig_cov.txt -P
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tcm_cont_cov -f Genomes/REFS/Tcm_b3v08.fasta -m 1000 -o Tcm -e 30aDR_contig_cov.txt -P
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tpa_cont_cov -f Genomes/REFS/Tpa_b3v08.fasta -m 1000 -o Tpa -e 30aDR_contig_cov.txt -P
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tps_cont_cov -f Genomes/REFS/Tps_b3v08.fasta -m 1000 -o Tps -e 30aDR_contig_cov.txt -P
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tge_cont_cov -f Genomes/REFS/Tge_b3v08.fasta -m 1000 -o Tge -e 30aDR_contig_cov.txt -P
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tms_cont_cov -f Genomes/REFS/Tms_b3v08.fasta -m 1000 -o Tms -e 30aDR_contig_cov.txt -P
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tdi_cont_cov -f Genomes/REFS/Tdi_b3v08.fasta -m 1000 -o Tdi -e 30aDR_contig_cov.txt -P
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tsi_cont_cov -f Genomes/REFS/Tsi_b3v08.fasta -m 1000 -o Tsi -e 30aDR_contig_cov.txt -P
python3 sex_chr_scripts/Timema_cov_tidier.py -i mapping_v8/mappingcoverage_ests_BWA_out_aDR/Tte_cont_cov -f Genomes/REFS/Tte_b3v08.fasta -m 1000 -o Tte -e 30aDR_contig_cov.txt -P

mkdir mapping_v8/v8aDR_cov_contig
mv *_contig_cov.txt mapping_v8/v8aDR_cov_contig

### Add linkage groups

for i in ./mapping_v8/v8aDR_cov_contig/*_contig_cov.txt; do
	sp=`echo $i | sed 's/.*\///' | sed 's/_.*//'`
	basename=`echo $i | sed 's/.txt//'`
	link_file=`echo "./linkage_map/"$sp"_scf_block_alignment.tsv"`
	
	echo $i
    echo $basename
    echo $link_file
	python3 ./sex_chr_scripts/class_scafs_to_lg.py -c $i -l $link_file -o $basename"_KF2"
done


##########################################################################################################################################################
##########################################################################################################################################################
### calc heterozygosity


###############################################################################################################################
### First indel realignment ##

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/GenomeAnalysisTK/3.7
module load UHTS/Analysis/picard-tools/2.9.0
module load UHTS/Analysis/samtools/1.3

## index and dict of fasta

for f in Genomes/REFS/*_b3v08.fasta ; do 
    dict_name=`echo $f | sed 's/.fasta/.dict/'`
    #samtools faidx $f
    picard-tools CreateSequenceDictionary R= $f O= $dict_name
done

# indel realignment 

for i in  ./mapping_v8/BWA_out/mapped_as_paired_merged/Tte_F_ReSeq_Te11*_mapqfilt_30aDR_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30aDR_sorted.bam/_mapqfilt_30aDRra_sorted.bam/'`
	flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
	interval_file=`echo $outbam | sed 's/.bam/.intervals/'`
	sp=`echo $i | sed 's/.*\///' | sed 's/_.*//'`
	ref_fa=`echo "Genomes/REFS/"$sp"_b3v08.fasta"`
	
	echo $i
	echo $outbam
	echo $flagstat_out_bam
	echo $interval_file
	echo $sp
	echo $ref_fa
	echo ""

	# index bam
	samtools index $i
	
	## make target intervals list
	GenomeAnalysisTK -T RealignerTargetCreator -R $ref_fa -I $i -o $interval_file
	
	## realign
	GenomeAnalysisTK -T IndelRealigner -R $ref_fa --targetIntervals $interval_file -I $i -o $outbam
	samtools flagstat $outbam > $flagstat_out_bam	

done

rm ./mapping_v8/BWA_out/mapped_as_paired_merged/*.intervals



#################################################################################
## ANGSD by scaffold  http://popgen.dk/angsd/index.php/Heterozygosity
## Note it does not use readgroups! One bam = one sample
# v8aDRra_angsD_out = with dup reads removed + supp reads removed + PCR dups removed + indels realigned 


module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/samtools/1.3

cd /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/

# running with 20 GB RAM - one sample at a time 
module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/ANGSD/0.921
module load UHTS/Analysis/samtools/1.3
module load UHTS/Aligner/bwa/0.7.15

# mkdir -p mapping_v8/v8_angsD_out/mapped_as_single
# mkdir -p mapping_v8/v8_angsD_out/mapped_as_paired

version="v8aDRra"
out_dir=`echo "mapping_v8/"$version"_angsD_out/mapped_as_paired"`
mkdir -p $out_dir

### index BAMs

for i in ./mapping_v8/BWA_out/mapped_as_paired_merged/*_sorted.bam; do
	in_name=`echo $i`
	echo $in_name
	samtools index $in_name
done

# first
# Need to pick a sensible read MAX depth. Maybe calc mean cov and go 2x higher?

## calculated with plot_genome_cov.R (above)
## WILL use 2x med cov after 0s filtered out.

cd /Users/dparker/Documents/University/Lausanne/Sex_chromosomes/mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/for_plotting
cat *_BWA_mapqfilt_30aDR_coverage_1000_cov.txtcovest.csv | grep -v 'sample' > All_v8_BWA_mapqfilt_30aDR_coverage_1000_cov.txtcovest.csv

## add to here (axiom)

mkdir -p mapping_v8/v8_angsD_out/
mkdir -p mapping_v8/v8aDR_angsD_out/
## for v8aDRra - just use same file as v8aDR
cp mapping_v8/v8aDR_angsD_out/All_v8_BWA_mapqfilt_30aDR_coverage_1000_cov.txtcovest.csv mapping_v8/v8aDRra_angsD_out/
cd mapping_v8/v8aDRra_angsD_outv8aDRra_angsD_out/
mv All_v8_BWA_mapqfilt_30aDR_coverage_1000_cov.txtcovest.csv All_v8_BWA_mapqfilt_30aDRra_coverage_1000_cov.txtcovest.csv
sed -i 's/30aDR/30aDRra/g' All_v8_BWA_mapqfilt_30aDRra_coverage_1000_cov.txtcovest.csv
 
### scaf by scaf for scafs >= 1000 bp

### first get scafs I want.

for i in ./Genomes/REFS/*fasta; do
	out_fa=`echo $i | sed 's/.fasta/_1000.fasta/'`
	python3 ~/fasta_tools/fasta_select_by_len.py $i max 1000 $out_fa
done

for i in ./Genomes/REFS/*_1000.fasta; do
	out_file=`echo $i | sed 's/.fasta/.names/'`
	grep ">" $i > $out_file
done

#### index fastas 

for i in ./Genomes/REFS/*_1000.fasta; do
	samtools faidx $i
done


### run angsD

version_30=`echo $version | sed 's/v8/30/'`

for i in "./mapping_v8/BWA_out/mapped_as_paired_merged/"*$version_30"_sorted.bam" ; do
	
sp=`echo $i | sed 's/.*\///' | sed  's/_.*//' `
wanted_seq_name=`echo "Genomes/REFS/"$sp"_b3v08_1000.names"`
genome_name=`echo "Genomes/REFS/"$sp"_b3v08_1000.fasta"`
out_file_name_a=`echo $i | sed 's/.*\///' | sed 's/_sorted.bam//'`
	
echo $i
echo $wanted_seq_name
echo $genome_name

echo -e "scaf\tN_homo\tN_hetero" > $out_dir"/"$out_file_name_a"_est.ml"	

export i
export version
	
max_cov=`python3 << END
import os
curr_v = os.environ['version']
curr_v_30 = curr_v.replace("v8", "30")
cov_file = open("mapping_v8/" + curr_v + "_angsD_out/All_v8_BWA_mapqfilt_" + curr_v_30 + "_coverage_1000_cov.txtcovest.csv")
cov_dict = {}
        
for line in cov_file:
    line = line.strip().replace('"','').split(",")
    cov_filename = line[0].split("_BWA_mapqfilt_30")[0]
    double_med_cov_after_filt = line[4]
    #print(line)
    #print(cov_filename)
    #print(double_med_cov_after_filt)
    cov_dict[cov_filename] = double_med_cov_after_filt
            
filenames = os.environ['i']
#print(filenames)
sample_name = filenames.strip("/").split("/")[-1].split("_BWA_mapqfilt_30")[0]
#print(sample_name)

curr_samp_cov = cov_dict.get(sample_name)
print(curr_samp_cov )
        
END`

echo $max_cov

while read line; do
		
	want_seq=`echo $line | sed 's/>//'`
	echo $want_seq
	want_bam=`echo $i`
	bam_name=` echo $want_bam | sed 's/.*\///' | sed 's/_sorted.bam//' `
	
	echo $bam_name
	# angsd
	angsd -i $want_bam \
	-anc $genome_name -P 8 -r $want_seq \
	-doSaf 1 -gl 1 -minQ 20 -minMapQ 40 -fold 1 -doCounts 1 -setMinDepth 5 -setMaxDepth $max_cov \
	-out $out_dir"/"$bam_name"_"$want_seq
		
	# get hetero
	realSFS $out_dir"/"$bam_name"_"$want_seq".saf.idx" > $out_dir"/"$bam_name"_"$want_seq"_est.ml"	
	
	wait
	
	rm $out_dir"/"$bam_name"_"$want_seq".arg"
	rm $out_dir"/"$bam_name"_"$want_seq".saf.gz"
	rm $out_dir"/"$bam_name"_"$want_seq".saf.idx"
	rm $out_dir"/"$bam_name"_"$want_seq".saf.pos.gz"

	echo $want_seq > $out_dir"/"$bam_name"_"$want_seq"_scaf.temp"
	paste $out_dir"/"$bam_name"_"$want_seq"_scaf.temp" $out_dir"/"$bam_name"_"$want_seq"_est.ml" > $out_dir"/"$bam_name"_"$want_seq"_est.ml2"

	cat $out_dir"/"$bam_name"_"$want_seq"_est.ml2" >> $out_dir"/"$out_file_name_a"_est.ml"	
	
	rm $out_dir"/"$bam_name"_"$want_seq"_est.ml"
	rm $out_dir"/"$bam_name"_"$want_seq"_est.ml2"
	rm $out_dir"/"$bam_name"_"$want_seq"_scaf.temp"
	
done <$wanted_seq_name

echo ""

done

############################################################################################
# join together with cov and linkage info


for sp in Tbi Tce Tcm Tpa Tps; do 
	echo $sp
	python3 add_hetero_info.py \
	-i /Users/dparker/Documents/University/Lausanne/Sex_chromosomes/v8aDRra_angsD_out/mapped_as_paired/ \
	-s $sp \
	-c "/Users/dparker/Documents/University/Lausanne/Sex_chromosomes/v8aDR_cov_contig/"$sp"_pairedcov_minlen=1000_contig_cov_KF2_wLGinfopos.csv"
done

#Samples:
#['Tbi_F_CC86B', 'Tbi_F_CC86C', 'Tbi_F_CC87B', 'Tbi_F_CC87C', 'Tbi_F_CC88B', 'Tbi_M_13_Tbi', 'Tbi_M_14_Tbi', 'Tbi_M_15_Tbi', 'Tbi_M_16_Tbi']
#
#Total number of scafs: 53814
#Total number where ansgD did not calc ests for ANY sites for at least 1 sample: 3348
#Number of scafs with missing angsD ests (i.e. in cov file but not angsD file): 0

#Samples:
#['Tce_F_CC22B', 'Tce_F_CC22C', 'Tce_F_CC24B', 'Tce_F_CC24C', 'Tce_F_CC25B', 'Tce_M_05_HM15', 'Tce_M_06_HM16', 'Tce_M_07_HM33', 'Tce_M_08_HM61']
#
#Total number of scafs: 47657
#Total number where ansgD did not calc ests for ANY sites for at least 1 sample: 6610
#Number of scafs with missing angsD ests (i.e. in cov file but not angsD file): 0

#Samples:
#['Tcm_F_HM217', 'Tcm_F_HM218', 'Tcm_F_HM219', 'Tcm_F_HM220', 'Tcm_F_HM221', 'Tcm_M_01_HM148', 'Tcm_M_02_HM149', 'Tcm_M_03_HM150', 'Tcm_M_04_HM151']
#
#Total number of scafs: 71329
#Total number where ansgD did not calc ests for ANY sites for at least 1 sample: 10371
#Number of scafs with missing angsD ests (i.e. in cov file but not angsD file): 0

#Samples:
#['Tpa_F_H54', 'Tpa_F_H56', 'Tpa_F_PA_CD', 'Tpa_F_PA_E', 'Tpa_F_Pa_AB', 'Tpa_M_09_Tpa', 'Tpa_M_10_Tpa', 'Tpa_M_11_Tpa', 'Tpa_M_12_Tpa']
#
#Total number of scafs: 155674
#Total number where ansgD did not calc ests for ANY sites for at least 1 sample: 19451
#Number of scafs with missing angsD ests (i.e. in cov file but not angsD file): 0

#Samples:
#['Tps_F_ReSeq_Ps08', 'Tps_F_ReSeq_Ps12', 'Tps_F_ReSeq_Ps14', 'Tps_F_ReSeq_Ps16', 'Tps_F_ReSeq_Ps18', 'Tps_M_17_HM99', 'Tps_M_18_HM100', 'Tps_M_19_HM101', 'Tps_M_20_15255']
#
#Total number of scafs: 90454
#Total number where ansgD did not calc ests for ANY sites for at least 1 sample: 12117
#Number of scafs with missing angsD ests (i.e. in cov file but not angsD file): 0



################################################################################################################
#### OUTPUT is stored here: data/coverage_and_heterozygosity 


