## Nosil_genome_read_mapping.sh

# get REF
#Using the 1.3 version (PRJNA417530)
#(https://drive.google.com/drive/folders/1Fb9nhKD5a_pE4vr-eFzW6YcVra10RqsX (the 1.3 version))

### prep ref

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.3
bwa index tcristinae_draft_1.3c2.fasta

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### map reads as paired reads with BWA
### removes multi + supp reads

cd     /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/
mkdir ./mapping_Nosil_Tce/BWA_out/
mkdir ./mapping_Nosil_Tce/BWA_out/mapped_as_paired
mkdir ./mapping_Nosil_Tce/BWA_out/flagstat_out_paired

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.3

read_dir="/scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/READS/trimmed_reads_by_RG"
ref_dir="/scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/mapping_Nosil_Tce/REFS"
map_out_dir="./mapping_Nosil_Tce/BWA_out/mapped_as_paired"
flag_out_dir="./mapping_Nosil_Tce/BWA_out/flagstat_out_paired"
mapqfilt="30"

for i in $read_dir/Tce_F*_R1_qtrimmed.fq.gz; do
	read_file_name_R1=`echo $i`
	read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
        base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
        base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
	badnames=`echo $base_read_name"_badnames.txt"`
        infile=`echo $i`
	sp=`echo $base_read_name | sed 's/_.*//'`
        ref_fa=`echo $ref_dir"/tcristinae_draft_1.3c2.fasta"`
        outsam=`echo $map_out_dir"/"$base_read_name"_to_"$sp"_Nosil_pe_BWA.sam"`
	outbam=`echo $map_out_dir"/"$base_read_name"_to_"$sp"_Nosil_pe_BWA_mapqfilt_"$mapqfilt"a.bam"`
	outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$sp"_Nosil_pe_BWA_mapqfilt_"$mapqfilt"a_sorted.bam"`
	flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$sp"_Nosil_pe_BWA_flagstat_out.txt"`
	flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$sp"_Nosil_pe_BWA_mapqfilt_"$mapqfilt"a_flagstat_out.txt"`
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
	samtools view -S  $outsam | fgrep -e 'XA:Z:' -e 'SA:Z:' | cut -f 1 > $badnames
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
##### merge bams by samp

cd /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker
mkdir /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged

samtools merge mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/Tce_F_CC22B_to_Tce_Nosil_pe_BWA_mapqfilt_30a_sorted.bam mapping_Nosil_Tce/BWA_out/mapped_as_paired/Tce_F_CC22B*_mapqfilt_30a_sorted.bam
samtools merge mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/Tce_F_CC22C_to_Tce_Nosil_pe_BWA_mapqfilt_30a_sorted.bam mapping_Nosil_Tce/BWA_out/mapped_as_paired/Tce_F_CC22C*_mapqfilt_30a_sorted.bam
samtools merge mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/Tce_F_CC24B_to_Tce_Nosil_pe_BWA_mapqfilt_30a_sorted.bam mapping_Nosil_Tce/BWA_out/mapped_as_paired/Tce_F_CC24B*_mapqfilt_30a_sorted.bam
samtools merge mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/Tce_F_CC24C_to_Tce_Nosil_pe_BWA_mapqfilt_30a_sorted.bam mapping_Nosil_Tce/BWA_out/mapped_as_paired/Tce_F_CC24C*_mapqfilt_30a_sorted.bam
samtools merge mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/Tce_F_CC25B_to_Tce_Nosil_pe_BWA_mapqfilt_30a_sorted.bam mapping_Nosil_Tce/BWA_out/mapped_as_paired/Tce_F_CC25B*_mapqfilt_30a_sorted.bam

samtools merge mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/Tce_M_05_HM15_to_Tce_Nosil_pe_BWA_mapqfilt_30a_sorted.bam mapping_Nosil_Tce/BWA_out/mapped_as_paired/Tce_M_05_HM15*_mapqfilt_30a_sorted.bam
samtools merge mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/Tce_M_06_HM16_to_Tce_Nosil_pe_BWA_mapqfilt_30a_sorted.bam mapping_Nosil_Tce/BWA_out/mapped_as_paired/Tce_M_06_HM16*_mapqfilt_30a_sorted.bam
samtools merge mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/Tce_M_07_HM33_to_Tce_Nosil_pe_BWA_mapqfilt_30a_sorted.bam mapping_Nosil_Tce/BWA_out/mapped_as_paired/Tce_M_07_HM33*_mapqfilt_30a_sorted.bam
samtools merge mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/Tce_M_08_HM61_to_Tce_Nosil_pe_BWA_mapqfilt_30a_sorted.bam mapping_Nosil_Tce/BWA_out/mapped_as_paired/Tce_M_08_HM61*_mapqfilt_30a_sorted.bam



#########################################################################################################################
##### remove PCR duplicates  (Not just mark, as I don't think angsD pays attention to this) (actually it should, but I also need to calc cov.)
###  60GB RAM 

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/picard-tools/2.9.0 
module load UHTS/Analysis/samtools/1.3

for i in  ./mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/*_mapqfilt_30a_sorted.bam; do
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

rm ./mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/*_mapqfilt_30a_sorted.bam 
rm ./mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/*_metric.txt


############################################################################################
###### calc cov per scaf

cd    /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.3
module load UHTS/Analysis/BEDTools/2.26.0

full_ref_dir="/scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/mapping_Nosil_Tce/REFS"

for i in ./mapping_Nosil_Tce/BWA_out/mapped_as_paired_merged/*_sorted.bam; do
    sp=`echo $i | sed 's/.*\///' | sed 's/_.*//'`
	inname=`echo $i`
    basename=`echo $i | sed 's/_sorted.bam//'`
    out_file=`echo $basename"_coverage.out"`
    ref_fa=`echo $ref_dir"/tcristinae_draft_1.3c2.fasta"`
    echo $sp
	echo $inname
    echo $basename
    echo $out_file
    echo $ref_fa
	echo ""
   
    genomeCoverageBed -ibam $inname -g $ref_fa > $out_file
	
done

####
mkdir ./mapping_Nosil_Tce/mappingcoverage_ests_BWA_out_merged_perscaf_aDR
mv mapping_Nosil_Tce/BWA_out/*_merged*/*_coverage.out ./mapping_Nosil_Tce/mappingcoverage_ests_BWA_out_merged_perscaf_aDR

### Calc cov PER SCAF for coverage analyses 

for i in ./mapping_Nosil_Tce/mappingcoverage_ests_BWA_out_merged_perscaf_aDR/*_coverage.out; do
	in_name=`echo $i`
    basename=`echo $i | sed 's/_coverage.out//'`	
	echo $in_name
    echo $basename
	python ~/Gen_BioInf/genomeCoverageBed_tidier.py $in_name max $basename
done




## produces Tce_Nosil_30aDR_cov.csv
## HERE: data/Nosil_mapping/Tce_Nosil_30aDR_cov.csv


## then Male_and_female_coverage_Nosil.R


# This was used as the basis to correct the linkage group assignments by removing X-linked scaffolds from autosomal linkage groups 1-12.
# These scaffolds were then assigned to a new, unordered collection of scaffolds from the X chromosome, together with X-linked scaffolds from linkage group 13 and from those not assigned to any linkage group in Nosil et al. (2018).
# Scaffolds from linkage group 13 that were not X-linked were assigned to linkage group NA.
# We classified scaffolds in the Nosil et al. (2018) assembly as X-linked if most of a scaffold was covered by aligned scaffolds from our assembly assigned to the X rather than autosomes
# (i.e. if aligned scaffolds assigned as X covered more than twice as many bases as those assigned as autosomal for a particular scaffold, it was classed as X-linked).


# this output provided in data/linkage_groups.tar.gz for convenience
