# 1_read_prep.sh

########################################################################################################################################
## download reads 
################################################################################################################################################

## NCBI Female genome reads bioproject accession number: PRJNA670663; Male genome reads bioproject accession number: PRJNA725673


########## Remove reads that failed Casava from files

#### explination - so reads that failed basic tests at the seq centre were not removed
##### but insteads just flagged with a 'Y' for failed and an 'N' for passed: e.g.

# @WINDU:16:C5UFBANXX:7:1310:11475:72909 1:N:0:GAGTGG ## passed
# @WINDU:16:C5UFBANXX:7:1310:11694:72820 1:Y:0:GAGTGG ## failed

# also seq company breaks reads up into parts of a file so they need cat-ing

### used something like this to get reads that passed - note also need to remove weird -- lines that appear from grepping
#zcat LIB01.fastq | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" > LIB01-filtered.fastq 2> LIB01-filtered.err 

### note if a read is flagged as failed it is failed in both R1 and R2:

#@HISEQ:262:CA3BFANXX:1:1101:1828:2088 1:Y:0:ATCACG
#@HISEQ:262:CA3BFANXX:1:1101:1828:2088 2:Y:0:ATCACG
#
#@HISEQ:262:CA3BFANXX:1:1101:2300:2213 1:Y:0:ATCACG
#@HISEQ:262:CA3BFANXX:1:1101:2300:2213 2:Y:0:ATCACG
#
#@HISEQ:262:CA3BFANXX:1:1101:2696:2087 1:Y:0:ATCACG
#@HISEQ:262:CA3BFANXX:1:1101:2696:2087 2:Y:0:ATCACG


### first let's rename the files to make life easier

mkdir temp_reads_for_cc

for r in 1 2 3 4 5 6 7 8 9 10 11; do
	for file in "raw_reads/round_"$r"/"*.gz; do
		echo $file
		i=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/_L._R._.*//' `
		lane_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/_R._.*//' | sed 's/.*_//'`
		R_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._//' | sed 's/_.*//'`
		file_number=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R._//' | sed 's/_.*//' `

		echo $i
		echo $lane_N	
		echo $R_N
		echo $file_number

		sp=""
        if [ $i = "01_HM148" ]; then
                sp="Tcm_M"
        elif [ $i = "02_HM149" ]; then
               sp="Tcm_M"
        elif [ $i = "03_HM150" ]; then
               sp="Tcm_M"
        elif [ $i = "04_HM151" ]; then
               sp="Tcm_M"

        elif [ $i = "13_Tbi" ]; then
               sp="Tbi_M"
        elif [ $i = "14_Tbi" ]; then
               sp="Tbi_M"
        elif [ $i = "15_Tbi" ]; then
               sp="Tbi_M"
        elif [ $i = "16_Tbi" ]; then
               sp="Tbi_M"

        elif [ $i = "05_HM15" ]; then
               sp="Tce_M"
        elif [ $i = "06_HM16" ]; then
               sp="Tce_M"
        elif [ $i = "07_HM33" ]; then
               sp="Tce_M"
        elif [ $i = "08_HM61" ]; then
               sp="Tce_M"
                                      		
        elif [ $i = "21_15232" ]; then
               sp="Tch_M"
        elif [ $i = "22_15233" ]; then
               sp="Tch_M"
        elif [ $i = "23_15234" ]; then
               sp="Tch_M"
        elif [ $i = "24_15235" ]; then
               sp="Tch_M"

        elif [ $i = "10091" ]; then
               sp="Tge_M"

        elif [ $i = "25_15055" ]; then
               sp="Tms_M"
        elif [ $i = "26_15056" ]; then
               sp="Tms_M"
        elif [ $i = "27_15057" ]; then
               sp="Tms_M"
        elif [ $i = "34_16002" ]; then
               sp="Tms_M"
	
        elif [ $i = "09_Tpa" ]; then
               sp="Tpa_M"
        elif [ $i = "10_Tpa" ]; then
               sp="Tpa_M"		   
        elif [ $i = "11_Tpa" ]; then
               sp="Tpa_M"		   
        elif [ $i = "12_Tpa" ]; then
               sp="Tpa_M"		   

        elif [ $i = "17_HM99" ]; then
               sp="Tps_M"	
        elif [ $i = "18_HM100" ]; then
               sp="Tps_M"
        elif [ $i = "19_HM101" ]; then
               sp="Tps_M"
        elif [ $i = "20_15255" ]; then
               sp="Tps_M"	

        elif [ $i = "15244" ]; then
               sp="Tch_F"	
        elif [ $i = "15245" ]; then
               sp="Tch_F"
        elif [ $i = "15246" ]; then
               sp="Tch_F"
        elif [ $i = "15247" ]; then
               sp="Tch_F"	   
	   
        else
                sp="nope"
        fi
		
	new_name=`echo $sp"_"$i"_"$R_N"_G"$r$lane_N"_"$file_number".fq.gz"` 
	echo $new_name
	
	cp $file "temp_reads_for_cc/"$new_name

	done
done		

#### cat and clean by read group 
mkdir cat_clean_by_read_group

for i in temp_reads_for_cc/*_R1_* ; do
	foo_R1=`echo $i`
	foo_R2=`echo $i | sed 's/_R1_/_R2_/g'`
	base_out_name=`echo $i | sed 's/.*\///' | sed 's/R1.*//' `
	read_group=`echo $i | sed 's/.*\///' | sed 's/.*R1_//' | sed 's/_.*//'`	
	
	out_R1=`echo "cat_clean_by_read_group/"$base_out_name$read_group"_R1_CC.fq"  `
	out_R2=`echo "cat_clean_by_read_group/"$base_out_name$read_group"_R2_CC.fq"  `
	echo $foo_R1
	echo $foo_R2
	echo $base_out_name
	echo $read_group
	echo $out_R1
	echo $out_R2
	echo ""
	zcat $foo_R1  | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> $out_R1
	zcat $foo_R2  | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> $out_R2
done

### Zip
sampleList=( cat_clean_by_read_group/*.fq )
echo "${sampleList[@]}"
sampleName=${sampleList[(($SLURM_ARRAY_TASK_ID-1))]}
echo "### Sample for the current run is: $sampleName"
gzip $sampleName

### tidyup 

rm -r temp_reads_for_cc/


### females
mkdir temp_reads_for_cc_reseq_females

for file in ./raw_reads_reseq_females/*.fastq.gz; do

	echo $file
	i=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/_L._R._.*//' `
	lane_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/_R._.*//' | sed 's/.*_//'`
	R_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R*//' | sed 's/_.*//'`	
	file_number=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R._//' | sed 's/_.*//' `
	Group_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R._..._//' | sed 's/_//g'  `
	echo $i
	echo $lane_N	
	echo $R_N
	echo $file_number
	echo $Group_N

	sp=""
	if [ $i = "CC86B" ]; then
			sp="Tbi_F"
	elif [ $i = "CC86C" ]; then
		   sp="Tbi_F"
	elif [ $i = "CC87B" ]; then
		   sp="Tbi_F"
	elif [ $i = "CC87C" ]; then
		   sp="Tbi_F"
	elif [ $i = "CC88B" ]; then
		   sp="Tbi_F"

	elif [ $i = "CC22B" ]; then
		   sp="Tce_F"
	elif [ $i = "CC22C" ]; then
		   sp="Tce_F"
	elif [ $i = "CC24B" ]; then
		   sp="Tce_F"
	elif [ $i = "CC24C" ]; then
		   sp="Tce_F"
	elif [ $i = "CC25B" ]; then
		   sp="Tce_F"

	elif [ $i = "HM217" ]; then
		   sp="Tcm_F"
	elif [ $i = "HM218" ]; then
		   sp="Tcm_F"
	elif [ $i = "HM219" ]; then
		   sp="Tcm_F"
	elif [ $i = "HM220" ]; then
		   sp="Tcm_F"
	elif [ $i = "HM221" ]; then
		   sp="Tcm_F"

	elif [ $i = "Pa_AB" ]; then
		   sp="Tpa_F"
	elif [ $i = "PA_CD" ]; then
		   sp="Tpa_F"
	elif [ $i = "PA_E" ]; then
		   sp="Tpa_F"
	elif [ $i = "H54" ]; then
		   sp="Tpa_F"
	elif [ $i = "H56" ]; then
		   sp="Tpa_F"

	elif [ $i = "ReSeq_Ms01" ]; then
		   sp="Tms_F"
	elif [ $i = "ReSeq_Ms02" ]; then
		   sp="Tms_F"
	elif [ $i = "ReSeq_Ms03" ]; then
		   sp="Tms_F"
	elif [ $i = "MS_Alp03b" ]; then
		   sp="Tms_F"
	elif [ $i = "MS_Alp04b" ]; then
		   sp="Tms_F"

	elif [ $i = "ReSeq_Ps14" ]; then
		   sp="Tps_F"
	elif [ $i = "ReSeq_Ps16" ]; then
		   sp="Tps_F"
	elif [ $i = "ReSeq_Ps18" ]; then
		   sp="Tps_F"
	elif [ $i = "ReSeq_Ps08" ]; then
		   sp="Tps_F"
	elif [ $i = "ReSeq_Ps12" ]; then
		   sp="Tps_F"

	elif [ $i = "CC59_A" ]; then
		   sp="Tge_F"
	elif [ $i = "CC59_C" ]; then
		   sp="Tge_F"
	elif [ $i = "CC65_B" ]; then
		   sp="Tge_F"
	elif [ $i = "CC66_A" ]; then
		   sp="Tge_F"
	elif [ $i = "CC67_A" ]; then
		   sp="Tge_F"

	elif [ $i = "ReSeq_Di02" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di04" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di06" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di08" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di10" ]; then
		   sp="Tdi_F"

	elif [ $i = "ReSeq_S14" ]; then
		   sp="Tsi_F"
	elif [ $i = "ReSeq_Si01" ]; then
		   sp="Tsi_F"
	elif [ $i = "ReSeq_Si03" ]; then
		   sp="Tsi_F"
	elif [ $i = "ReSeq_Si16" ]; then
		   sp="Tsi_F"
	elif [ $i = "ReSeq_Si18" ]; then
		   sp="Tsi_F"

	elif [ $i = "ReSeq_Te07" ]; then
		   sp="Tte_F"
	elif [ $i = "ReSeq_Te08" ]; then
		   sp="Tte_F"
	elif [ $i = "ReSeq_Te09" ]; then
		   sp="Tte_F"
	elif [ $i = "ReSeq_Te10" ]; then
		   sp="Tte_F"
	elif [ $i = "ReSeq_Te11" ]; then
		   sp="Tte_F"
		   
	else
			foo5="nope"
	fi

	new_name=`echo $sp"_"$i"_R"$R_N"_G"$Group_N$lane_N"_"$file_number".fq.gz"` 
	echo $new_name
	echo ""
	cp $file "temp_reads_for_cc_reseq_females/"$new_name
done


#### cat and clean by read group 
mkdir cat_clean_by_read_group

for i in temp_reads_for_cc_reseq_females/*_R1_* ; do
	foo_R1=`echo $i`
	foo_R2=`echo $i | sed 's/_R1_/_R2_/g'`
	base_out_name=`echo $i | sed 's/.*\///' | sed 's/R1.*//' `
	read_group=`echo $i | sed 's/.*\///' | sed 's/.*R1_//' | sed 's/_.*//'`	
	
	out_R1=`echo "cat_clean_by_read_group/"$base_out_name$read_group"_R1_CC.fq"  `
	out_R2=`echo "cat_clean_by_read_group/"$base_out_name$read_group"_R2_CC.fq"  `
	echo $foo_R1
	echo $foo_R2
	echo $base_out_name
	echo $read_group
	echo $out_R1
	echo $out_R2
	echo ""
	zcat $foo_R1  | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> $out_R1
	zcat $foo_R2  | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v "^--$" >> $out_R2
done


### zip

sampleList=( cat_clean_by_read_group/*.fq )
echo "${sampleList[@]}"
sampleName=${sampleList[(($SLURM_ARRAY_TASK_ID-1))]}
echo "### Sample for the current run is: $sampleName"

gzip $sampleName
rm -r temp_reads_for_cc_reseq_females 


####################################
### add all cat clean reads into
# cat_clean_by_read_group


#####################################################################################################
#####################################################################################################
### read trimming

mkdir trimmed_reads_by_RG
cd trimmed_reads_by_RG

cp data/AllIllumina-PEadapters.fa .

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/trimmomatic/0.36

for i in READS/cat_clean_by_read_group/*_R1_CC.fq.gz ; do
        foo1=`echo $i`
		basename=`echo $foo1 | sed 's/_R1_CC.fq.gz.*//' | sed 's/.*\///'`
        infileR1=`echo $foo1`
        infileR2=`echo $foo1 | sed 's/_R1_CC.fq.gz/_R2_CC.fq.gz/'`
        outfileR1=`echo "./"$basename"_R1_qtrimmed.fq"`
        outfileR2=`echo "./"$basename"_R2_qtrimmed.fq"`
        outfileR1_UP=`echo "./"$basename"_R1_qtrimmed_UNPAIRED.fq"`
        outfileR2_UP=`echo "./"$basename"_R2_qtrimmed_UNPAIRED.fq"`
        
        echo $infileR1
        echo $infileR2
        echo $outfileR1
        echo $outfileR1_UP
        echo $outfileR2
        echo $outfileR2_UP
		
        trimmomatic PE -threads 30 $infileR1 $infileR2 $outfileR1 $outfileR1_UP $outfileR2 $outfileR2_UP ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:90
done

