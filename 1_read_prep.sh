# 1_read_prep.sh

########################################################################################################################################
## download reads 
################################################################################################################################################


####################################
### add all cat clean reads into
# /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/READS/cat_clean_by_read_group

### count reads

cd /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/READS
~/DJP_shell_sc/count_reads_in_fq.gz_files.sh cat_clean_by_read_group/ read_counts/SChr_reads_rounds_1-11_cat_clean_by_read_group_all_asex.txt


#####################################################################################################
#####################################################################################################
### read trimming

cd /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/READS
mkdir trimmed_reads_by_RG
cd trimmed_reads_by_RG

cp ~/Gen_BioInf/Trimmomatic_adapters/AllIllumina-PEadapters.fa .

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/trimmomatic/0.36

for i in /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/READS/cat_clean_by_read_group/*_R1_CC.fq.gz ; do
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




