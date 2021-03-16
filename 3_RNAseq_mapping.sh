#### 3_RNAseq_mapping.sh

### Map RNA-seq reads to genome

## STRAT
# I will map with STAR - get read counts with HTseq

cd /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker
mkdir mapping_RNAseq_v8
mkdir mapping_RNAseq_v8/REFS
mkdir mapping_RNAseq_v8/gffs

### get ref genomes - put each in a sep dir for making STAR indexes in

mkdir mapping_RNAseq_v8/REFS/Tbi
mkdir mapping_RNAseq_v8/REFS/Tce
mkdir mapping_RNAseq_v8/REFS/Tcm
mkdir mapping_RNAseq_v8/REFS/Tpa
mkdir mapping_RNAseq_v8/REFS/Tps

cp Genomes/REFS/Tbi_b3v08.fasta  mapping_RNAseq_v8/REFS/Tbi
cp Genomes/REFS/Tce_b3v08.fasta  mapping_RNAseq_v8/REFS/Tce
cp Genomes/REFS/Tcm_b3v08.fasta  mapping_RNAseq_v8/REFS/Tcm
cp Genomes/REFS/Tpa_b3v08.fasta  mapping_RNAseq_v8/REFS/Tpa
cp Genomes/REFS/Tps_b3v08.fasta  mapping_RNAseq_v8/REFS/Tps

### get gffs

cp  Genomes/gffs/*.gff mapping_RNAseq_v8/gffs

## Use STAR

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/STAR/2.6.0c 

### index genome

#2.2.6 Genome with a large number of references.
#If you are using a genome with a large (>5,000) number of references (chrosomes/scaffolds), you may
#need to reduce the --genomeChrBinNbits to reduce RAM consumption. The following scaling is
#recomended: --genomeChrBinNbits = min(18, log2(GenomeLength/NumberOfReferences)). For
#example, for 3 gigaBase genome with 100,000 chromosomes/scaffolds, this is equal to 15.

# ALSO
# it will fail if not allocated enough RAM, but it will tell you how much it needs, which you can the specify with
#  --limitGenomeGenerateRAM 
# need to use --sjdbGTFtagExonParentTranscript Parent if using a gff rather than gtf 

STAR --runThreadN 12 --runMode genomeGenerate \
--genomeDir mapping_RNAseq_v8/REFS/Tbi \
--genomeFastaFiles mapping_RNAseq_v8/REFS/Tbi/Tbi_b3v08.fasta \
--sjdbGTFfile mapping_RNAseq_v8/gffs/Tbi_b3v08.max_arth_b2g_droso_b2g.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 --genomeChrBinNbits 15  --limitGenomeGenerateRAM 79000000000

STAR --runThreadN 12 --runMode genomeGenerate \
--genomeDir mapping_RNAseq_v8/REFS/Tce \
--genomeFastaFiles mapping_RNAseq_v8/REFS/Tce/Tce_b3v08.fasta \
--sjdbGTFfile mapping_RNAseq_v8/gffs/Tce_b3v08.max_arth_b2g_droso_b2g.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 --genomeChrBinNbits 15  --limitGenomeGenerateRAM 79000000000

STAR --runThreadN 12 --runMode genomeGenerate \
--genomeDir mapping_RNAseq_v8/REFS/Tcm \
--genomeFastaFiles mapping_RNAseq_v8/REFS/Tcm/Tcm_b3v08.fasta \
--sjdbGTFfile mapping_RNAseq_v8/gffs/Tcm_b3v08.max_arth_b2g_droso_b2g.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 --genomeChrBinNbits 15  --limitGenomeGenerateRAM 79000000000

STAR --runThreadN 12 --runMode genomeGenerate \
--genomeDir mapping_RNAseq_v8/REFS/Tpa \
--genomeFastaFiles mapping_RNAseq_v8/REFS/Tpa/Tpa_b3v08.fasta \
--sjdbGTFfile mapping_RNAseq_v8/gffs/Tpa_b3v08.max_arth_b2g_droso_b2g.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 --genomeChrBinNbits 15  --limitGenomeGenerateRAM 79000000000

STAR --runThreadN 12 --runMode genomeGenerate \
--genomeDir mapping_RNAseq_v8/REFS/Tps \
--genomeFastaFiles mapping_RNAseq_v8/REFS/Tps/Tps_b3v08.fasta \
--sjdbGTFfile mapping_RNAseq_v8/gffs/Tps_b3v08.max_arth_b2g_droso_b2g.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 --genomeChrBinNbits 15  --limitGenomeGenerateRAM 79000000000


########################################################################################################################
##### Download RNAseq Reads

########################################################################################################################
##### Trim RNAseq Reads

# ADD

#### map reads
#### MEMORY - mapping failed with 20GB - upped to 60GB

### Mapping to sp ref

mkdir mapping_RNAseq_v8/STAR_out
	
for i in /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/READS/RNA_seq/Tissue_specific/Qual_Trimmed_1a/*.fq.gz; do
        cd /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/mapping_RNAseq_v8/STAR_out/
		foo2_filename=`echo "$i"`
                foo3_base=`echo "$i" | sed 's/AandQtrimmed.fq.gz.*//' | sed 's/_Re1_/_Md_Re1/' | sed 's/_Re2_/_Md_Re2/' | sed 's/_Re3_/_Md_Re3/' | sed 's/.*\///' | sed 's/Tte_AF_GN_Md/Tte_AF_GN_Vi/' | sed 's/Tte_AF_HD_Md/Tte_AF_HD_Vi/' | sed 's/Tte_AF_LG_Md/Tte_AF_LG_Vi/' | sed 's/Tdi_AF_GN_Md/Tdi_AF_GN_Vi/' | sed 's/Tdi_AF_HD_Md/Tdi_AF_HD_Vi/' | sed 's/Tdi_AF_LG_Md/Tdi_AF_LG_Vi/' | sed 's/Tsi_AF_GN_Md/Tsi_AF_GN_Vi/' | sed 's/Tsi_AF_HD_Md/Tsi_AF_HD_Vi/' | sed 's/Tsi_AF_LG_Md/Tsi_AF_LG_Vi/' | sed 's/Tge_AF_GN_Md/Tge_AF_GN_Vi/' | sed 's/Tge_AF_HD_Md/Tge_AF_HD_Vi/' | sed 's/Tge_AF_LG_Md/Tge_AF_LG_Vi/' | sed 's/Tms_AF_GN_Md/Tms_AF_GN_Vi/' | sed 's/Tms_AF_HD_Md/Tms_AF_HD_Vi/' | sed 's/Tms_AF_LG_Md/Tms_AF_LG_Vi/'`
                sp=`echo "$i" | sed 's/.fq//' | sed 's/.*\///' | sed 's/_.*//'`
                echo $foo2_filename
                echo $foo3_base
                 echo $sp
		
		mkdir $foo3_base
		cd $foo3_base

		ref_name=`echo "/scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/mapping_RNAseq_v8/REFS/"$sp`
		echo $ref_name
		echo ""
		
		STAR --runMode alignReads --outSAMtype BAM Unsorted --readFilesCommand zcat \
		--genomeDir $ref_name \
		--outFileNamePrefix $foo3_base \
		--readFilesIn $foo2_filename	
done

### then I need to sort the bams for HTseq by name

mkdir mapping_RNAseq_v8/STAR_sorted_bams

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/samtools/1.3

for f in mapping_RNAseq_v8/STAR_out/*/*Aligned.out.bam; do
	sample_name=`echo $f | sed 's/.*\///' | sed 's/Aligned.out.bam//'`
	out_name=`echo "mapping_RNAseq_v8/STAR_sorted_bams/"$sample_name"_sorted.bam"`
	
	echo $f
	echo $sample_name
	echo $out_name

	samtools sort -n -o $out_name $f
done


### tidy up - keep final log files, remove STAR out bams

cd mapping_RNAseq_v8/
mkdir STAR_outlogs
mv STAR_out/*/*Log.final.out STAR_outlogs/
rm -r STAR_out/

### need to convert gff
# use Maker_gff_to_HTseq_gff.py

for i in Genomes/gffs/*.gff; do
	echo $i
	out_pre=`echo $i | sed 's/_droso_b2g.gff/_droso_b2g/'`
	echo $out_pre
	python3 ./accessory_scripts/Maker_gff_to_HTseq_gff.py -i $i -o $out_pre
done

## run HTseq-count
#  Part of the 'HTSeq' framework, version 0.9.1.

mkdir mapping_RNAseq_v8/HTseq_out

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/HTSeq/0.9.1

### exon (proper way) (reads in exons)

for f in mapping_RNAseq_v8/STAR_sorted_bams/*.bam ; do
        sample_name=`echo $f | sed 's/.*\///' | sed 's/_sorted.bam//'`
        out_name=`echo "mapping_RNAseq_v8/HTseq_out/"$sample_name".counts"`
        sp=`echo $f | sed 's/.*\///' | sed 's/_sorted.bam//' | sed 's/_.*//'`
        gff=`echo "Genomes/gffs/"$sp"_b3v08.max_arth_b2g_droso_b2g_forHTSeq.gff"`

        echo $f
        echo $sample_name
        echo $sp
        echo $out_name
        echo $gff       
        echo ""

	python2.7 -m HTSeq.scripts.count --order=name --type=exon --idattr=gene_id --stranded=reverse --format=bam $f $gff > $out_name
		
done


### clean up

rm -r STAR_sorted_bams/


##################################################################################################################################
###### stick read counts together

cd mapping_RNAseq_v8/

mkdir HTseq_out_ALL
mkdir HTseq_out_ALL/Tbi
mkdir HTseq_out_ALL/Tce
mkdir HTseq_out_ALL/Tcm
mkdir HTseq_out_ALL/Tpa
mkdir HTseq_out_ALL/Tps

cp HTseq_out/Tbi* HTseq_out_ALL/Tbi
cp HTseq_out/Tce* HTseq_out_ALL/Tce
cp HTseq_out/Tcm* HTseq_out_ALL/Tcm
cp HTseq_out/Tpa* HTseq_out_ALL/Tpa
cp HTseq_out/Tps* HTseq_out_ALL/Tps

cd data/output/HTseq_out_ALL

python3 ../../../accessory_scripts/HTSeq_to_edgeR.py -i Tbi -o Tbi_WBHDLGRT_v8
python3 ../../../accessory_scripts/HTSeq_to_edgeR.py -i Tce -o Tce_WBHDLGRT_v8
python3 ../../../accessory_scripts/HTSeq_to_edgeR.py -i Tcm -o Tcm_WBHDLGRT_v8
python3 ../../../accessory_scripts/HTSeq_to_edgeR.py -i Tpa -o Tpa_WBHDLGRT_v8
python3 ../../../accessory_scripts/HTSeq_to_edgeR.py -i Tps -o Tps_WBHDLGRT_v8


###################################################################################################################################
####### get total exon len by gene


mkdir data/output/gene_lens
cd    data/output/gene_lens

for f in ../../Genomes/gffs/*_b2g_forHTSeq.gff ; do
	sp=`echo $f | sed 's/.*\///' | sed 's/_.*//'`
	echo $sp
	python3 ../../../accessory_scripts/gff_feature_lengths.py -i $f -o $sp
done


	
###################################################################################################################################
### get rel position of genes in Nosil's linkage map

### select largest alignment block. most have one large block and many little ones.
### work out rel position of each gene.
# NOTE I match the midpoint of the largest alignment block to midpoint of the part of UNIL scaf it aligned to.

mkdir data/output/LG_pos

python3 accessory_scripts/class_genes_to_lg.py -l data/linkage_groups/Tbi_scf_block_alignment.tsv \
                                               -g Tbi_b3v08.max_arth_b2g_droso_b2g.gff \
                                               -o data/output/LG_pos/Tbi

python3 accessory_scripts/class_genes_to_lg.py -l data/linkage_groups/Tce_scf_block_alignment.tsv \
                                               -g Tce_b3v08.max_arth_b2g_droso_b2g.gff \
                                               -o data/output/LG_pos/Tce

python3 accessory_scripts/class_genes_to_lg.py -l data/linkage_groups/Tcm_scf_block_alignment.tsv \
                                               -g Tcm_b3v08.max_arth_b2g_droso_b2g.gff \
                                               -o data/output/LG_pos/Tcm

python3 accessory_scripts/class_genes_to_lg.py -l data/linkage_groups/Tpa_scf_block_alignment.tsv \
                                               -g Tpa_b3v08.max_arth_b2g_droso_b2g.gff \
                                               -o data/output/LG_pos/Tpa

python3 accessory_scripts/class_genes_to_lg.py -l data/linkage_groups/Tps_scf_block_alignment.tsv \
                                               -g Tps_b3v08.max_arth_b2g_droso_b2g.gff \
                                               -o data/output/LG_pos/Tps


###################################################################################################################################
### bring coverage, orthologs, and counts together
##### STORED output in data/counts for convenience

mkdir data/counts

for sp in "Tbi" "Tce" "Tcm" "Tpa" "Tps"; do
echo $sp
for scaf_len in 1000 5000; do
echo $scaf_len        

python3 accessory_scripts/sex_chr_cov_readcounts_tidier.py \
-c "data/output/MF_cov/"$sp"_MFcov_filt_"$scaf_len".csv"  \
-g "data/gffs/"$sp"_b3v08.max_arth_b2g_droso_b2g_forHTSeq.gff" \
-l "data/output/gene_lens/"$sp"_exon_by_gene_id_lengths.csv" \
-L "data/output/LG_pos/"$sp"_gene_rel_pos_lg.csv" \
-r "data/output/HTseq_out_ALL/"$sp"_WBHDLGRT_v8_H2E.counts.csv" \
-y data/TBITCETCMTPATPS_HOG_matrix.txt \
-o "data/counts/"$sp"_"$scaf_len

done
done



### Expression_analyses.R


cat FT_MB_FB_X_A_*FPKM* > FT_MB_FB_X_A_ALL_FPKM.csv
cat FT_MB_FB_X_A_*TPM* > FT_MB_FB_X_A_ALL_TPM.csv

