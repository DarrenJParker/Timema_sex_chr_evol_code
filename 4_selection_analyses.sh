# 4_selection_analyses.sh


## get GC and CUB
### get sequences from alignments (data/selection/Selectome_v07_Timema-nt_masked.zip)
## Unzip first

### remove masked codons (nnn) and gaps (-) and make seperate file for each 

python3 accessory_scripts/alignments_to_seq.py -i data/selection/Selectome_v07_Timema-nt_masked/ -o data/output/sel_out/Ortho

### get genomes

# version 8 genomes (bioproject accession number: also here: https://doi.org/10.5281/zenodo.5636226).
# download and add here:data/genomes
mkdir data/genomes
cd data/genomes
wget https://zenodo.org/record/5636226/files/Tbi_b3v08.fasta
wget https://zenodo.org/record/5636226/files/Tce_b3v08.fasta
wget https://zenodo.org/record/5636226/files/Tcm_b3v08.fasta
wget https://zenodo.org/record/5636226/files/Tpa_b3v08.fasta
wget https://zenodo.org/record/5636226/files/Tps_b3v08.fasta
cd ../../

### remove contigs < 1000 BP
for f in data/genomes/*fasta; do
out_f=`echo $f | sed 's/.fasta/_1000.fasta/'`
echo $f
echo $out_f
python3 accessory_scripts/fasta_select_by_len.py $f max 1000 $out_f
done

rm data/genomes/*_b3v08.fasta


### get GC for genomes with GC_CUB.R


####################################
## calc GC for all contigs with GC_CUB.R

## add to chr info

python3 accessory_scripts/add_GC_contigs.py -g data/output/sel_out/GC_g.csv -c data/output/MF_cov/Tbi_MFcov_filt_1000.csv
python3 accessory_scripts/add_GC_contigs.py -g data/output/sel_out/GC_g.csv -c data/output/MF_cov/Tce_MFcov_filt_1000.csv
python3 accessory_scripts/add_GC_contigs.py -g data/output/sel_out/GC_g.csv -c data/output/MF_cov/Tcm_MFcov_filt_1000.csv
python3 accessory_scripts/add_GC_contigs.py -g data/output/sel_out/GC_g.csv -c data/output/MF_cov/Tpa_MFcov_filt_1000.csv
python3 accessory_scripts/add_GC_contigs.py -g data/output/sel_out/GC_g.csv -c data/output/MF_cov/Tps_MFcov_filt_1000.csv 


### then GC_genomic.R to plot


### Then selection_on_the_X.R to look at selection on the X




#########################
### M8a

mkdir data/output/M8a
## need just sexual species - keeping only HOGs with >= 4 sexual species
### get sequences from alignments (data/selection/Selectome_v07_Timema-nt_masked.zip)
## Unzip first and add to data/output/M8a

cd data/output/M8a

for f in Selectome_v07_Timema-nt_masked/*fas; do
python3  ../../../accessory_scripts/filter_alignment.py -f $f
done

mkdir Selectome_v07_Timema-nt_masked_sexualspeciesonly
mv Selectome_v07_Timema-nt_masked/*_sexonly.fas Selectome_v07_Timema-nt_masked_sexualspeciesonly/

ls -l Selectome_v07_Timema-nt_masked/*fas | wc -l
#7157
ls -l Selectome_v07_Timema-nt_masked_sexualspeciesonly/*fas | wc -l
#6296










