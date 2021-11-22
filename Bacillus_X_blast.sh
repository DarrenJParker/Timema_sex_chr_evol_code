### blast to bacillus.

## genome Using the Bacillus rossius genom from T. schwander

## Longest isoforms of Timema X orthologs
## data/Timema_X_soft_HOG_longest.fa


### blast

makeblastdb -in Brsri_v3.fasta -out Brsri_v3_db -dbtype nucl
tblastn -query X_soft_HOG_longest_aa.fa -db Brsri_v3_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out X_soft_HOG_longest_aa_to_Brsri_v3_blast_out.csv
