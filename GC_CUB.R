# GC_CUB

### libs

library(ggplot2)
library(stringr)
library(modeest)
library(plyr)
library(Rtsne)
library(cowplot)
library(edgeR)
library(hash)
library(pheatmap)
library(seqinr)
library(coRdon)

#### get data
## genes
Tbi_seqs <- read.fasta(file="data/output/sel_out/Ortho_Tbi.fa")
Tbi_seq_names <- as.list(getName(Tbi_seqs ))
Tce_seqs <- read.fasta(file="data/output/sel_out/Ortho_Tce.fa")
Tce_seq_names <- as.list(getName(Tce_seqs ))
Tcm_seqs <- read.fasta(file="data/output/sel_out/Ortho_Tcm.fa")
Tcm_seq_names <- as.list(getName(Tcm_seqs ))
Tpa_seqs <- read.fasta(file="data/output/sel_out/Ortho_Tpa.fa")
Tpa_seq_names <- as.list(getName(Tpa_seqs ))
Tps_seqs <- read.fasta(file="data/output/sel_out/Ortho_Tps.fa")
Tps_seq_names <- as.list(getName(Tps_seqs ))

#### calc GC

get_GC <- function(sp){
  GC_all  <- c()
  GC_1_all <- c()
  GC_2_all <- c()
  GC_3_all <- c()
  
  N_genes = 0
  g_name = c()
  for(s in eval(parse(text=paste(sp, '_seq_names',sep='')))){
    #print(s)
    N_genes = N_genes + 1
    g_name <- c(g_name, s)
    GC_n   = GC(eval(parse(text=paste(sp, '_seqs','$',s,sep=''))), forceToLower = TRUE, exact = TRUE)
    GC_1_n = GC1(eval(parse(text=paste(sp, '_seqs','$',s,sep=''))), forceToLower = TRUE, exact = TRUE)
    GC_2_n = GC2(eval(parse(text=paste(sp, '_seqs','$',s,sep=''))), forceToLower = TRUE, exact = TRUE)
    GC_3_n = GC3(eval(parse(text=paste(sp, '_seqs','$',s,sep=''))), forceToLower = TRUE, exact = TRUE)
  
    #print(A1)
  
    GC_all <- c(GC_all, GC_n)
    GC_1_all <- c(GC_1_all, GC_1_n)
    GC_2_all <- c(GC_2_all, GC_2_n)
    GC_3_all <- c(GC_3_all, GC_3_n)
  }
  print(N_genes)
  out_df_GC <- cbind(g_name, GC_all, GC_1_all, GC_2_all, GC_3_all)
  return(out_df_GC)
}

Tbi_GC <- get_GC("Tbi")
Tce_GC <- get_GC("Tce")
Tcm_GC <- get_GC("Tcm")
Tpa_GC <- get_GC("Tpa")
Tps_GC <- get_GC("Tps")

## output
write.csv(as.data.frame(rbind(Tbi_GC, Tce_GC, Tcm_GC, Tpa_GC, Tps_GC)), "data/output/sel_out/GC.csv", row.names = F, quote = F)



## contigs ## just the contigs >= 1000 bp 


#### calc genome GC

get_g_GC <- function(sp){
  GC_all  <- c()
  N_genes = 0
  g_name = c()
  for(s in eval(parse(text=paste(sp, '_g_seq_names',sep='')))){
    #print(s)
    N_genes = N_genes + 1
    g_name <- c(g_name, s)
    GC_n   = GC(eval(parse(text=paste(sp, '_g_seqs','$"',s ,'"',sep=''))), forceToLower = TRUE, exact = TRUE)
    
    #print(A1)
    
    GC_all <- c(GC_all, GC_n)
  }
  print(N_genes)
  out_df_GC <- cbind(g_name, GC_all)
  return(out_df_GC)
}

Tbi_g_seqs <- read.fasta(file="data/genomes/Tbi_b3v08_1000.fasta")
Tbi_g_seq_names <- as.list(getName(Tbi_g_seqs ))
Tbi_g_seqs$"4_Tbi_b3v08_scaf000001"
Tbi_g_GC <- get_g_GC("Tbi")
rm(Tbi_g_seqs)
rm(Tbi_g_seq_names)
write.csv(as.data.frame(rbind(Tbi_g_GC)), "data/output/sel_out/Tbi_g_GC.csv", row.names = F, quote = F)



Tce_g_seqs <- read.fasta(file="data/genomes/Tce_b3v08_1000.fasta")
Tce_g_seq_names <- as.list(getName(Tce_g_seqs ))
Tce_g_GC <- get_g_GC("Tce")
rm(Tce_g_seqs)
rm(Tce_g_seq_names)
write.csv(as.data.frame(rbind(Tce_g_GC)), "data/output/sel_out/Tce_g_GC.csv", row.names = F, quote = F)




Tcm_g_seqs <- read.fasta(file="data/genomes/Tcm_b3v08_1000.fasta")
Tcm_g_seq_names <- as.list(getName(Tcm_g_seqs ))
Tcm_g_GC <- get_g_GC("Tcm")
rm(Tcm_g_seqs)
rm(Tcm_g_seq_names)
write.csv(as.data.frame(rbind(Tcm_g_GC)), "data/output/sel_out/Tcm_g_GC.csv", row.names = F, quote = F)



Tpa_g_seqs <- read.fasta(file="data/genomes/Tpa_b3v08_1000.fasta")
Tpa_g_seq_names <- as.list(getName(Tpa_g_seqs ))
Tpa_g_GC <- get_g_GC("Tpa")
rm(Tpa_g_seqs)
rm(Tpa_g_seq_names)
write.csv(as.data.frame(rbind(Tpa_g_GC)), "data/output/sel_out/Tpa_g_GC.csv", row.names = F, quote = F)



Tps_g_seqs <- read.fasta(file="data/genomes/Tps_b3v08_1000.fasta")
Tps_g_seq_names <- as.list(getName(Tps_g_seqs ))
Tps_g_GC <- get_g_GC("Tps")
rm(Tps_g_seqs)
rm(Tps_g_seq_names)
write.csv(as.data.frame(rbind(Tps_g_GC)), "data/output/sel_out/Tps_g_GC.csv", row.names = F, quote = F)


## output
write.csv(as.data.frame(rbind(Tbi_g_GC, Tce_g_GC, Tcm_g_GC, Tpa_g_GC, Tps_g_GC)), "data/output/sel_out/GC_g.csv", row.names = F, quote = F)




#### calc CUB

Tbi_seqset <- readSet(file="data/output/sel_out/Ortho_Tbi.fa")
Tce_seqset <- readSet(file="data/output/sel_out/Ortho_Tce.fa")
Tcm_seqset <- readSet(file="data/output/sel_out/Ortho_Tcm.fa")
Tpa_seqset <- readSet(file="data/output/sel_out/Ortho_Tpa.fa")
Tps_seqset <- readSet(file="data/output/sel_out/Ortho_Tps.fa")

# Tte_seqset <- readSet(file="data/output/sel_out/Ortho_Tte.fa")
# Tms_seqset <- readSet(file="data/output/sel_out/Ortho_Tms.fa")
# Tsi_seqset <- readSet(file="data/output/sel_out/Ortho_Tsi.fa")
# Tge_seqset <- readSet(file="data/output/sel_out/Ortho_Tge.fa")
# Tdi_seqset <- readSet(file="data/output/sel_out/Ortho_Tdi.fa")


## ENC
# ENC (Wright 1990) ENC is a non-directional measure of CUB (Subramanian and Sarkar 2015), 
# and its values can range from 20 (high CUB) to 61 (no CUB).

get_CUB <- function(seq_set){
  ct_table <- codonTable(seq_set)
  out_df <- as.data.frame(cbind(
    getID(ct_table),
    ENC(ct_table)))
  
  colnames(out_df) <- c("Gene", "ENC")
  return(out_df)

}

Tbi_ENC <- get_CUB (Tbi_seqset)
Tce_ENC <- get_CUB (Tce_seqset)
Tcm_ENC <- get_CUB (Tcm_seqset)
Tpa_ENC <- get_CUB (Tpa_seqset)
Tps_ENC <- get_CUB (Tps_seqset)




### are pref codons GC biased?

Tbi_ct <- codonTable(Tbi_seqset )
Tce_ct <- codonTable(Tce_seqset )
Tcm_ct <- codonTable(Tcm_seqset )
Tpa_ct <- codonTable(Tpa_seqset )
Tps_ct <- codonTable(Tps_seqset )

get_all_codon_counts <- function(ct){
  cc <- as.data.frame(ct@counts)
  #print(cc)
  
  aa_list <- as.list(Biostrings::getGeneticCode(id_or_name2="1"))
  
  out_df = c()
  for(c in names(as.list(Biostrings::getGeneticCode(id_or_name2="1")))){
    #print(c)
    c_count <- sum(eval(parse(text=paste("cc",'$',c,sep=''))))
    AA <- eval(parse(text=paste("aa_list",'$"',c, '"', sep='')))
    #print(c_count)
    #print(AA)
    gc <- GC(c(substr(c,1,1), substr(c,2,2), substr(c,3,3)))
    #print(gc)
    out_df <- rbind(out_df, c(c, AA, c_count, gc))
  }
  
  out_df <- as.data.frame(out_df)
  colnames(out_df) <- c("codon", "AA", "count", "GC")
  out_df$count <- as.numeric(out_df$count)
  out_df$GC <- as.numeric(out_df$GC)
  
  ##########################################################################
  ### are pref codon GC-biased?
  
  aa_to_use <-  c("F", "L", "S", "Y", "C", "P", "H", "Q", "R", "I", "T", "N", "K", "V", "A", "D", "E", "G") ## not stops, or W or M (as have only 1 codon)
  N_aa_pref_GC_highest = 0
  for(a in aa_to_use){
    aa_s <- subset(out_df, out_df$AA == a)
    print(aa_s)
    GC_of_pref_codon <- subset(aa_s, aa_s$count == max(aa_s$count))$GC
    highest_GC_pos   <- max(aa_s$GC)
    print(GC_of_pref_codon )   
    print(highest_GC_pos)
    
    pref_GC <- ifelse(GC_of_pref_codon >= highest_GC_pos, 1, 0)
    print(pref_GC)
    N_aa_pref_GC_highest =  N_aa_pref_GC_highest + pref_GC
    
  }
  
  print("N aa with syn codons where a codon with the highest GC is also the prefered codon")
  print(N_aa_pref_GC_highest)
  print("Total N aa with syn codons")
  print(length(aa_to_use))


  aa_to_use <-  c("P", "T", "V", "A", "G") ## not stops, or W or M (as have only 1 codon)
  N_aa_pref_GC_highest = 0
  for(a in aa_to_use){
    aa_s <- subset(out_df, out_df$AA == a)
    print(aa_s)
    GC_of_pref_codon <- subset(aa_s, aa_s$count == max(aa_s$count))$GC
    highest_GC_pos   <- max(aa_s$GC)
    print(GC_of_pref_codon )   
    print(highest_GC_pos)
    
    pref_GC <- ifelse(GC_of_pref_codon >= highest_GC_pos, 1, 0)
    print(pref_GC)
    N_aa_pref_GC_highest =  N_aa_pref_GC_highest + pref_GC
    
  }
  
  print("N aa with syn codons where a codon with the highest GC is also the prefered codon")
  print(N_aa_pref_GC_highest)
  print("Total N aa with syn codons")
  print(length(aa_to_use))
  
    
  return(out_df)
}

Tbi_aa <- get_all_codon_counts(Tbi_ct)
Tce_aa <- get_all_codon_counts(Tce_ct)
Tcm_aa <- get_all_codon_counts(Tcm_ct)
Tpa_aa <- get_all_codon_counts(Tpa_ct)
Tps_aa <- get_all_codon_counts(Tps_ct)



## output
write.csv(as.data.frame(rbind(Tbi_ENC, Tce_ENC, Tcm_ENC, Tpa_ENC, Tps_ENC)), "data/output/sel_out/ENC.csv", row.names = F, quote = F)

### 
writeLines(capture.output(sessionInfo()), "GC_CUB_sess_info.txt")

