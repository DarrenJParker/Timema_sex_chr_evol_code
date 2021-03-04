# Expression_analyses.R

### libs

library(ggplot2)
library(stringr)
library(modeest)
library(plyr)
library(Rtsne)
library(cowplot)
library(edgeR)
library(hash)



### read data

setwd("data/counts")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### data

dat1_Tbi_raw <- read.table ("Tbi_1000_chr_info_and_counts.csv", header = T, sep = ',')
dat1_Tce_raw <- read.table ("Tce_1000_chr_info_and_counts.csv", header = T, sep = ',')
dat1_Tcm_raw <- read.table ("Tcm_1000_chr_info_and_counts.csv", header = T, sep = ',')
dat1_Tpa_raw <- read.table ("Tpa_1000_chr_info_and_counts.csv", header = T, sep = ',')
dat1_Tps_raw <- read.table ("Tps_1000_chr_info_and_counts.csv", header = T, sep = ',')

dat1_Tbi_raw$lg <- str_split_fixed(as.character(dat1_Tbi_raw$lg_pos), "_",2)[,1]
dat1_Tce_raw$lg <- str_split_fixed(as.character(dat1_Tce_raw$lg_pos), "_",2)[,1]
dat1_Tcm_raw$lg <- str_split_fixed(as.character(dat1_Tcm_raw$lg_pos), "_",2)[,1]
dat1_Tpa_raw$lg <- str_split_fixed(as.character(dat1_Tpa_raw$lg_pos), "_",2)[,1]
dat1_Tps_raw$lg <- str_split_fixed(as.character(dat1_Tps_raw$lg_pos), "_",2)[,1]

#### change wd to output folder
dir.create("../Exp_out")
setwd("../Exp_out")


#########################################################################################################################
## functions

filt_and_norm_maj_FPKM <- function(y,cpm_cut,cut_in_Nsams,gene_lens){
	
	cat("\nNumber number of genes / samples in orig data\n")
	print(dim(y)) ### number of genes / samples
	#print(head(rpkm(y, gene.length=gene_lens)))
	keep <- 
	  rowSums(rpkm(y[,1:3], gene.length=gene_lens,log=FALSE)> cpm_cut) >= cut_in_Nsams & 
      rowSums(rpkm(y[,4:6], gene.length=gene_lens,log=FALSE)> cpm_cut) >= cut_in_Nsams 	

	y <- y[keep,]
	
	y$samples$lib.size <- colSums(y$counts) 
	y <- calcNormFactors(y)
	cat("\nLib norm factors\n")
	print(y$samples)	
	cat("\nNumber number of genes / samples after filtering\n")
	print(dim(y))
	
	return(y)
}


filt_and_norm_maj_TPM <- function(y,cpm_cut,cut_in_Nsams,gene_lens){
	cat("\nNumber number of genes / samples in orig data\n")
	print(dim(y)) ### number of genes / samples

	rkpm_df <- rpkm(y, gene.length=gene_lens, log=FALSE)
	TPM_df  <- c()
	
	for(i in seq(1: length(colnames(rkpm_df)))){
		
		TPM <- 10^6 * (rkpm_df[,i] / sum(rkpm_df[,i]))
		TPM_df <- cbind(TPM_df, TPM)
	}
	
	# print(head(rkpm_df))
	# print(head(TPM_df))	
	
	keep <- 
	  rowSums(TPM_df[,1:3] > cpm_cut) >= cut_in_Nsams & 
      rowSums(TPM_df[,4:6] > cpm_cut) >= cut_in_Nsams 		

	y <- y[keep,]
	
	y$samples$lib.size <- colSums(y$counts) 
	y <- calcNormFactors(y)
	cat("\nLib norm factors\n")
	print(y$samples)	
	cat("\nNumber number of genes / samples after filtering\n")
	print(dim(y))
	
	return(y)

}


####################
## exper analyses 
###################################################################################
## filter out genes I did not classify with coverage

dat1_Tbi <- subset(dat1_Tbi_raw, !is.na(dat1_Tbi_raw$chr_hard))
dat1_Tce <- subset(dat1_Tce_raw, !is.na(dat1_Tce_raw$chr_hard))
dat1_Tcm <- subset(dat1_Tcm_raw, !is.na(dat1_Tcm_raw$chr_hard))
dat1_Tpa <- subset(dat1_Tpa_raw, !is.na(dat1_Tpa_raw$chr_hard))
dat1_Tps <- subset(dat1_Tps_raw, !is.na(dat1_Tps_raw$chr_hard))

## drops a few genes
length(dat1_Tbi_raw[,1]) - length(dat1_Tbi[,1])
length(dat1_Tce_raw[,1]) - length(dat1_Tce[,1])
length(dat1_Tcm_raw[,1]) - length(dat1_Tcm[,1])
length(dat1_Tpa_raw[,1]) - length(dat1_Tpa[,1])
length(dat1_Tps_raw[,1]) - length(dat1_Tps[,1])


#####################################
### read into DGE obj

## RT
y_Tbi_RT_UF <- DGEList(counts=dat1_Tbi[,c(
"Tbi_SF_RT_Md_Re1", "Tbi_SF_RT_Md_Re2", "Tbi_SF_RT_Md_Re3",
"Tbi_SM_RT_Md_Re1", "Tbi_SM_RT_Md_Re2", "Tbi_SM_RT_Md_Re3"
)], genes=dat1_Tbi$gene_id)

y_Tce_RT_UF <- DGEList(counts=dat1_Tce[,c(
"Tce_SF_RT_Md_Re1", "Tce_SF_RT_Md_Re2", "Tce_SF_RT_Md_Re3",
"Tce_SM_RT_Md_Re1", "Tce_SM_RT_Md_Re2", "Tce_SM_RT_Md_Re3"
)], genes=dat1_Tce$gene_id)

y_Tcm_RT_UF <- DGEList(counts=dat1_Tcm[,c(
"Tcm_SF_RT_Md_Re1", "Tcm_SF_RT_Md_Re2", "Tcm_SF_RT_Md_Re3",
"Tcm_SM_RT_Md_Re1", "Tcm_SM_RT_Md_Re2", "Tcm_SM_RT_Md_Re3"
)], genes=dat1_Tcm$gene_id)

y_Tpa_RT_UF <- DGEList(counts=dat1_Tpa[,c(
"Tpa_SF_RT_Md_Re1", "Tpa_SF_RT_Md_Re2", "Tpa_SF_RT_Md_Re3",
"Tpa_SM_RT_Md_Re1", "Tpa_SM_RT_Md_Re2", "Tpa_SM_RT_Md_Re3"
)], genes=dat1_Tpa$gene_id)

y_Tps_RT_UF <- DGEList(counts=dat1_Tps[,c(
"Tps_SF_RT_Md_Re1", "Tps_SF_RT_Md_Re2", "Tps_SF_RT_Md_Re3",
"Tps_SM_RT_Md_Re1", "Tps_SM_RT_Md_Re2", "Tps_SM_RT_Md_Re3"
)], genes=dat1_Tps$gene_id)


# HD
y_Tbi_HD_UF <- DGEList(counts=dat1_Tbi[,c(
"Tbi_SF_HD_Md_Re1", "Tbi_SF_HD_Md_Re2", "Tbi_SF_HD_Md_Re3",
"Tbi_SM_HD_Md_Re1", "Tbi_SM_HD_Md_Re2", "Tbi_SM_HD_Md_Re3"
)], genes=dat1_Tbi$gene_id)

y_Tce_HD_UF <- DGEList(counts=dat1_Tce[,c(
"Tce_SF_HD_Md_Re1", "Tce_SF_HD_Md_Re2", "Tce_SF_HD_Md_Re3",
"Tce_SM_HD_Md_Re1", "Tce_SM_HD_Md_Re2", "Tce_SM_HD_Md_Re3"
)], genes=dat1_Tce$gene_id)

y_Tcm_HD_UF <- DGEList(counts=dat1_Tcm[,c(
"Tcm_SF_HD_Md_Re1", "Tcm_SF_HD_Md_Re2", "Tcm_SF_HD_Md_Re3",
"Tcm_SM_HD_Md_Re1", "Tcm_SM_HD_Md_Re2", "Tcm_SM_HD_Md_Re3"
)], genes=dat1_Tcm$gene_id)

y_Tpa_HD_UF <- DGEList(counts=dat1_Tpa[,c(
"Tpa_SF_HD_Md_Re1", "Tpa_SF_HD_Md_Re2", "Tpa_SF_HD_Md_Re3",
"Tpa_SM_HD_Md_Re1", "Tpa_SM_HD_Md_Re2", "Tpa_SM_HD_Md_Re3"
)], genes=dat1_Tpa$gene_id)

y_Tps_HD_UF <- DGEList(counts=dat1_Tps[,c(
"Tps_SF_HD_Md_Re1", "Tps_SF_HD_Md_Re2", "Tps_SF_HD_Md_Re3",
"Tps_SM_HD_Md_Re1", "Tps_SM_HD_Md_Re2", "Tps_SM_HD_Md_Re3"
)], genes=dat1_Tps$gene_id)


# LG
y_Tbi_LG_UF <- DGEList(counts=dat1_Tbi[,c(
"Tbi_SF_LG_Md_Re1", "Tbi_SF_LG_Md_Re2", "Tbi_SF_LG_Md_Re3",
"Tbi_SM_LG_Md_Re1", "Tbi_SM_LG_Md_Re2", "Tbi_SM_LG_Md_Re3"
)], genes=dat1_Tbi$gene_id)


y_Tce_LG_UF <- DGEList(counts=dat1_Tce[,c(
"Tce_SF_LG_Md_Re1", "Tce_SF_LG_Md_Re2", "Tce_SF_LG_Md_Re3",
"Tce_SM_LG_Md_Re1", "Tce_SM_LG_Md_Re2", "Tce_SM_LG_Md_Re3"
)], genes=dat1_Tce$gene_id)

y_Tcm_LG_UF <- DGEList(counts=dat1_Tcm[,c(
"Tcm_SF_LG_Md_Re1", "Tcm_SF_LG_Md_Re2", "Tcm_SF_LG_Md_Re3",
"Tcm_SM_LG_Md_Re1", "Tcm_SM_LG_Md_Re2", "Tcm_SM_LG_Md_Re3"
)], genes=dat1_Tcm$gene_id)

y_Tpa_LG_UF <- DGEList(counts=dat1_Tpa[,c(
"Tpa_SF_LG_Md_Re1", "Tpa_SF_LG_Md_Re2", "Tpa_SF_LG_Md_Re3",
"Tpa_SM_LG_Md_Re1", "Tpa_SM_LG_Md_Re2", "Tpa_SM_LG_Md_Re3"
)], genes=dat1_Tpa$gene_id)

y_Tps_LG_UF <- DGEList(counts=dat1_Tps[,c(
"Tps_SF_LG_Md_Re1", "Tps_SF_LG_Md_Re2", "Tps_SF_LG_Md_Re3",
"Tps_SM_LG_Md_Re1", "Tps_SM_LG_Md_Re2", "Tps_SM_LG_Md_Re3"
)], genes=dat1_Tps$gene_id)




#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter low expressed gene by FPKM or TPM, then TMM normalisation 

# FPKM 

FPKM_filt_value = 2

y_Tbi_RT_F_FPKM <- filt_and_norm_maj_FPKM(y_Tbi_RT_UF,FPKM_filt_value,2,dat1_Tbi$tot_exon_len)
y_Tce_RT_F_FPKM <- filt_and_norm_maj_FPKM(y_Tce_RT_UF,FPKM_filt_value,2,dat1_Tce$tot_exon_len)
y_Tcm_RT_F_FPKM <- filt_and_norm_maj_FPKM(y_Tcm_RT_UF,FPKM_filt_value,2,dat1_Tcm$tot_exon_len)
y_Tpa_RT_F_FPKM <- filt_and_norm_maj_FPKM(y_Tpa_RT_UF,FPKM_filt_value,2,dat1_Tpa$tot_exon_len)
y_Tps_RT_F_FPKM <- filt_and_norm_maj_FPKM(y_Tps_RT_UF,FPKM_filt_value,2,dat1_Tps$tot_exon_len)

y_Tbi_HD_F_FPKM <- filt_and_norm_maj_FPKM(y_Tbi_HD_UF,FPKM_filt_value,2,dat1_Tbi$tot_exon_len)
y_Tce_HD_F_FPKM <- filt_and_norm_maj_FPKM(y_Tce_HD_UF,FPKM_filt_value,2,dat1_Tce$tot_exon_len)
y_Tcm_HD_F_FPKM <- filt_and_norm_maj_FPKM(y_Tcm_HD_UF,FPKM_filt_value,2,dat1_Tcm$tot_exon_len)
y_Tpa_HD_F_FPKM <- filt_and_norm_maj_FPKM(y_Tpa_HD_UF,FPKM_filt_value,2,dat1_Tpa$tot_exon_len)
y_Tps_HD_F_FPKM <- filt_and_norm_maj_FPKM(y_Tps_HD_UF,FPKM_filt_value,2,dat1_Tps$tot_exon_len)

y_Tbi_LG_F_FPKM <- filt_and_norm_maj_FPKM(y_Tbi_LG_UF,FPKM_filt_value,2,dat1_Tbi$tot_exon_len)
y_Tce_LG_F_FPKM <- filt_and_norm_maj_FPKM(y_Tce_LG_UF,FPKM_filt_value,2,dat1_Tce$tot_exon_len)
y_Tcm_LG_F_FPKM <- filt_and_norm_maj_FPKM(y_Tcm_LG_UF,FPKM_filt_value,2,dat1_Tcm$tot_exon_len)
y_Tpa_LG_F_FPKM <- filt_and_norm_maj_FPKM(y_Tpa_LG_UF,FPKM_filt_value,2,dat1_Tpa$tot_exon_len)
y_Tps_LG_F_FPKM <- filt_and_norm_maj_FPKM(y_Tps_LG_UF,FPKM_filt_value,2,dat1_Tps$tot_exon_len)


# TPM 

TPM_filt_value = 2

y_Tbi_RT_F_TPM <- filt_and_norm_maj_TPM(y_Tbi_RT_UF,TPM_filt_value,2,dat1_Tbi$tot_exon_len)
y_Tce_RT_F_TPM <- filt_and_norm_maj_TPM(y_Tce_RT_UF,TPM_filt_value,2,dat1_Tce$tot_exon_len)
y_Tcm_RT_F_TPM <- filt_and_norm_maj_TPM(y_Tcm_RT_UF,TPM_filt_value,2,dat1_Tcm$tot_exon_len)
y_Tpa_RT_F_TPM <- filt_and_norm_maj_TPM(y_Tpa_RT_UF,TPM_filt_value,2,dat1_Tpa$tot_exon_len)
y_Tps_RT_F_TPM <- filt_and_norm_maj_TPM(y_Tps_RT_UF,TPM_filt_value,2,dat1_Tps$tot_exon_len)

y_Tbi_HD_F_TPM <- filt_and_norm_maj_TPM(y_Tbi_HD_UF,TPM_filt_value,2,dat1_Tbi$tot_exon_len)
y_Tce_HD_F_TPM <- filt_and_norm_maj_TPM(y_Tce_HD_UF,TPM_filt_value,2,dat1_Tce$tot_exon_len)
y_Tcm_HD_F_TPM <- filt_and_norm_maj_TPM(y_Tcm_HD_UF,TPM_filt_value,2,dat1_Tcm$tot_exon_len)
y_Tpa_HD_F_TPM <- filt_and_norm_maj_TPM(y_Tpa_HD_UF,TPM_filt_value,2,dat1_Tpa$tot_exon_len)
y_Tps_HD_F_TPM <- filt_and_norm_maj_TPM(y_Tps_HD_UF,TPM_filt_value,2,dat1_Tps$tot_exon_len)

y_Tbi_LG_F_TPM <- filt_and_norm_maj_TPM(y_Tbi_LG_UF,TPM_filt_value,2,dat1_Tbi$tot_exon_len)
y_Tce_LG_F_TPM <- filt_and_norm_maj_TPM(y_Tce_LG_UF,TPM_filt_value,2,dat1_Tce$tot_exon_len)
y_Tcm_LG_F_TPM <- filt_and_norm_maj_TPM(y_Tcm_LG_UF,TPM_filt_value,2,dat1_Tcm$tot_exon_len)
y_Tpa_LG_F_TPM <- filt_and_norm_maj_TPM(y_Tpa_LG_UF,TPM_filt_value,2,dat1_Tpa$tot_exon_len)
y_Tps_LG_F_TPM <- filt_and_norm_maj_TPM(y_Tps_LG_UF,TPM_filt_value,2,dat1_Tps$tot_exon_len)



#################################################################################################################################################################
### get FPKM and TPMs


get_FPKM <- function(y_F, full_df, female_prefix, male_prefix){
	y_F_genenames_df = as.data.frame(y_F$genes)
	y_F_genenames_df$genes <- as.character(y_F_genenames_df$genes)
	subset_df <- subset(full_df, full_df$gene_id %in% y_F_genenames_df$genes)
	FPKM_df   <- rpkm(y_F, gene.length=subset_df$tot_exon_len, normalized.lib.sizes=TRUE, log=FALSE)

	col_out = c()
	for(i in colnames(FPKM_df)){
		col_out <- c(col_out, paste(i, "FPKM", sep = "_"))
	}
	colnames(FPKM_df) <- col_out

	out_df   <- cbind(subset_df[,1:8], subset_df[,(length(subset_df) - 4): length(subset_df)], FPKM_df)
	col_out_2 <- colnames(out_df )
	
	#### calc mean FPKM for males and females
	
	subData_F <- out_df[,grepl(female_prefix,colnames(out_df))]
	print("cols for females")
	print(colnames(subData_F))
	
	subData_M <- out_df[,grepl(male_prefix,colnames(out_df))]
	print("cols for males")
	print(colnames(subData_M))
		
	out_df$female_FPKM_mean <-rowMeans(subData_F)
	out_df$male_FPKM_mean   <-rowMeans(subData_M)	
	out_df$FPKM_mean <- (out_df$female_FPKM_mean + out_df$male_FPKM_mean) / 2	
	colnames(out_df) <- c(col_out_2, paste(female_prefix, "meanFPKM", sep = "_"), paste(male_prefix, "meanFPKM", sep = "_"),  paste(str_replace(male_prefix, "_SM", ""), "meanFPKM", sep = "_"))
	return(out_df )		
	
}

get_TPM <- function(y_F, full_df, female_prefix, male_prefix){
	y_F_genenames_df = as.data.frame(y_F$genes)
	y_F_genenames_df$genes <- as.character(y_F_genenames_df$genes)
	subset_df <- subset(full_df, full_df$gene_id %in% y_F_genenames_df$genes)
	FPKM_df   <- rpkm(y_F, gene.length=subset_df$tot_exon_len, normalized.lib.sizes=TRUE, log=FALSE)
	TPM_df  <- c()
	for(i in seq(1: length(colnames(FPKM_df)))){
		
		TPM <- 10^6 * (FPKM_df[,i] / sum(FPKM_df[,i]))
		TPM_df <- cbind(TPM_df, TPM)
	}
	
	col_out = c()
	for(i in colnames(FPKM_df)){
		col_out <- c(col_out, paste(i, "TPM", sep = "_"))
	}
	
	colnames(TPM_df) <- col_out
	out_df   <- cbind(subset_df[,1:8], subset_df[,(length(subset_df) - 4): length(subset_df)], TPM_df)
	col_out_2 <- colnames(out_df )
	
	#### calc mean TPM for males and females
	
	subData_F <- out_df[,grepl(female_prefix,colnames(out_df))]
	print("cols for females")
	print(colnames(subData_F))
	
	subData_M <- out_df[,grepl(male_prefix,colnames(out_df))]
	print("cols for males")
	print(colnames(subData_M))
		
	out_df$female_TPM_mean <-rowMeans(subData_F)
	out_df$male_TPM_mean   <-rowMeans(subData_M)	
	out_df$TPM_mean <- (out_df$female_TPM_mean + out_df$male_TPM_mean) / 2
	
	colnames(out_df) <- c(col_out_2, paste(female_prefix, "meanTPM", sep = "_"), paste(male_prefix, "meanTPM", sep = "_"),  paste(str_replace(male_prefix, "_SM", ""), "meanTPM", sep = "_"))
	return(out_df )	

}


Tbi_RT_F_FPKM <- get_FPKM(y_Tbi_RT_F_FPKM, dat1_Tbi, "Tbi_SF_RT", "Tbi_SM_RT")
Tce_RT_F_FPKM <- get_FPKM(y_Tce_RT_F_FPKM, dat1_Tce, "Tce_SF_RT", "Tce_SM_RT")
Tcm_RT_F_FPKM <- get_FPKM(y_Tcm_RT_F_FPKM, dat1_Tcm, "Tcm_SF_RT", "Tcm_SM_RT")
Tpa_RT_F_FPKM <- get_FPKM(y_Tpa_RT_F_FPKM, dat1_Tpa, "Tpa_SF_RT", "Tpa_SM_RT")
Tps_RT_F_FPKM <- get_FPKM(y_Tps_RT_F_FPKM, dat1_Tps, "Tps_SF_RT", "Tps_SM_RT")

Tbi_HD_F_FPKM <- get_FPKM(y_Tbi_HD_F_FPKM, dat1_Tbi, "Tbi_SF_HD", "Tbi_SM_HD")
Tce_HD_F_FPKM <- get_FPKM(y_Tce_HD_F_FPKM, dat1_Tce, "Tce_SF_HD", "Tce_SM_HD")
Tcm_HD_F_FPKM <- get_FPKM(y_Tcm_HD_F_FPKM, dat1_Tcm, "Tcm_SF_HD", "Tcm_SM_HD")
Tpa_HD_F_FPKM <- get_FPKM(y_Tpa_HD_F_FPKM, dat1_Tpa, "Tpa_SF_HD", "Tpa_SM_HD")
Tps_HD_F_FPKM <- get_FPKM(y_Tps_HD_F_FPKM, dat1_Tps, "Tps_SF_HD", "Tps_SM_HD")

Tbi_LG_F_FPKM <- get_FPKM(y_Tbi_LG_F_FPKM, dat1_Tbi, "Tbi_SF_LG", "Tbi_SM_LG")
Tce_LG_F_FPKM <- get_FPKM(y_Tce_LG_F_FPKM, dat1_Tce, "Tce_SF_LG", "Tce_SM_LG")
Tcm_LG_F_FPKM <- get_FPKM(y_Tcm_LG_F_FPKM, dat1_Tcm, "Tcm_SF_LG", "Tcm_SM_LG")
Tpa_LG_F_FPKM <- get_FPKM(y_Tpa_LG_F_FPKM, dat1_Tpa, "Tpa_SF_LG", "Tpa_SM_LG")
Tps_LG_F_FPKM <- get_FPKM(y_Tps_LG_F_FPKM, dat1_Tps, "Tps_SF_LG", "Tps_SM_LG")

Tbi_RT_F_TPM <- get_TPM(y_Tbi_RT_F_TPM, dat1_Tbi, "Tbi_SF_RT", "Tbi_SM_RT")
Tce_RT_F_TPM <- get_TPM(y_Tce_RT_F_TPM, dat1_Tce, "Tce_SF_RT", "Tce_SM_RT")
Tcm_RT_F_TPM <- get_TPM(y_Tcm_RT_F_TPM, dat1_Tcm, "Tcm_SF_RT", "Tcm_SM_RT")
Tpa_RT_F_TPM <- get_TPM(y_Tpa_RT_F_TPM, dat1_Tpa, "Tpa_SF_RT", "Tpa_SM_RT")
Tps_RT_F_TPM <- get_TPM(y_Tps_RT_F_TPM, dat1_Tps, "Tps_SF_RT", "Tps_SM_RT")

Tbi_HD_F_TPM <- get_TPM(y_Tbi_HD_F_TPM, dat1_Tbi, "Tbi_SF_HD", "Tbi_SM_HD")
Tce_HD_F_TPM <- get_TPM(y_Tce_HD_F_TPM, dat1_Tce, "Tce_SF_HD", "Tce_SM_HD")
Tcm_HD_F_TPM <- get_TPM(y_Tcm_HD_F_TPM, dat1_Tcm, "Tcm_SF_HD", "Tcm_SM_HD")
Tpa_HD_F_TPM <- get_TPM(y_Tpa_HD_F_TPM, dat1_Tpa, "Tpa_SF_HD", "Tpa_SM_HD")
Tps_HD_F_TPM <- get_TPM(y_Tps_HD_F_TPM, dat1_Tps, "Tps_SF_HD", "Tps_SM_HD")

Tbi_LG_F_TPM <- get_TPM(y_Tbi_LG_F_TPM, dat1_Tbi, "Tbi_SF_LG", "Tbi_SM_LG")
Tce_LG_F_TPM <- get_TPM(y_Tce_LG_F_TPM, dat1_Tce, "Tce_SF_LG", "Tce_SM_LG")
Tcm_LG_F_TPM <- get_TPM(y_Tcm_LG_F_TPM, dat1_Tcm, "Tcm_SF_LG", "Tcm_SM_LG")
Tpa_LG_F_TPM <- get_TPM(y_Tpa_LG_F_TPM, dat1_Tpa, "Tpa_SF_LG", "Tpa_SM_LG")
Tps_LG_F_TPM <- get_TPM(y_Tps_LG_F_TPM, dat1_Tps, "Tps_SF_LG", "Tps_SM_LG")


head(Tbi_WB_F_TPM )

###################################################################################
### Dosage comp

make_long_table_avEXP = function(dat_Tbi,dat_Tce,dat_Tcm,dat_Tpa,dat_Tps,tiss, exp_type){
	
	sp_u = "Tbi"
	df_t1_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tbi,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tbi,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tbi,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tbi,'$chr_hard',sep='')))))) 
	df_t1_F$sp      = rep(sp_u, length(df_t1_F[,1]))
	df_t1_F$sex     = rep("F", length(df_t1_F[,1]))	

	df_t1_M = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tbi,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tbi,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tbi,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tbi,'$chr_hard',sep='')))))) 
	df_t1_M$sp      = rep(sp_u, length(df_t1_M[,1]))
	df_t1_M$sex     = rep("M", length(df_t1_M[,1]))		

	sp_u = "Tce"
	df_t2_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tce,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tce,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tce,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tce,'$chr_hard',sep='')))))) 
	df_t2_F$sp      = rep(sp_u, length(df_t2_F[,1]))
	df_t2_F$sex     = rep("F", length(df_t2_F[,1]))	

	df_t2_M = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tce,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tce,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tce,'$chr_soft',sep='')))), 
		as.character(eval(parse(text=paste(dat_Tce,'$chr_hard',sep='')))))) 
	df_t2_M$sp      = rep(sp_u, length(df_t2_M[,1]))
	df_t2_M$sex     = rep("M", length(df_t2_M[,1]))		

	sp_u = "Tcm"
	df_t3_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tcm,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tcm,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tcm,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tcm,'$chr_hard',sep='')))))) 
	df_t3_F$sp      = rep(sp_u, length(df_t3_F[,1]))
	df_t3_F$sex     = rep("F", length(df_t3_F[,1]))	

	df_t3_M = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tcm,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tcm,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tcm,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tcm,'$chr_hard',sep='')))))) 
	df_t3_M$sp      = rep(sp_u, length(df_t3_M[,1]))
	df_t3_M$sex     = rep("M", length(df_t3_M[,1]))		

	sp_u = "Tpa"
	df_t4_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tpa,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tpa,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tpa,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tpa,'$chr_hard',sep='')))))) 
	df_t4_F$sp      = rep(sp_u, length(df_t4_F[,1]))
	df_t4_F$sex     = rep("F", length(df_t4_F[,1]))	

	df_t4_M = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tpa,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tpa,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tpa,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tpa,'$chr_hard',sep='')))))) 
	df_t4_M$sp      = rep(sp_u, length(df_t4_M[,1]))
	df_t4_M$sex     = rep("M", length(df_t4_M[,1]))		

	sp_u = "Tps"
	df_t5_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tps,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tps,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tps,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tps,'$chr_hard',sep='')))))) 
	df_t5_F$sp      = rep(sp_u, length(df_t5_F[,1]))
	df_t5_F$sex     = rep("F", length(df_t5_F[,1]))	

	df_t5_M = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tps,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tps,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tps,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tps,'$chr_hard',sep='')))))) 
	df_t5_M$sp      = rep(sp_u, length(df_t5_M[,1]))
	df_t5_M$sex     = rep("M", length(df_t5_M[,1]))		


	df_t2 <- rbind(df_t1_F, df_t1_M, df_t2_F, df_t2_M, df_t3_F, df_t3_M, df_t4_F, df_t4_M, df_t5_F, df_t5_M)
	df_t2$V2 <- as.numeric(as.character(df_t2$V2))
	colnames(df_t2) <- c("gene_id", paste("Av", exp_type, sep=""), "chr_soft", "chr_hard", "sp", "sex") 
	
	df_t = na.omit(df_t2)
	
	
	df_t$tiss <- as.factor(rep(tiss, length(df_t[,1])))
	df_t$sex <- as.factor(df_t$sex)
	df_t$sp  <- as.factor(df_t$sp)
	df_t$group_soft <- as.factor(paste(df_t$sp, df_t$sex, df_t$chr_soft, sep = "_"))
	df_t$sp_chr_soft <- as.factor(paste(df_t$sp, df_t$chr_soft, sep = "_"))	
	df_t$sex_chr_soft <- as.factor(paste(df_t$sex, df_t$chr_soft, sep = "_"))	
	df_t$group_soft_ord <- ordered(df_t$group_soft, levels = c(
	"Tbi_M_A", "Tbi_M_X", "Tbi_F_A", "Tbi_F_X",
	"Tce_M_A", "Tce_M_X", "Tce_F_A", "Tce_F_X",
	"Tps_M_A", "Tps_M_X", "Tps_F_A", "Tps_F_X",
	"Tcm_M_A", "Tcm_M_X", "Tcm_F_A", "Tcm_F_X",
	"Tpa_M_A", "Tpa_M_X", "Tpa_F_A", "Tpa_F_X"	
	 ))
	df_t$group_hard <- as.factor(paste(df_t$sp, df_t$sex, df_t$chr_hard, sep = "_"))
	df_t$sp_chr_hard <- as.factor(paste(df_t$sp, df_t$chr_hard, sep = "_"))	
	df_t$sex_chr_hard <- as.factor(paste(df_t$sex, df_t$chr_hard, sep = "_"))	
	df_t$group_hard_ord <- ordered(df_t$group_hard, levels = c(
	"Tbi_M_A", "Tbi_M_X", "Tbi_F_A", "Tbi_F_X",
	"Tce_M_A", "Tce_M_X", "Tce_F_A", "Tce_F_X",
	"Tps_M_A", "Tps_M_X", "Tps_F_A", "Tps_F_X",
	"Tcm_M_A", "Tcm_M_X", "Tcm_F_A", "Tcm_F_X",
	"Tpa_M_A", "Tpa_M_X", "Tpa_F_A", "Tpa_F_X"	
	 ))
	return(df_t)
}

All_RT_F_FPKM_long <- make_long_table_avEXP("Tbi_RT_F_FPKM", "Tce_RT_F_FPKM", "Tcm_RT_F_FPKM", "Tpa_RT_F_FPKM", "Tps_RT_F_FPKM",  "RT", "FPKM")
All_HD_F_FPKM_long <- make_long_table_avEXP("Tbi_HD_F_FPKM", "Tce_HD_F_FPKM", "Tcm_HD_F_FPKM", "Tpa_HD_F_FPKM", "Tps_HD_F_FPKM",  "HD", "FPKM")
All_LG_F_FPKM_long <- make_long_table_avEXP("Tbi_LG_F_FPKM", "Tce_LG_F_FPKM", "Tcm_LG_F_FPKM", "Tpa_LG_F_FPKM", "Tps_LG_F_FPKM",  "LG", "FPKM")

All_RT_F_TPM_long <- make_long_table_avEXP("Tbi_RT_F_TPM", "Tce_RT_F_TPM", "Tcm_RT_F_TPM", "Tpa_RT_F_TPM", "Tps_RT_F_TPM",  "RT", "TPM")
All_HD_F_TPM_long <- make_long_table_avEXP("Tbi_HD_F_TPM", "Tce_HD_F_TPM", "Tcm_HD_F_TPM", "Tpa_HD_F_TPM", "Tps_HD_F_TPM",  "HD", "TPM")
All_LG_F_TPM_long <- make_long_table_avEXP("Tbi_LG_F_TPM", "Tce_LG_F_TPM", "Tcm_LG_F_TPM", "Tpa_LG_F_TPM", "Tps_LG_F_TPM",  "LG", "TPM")

###############################################################################################################################
##### box plots

plot_avFPKM <- function(df, max_y, tissue){
	P1_hard <- ggplot(df, aes(sp, AvFPKM)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(sex_chr_hard)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(0,max_y)) +
		ylab ("FPKM") +
		xlab ("") + 
		scale_fill_manual(values=c("grey", "firebrick2","white", "royalblue2"), labels=c("Female A", "Female X", "Male A", "Male X")) +
		ggtitle(paste(tissue, "|", "hard class" ))

	P1_soft <- ggplot(df, aes(sp, AvFPKM)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(sex_chr_soft)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(0,max_y)) +
		ylab ("FPKM") +
		xlab ("") + 
		scale_fill_manual(values=c("grey", "firebrick2","white", "royalblue2"), labels=c("Female A", "Female X", "Male A", "Male X")) +
		ggtitle(paste(tissue, "|", "soft class" ))

	outlist = list("P1_hard" = P1_hard, "P1_soft" = P1_soft) 
	return(outlist)	
}

plot_avTPM <- function(df, max_y, tissue){
	P1_hard <- ggplot(df, aes(sp, AvTPM)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(sex_chr_hard)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(0,max_y)) +
		ylab ("TPM") +
		xlab ("") + 
		scale_fill_manual(values=c("grey", "firebrick2","white", "royalblue2"), labels=c("Female A", "Female X", "Male A", "Male X")) +
		ggtitle(paste(tissue, "|", "hard class" ))

	P1_soft <- ggplot(df, aes(sp, AvTPM)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(sex_chr_soft)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(0,max_y)) +
		ylab ("TPM") +
		xlab ("") + 
		scale_fill_manual(values=c("grey", "firebrick2","white", "royalblue2"), labels=c("Female A", "Female X", "Male A", "Male X")) +
		ggtitle(paste(tissue, "|", "soft class" ))

	outlist = list("P1_hard" = P1_hard, "P1_soft" = P1_soft) 
	return(outlist)	
}


RT_avFPKM <- plot_avFPKM(All_RT_F_FPKM_long, 200, "RT")
HD_avFPKM <- plot_avFPKM(All_HD_F_FPKM_long, 80,  "HD")
LG_avFPKM <- plot_avFPKM(All_LG_F_FPKM_long, 80,  "LG")

RT_avTPM <- plot_avTPM(All_RT_F_TPM_long, 200, "RT")
HD_avTPM <- plot_avTPM(All_HD_F_TPM_long, 80,  "HD")
LG_avTPM <- plot_avTPM(All_LG_F_TPM_long, 80,  "LG")


### output

pdf("RTHDLG_avFPKM_soft.pdf", width = 9, height = 9)
plot_grid(
HD_avFPKM$P1_soft, 
LG_avFPKM$P1_soft, 
RT_avFPKM$P1_soft, ncol = 1)
dev.off()
getwd() ## where has my plot gone....

pdf("RTHDLG_avTPM_soft.pdf", width = 9, height = 9)
plot_grid( 
HD_avTPM$P1_soft, 
LG_avTPM$P1_soft, 
RT_avTPM$P1_soft, ncol = 1)
dev.off()
getwd() ## where has my plot gone....


pdf("RTHDLG_avFPKM_hard.pdf", width = 9, height = 9)
plot_grid(
HD_avFPKM$P1_hard, 
LG_avFPKM$P1_hard,
RT_avFPKM$P1_hard,  ncol = 1)
dev.off()
getwd() ## where has my plot gone....

pdf("RTHDLG_avTPM_hard.pdf", width = 9, height = 9)
plot_grid(
HD_avTPM$P1_hard, 
LG_avTPM$P1_hard, 
RT_avTPM$P1_hard, ncol = 1)
dev.off()
getwd() ## where has my plot gone....



################################################################################################
#### wilcoxon tests on expression



wilcox_exp <- function(chr_type, exp_type){
	sp_want   <- c("Tbi", "Tce", "Tcm", "Tpa", "Tps")
	tiss_want <- c("RT", "HD", "LG")

	XA_wilcox_SFSM_df <- c()
	for(sp in sp_want){
		for(tiss in tiss_want){
	
		curr_df <- eval(parse(text=paste(sp , '_', tiss, '_F_', exp_type ,sep='')))
		X_df <- subset(curr_df, eval(parse(text=paste('curr_df' , '$' , 'chr_', chr_type, sep = ''))) == "X")
		A_df <- subset(curr_df, eval(parse(text=paste('curr_df' , '$' , 'chr_', chr_type, sep = ''))) == "A")
	
		SF_XA_wx <- wilcox.test(eval(parse(text=paste('X_df', '$' , sp, '_SF_', tiss, '_mean', exp_type ,sep = ''))),   eval(parse(text=paste('A_df' , '$' , sp, '_SF_', tiss, '_mean', exp_type,sep=''))), paired = FALSE)	
		SM_XA_wx <- wilcox.test(eval(parse(text=paste('X_df', '$' , sp, '_SM_', tiss, '_mean', exp_type ,sep = ''))),   eval(parse(text=paste('A_df' , '$' , sp, '_SM_', tiss, '_mean', exp_type,sep=''))), paired = FALSE)
	
		XA_wilcox_SFSM_df <- as.data.frame(rbind(
		XA_wilcox_SFSM_df,
		c(sp, tiss, "SF", chr_type, exp_type, SF_XA_wx$p.value),
		c(sp, tiss, "SM", chr_type, exp_type, SM_XA_wx$p.value)		
		))
		}
	}
	colnames(XA_wilcox_SFSM_df) <- c("sp", "tiss", "Sex", "chr_class", "exp_type", "wilcox_p")
	XA_wilcox_SFSM_df$FDR <- p.adjust(as.numeric(as.character(XA_wilcox_SFSM_df$wilcox_p)), method = "BH")
	return(XA_wilcox_SFSM_df )
}



#### output

write.csv(wilcox_exp("hard", "FPKM"), "XA_wilcox_SFSM_hard_FPKM_df.csv", row.names = F)
write.csv(wilcox_exp("soft", "FPKM"), "XA_wilcox_SFSM_soft_FPKM_df.csv", row.names = F)
write.csv(wilcox_exp("hard", "TPM"), "XA_wilcox_SFSM_hard_TPM_df.csv", row.names = F)
write.csv(wilcox_exp("soft", "TPM"), "XA_wilcox_SFSM_soft_TPM_df.csv", row.names = F)



###############################################################################################################################
##### log2(M/F)

make_long_table_log2MF = function(dat_Tbi,dat_Tce,dat_Tcm,dat_Tpa,dat_Tps,tiss, exp_type ){
	
	sp_u = "Tbi"
	df_t1_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tbi,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tbi,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
		eval(parse(text=paste(dat_Tbi,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tbi,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tbi,'$chr_hard',sep='')))),
		as.character(eval(parse(text=paste(dat_Tbi,'$lg',sep='')))),
		as.character(eval(parse(text=paste(dat_Tbi,'$multi_lg',sep=''))))))  
	df_t1_F$sp      = rep(sp_u, length(df_t1_F[,1]))

	sp_u = "Tce"
	df_t2_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tce,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tce,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
 		eval(parse(text=paste(dat_Tce,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tce,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tce,'$chr_hard',sep='')))),
		as.character(eval(parse(text=paste(dat_Tce,'$lg',sep='')))),
		as.character(eval(parse(text=paste(dat_Tce,'$multi_lg',sep='')))))) 
	df_t2_F$sp      = rep(sp_u, length(df_t2_F[,1]))

	sp_u = "Tcm"
	df_t3_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tcm,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tcm,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
		eval(parse(text=paste(dat_Tcm,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tcm,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tcm,'$chr_hard',sep='')))),
		as.character(eval(parse(text=paste(dat_Tcm,'$lg',sep='')))),
		as.character(eval(parse(text=paste(dat_Tcm,'$multi_lg',sep='')))))) 
	df_t3_F$sp      = rep(sp_u, length(df_t3_F[,1]))
	
	sp_u = "Tpa"
	df_t4_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tpa,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tpa,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
 		eval(parse(text=paste(dat_Tpa,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tpa,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tpa,'$chr_hard',sep='')))),
		as.character(eval(parse(text=paste(dat_Tpa,'$lg',sep='')))),
		as.character(eval(parse(text=paste(dat_Tpa,'$multi_lg',sep='')))))) 
	df_t4_F$sp      = rep(sp_u, length(df_t4_F[,1]))

	sp_u = "Tps"
	df_t5_F = as.data.frame(cbind(
		as.character(eval(parse(text=paste(dat_Tps,'$gene_id',sep='')))), 
		eval(parse(text=paste(dat_Tps,'$', sp_u, '_SF_' ,tiss,'_mean', exp_type,sep=''))), 
		eval(parse(text=paste(dat_Tps,'$', sp_u, '_SM_' ,tiss,'_mean', exp_type,sep=''))), 
		as.character(eval(parse(text=paste(dat_Tps,'$chr_soft',sep='')))),
		as.character(eval(parse(text=paste(dat_Tps,'$chr_hard',sep='')))),
		as.character(eval(parse(text=paste(dat_Tps,'$lg',sep='')))),
		as.character(eval(parse(text=paste(dat_Tps,'$multi_lg',sep='')))))) 
	df_t5_F$sp      = rep(sp_u, length(df_t5_F[,1]))


	df_t2 <- rbind(df_t1_F, df_t2_F, df_t3_F, df_t4_F, df_t5_F)
	
	df_t2$V2 <- as.numeric(as.character(df_t2$V2)) ## female
	df_t2$V3 <- as.numeric(as.character(df_t2$V3))
	df_t2$log2MF <- log2(df_t2$V3 / df_t2$V2)

	colnames(df_t2) <- c("gene_id", paste("Av", exp_type, "_Female", sep=""), paste("Av", exp_type, "_Male", sep=""), "chr_soft", "chr_hard", "lg", "multi_lg", "sp", paste("Av", exp_type, "_log2MF", sep="")) 
	
	#df_t = na.omit(df_t2)
	df_t <-  df_t2
	
	df_t$tiss <- as.factor(rep(tiss, length(df_t[,1])))
	df_t$sp  <- as.factor(df_t$sp)
	df_t$group_soft <- as.factor(paste(df_t$sp, df_t$chr_soft, sep = "_"))
	df_t$group_soft_ord <- ordered(df_t$group_soft, levels = c(
	"Tbi_A", "Tbi_X",
	"Tce_A", "Tce_X",
	"Tps_A", "Tps_X",
	"Tcm_A", "Tcm_X",
	"Tpa_A", "Tpa_X"	
	 ))
	df_t$group_hard <- as.factor(paste(df_t$sp, df_t$chr_hard, sep = "_"))
	df_t$group_hard_ord <- ordered(df_t$group_hard, levels = c(
	"Tbi_A", "Tbi_X",
	"Tce_A", "Tce_X",
	"Tps_A", "Tps_X",
	"Tcm_A", "Tcm_X",
	"Tpa_A", "Tpa_X"	
	 ))
	return(df_t)
}


All_RT_F_FPKM_log2MF_long <- make_long_table_log2MF("Tbi_RT_F_FPKM", "Tce_RT_F_FPKM", "Tcm_RT_F_FPKM", "Tpa_RT_F_FPKM", "Tps_RT_F_FPKM",  "RT", "FPKM")
All_HD_F_FPKM_log2MF_long <- make_long_table_log2MF("Tbi_HD_F_FPKM", "Tce_HD_F_FPKM", "Tcm_HD_F_FPKM", "Tpa_HD_F_FPKM", "Tps_HD_F_FPKM",  "HD", "FPKM")
All_LG_F_FPKM_log2MF_long <- make_long_table_log2MF("Tbi_LG_F_FPKM", "Tce_LG_F_FPKM", "Tcm_LG_F_FPKM", "Tpa_LG_F_FPKM", "Tps_LG_F_FPKM",  "LG", "FPKM")

All_RT_F_TPM_log2MF_long <- make_long_table_log2MF("Tbi_RT_F_TPM", "Tce_RT_F_TPM", "Tcm_RT_F_TPM", "Tpa_RT_F_TPM", "Tps_RT_F_TPM",  "RT", "TPM")
All_HD_F_TPM_log2MF_long <- make_long_table_log2MF("Tbi_HD_F_TPM", "Tce_HD_F_TPM", "Tcm_HD_F_TPM", "Tpa_HD_F_TPM", "Tps_HD_F_TPM",  "HD", "TPM")
All_LG_F_TPM_log2MF_long <- make_long_table_log2MF("Tbi_LG_F_TPM", "Tce_LG_F_TPM", "Tcm_LG_F_TPM", "Tpa_LG_F_TPM", "Tps_LG_F_TPM",  "LG", "TPM")


###############################################################################################################################
##### boxplots


plot_FPKM_log2MF <- function(df, tissue){
	
	P1_hard <- ggplot(df, aes(sp, AvFPKM_log2MF)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(chr_hard)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(-4,4)) +
		ylab ("log2(Male FPKM / Female FPKM)") +
		xlab ("Species") + 
		scale_fill_manual(values=c("grey", "darkorange"))  + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
		ggtitle(paste(tissue, "|", "hard class" ))

	P1_soft <- ggplot(df, aes(sp, AvFPKM_log2MF)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(chr_soft)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(-4,4)) +
		ylab ("log2(Male FPKM / Female FPKM)") +
		xlab ("Species") + 
		scale_fill_manual(values=c("grey", "darkorange"))  + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
		ggtitle(paste(tissue, "|", "soft class" ))
	
	outlist = list("P1_hard" = P1_hard, "P1_soft" = P1_soft) 
	return(outlist)	
}



plot_TPM_log2MF <- function(df, tissue){
	
	P1_hard <- ggplot(df, aes(sp, AvTPM_log2MF)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(chr_hard)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(-4,4)) +
		ylab ("log2(Male TPM / Female TPM)") +
		xlab ("Species") + 
		scale_fill_manual(values=c("grey", "darkorange"))  + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
		ggtitle(paste(tissue, "|", "hard class" ))

	P1_soft <- ggplot(df, aes(sp, AvTPM_log2MF)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(chr_soft)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(-4,4)) +
		ylab ("log2(Male TPM / Female TPM)") +
		xlab ("Species") + 
		scale_fill_manual(values=c("grey", "darkorange"))  + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
		ggtitle(paste(tissue, "|", "soft class" ))
	
	outlist = list("P1_hard" = P1_hard, "P1_soft" = P1_soft) 
	return(outlist)	
}


RT_FPKM_log2MF <- plot_FPKM_log2MF(All_RT_F_FPKM_log2MF_long, "RT")
HD_FPKM_log2MF <- plot_FPKM_log2MF(All_HD_F_FPKM_log2MF_long, "HD")
LG_FPKM_log2MF <- plot_FPKM_log2MF(All_LG_F_FPKM_log2MF_long, "LG")

RT_TPM_log2MF <- plot_TPM_log2MF(All_RT_F_TPM_log2MF_long, "RT")
HD_TPM_log2MF <- plot_TPM_log2MF(All_HD_F_TPM_log2MF_long, "HD")
LG_TPM_log2MF <- plot_TPM_log2MF(All_LG_F_TPM_log2MF_long, "LG")



######## by LG


plot_FPKM_log2MF_LG <- function(df, tissue){
	df_sub <- subset(df, df$multi_lg == "NO")
	df_sub$lg_ord <- ordered(df_sub$lg, levels = c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8", "lg9", "lg10", "lg11", "lg12", "lgX"))
	
	print(head(df_sub))

	P1 <- ggplot(df_sub, aes(sp, AvFPKM_log2MF)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(lg_ord)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(-4,4)) +
		ylab ("log2(Male FPKM / Female FPKM)") +
		xlab ("Species") + 
		scale_fill_manual(values=c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","darkorange"))  + 
		geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
		ggtitle(tissue)
		
	return(P1)
}

plot_TPM_log2MF_LG <- function(df, tissue){
	df_sub <- subset(df, df$multi_lg == "NO")
	df_sub$lg_ord <- ordered(df_sub$lg, levels = c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8", "lg9", "lg10", "lg11", "lg12", "lgX"))
	
	print(head(df_sub))

	P1 <- ggplot(df_sub, aes(sp, AvTPM_log2MF)) + 
		theme_classic() +
		geom_boxplot(aes(fill = factor(lg_ord)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
		coord_cartesian(ylim=c(-4,4)) +
		ylab ("log2(Male TPM / Female TPM)") +
		xlab ("Species") + 
		scale_fill_manual(values=c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","darkorange"))  + 
		geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
		ggtitle(tissue)
		
	return(P1)
}

### output

pdf("RTHDLG_FPKM_log2MF_soft.pdf", width = 6, height = 9)
plot_grid( 
HD_FPKM_log2MF$P1_soft, 
LG_FPKM_log2MF$P1_soft,
RT_FPKM_log2MF$P1_soft, ncol = 1)
dev.off()

pdf("RTHDLG_FPKM_log2MF_hard.pdf", width = 6, height = 9)
plot_grid( 
HD_FPKM_log2MF$P1_hard, 
LG_FPKM_log2MF$P1_hard,
RT_FPKM_log2MF$P1_hard, ncol = 1)
dev.off()

pdf("RTHDLG_TPM_log2MF_soft.pdf", width = 6, height = 9)
plot_grid(
HD_TPM_log2MF$P1_soft, 
LG_TPM_log2MF$P1_soft,
RT_TPM_log2MF$P1_soft, ncol = 1)
dev.off()

pdf("RTHDLG_TPM_log2MF_hard.pdf", width = 6, height = 9)
plot_grid(
HD_TPM_log2MF$P1_hard, 
LG_TPM_log2MF$P1_hard,
RT_TPM_log2MF$P1_hard, ncol = 1)
dev.off()

pdf("RTHDLG_FPKM_log2MF_LG.pdf", width = 9, height =  13.5)
plot_grid(
plot_FPKM_log2MF_LG(All_HD_F_FPKM_log2MF_long, "HD"),
plot_FPKM_log2MF_LG(All_LG_F_FPKM_log2MF_long, "LG"),
plot_FPKM_log2MF_LG(All_RT_F_FPKM_log2MF_long, "RT"), ncol = 1)
dev.off()

pdf("RTHDLG_TPM_log2MF_LG.pdf", width = 9, height =  13.5)
plot_grid(
plot_TPM_log2MF_LG(All_HD_F_TPM_log2MF_long, "HD"),
plot_TPM_log2MF_LG(All_LG_F_TPM_log2MF_long, "LG"),
plot_TPM_log2MF_LG(All_RT_F_TPM_log2MF_long, "RT"), ncol = 1)
dev.off()




################################################################################################
#### wilcoxon tests on log2MF

wilcox_log2MF <- function(chr_type, exp_type){
	sp_want   <- c("Tbi", "Tce", "Tcm", "Tpa", "Tps")
	tiss_want <- c("RT", "HD", "LG")

	XA_wilcox_log2M_df <- c()
	for(s in sp_want){
		for(tiss in tiss_want){
	
		curr_df <- eval(parse(text=paste('All', '_', tiss, '_F_', exp_type, '_log2MF_long',sep='')))
		
		X_df <- subset(curr_df, eval(parse(text=paste('curr_df' , '$' , 'group_', chr_type, sep = ''))) == paste(s, "_X", sep = ""))
		A_df <- subset(curr_df, eval(parse(text=paste('curr_df' , '$' , 'group_', chr_type, sep = ''))) == paste(s, "_A", sep = ""))
		
		log2MF_XA_wx <- wilcox.test(eval(parse(text=paste('X_df', '$' , 'Av', exp_type, '_log2MF',sep = ''))),   eval(parse(text=paste('A_df', '$' , 'Av', exp_type, '_log2MF',sep = ''))), paired = FALSE)	

		XA_wilcox_log2M_df <- as.data.frame(rbind(
		XA_wilcox_log2M_df,
		c(s, tiss, chr_type, exp_type, log2MF_XA_wx$p.value)
		))
		}
	}
	colnames(XA_wilcox_log2M_df) <- c("sp", "tiss", "chr_class", "exp_type", "wilcox_p")
	XA_wilcox_log2M_df$FDR <- p.adjust(as.numeric(as.character(XA_wilcox_log2M_df$wilcox_p)), method = "BH")
	return(XA_wilcox_log2M_df)
}

#### output

write.csv(wilcox_log2MF("hard", "FPKM"), "XA_wilcox_log2MF_hard_FPKM_df.csv", row.names = F)
write.csv(wilcox_log2MF("soft", "FPKM"), "XA_wilcox_log2MF_soft_FPKM_df.csv", row.names = F)
write.csv(wilcox_log2MF("hard", "TPM"), "XA_wilcox_log2MF_hard_TPM_df.csv", row.names = F)
write.csv(wilcox_log2MF("soft", "TPM"), "XA_wilcox_log2MF_soft_TPM_df.csv", row.names = F)





###############################################################################################################################
##### hists

# species that lack complete dosage compensation, which often show a major peak of genes with strongly reduced expression in the heterogametic sex


FPKM_MF_hists <- function(df1, sp1, tiss1){
	df2 <- subset(df1, df1$sp == sp1)
	
	P1_hard <- ggplot(df2, aes(x=AvFPKM_log2MF, fill=chr_hard)) +
 	      theme_bw() +
 	      geom_density(data=subset(df2,chr_hard == "A"), color="darkblue", fill="blue", alpha = 0, size = 1) 

	P2_hard <- ggplot(df2, aes(x=AvFPKM_log2MF, fill=chr_hard)) +
 	   	theme_bw() +
 	   	geom_density(data=subset(df2,chr_hard == "X"), color="darkred", fill="blue", alpha = 0,  size = 1) 

	
	P1_hard_max <- max(ggplot_build(P1_hard)$data[[1]]$density)
	P2_hard_max <- max(ggplot_build(P2_hard)$data[[1]]$density)
	
	print(P1_hard_max)
	print(P2_hard_max)
	
	Padj_hard <- ifelse(P1_hard_max >= P2_hard_max, P1_hard_max, P2_hard_max)
	
	print(Padj_hard)
	
	Padj_hard1 = Padj_hard * 1.05
	
	P3_hard <- ggplot(df2, aes(x=AvFPKM_log2MF, fill=chr_hard)) +
 	   theme_bw() +
 	   geom_density(data=subset(df2,chr_hard == "A"), color="black", fill="blue", alpha = 0, size = 1) + 
 	   geom_density(data=subset(df2,chr_hard == "X"), color="darkorange", fill="blue", alpha = 0,  size = 1) + 
 	   scale_y_continuous(expand = c(0,0), limits = c(0, Padj_hard1)) + 
 	   scale_x_continuous(limits = c(-6, 6)) +  
 	   geom_vline(xintercept = 0,  linetype = "longdash") +
	    geom_hline(yintercept = 0) +
 	   xlab("log2(Male FPKM / Female FPKM)") + 
		ggtitle(paste(sp1,tiss1, "hard class"))


	P1_soft <- ggplot(df2, aes(x=AvFPKM_log2MF, fill=chr_soft)) +
 	      theme_bw() +
 	      geom_density(data=subset(df2,chr_soft == "A"), color="darkblue", fill="blue", alpha = 0, size = 1) 

	P2_soft <- ggplot(df2, aes(x=AvFPKM_log2MF, fill=chr_soft)) +
 	   	theme_bw() +
 	   	geom_density(data=subset(df2,chr_soft == "X"), color="darkred", fill="blue", alpha = 0,  size = 1) 

	
	P1_soft_max <- max(ggplot_build(P1_soft)$data[[1]]$density)
	P2_soft_max <- max(ggplot_build(P2_soft)$data[[1]]$density)
	
	print(P1_soft_max)
	print(P2_soft_max)
	
	Padj_soft <- ifelse(P1_soft_max >= P2_soft_max, P1_soft_max, P2_soft_max)
	
	print(Padj_soft)

	
	Padj_soft1 = Padj_soft * 1.05
	
	P3_soft <- ggplot(df2, aes(x=AvFPKM_log2MF, fill=chr_soft)) +
 	   theme_bw() +
 	   geom_density(data=subset(df2,chr_soft == "A"), color="black", fill="blue", alpha = 0, size = 1) + 
 	   geom_density(data=subset(df2,chr_soft == "X"), color="darkorange", fill="blue", alpha = 0,  size = 1) + 
 	   scale_y_continuous(expand = c(0,0), limits = c(0, Padj_soft1)) + 
 	   scale_x_continuous(limits = c(-6, 6)) +  
 	   geom_vline(xintercept = 0,  linetype = "longdash") +
	    geom_hline(yintercept = 0) +
 	   xlab("log2(Male FPKM / Female FPKM)") + 
		ggtitle(paste(sp1,tiss1,"soft class"))

	outlist = list("P3_hard" = P3_hard, "P3_soft" = P3_soft) 
	return(outlist)	
		
}



FPKM_log2MF_hist_Tbi_RT <- FPKM_MF_hists(All_RT_F_FPKM_log2MF_long,"Tbi","RT")
FPKM_log2MF_hist_Tbi_HD <- FPKM_MF_hists(All_HD_F_FPKM_log2MF_long,"Tbi","HD")
FPKM_log2MF_hist_Tbi_LG <- FPKM_MF_hists(All_LG_F_FPKM_log2MF_long,"Tbi","LG")

FPKM_log2MF_hist_Tce_RT <- FPKM_MF_hists(All_RT_F_FPKM_log2MF_long,"Tce","RT")
FPKM_log2MF_hist_Tce_HD <- FPKM_MF_hists(All_HD_F_FPKM_log2MF_long,"Tce","HD")
FPKM_log2MF_hist_Tce_LG <- FPKM_MF_hists(All_LG_F_FPKM_log2MF_long,"Tce","LG")

FPKM_log2MF_hist_Tcm_RT <- FPKM_MF_hists(All_RT_F_FPKM_log2MF_long,"Tcm","RT")
FPKM_log2MF_hist_Tcm_HD <- FPKM_MF_hists(All_HD_F_FPKM_log2MF_long,"Tcm","HD")
FPKM_log2MF_hist_Tcm_LG <- FPKM_MF_hists(All_LG_F_FPKM_log2MF_long,"Tcm","LG")

FPKM_log2MF_hist_Tpa_RT <- FPKM_MF_hists(All_RT_F_FPKM_log2MF_long,"Tpa","RT")
FPKM_log2MF_hist_Tpa_HD <- FPKM_MF_hists(All_HD_F_FPKM_log2MF_long,"Tpa","HD")
FPKM_log2MF_hist_Tpa_LG <- FPKM_MF_hists(All_LG_F_FPKM_log2MF_long,"Tpa","LG")

FPKM_log2MF_hist_Tps_RT <- FPKM_MF_hists(All_RT_F_FPKM_log2MF_long,"Tps","RT")
FPKM_log2MF_hist_Tps_HD <- FPKM_MF_hists(All_HD_F_FPKM_log2MF_long,"Tps","HD")
FPKM_log2MF_hist_Tps_LG <- FPKM_MF_hists(All_LG_F_FPKM_log2MF_long,"Tps","LG")


FPKM_log2MF_hist_soft_all <- 
plot_grid(
FPKM_log2MF_hist_Tbi_HD$P3_soft, FPKM_log2MF_hist_Tce_HD$P3_soft, FPKM_log2MF_hist_Tcm_HD$P3_soft, FPKM_log2MF_hist_Tpa_HD$P3_soft, FPKM_log2MF_hist_Tps_HD$P3_soft,
FPKM_log2MF_hist_Tbi_LG$P3_soft, FPKM_log2MF_hist_Tce_LG$P3_soft, FPKM_log2MF_hist_Tcm_LG$P3_soft, FPKM_log2MF_hist_Tpa_LG$P3_soft, FPKM_log2MF_hist_Tps_LG$P3_soft,
FPKM_log2MF_hist_Tbi_RT$P3_soft, FPKM_log2MF_hist_Tce_RT$P3_soft, FPKM_log2MF_hist_Tcm_RT$P3_soft, FPKM_log2MF_hist_Tpa_RT$P3_soft, FPKM_log2MF_hist_Tps_RT$P3_soft,
ncol = 5, nrow = 3)

pdf("FPKM_log2MF_hist_soft_all.pdf", width = 14, height = 9)
FPKM_log2MF_hist_soft_all
dev.off()
getwd() ## where has my plot gone....

FPKM_log2MF_hist_soft_all_2 <- 
plot_grid(
FPKM_log2MF_hist_Tbi_HD$P3_soft, FPKM_log2MF_hist_Tbi_LG$P3_soft, FPKM_log2MF_hist_Tbi_RT$P3_soft, 
FPKM_log2MF_hist_Tce_HD$P3_soft, FPKM_log2MF_hist_Tce_LG$P3_soft, FPKM_log2MF_hist_Tce_RT$P3_soft,
FPKM_log2MF_hist_Tcm_HD$P3_soft, FPKM_log2MF_hist_Tcm_LG$P3_soft, FPKM_log2MF_hist_Tcm_RT$P3_soft, 
FPKM_log2MF_hist_Tpa_HD$P3_soft, FPKM_log2MF_hist_Tpa_LG$P3_soft, FPKM_log2MF_hist_Tpa_RT$P3_soft, 
FPKM_log2MF_hist_Tps_HD$P3_soft, FPKM_log2MF_hist_Tps_LG$P3_soft, FPKM_log2MF_hist_Tps_RT$P3_soft, ncol = 3, nrow = 5)

pdf("FPKM_log2MF_hist_soft_all_2.pdf", width = 9, height = 14)
FPKM_log2MF_hist_soft_all_2
dev.off()
getwd() ## where has my plot gone....

FPKM_log2MF_hist_hard_all <- 
plot_grid(
FPKM_log2MF_hist_Tbi_HD$P3_hard, FPKM_log2MF_hist_Tce_HD$P3_hard, FPKM_log2MF_hist_Tcm_HD$P3_hard, FPKM_log2MF_hist_Tpa_HD$P3_hard, FPKM_log2MF_hist_Tps_HD$P3_hard,
FPKM_log2MF_hist_Tbi_LG$P3_hard, FPKM_log2MF_hist_Tce_LG$P3_hard, FPKM_log2MF_hist_Tcm_LG$P3_hard, FPKM_log2MF_hist_Tpa_LG$P3_hard, FPKM_log2MF_hist_Tps_LG$P3_hard,
FPKM_log2MF_hist_Tbi_RT$P3_hard, FPKM_log2MF_hist_Tce_RT$P3_hard, FPKM_log2MF_hist_Tcm_RT$P3_hard, FPKM_log2MF_hist_Tpa_RT$P3_hard, FPKM_log2MF_hist_Tps_RT$P3_hard,
ncol = 5, nrow = 3)

pdf("FPKM_log2MF_hist_hard_all.pdf", width = 14, height = 9)
FPKM_log2MF_hist_hard_all
dev.off()
getwd() ## where has my plot gone....

FPKM_log2MF_hist_hard_all_2 <- 
plot_grid(
FPKM_log2MF_hist_Tbi_HD$P3_hard, FPKM_log2MF_hist_Tbi_LG$P3_hard, FPKM_log2MF_hist_Tbi_RT$P3_hard, 
FPKM_log2MF_hist_Tce_HD$P3_hard, FPKM_log2MF_hist_Tce_LG$P3_hard, FPKM_log2MF_hist_Tce_RT$P3_hard,
FPKM_log2MF_hist_Tcm_HD$P3_hard, FPKM_log2MF_hist_Tcm_LG$P3_hard, FPKM_log2MF_hist_Tcm_RT$P3_hard, 
FPKM_log2MF_hist_Tpa_HD$P3_hard, FPKM_log2MF_hist_Tpa_LG$P3_hard, FPKM_log2MF_hist_Tpa_RT$P3_hard, 
FPKM_log2MF_hist_Tps_HD$P3_hard, FPKM_log2MF_hist_Tps_LG$P3_hard, FPKM_log2MF_hist_Tps_RT$P3_hard, ncol = 3, nrow = 5)

pdf("FPKM_log2MF_hist_hard_all_2.pdf", width = 9, height = 14)
FPKM_log2MF_hist_hard_all_2
dev.off()
getwd() ## where has my plot gone....


TPM_MF_hists <- function(df1, sp1, tiss1){
	df2 <- subset(df1, df1$sp == sp1)
	
	P1_hard <- ggplot(df2, aes(x=AvTPM_log2MF, fill=chr_hard)) +
 	      theme_bw() +
 	      geom_density(data=subset(df2,chr_hard == "A"), color="darkblue", fill="blue", alpha = 0, size = 1) 

	P2_hard <- ggplot(df2, aes(x=AvTPM_log2MF, fill=chr_hard)) +
 	   	theme_bw() +
 	   	geom_density(data=subset(df2,chr_hard == "X"), color="darkred", fill="blue", alpha = 0,  size = 1) 

	
	P1_hard_max <- max(ggplot_build(P1_hard)$data[[1]]$density)
	P2_hard_max <- max(ggplot_build(P2_hard)$data[[1]]$density)
	
	print(P1_hard_max)
	print(P2_hard_max)
	
	Padj_hard <- ifelse(P1_hard_max >= P2_hard_max, P1_hard_max, P2_hard_max)
	
	print(Padj_hard)
	
	Padj_hard1 = Padj_hard * 1.05
	
	P3_hard <- ggplot(df2, aes(x=AvTPM_log2MF, fill=chr_hard)) +
 	   theme_bw() +
 	   geom_density(data=subset(df2,chr_hard == "A"), color="black", fill="blue", alpha = 0, size = 1) + 
 	   geom_density(data=subset(df2,chr_hard == "X"), color="darkorange", fill="blue", alpha = 0,  size = 1) + 
 	   scale_y_continuous(expand = c(0,0), limits = c(0, Padj_hard1)) + 
 	   scale_x_continuous(limits = c(-6, 6)) +  
 	   geom_vline(xintercept = 0,  linetype = "longdash") +
	    geom_hline(yintercept = 0) +
 	   xlab("log2(Male TPM / Female TPM)") + 
		ggtitle(paste(sp1,tiss1, "hard class"))


	P1_soft <- ggplot(df2, aes(x=AvTPM_log2MF, fill=chr_soft)) +
 	      theme_bw() +
 	      geom_density(data=subset(df2,chr_soft == "A"), color="darkblue", fill="blue", alpha = 0, size = 1) 

	P2_soft <- ggplot(df2, aes(x=AvTPM_log2MF, fill=chr_soft)) +
 	   	theme_bw() +
 	   	geom_density(data=subset(df2,chr_soft == "X"), color="darkred", fill="blue", alpha = 0,  size = 1) 

	
	P1_soft_max <- max(ggplot_build(P1_soft)$data[[1]]$density)
	P2_soft_max <- max(ggplot_build(P2_soft)$data[[1]]$density)
	
	print(P1_soft_max)
	print(P2_soft_max)
	
	Padj_soft <- ifelse(P1_soft_max >= P2_soft_max, P1_soft_max, P2_soft_max)
	
	print(Padj_soft)

	
	Padj_soft1 = Padj_soft * 1.05
	
	P3_soft <- ggplot(df2, aes(x=AvTPM_log2MF, fill=chr_soft)) +
 	   theme_bw() +
 	   geom_density(data=subset(df2,chr_soft == "A"), color="black", fill="blue", alpha = 0, size = 1) + 
 	   geom_density(data=subset(df2,chr_soft == "X"), color="darkorange", fill="blue", alpha = 0,  size = 1) + 
 	   scale_y_continuous(expand = c(0,0), limits = c(0, Padj_soft1)) + 
 	   scale_x_continuous(limits = c(-6, 6)) +  
 	   geom_vline(xintercept = 0,  linetype = "longdash") +
	    geom_hline(yintercept = 0) +
 	   xlab("log2(Male TPM / Female TPM)") + 
		ggtitle(paste(sp1,tiss1,"soft class"))

	outlist = list("P3_hard" = P3_hard, "P3_soft" = P3_soft) 
	return(outlist)	
		
}


TPM_log2MF_hist_Tbi_RT <- TPM_MF_hists(All_RT_F_TPM_log2MF_long,"Tbi","RT")
TPM_log2MF_hist_Tbi_HD <- TPM_MF_hists(All_HD_F_TPM_log2MF_long,"Tbi","HD")
TPM_log2MF_hist_Tbi_LG <- TPM_MF_hists(All_LG_F_TPM_log2MF_long,"Tbi","LG")

TPM_log2MF_hist_Tce_RT <- TPM_MF_hists(All_RT_F_TPM_log2MF_long,"Tce","RT")
TPM_log2MF_hist_Tce_HD <- TPM_MF_hists(All_HD_F_TPM_log2MF_long,"Tce","HD")
TPM_log2MF_hist_Tce_LG <- TPM_MF_hists(All_LG_F_TPM_log2MF_long,"Tce","LG")
)
TPM_log2MF_hist_Tcm_RT <- TPM_MF_hists(All_RT_F_TPM_log2MF_long,"Tcm","RT")
TPM_log2MF_hist_Tcm_HD <- TPM_MF_hists(All_HD_F_TPM_log2MF_long,"Tcm","HD")
TPM_log2MF_hist_Tcm_LG <- TPM_MF_hists(All_LG_F_TPM_log2MF_long,"Tcm","LG")

TPM_log2MF_hist_Tpa_RT <- TPM_MF_hists(All_RT_F_TPM_log2MF_long,"Tpa","RT")
TPM_log2MF_hist_Tpa_HD <- TPM_MF_hists(All_HD_F_TPM_log2MF_long,"Tpa","HD")
TPM_log2MF_hist_Tpa_LG <- TPM_MF_hists(All_LG_F_TPM_log2MF_long,"Tpa","LG")

TPM_log2MF_hist_Tps_RT <- TPM_MF_hists(All_RT_F_TPM_log2MF_long,"Tps","RT")
TPM_log2MF_hist_Tps_HD <- TPM_MF_hists(All_HD_F_TPM_log2MF_long,"Tps","HD")
TPM_log2MF_hist_Tps_LG <- TPM_MF_hists(All_LG_F_TPM_log2MF_long,"Tps","LG")


TPM_log2MF_hist_soft_all <- 
plot_grid(
TPM_log2MF_hist_Tbi_HD$P3_soft, TPM_log2MF_hist_Tce_HD$P3_soft, TPM_log2MF_hist_Tcm_HD$P3_soft, TPM_log2MF_hist_Tpa_HD$P3_soft, TPM_log2MF_hist_Tps_HD$P3_soft,
TPM_log2MF_hist_Tbi_LG$P3_soft, TPM_log2MF_hist_Tce_LG$P3_soft, TPM_log2MF_hist_Tcm_LG$P3_soft, TPM_log2MF_hist_Tpa_LG$P3_soft, TPM_log2MF_hist_Tps_LG$P3_soft,
TPM_log2MF_hist_Tbi_RT$P3_soft, TPM_log2MF_hist_Tce_RT$P3_soft, TPM_log2MF_hist_Tcm_RT$P3_soft, TPM_log2MF_hist_Tpa_RT$P3_soft, TPM_log2MF_hist_Tps_RT$P3_soft,
ncol = 5, nrow = 3)

pdf("TPM_log2MF_hist_soft_all.pdf", width = 14, height = 9)
TPM_log2MF_hist_soft_all
dev.off()
getwd() ## where has my plot gone....

TPM_log2MF_hist_soft_all_2 <- 
plot_grid(
TPM_log2MF_hist_Tbi_HD$P3_soft, TPM_log2MF_hist_Tbi_LG$P3_soft, TPM_log2MF_hist_Tbi_RT$P3_soft, 
TPM_log2MF_hist_Tce_HD$P3_soft, TPM_log2MF_hist_Tce_LG$P3_soft, TPM_log2MF_hist_Tce_RT$P3_soft,
TPM_log2MF_hist_Tcm_HD$P3_soft, TPM_log2MF_hist_Tcm_LG$P3_soft, TPM_log2MF_hist_Tcm_RT$P3_soft, 
TPM_log2MF_hist_Tpa_HD$P3_soft, TPM_log2MF_hist_Tpa_LG$P3_soft, TPM_log2MF_hist_Tpa_RT$P3_soft, 
TPM_log2MF_hist_Tps_HD$P3_soft, TPM_log2MF_hist_Tps_LG$P3_soft, TPM_log2MF_hist_Tps_RT$P3_soft, ncol = 3, nrow = 5)

pdf("TPM_log2MF_hist_soft_all_2.pdf", width = 9, height = 14)
TPM_log2MF_hist_soft_all_2
dev.off()
getwd() ## where has my plot gone....

TPM_log2MF_hist_hard_all <- 
plot_grid(
TPM_log2MF_hist_Tbi_HD$P3_hard, TPM_log2MF_hist_Tce_HD$P3_hard, TPM_log2MF_hist_Tcm_HD$P3_hard, TPM_log2MF_hist_Tpa_HD$P3_hard, TPM_log2MF_hist_Tps_HD$P3_hard,
TPM_log2MF_hist_Tbi_LG$P3_hard, TPM_log2MF_hist_Tce_LG$P3_hard, TPM_log2MF_hist_Tcm_LG$P3_hard, TPM_log2MF_hist_Tpa_LG$P3_hard, TPM_log2MF_hist_Tps_LG$P3_hard,
TPM_log2MF_hist_Tbi_RT$P3_hard, TPM_log2MF_hist_Tce_RT$P3_hard, TPM_log2MF_hist_Tcm_RT$P3_hard, TPM_log2MF_hist_Tpa_RT$P3_hard, TPM_log2MF_hist_Tps_RT$P3_hard,
ncol = 5, nrow = 3)

pdf("TPM_log2MF_hist_hard_all.pdf", width = 14, height = 9)
TPM_log2MF_hist_hard_all
dev.off()
getwd() ## where has my plot gone....

TPM_log2MF_hist_hard_all_2 <- 
plot_grid(
TPM_log2MF_hist_Tbi_HD$P3_hard, TPM_log2MF_hist_Tbi_LG$P3_hard, TPM_log2MF_hist_Tbi_RT$P3_hard, 
TPM_log2MF_hist_Tce_HD$P3_hard, TPM_log2MF_hist_Tce_LG$P3_hard, TPM_log2MF_hist_Tce_RT$P3_hard,
TPM_log2MF_hist_Tcm_HD$P3_hard, TPM_log2MF_hist_Tcm_LG$P3_hard, TPM_log2MF_hist_Tcm_RT$P3_hard, 
TPM_log2MF_hist_Tpa_HD$P3_hard, TPM_log2MF_hist_Tpa_LG$P3_hard, TPM_log2MF_hist_Tpa_RT$P3_hard, 
TPM_log2MF_hist_Tps_HD$P3_hard, TPM_log2MF_hist_Tps_LG$P3_hard, TPM_log2MF_hist_Tps_RT$P3_hard, ncol = 3, nrow = 5)

pdf("TPM_log2MF_hist_hard_all_2.pdf", width = 9, height = 14)
TPM_log2MF_hist_hard_all_2
dev.off()
getwd() ## where has my plot gone....


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### GET SB genes

### get DE_gene_table (code)

## returns full table and sig DE genes as a vector
get_DE_genes <- function(fita,FDRa){
	TT1 = topTags(fita, n =3000000000)
	TT2 = TT1$table
	temp_sig <- subset(TT2, TT2$FDR <= FDRa)	
	sig_genes = temp_sig$genes
	sig_logFC = temp_sig$logFC
	
	N_sig_genes <- length(sig_genes)
	cat("Number of sig genes: ", N_sig_genes )
	
	r_list <- list("table" = TT2, "S_gene_list" = sig_genes, "S_logFC_list" = sig_logFC )
	return(r_list)
}

####### design matrix

## using 0M for male and 1F for female samples

sex_SB = factor(c(
"1F","1F","1F","0M","0M","0M"
))

data.frame(Sample=colnames(y_Tbi_WB_F ),sex_SB)


call_SB <- function(y, sig_level){
	
	### design matrix
	design <- model.matrix(~sex_SB) 
	rownames(design) <- colnames(y)
	print(design)

	### Est dispersion 
	y <- estimateDisp(y, design)
	
	### Fit model
	fit_y <- glmFit(y, design,robust=TRUE)
	
	### Test
	print(colnames(fit_y))
	fit_c <-  glmLRT(fit_y,coef=2)
	
	### get_SB gene
	TTT_SB <- get_DE_genes(fit_c, sig_level)
	
	return(TTT_SB)
}


call_SB <- function(y, sig_level){
	
	### design matrix
	design <- model.matrix(~sex_SB) 
	rownames(design) <- colnames(y)
	print(design)

	### Est dispersion 
	y <- estimateDisp(y, design)
	
	### Fit model
	fit_y <- glmFit(y, design,robust=TRUE)
	
	### Test
	print(colnames(fit_y))
	fit_c <-  glmLRT(fit_y,coef=2)
	
	### get_SB gene
	TTT_SB <- get_DE_genes(fit_c, sig_level)
	
	return(TTT_SB)
}

TTT_RT_Tbi_sex_bias_FPKM <- call_SB(y_Tbi_RT_F_FPKM, 0.05)
TTT_RT_Tce_sex_bias_FPKM <- call_SB(y_Tce_RT_F_FPKM, 0.05)
TTT_RT_Tcm_sex_bias_FPKM <- call_SB(y_Tcm_RT_F_FPKM, 0.05)
TTT_RT_Tpa_sex_bias_FPKM <- call_SB(y_Tpa_RT_F_FPKM, 0.05)
TTT_RT_Tps_sex_bias_FPKM <- call_SB(y_Tps_RT_F_FPKM, 0.05)

TTT_HD_Tbi_sex_bias_FPKM <- call_SB(y_Tbi_HD_F_FPKM, 0.05)
TTT_HD_Tce_sex_bias_FPKM <- call_SB(y_Tce_HD_F_FPKM, 0.05)
TTT_HD_Tcm_sex_bias_FPKM <- call_SB(y_Tcm_HD_F_FPKM, 0.05)
TTT_HD_Tpa_sex_bias_FPKM <- call_SB(y_Tpa_HD_F_FPKM, 0.05)
TTT_HD_Tps_sex_bias_FPKM <- call_SB(y_Tps_HD_F_FPKM, 0.05)

TTT_LG_Tbi_sex_bias_FPKM <- call_SB(y_Tbi_LG_F_FPKM, 0.05)
TTT_LG_Tce_sex_bias_FPKM <- call_SB(y_Tce_LG_F_FPKM, 0.05)
TTT_LG_Tcm_sex_bias_FPKM <- call_SB(y_Tcm_LG_F_FPKM, 0.05)
TTT_LG_Tpa_sex_bias_FPKM <- call_SB(y_Tpa_LG_F_FPKM, 0.05)
TTT_LG_Tps_sex_bias_FPKM <- call_SB(y_Tps_LG_F_FPKM, 0.05)


N_SB_genes_FPKM <- as.data.frame(cbind(
rep(c("Tbi", "Tce", "Tcm", "Tpa", "Tps"), 3),
c(rep("RT", 5), rep("HD", 5), rep("LG", 5)),
c(

length(TTT_RT_Tbi_sex_bias_FPKM$S_gene_list),
length(TTT_RT_Tce_sex_bias_FPKM$S_gene_list),
length(TTT_RT_Tcm_sex_bias_FPKM$S_gene_list),
length(TTT_RT_Tpa_sex_bias_FPKM$S_gene_list),
length(TTT_RT_Tps_sex_bias_FPKM$S_gene_list),

length(TTT_HD_Tbi_sex_bias_FPKM$S_gene_list),
length(TTT_HD_Tce_sex_bias_FPKM$S_gene_list),
length(TTT_HD_Tcm_sex_bias_FPKM$S_gene_list),
length(TTT_HD_Tpa_sex_bias_FPKM$S_gene_list),
length(TTT_HD_Tps_sex_bias_FPKM$S_gene_list),

length(TTT_LG_Tbi_sex_bias_FPKM$S_gene_list),
length(TTT_LG_Tce_sex_bias_FPKM$S_gene_list),
length(TTT_LG_Tcm_sex_bias_FPKM$S_gene_list),
length(TTT_LG_Tpa_sex_bias_FPKM$S_gene_list),
length(TTT_LG_Tps_sex_bias_FPKM$S_gene_list))))

colnames(N_SB_genes_FPKM ) <- c("sp", "tiss", "N_SB")



TTT_RT_Tbi_sex_bias_TPM <- call_SB(y_Tbi_RT_F_TPM, 0.05)
TTT_RT_Tce_sex_bias_TPM <- call_SB(y_Tce_RT_F_TPM, 0.05)
TTT_RT_Tcm_sex_bias_TPM <- call_SB(y_Tcm_RT_F_TPM, 0.05)
TTT_RT_Tpa_sex_bias_TPM <- call_SB(y_Tpa_RT_F_TPM, 0.05)
TTT_RT_Tps_sex_bias_TPM <- call_SB(y_Tps_RT_F_TPM, 0.05)

TTT_HD_Tbi_sex_bias_TPM <- call_SB(y_Tbi_HD_F_TPM, 0.05)
TTT_HD_Tce_sex_bias_TPM <- call_SB(y_Tce_HD_F_TPM, 0.05)
TTT_HD_Tcm_sex_bias_TPM <- call_SB(y_Tcm_HD_F_TPM, 0.05)
TTT_HD_Tpa_sex_bias_TPM <- call_SB(y_Tpa_HD_F_TPM, 0.05)
TTT_HD_Tps_sex_bias_TPM <- call_SB(y_Tps_HD_F_TPM, 0.05)

TTT_LG_Tbi_sex_bias_TPM <- call_SB(y_Tbi_LG_F_TPM, 0.05)
TTT_LG_Tce_sex_bias_TPM <- call_SB(y_Tce_LG_F_TPM, 0.05)
TTT_LG_Tcm_sex_bias_TPM <- call_SB(y_Tcm_LG_F_TPM, 0.05)
TTT_LG_Tpa_sex_bias_TPM <- call_SB(y_Tpa_LG_F_TPM, 0.05)
TTT_LG_Tps_sex_bias_TPM <- call_SB(y_Tps_LG_F_TPM, 0.05)


N_SB_genes_TPM <- as.data.frame(cbind(
rep(c("Tbi", "Tce", "Tcm", "Tpa", "Tps"), 3),
c(rep("RT", 5), rep("HD", 5), rep("LG", 5)),
c(

length(TTT_RT_Tbi_sex_bias_TPM$S_gene_list),
length(TTT_RT_Tce_sex_bias_TPM$S_gene_list),
length(TTT_RT_Tcm_sex_bias_TPM$S_gene_list),
length(TTT_RT_Tpa_sex_bias_TPM$S_gene_list),
length(TTT_RT_Tps_sex_bias_TPM$S_gene_list),

length(TTT_HD_Tbi_sex_bias_TPM$S_gene_list),
length(TTT_HD_Tce_sex_bias_TPM$S_gene_list),
length(TTT_HD_Tcm_sex_bias_TPM$S_gene_list),
length(TTT_HD_Tpa_sex_bias_TPM$S_gene_list),
length(TTT_HD_Tps_sex_bias_TPM$S_gene_list),

length(TTT_LG_Tbi_sex_bias_TPM$S_gene_list),
length(TTT_LG_Tce_sex_bias_TPM$S_gene_list),
length(TTT_LG_Tcm_sex_bias_TPM$S_gene_list),
length(TTT_LG_Tpa_sex_bias_TPM$S_gene_list),
length(TTT_LG_Tps_sex_bias_TPM$S_gene_list))))

colnames(N_SB_genes_TPM ) <- c("sp", "tiss", "N_SB")


## set  FC_th to -1 for no FC thresh
plot_smear <- function(expr_df, expr_df_name, exp_name, DE_genes_tab, FDR_th, tit_txt, exp_type, max_FC, FC_th){

	DE_genes_df <- DE_genes_tab

	### get exp vals
	gene_expr_dict <- hash()

	for(i in seq(1:length(expr_df[,1]))){
		gene_n <- expr_df$gene_id[i]
		expr_n <- eval(parse(text=paste(expr_df_name,'$',exp_name,sep='')))[i]
		gene_expr_dict[[gene_n]] <- expr_n
		}

	want_exp <- c()
	for(i in seq(1:length(DE_genes_df[,1]))){
		gene_n <- DE_genes_df$genes[i]
		expr_n <- gene_expr_dict[[gene_n]]
		#### returns a NULL if key not found. Turn into NA
        if(length(expr_n) == 0){
			expr_n = NA
			}
		want_exp <- c(want_exp, expr_n)
		}
	
	DE_genes_df$log_expr <- log2(want_exp)
	DE_genes_df$logFC_abs <- sqrt(DE_genes_df$logFC * DE_genes_df$logFC)

	DE_genes_df$sig <- ifelse(DE_genes_df$FDR <= FDR_th & DE_genes_df$logFC_abs >= FC_th, "sig", "non-sig")
	
	print(head(DE_genes_df))
	
	DE_genes_df_2 <- 
	rbind(
	subset(DE_genes_df, DE_genes_df$sig == "non-sig"),
	subset(DE_genes_df, DE_genes_df$sig == "sig"))

	min_FC <- max_FC * -1
	
	P1 <- ggplot(DE_genes_df_2 , aes(x=log_expr,logFC)) + 
		geom_point(aes(colour = factor(sig)), size = 0.5) +
		theme_bw() +
		scale_color_manual(values=c( "gray30", "red3")) + 
		theme(legend.position="none") +
		coord_cartesian(ylim = c(min_FC,max_FC)) +
		ggtitle(tit_txt) + xlab(paste("log2", exp_type)) 
	
	return(P1)	
}

plot_smear_XA <- function(expr_df, expr_df_name, exp_name, DE_genes_tab, FDR_th, tit_txt, exp_type, max_FC, chr_type, FC_th){

	DE_genes_df <- DE_genes_tab

	### get exp vals
	gene_expr_dict     <- hash()
	gene_chr_dict <- hash()

	for(i in seq(1:length(expr_df[,1]))){
		gene_n <- expr_df$gene_id[i]
		expr_n <- eval(parse(text=paste(expr_df_name,'$',exp_name,sep='')))[i]
		chr_n <- eval(parse(text=paste(expr_df_name,'$chr_', chr_type,sep='')))[i]
		gene_expr_dict[[gene_n]] <- expr_n
		gene_chr_dict[[gene_n]] <- chr_n
		}


	want_exp <- c()
	chr_v <- c()
	for(i in seq(1:length(DE_genes_df[,1]))){
		gene_n <- DE_genes_df$genes[i]
		expr_n <- gene_expr_dict[[gene_n]]
		chr_n <- gene_chr_dict[[gene_n]]
		#### returns a NULL if key not found. Turn into NA
        if(length(expr_n) == 0){
			expr_n = NA
			}
		#### returns a NULL if key not found. Turn into NA
        if(length(chr_n) == 0){
			chr_n = NA
			}
		want_exp <- c(want_exp, expr_n)
		chr_v <- c(chr_v, chr_n)
		}
	
	DE_genes_df$log_expr <- log2(want_exp)
	DE_genes_df$logFC_abs <- sqrt(DE_genes_df$logFC * DE_genes_df$logFC)

	DE_genes_df$sig <- ifelse(DE_genes_df$FDR <= FDR_th & DE_genes_df$logFC_abs >= FC_th, "sig", "non-sig")
	DE_genes_df$chr_type  <- chr_v
	
	
	DE_genes_df_2 <- 
	rbind(
	subset(DE_genes_df, DE_genes_df$sig == "non-sig"),
	subset(DE_genes_df, DE_genes_df$sig == "sig"))
	

	min_FC <- max_FC * -1
	
	DE_genes_df_2_X <- subset(DE_genes_df_2, DE_genes_df_2$chr_type == "X")
	DE_genes_df_2_A <- subset(DE_genes_df_2, DE_genes_df_2$chr_type == "A")	
	
	P1_X <- ggplot(DE_genes_df_2_X , aes(x=log_expr,logFC)) + 
		geom_point(aes(colour = factor(sig)), size = 0.5) +
		theme_bw() +
		scale_color_manual(values=c( "gray30", "red3")) + 
		theme(legend.position="none") +
		coord_cartesian(ylim = c(min_FC,max_FC)) +
		ggtitle(paste(tit_txt,  chr_type, "X"))  + xlab(paste("log2", exp_type, chr_type)) 	

	P1_A <- ggplot(DE_genes_df_2_A , aes(x=log_expr,logFC)) + 
		geom_point(aes(colour = factor(sig)), size = 0.5) +
		theme_bw() +
		scale_color_manual(values=c( "gray30", "red3")) + 
		theme(legend.position="none") +
		coord_cartesian(ylim = c(min_FC,max_FC)) +
		ggtitle(paste(tit_txt,  chr_type, "A")) + xlab(paste("log2", exp_type)) 	
		
	P1_XA <- plot_grid(P1_A, P1_X)
	return(P1_XA)
	
}


RT_max_FC <- ceiling(max(abs(c(TTT_RT_Tbi_sex_bias_FPKM$table$logFC, TTT_RT_Tce_sex_bias_FPKM$table$logFC, TTT_RT_Tcm_sex_bias_FPKM$table$logFC,TTT_RT_Tpa_sex_bias_FPKM$table$logFC,TTT_RT_Tps_sex_bias_FPKM$table$logFC))))
SB_RT_FPKM_smear <- plot_grid(
plot_smear(Tbi_RT_F_FPKM, "Tbi_RT_F_FPKM", "Tbi_RT_meanFPKM", TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, "Tbi RT", "FPKM", RT_max_FC, -1),
plot_smear(Tce_RT_F_FPKM, "Tce_RT_F_FPKM", "Tce_RT_meanFPKM", TTT_RT_Tce_sex_bias_FPKM$table, 0.05, "Tce RT", "FPKM", RT_max_FC, -1),
plot_smear(Tcm_RT_F_FPKM, "Tcm_RT_F_FPKM", "Tcm_RT_meanFPKM", TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, "Tcm RT", "FPKM", RT_max_FC, -1),
plot_smear(Tpa_RT_F_FPKM, "Tpa_RT_F_FPKM", "Tpa_RT_meanFPKM", TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, "Tpa RT", "FPKM", RT_max_FC, -1),
plot_smear(Tps_RT_F_FPKM, "Tps_RT_F_FPKM", "Tps_RT_meanFPKM", TTT_RT_Tps_sex_bias_FPKM$table, 0.05, "Tps RT", "FPKM", RT_max_FC, -1), ncol = 2)

HD_max_FC <- ceiling(max(abs(c(TTT_HD_Tbi_sex_bias_FPKM$table$logFC, TTT_HD_Tce_sex_bias_FPKM$table$logFC, TTT_HD_Tcm_sex_bias_FPKM$table$logFC,TTT_HD_Tpa_sex_bias_FPKM$table$logFC,TTT_HD_Tps_sex_bias_FPKM$table$logFC))))
SB_HD_FPKM_smear <- plot_grid(
plot_smear(Tbi_HD_F_FPKM, "Tbi_HD_F_FPKM", "Tbi_HD_meanFPKM", TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, "Tbi HD", "FPKM", HD_max_FC, -1),
plot_smear(Tce_HD_F_FPKM, "Tce_HD_F_FPKM", "Tce_HD_meanFPKM", TTT_HD_Tce_sex_bias_FPKM$table, 0.05, "Tce HD", "FPKM", HD_max_FC, -1),
plot_smear(Tcm_HD_F_FPKM, "Tcm_HD_F_FPKM", "Tcm_HD_meanFPKM", TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, "Tcm HD", "FPKM", HD_max_FC, -1),
plot_smear(Tpa_HD_F_FPKM, "Tpa_HD_F_FPKM", "Tpa_HD_meanFPKM", TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, "Tpa HD", "FPKM", HD_max_FC, -1),
plot_smear(Tps_HD_F_FPKM, "Tps_HD_F_FPKM", "Tps_HD_meanFPKM", TTT_HD_Tps_sex_bias_FPKM$table, 0.05, "Tps HD", "FPKM", HD_max_FC, -1), ncol = 2)

LG_max_FC <- ceiling(max(abs(c(TTT_LG_Tbi_sex_bias_FPKM$table$logFC, TTT_LG_Tce_sex_bias_FPKM$table$logFC, TTT_LG_Tcm_sex_bias_FPKM$table$logFC,TTT_LG_Tpa_sex_bias_FPKM$table$logFC,TTT_LG_Tps_sex_bias_FPKM$table$logFC))))
SB_LG_FPKM_smear <- plot_grid(
plot_smear(Tbi_LG_F_FPKM, "Tbi_LG_F_FPKM", "Tbi_LG_meanFPKM", TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, "Tbi LG", "FPKM", LG_max_FC, -1),
plot_smear(Tce_LG_F_FPKM, "Tce_LG_F_FPKM", "Tce_LG_meanFPKM", TTT_LG_Tce_sex_bias_FPKM$table, 0.05, "Tce LG", "FPKM", LG_max_FC, -1),
plot_smear(Tcm_LG_F_FPKM, "Tcm_LG_F_FPKM", "Tcm_LG_meanFPKM", TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, "Tcm LG", "FPKM", LG_max_FC, -1),
plot_smear(Tpa_LG_F_FPKM, "Tpa_LG_F_FPKM", "Tpa_LG_meanFPKM", TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, "Tpa LG", "FPKM", LG_max_FC, -1),
plot_smear(Tps_LG_F_FPKM, "Tps_LG_F_FPKM", "Tps_LG_meanFPKM", TTT_LG_Tps_sex_bias_FPKM$table, 0.05, "Tps LG", "FPKM", LG_max_FC, -1), ncol = 2)

pdf("SB_RT_FPKM_smear.pdf",  width = 10, height = 15)
SB_RT_FPKM_smear
dev.off()
getwd() ## where has my plot gone....

pdf("SB_HD_FPKM_smear.pdf",  width = 10, height = 15)
SB_HD_FPKM_smear
dev.off()
getwd() ## where has my plot gone....

pdf("SB_LG_FPKM_smear.pdf",  width = 10, height = 15)
SB_LG_FPKM_smear
dev.off()
getwd() ## where has my plot gone....



RT_max_FC <- ceiling(max(abs(c(TTT_RT_Tbi_sex_bias_TPM$table$logFC, TTT_RT_Tce_sex_bias_TPM$table$logFC, TTT_RT_Tcm_sex_bias_TPM$table$logFC,TTT_RT_Tpa_sex_bias_TPM$table$logFC,TTT_RT_Tps_sex_bias_TPM$table$logFC))))
SB_RT_TPM_smear <- plot_grid(
plot_smear(Tbi_RT_F_TPM, "Tbi_RT_F_TPM", "Tbi_RT_meanTPM", TTT_RT_Tbi_sex_bias_TPM$table, 0.05, "Tbi RT", "TPM", RT_max_FC, -1),
plot_smear(Tce_RT_F_TPM, "Tce_RT_F_TPM", "Tce_RT_meanTPM", TTT_RT_Tce_sex_bias_TPM$table, 0.05, "Tce RT", "TPM", RT_max_FC, -1),
plot_smear(Tcm_RT_F_TPM, "Tcm_RT_F_TPM", "Tcm_RT_meanTPM", TTT_RT_Tcm_sex_bias_TPM$table, 0.05, "Tcm RT", "TPM", RT_max_FC, -1),
plot_smear(Tpa_RT_F_TPM, "Tpa_RT_F_TPM", "Tpa_RT_meanTPM", TTT_RT_Tpa_sex_bias_TPM$table, 0.05, "Tpa RT", "TPM", RT_max_FC, -1),
plot_smear(Tps_RT_F_TPM, "Tps_RT_F_TPM", "Tps_RT_meanTPM", TTT_RT_Tps_sex_bias_TPM$table, 0.05, "Tps RT", "TPM", RT_max_FC, -1), ncol = 2)

HD_max_FC <- ceiling(max(abs(c(TTT_HD_Tbi_sex_bias_TPM$table$logFC, TTT_HD_Tce_sex_bias_TPM$table$logFC, TTT_HD_Tcm_sex_bias_TPM$table$logFC,TTT_HD_Tpa_sex_bias_TPM$table$logFC,TTT_HD_Tps_sex_bias_TPM$table$logFC))))
SB_HD_TPM_smear <- plot_grid(
plot_smear(Tbi_HD_F_TPM, "Tbi_HD_F_TPM", "Tbi_HD_meanTPM", TTT_HD_Tbi_sex_bias_TPM$table, 0.05, "Tbi HD", "TPM", HD_max_FC, -1),
plot_smear(Tce_HD_F_TPM, "Tce_HD_F_TPM", "Tce_HD_meanTPM", TTT_HD_Tce_sex_bias_TPM$table, 0.05, "Tce HD", "TPM", HD_max_FC, -1),
plot_smear(Tcm_HD_F_TPM, "Tcm_HD_F_TPM", "Tcm_HD_meanTPM", TTT_HD_Tcm_sex_bias_TPM$table, 0.05, "Tcm HD", "TPM", HD_max_FC, -1),
plot_smear(Tpa_HD_F_TPM, "Tpa_HD_F_TPM", "Tpa_HD_meanTPM", TTT_HD_Tpa_sex_bias_TPM$table, 0.05, "Tpa HD", "TPM", HD_max_FC, -1),
plot_smear(Tps_HD_F_TPM, "Tps_HD_F_TPM", "Tps_HD_meanTPM", TTT_HD_Tps_sex_bias_TPM$table, 0.05, "Tps HD", "TPM", HD_max_FC, -1), ncol = 2)

LG_max_FC <- ceiling(max(abs(c(TTT_LG_Tbi_sex_bias_TPM$table$logFC, TTT_LG_Tce_sex_bias_TPM$table$logFC, TTT_LG_Tcm_sex_bias_TPM$table$logFC,TTT_LG_Tpa_sex_bias_TPM$table$logFC,TTT_LG_Tps_sex_bias_TPM$table$logFC))))
SB_LG_TPM_smear <- plot_grid(
plot_smear(Tbi_LG_F_TPM, "Tbi_LG_F_TPM", "Tbi_LG_meanTPM", TTT_LG_Tbi_sex_bias_TPM$table, 0.05, "Tbi LG", "TPM", LG_max_FC, -1),
plot_smear(Tce_LG_F_TPM, "Tce_LG_F_TPM", "Tce_LG_meanTPM", TTT_LG_Tce_sex_bias_TPM$table, 0.05, "Tce LG", "TPM", LG_max_FC, -1),
plot_smear(Tcm_LG_F_TPM, "Tcm_LG_F_TPM", "Tcm_LG_meanTPM", TTT_LG_Tcm_sex_bias_TPM$table, 0.05, "Tcm LG", "TPM", LG_max_FC, -1),
plot_smear(Tpa_LG_F_TPM, "Tpa_LG_F_TPM", "Tpa_LG_meanTPM", TTT_LG_Tpa_sex_bias_TPM$table, 0.05, "Tpa LG", "TPM", LG_max_FC, -1),
plot_smear(Tps_LG_F_TPM, "Tps_LG_F_TPM", "Tps_LG_meanTPM", TTT_LG_Tps_sex_bias_TPM$table, 0.05, "Tps LG", "TPM", LG_max_FC, -1), ncol = 2)

pdf("SB_RT_TPM_smear.pdf",  width = 10, height = 15)
SB_RT_TPM_smear
dev.off()
getwd() ## where has my plot gone....

pdf("SB_HD_TPM_smear.pdf",  width = 10, height = 15)
SB_HD_TPM_smear
dev.off()
getwd() ## where has my plot gone....

pdf("SB_LG_TPM_smear.pdf",  width = 10, height = 15)
SB_LG_TPM_smear
dev.off()
getwd() ## where has my plot gone....

				


### XA

SB_RT_FPKM_smear_XA_soft <- plot_grid(
plot_smear_XA(Tbi_RT_F_FPKM, "Tbi_RT_F_FPKM", "Tbi_RT_meanFPKM", TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, "Tbi RT", "FPKM", RT_max_FC, "soft", -1),
plot_smear_XA(Tce_RT_F_FPKM, "Tce_RT_F_FPKM", "Tce_RT_meanFPKM", TTT_RT_Tce_sex_bias_FPKM$table, 0.05, "Tce RT", "FPKM", RT_max_FC, "soft", -1),
plot_smear_XA(Tcm_RT_F_FPKM, "Tcm_RT_F_FPKM", "Tcm_RT_meanFPKM", TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, "Tcm RT", "FPKM", RT_max_FC, "soft", -1),
plot_smear_XA(Tpa_RT_F_FPKM, "Tpa_RT_F_FPKM", "Tpa_RT_meanFPKM", TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, "Tpa RT", "FPKM", RT_max_FC, "soft", -1),
plot_smear_XA(Tps_RT_F_FPKM, "Tps_RT_F_FPKM", "Tps_RT_meanFPKM", TTT_RT_Tps_sex_bias_FPKM$table, 0.05, "Tps RT", "FPKM", RT_max_FC, "soft", -1), ncol = 1)

SB_HD_FPKM_smear_XA_soft <- plot_grid(
plot_smear_XA(Tbi_HD_F_FPKM, "Tbi_HD_F_FPKM", "Tbi_HD_meanFPKM", TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, "Tbi HD", "FPKM", HD_max_FC, "soft", -1),
plot_smear_XA(Tce_HD_F_FPKM, "Tce_HD_F_FPKM", "Tce_HD_meanFPKM", TTT_HD_Tce_sex_bias_FPKM$table, 0.05, "Tce HD", "FPKM", HD_max_FC, "soft", -1),
plot_smear_XA(Tcm_HD_F_FPKM, "Tcm_HD_F_FPKM", "Tcm_HD_meanFPKM", TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, "Tcm HD", "FPKM", HD_max_FC, "soft", -1),
plot_smear_XA(Tpa_HD_F_FPKM, "Tpa_HD_F_FPKM", "Tpa_HD_meanFPKM", TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, "Tpa HD", "FPKM", HD_max_FC, "soft", -1),
plot_smear_XA(Tps_HD_F_FPKM, "Tps_HD_F_FPKM", "Tps_HD_meanFPKM", TTT_HD_Tps_sex_bias_FPKM$table, 0.05, "Tps HD", "FPKM", HD_max_FC, "soft", -1), ncol = 1)

SB_LG_FPKM_smear_XA_soft <- plot_grid(
plot_smear_XA(Tbi_LG_F_FPKM, "Tbi_LG_F_FPKM", "Tbi_LG_meanFPKM", TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, "Tbi LG", "FPKM", LG_max_FC, "soft", -1),
plot_smear_XA(Tce_LG_F_FPKM, "Tce_LG_F_FPKM", "Tce_LG_meanFPKM", TTT_LG_Tce_sex_bias_FPKM$table, 0.05, "Tce LG", "FPKM", LG_max_FC, "soft", -1),
plot_smear_XA(Tcm_LG_F_FPKM, "Tcm_LG_F_FPKM", "Tcm_LG_meanFPKM", TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, "Tcm LG", "FPKM", LG_max_FC, "soft", -1),
plot_smear_XA(Tpa_LG_F_FPKM, "Tpa_LG_F_FPKM", "Tpa_LG_meanFPKM", TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, "Tpa LG", "FPKM", LG_max_FC, "soft", -1),
plot_smear_XA(Tps_LG_F_FPKM, "Tps_LG_F_FPKM", "Tps_LG_meanFPKM", TTT_LG_Tps_sex_bias_FPKM$table, 0.05, "Tps LG", "FPKM", LG_max_FC, "soft", -1), ncol = 1)

pdf("SB_RT_FPKM_smear_XA_soft.pdf",  width = 9, height = 20)
SB_RT_FPKM_smear_XA_soft
dev.off()
getwd() ## where has my plot gone....

pdf("SB_HD_FPKM_smear_XA_soft.pdf",  width = 9, height = 20)
SB_HD_FPKM_smear_XA_soft
dev.off()
getwd() ## where has my plot gone....

pdf("SB_LG_FPKM_smear_XA_soft.pdf",  width = 9, height = 20)
SB_LG_FPKM_smear_XA_soft
dev.off()
getwd() ## where has my plot gone....


SB_RT_FPKM_smear_XA_hard <- plot_grid(
plot_smear_XA(Tbi_RT_F_FPKM, "Tbi_RT_F_FPKM", "Tbi_RT_meanFPKM", TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, "Tbi RT", "FPKM", RT_max_FC, "hard", -1),
plot_smear_XA(Tce_RT_F_FPKM, "Tce_RT_F_FPKM", "Tce_RT_meanFPKM", TTT_RT_Tce_sex_bias_FPKM$table, 0.05, "Tce RT", "FPKM", RT_max_FC, "hard", -1),
plot_smear_XA(Tcm_RT_F_FPKM, "Tcm_RT_F_FPKM", "Tcm_RT_meanFPKM", TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, "Tcm RT", "FPKM", RT_max_FC, "hard", -1),
plot_smear_XA(Tpa_RT_F_FPKM, "Tpa_RT_F_FPKM", "Tpa_RT_meanFPKM", TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, "Tpa RT", "FPKM", RT_max_FC, "hard", -1),
plot_smear_XA(Tps_RT_F_FPKM, "Tps_RT_F_FPKM", "Tps_RT_meanFPKM", TTT_RT_Tps_sex_bias_FPKM$table, 0.05, "Tps RT", "FPKM", RT_max_FC, "hard", -1), ncol = 1)

SB_HD_FPKM_smear_XA_hard <- plot_grid(
plot_smear_XA(Tbi_HD_F_FPKM, "Tbi_HD_F_FPKM", "Tbi_HD_meanFPKM", TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, "Tbi HD", "FPKM", HD_max_FC, "hard", -1),
plot_smear_XA(Tce_HD_F_FPKM, "Tce_HD_F_FPKM", "Tce_HD_meanFPKM", TTT_HD_Tce_sex_bias_FPKM$table, 0.05, "Tce HD", "FPKM", HD_max_FC, "hard", -1),
plot_smear_XA(Tcm_HD_F_FPKM, "Tcm_HD_F_FPKM", "Tcm_HD_meanFPKM", TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, "Tcm HD", "FPKM", HD_max_FC, "hard", -1),
plot_smear_XA(Tpa_HD_F_FPKM, "Tpa_HD_F_FPKM", "Tpa_HD_meanFPKM", TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, "Tpa HD", "FPKM", HD_max_FC, "hard", -1),
plot_smear_XA(Tps_HD_F_FPKM, "Tps_HD_F_FPKM", "Tps_HD_meanFPKM", TTT_HD_Tps_sex_bias_FPKM$table, 0.05, "Tps HD", "FPKM", HD_max_FC, "hard", -1), ncol = 1)

SB_LG_FPKM_smear_XA_hard <- plot_grid(
plot_smear_XA(Tbi_LG_F_FPKM, "Tbi_LG_F_FPKM", "Tbi_LG_meanFPKM", TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, "Tbi LG", "FPKM", LG_max_FC, "hard", -1),
plot_smear_XA(Tce_LG_F_FPKM, "Tce_LG_F_FPKM", "Tce_LG_meanFPKM", TTT_LG_Tce_sex_bias_FPKM$table, 0.05, "Tce LG", "FPKM", LG_max_FC, "hard", -1),
plot_smear_XA(Tcm_LG_F_FPKM, "Tcm_LG_F_FPKM", "Tcm_LG_meanFPKM", TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, "Tcm LG", "FPKM", LG_max_FC, "hard", -1),
plot_smear_XA(Tpa_LG_F_FPKM, "Tpa_LG_F_FPKM", "Tpa_LG_meanFPKM", TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, "Tpa LG", "FPKM", LG_max_FC, "hard", -1),
plot_smear_XA(Tps_LG_F_FPKM, "Tps_LG_F_FPKM", "Tps_LG_meanFPKM", TTT_LG_Tps_sex_bias_FPKM$table, 0.05, "Tps LG", "FPKM", LG_max_FC, "hard", -1), ncol = 1)

pdf("SB_RT_FPKM_smear_XA_hard.pdf",  width = 9, height = 20)
SB_RT_FPKM_smear_XA_hard
dev.off()

pdf("SB_HD_FPKM_smear_XA_hard.pdf",  width = 9, height = 20)
SB_HD_FPKM_smear_XA_hard
dev.off()

pdf("SB_LG_FPKM_smear_XA_hard.pdf",  width = 9, height = 20)
SB_LG_FPKM_smear_XA_hard
dev.off()


SB_RT_TPM_smear_XA_soft <- plot_grid(
plot_smear_XA(Tbi_RT_F_TPM, "Tbi_RT_F_TPM", "Tbi_RT_meanTPM", TTT_RT_Tbi_sex_bias_TPM$table, 0.05, "Tbi RT", "TPM", RT_max_FC, "soft", -1),
plot_smear_XA(Tce_RT_F_TPM, "Tce_RT_F_TPM", "Tce_RT_meanTPM", TTT_RT_Tce_sex_bias_TPM$table, 0.05, "Tce RT", "TPM", RT_max_FC, "soft", -1),
plot_smear_XA(Tcm_RT_F_TPM, "Tcm_RT_F_TPM", "Tcm_RT_meanTPM", TTT_RT_Tcm_sex_bias_TPM$table, 0.05, "Tcm RT", "TPM", RT_max_FC, "soft", -1),
plot_smear_XA(Tpa_RT_F_TPM, "Tpa_RT_F_TPM", "Tpa_RT_meanTPM", TTT_RT_Tpa_sex_bias_TPM$table, 0.05, "Tpa RT", "TPM", RT_max_FC, "soft", -1),
plot_smear_XA(Tps_RT_F_TPM, "Tps_RT_F_TPM", "Tps_RT_meanTPM", TTT_RT_Tps_sex_bias_TPM$table, 0.05, "Tps RT", "TPM", RT_max_FC, "soft", -1), ncol = 1)

SB_HD_TPM_smear_XA_soft <- plot_grid(
plot_smear_XA(Tbi_HD_F_TPM, "Tbi_HD_F_TPM", "Tbi_HD_meanTPM", TTT_HD_Tbi_sex_bias_TPM$table, 0.05, "Tbi HD", "TPM", HD_max_FC, "soft", -1),
plot_smear_XA(Tce_HD_F_TPM, "Tce_HD_F_TPM", "Tce_HD_meanTPM", TTT_HD_Tce_sex_bias_TPM$table, 0.05, "Tce HD", "TPM", HD_max_FC, "soft", -1),
plot_smear_XA(Tcm_HD_F_TPM, "Tcm_HD_F_TPM", "Tcm_HD_meanTPM", TTT_HD_Tcm_sex_bias_TPM$table, 0.05, "Tcm HD", "TPM", HD_max_FC, "soft", -1),
plot_smear_XA(Tpa_HD_F_TPM, "Tpa_HD_F_TPM", "Tpa_HD_meanTPM", TTT_HD_Tpa_sex_bias_TPM$table, 0.05, "Tpa HD", "TPM", HD_max_FC, "soft", -1),
plot_smear_XA(Tps_HD_F_TPM, "Tps_HD_F_TPM", "Tps_HD_meanTPM", TTT_HD_Tps_sex_bias_TPM$table, 0.05, "Tps HD", "TPM", HD_max_FC, "soft", -1), ncol = 1)

SB_LG_TPM_smear_XA_soft <- plot_grid(
plot_smear_XA(Tbi_LG_F_TPM, "Tbi_LG_F_TPM", "Tbi_LG_meanTPM", TTT_LG_Tbi_sex_bias_TPM$table, 0.05, "Tbi LG", "TPM", LG_max_FC, "soft", -1),
plot_smear_XA(Tce_LG_F_TPM, "Tce_LG_F_TPM", "Tce_LG_meanTPM", TTT_LG_Tce_sex_bias_TPM$table, 0.05, "Tce LG", "TPM", LG_max_FC, "soft", -1),
plot_smear_XA(Tcm_LG_F_TPM, "Tcm_LG_F_TPM", "Tcm_LG_meanTPM", TTT_LG_Tcm_sex_bias_TPM$table, 0.05, "Tcm LG", "TPM", LG_max_FC, "soft", -1),
plot_smear_XA(Tpa_LG_F_TPM, "Tpa_LG_F_TPM", "Tpa_LG_meanTPM", TTT_LG_Tpa_sex_bias_TPM$table, 0.05, "Tpa LG", "TPM", LG_max_FC, "soft", -1),
plot_smear_XA(Tps_LG_F_TPM, "Tps_LG_F_TPM", "Tps_LG_meanTPM", TTT_LG_Tps_sex_bias_TPM$table, 0.05, "Tps LG", "TPM", LG_max_FC, "soft", -1), ncol = 1)

pdf("SB_RT_TPM_smear_XA_soft.pdf",  width = 9, height = 20)
SB_RT_TPM_smear_XA_soft
dev.off()

pdf("SB_HD_TPM_smear_XA_soft.pdf",  width = 9, height = 20)
SB_HD_TPM_smear_XA_soft
dev.off()

pdf("SB_LG_TPM_smear_XA_soft.pdf",  width = 9, height = 20)
SB_LG_TPM_smear_XA_soft
dev.off()



SB_RT_TPM_smear_XA_hard <- plot_grid(
plot_smear_XA(Tbi_RT_F_TPM, "Tbi_RT_F_TPM", "Tbi_RT_meanTPM", TTT_RT_Tbi_sex_bias_TPM$table, 0.05, "Tbi RT", "TPM", RT_max_FC, "hard", -1),
plot_smear_XA(Tce_RT_F_TPM, "Tce_RT_F_TPM", "Tce_RT_meanTPM", TTT_RT_Tce_sex_bias_TPM$table, 0.05, "Tce RT", "TPM", RT_max_FC, "hard", -1),
plot_smear_XA(Tcm_RT_F_TPM, "Tcm_RT_F_TPM", "Tcm_RT_meanTPM", TTT_RT_Tcm_sex_bias_TPM$table, 0.05, "Tcm RT", "TPM", RT_max_FC, "hard", -1),
plot_smear_XA(Tpa_RT_F_TPM, "Tpa_RT_F_TPM", "Tpa_RT_meanTPM", TTT_RT_Tpa_sex_bias_TPM$table, 0.05, "Tpa RT", "TPM", RT_max_FC, "hard", -1),
plot_smear_XA(Tps_RT_F_TPM, "Tps_RT_F_TPM", "Tps_RT_meanTPM", TTT_RT_Tps_sex_bias_TPM$table, 0.05, "Tps RT", "TPM", RT_max_FC, "hard", -1), ncol = 1)

SB_HD_TPM_smear_XA_hard <- plot_grid(
plot_smear_XA(Tbi_HD_F_TPM, "Tbi_HD_F_TPM", "Tbi_HD_meanTPM", TTT_HD_Tbi_sex_bias_TPM$table, 0.05, "Tbi HD", "TPM", HD_max_FC, "hard", -1),
plot_smear_XA(Tce_HD_F_TPM, "Tce_HD_F_TPM", "Tce_HD_meanTPM", TTT_HD_Tce_sex_bias_TPM$table, 0.05, "Tce HD", "TPM", HD_max_FC, "hard", -1),
plot_smear_XA(Tcm_HD_F_TPM, "Tcm_HD_F_TPM", "Tcm_HD_meanTPM", TTT_HD_Tcm_sex_bias_TPM$table, 0.05, "Tcm HD", "TPM", HD_max_FC, "hard", -1),
plot_smear_XA(Tpa_HD_F_TPM, "Tpa_HD_F_TPM", "Tpa_HD_meanTPM", TTT_HD_Tpa_sex_bias_TPM$table, 0.05, "Tpa HD", "TPM", HD_max_FC, "hard", -1),
plot_smear_XA(Tps_HD_F_TPM, "Tps_HD_F_TPM", "Tps_HD_meanTPM", TTT_HD_Tps_sex_bias_TPM$table, 0.05, "Tps HD", "TPM", HD_max_FC, "hard", -1), ncol = 1)

SB_LG_TPM_smear_XA_hard <- plot_grid(
plot_smear_XA(Tbi_LG_F_TPM, "Tbi_LG_F_TPM", "Tbi_LG_meanTPM", TTT_LG_Tbi_sex_bias_TPM$table, 0.05, "Tbi LG", "TPM", LG_max_FC, "hard", -1),
plot_smear_XA(Tce_LG_F_TPM, "Tce_LG_F_TPM", "Tce_LG_meanTPM", TTT_LG_Tce_sex_bias_TPM$table, 0.05, "Tce LG", "TPM", LG_max_FC, "hard", -1),
plot_smear_XA(Tcm_LG_F_TPM, "Tcm_LG_F_TPM", "Tcm_LG_meanTPM", TTT_LG_Tcm_sex_bias_TPM$table, 0.05, "Tcm LG", "TPM", LG_max_FC, "hard", -1),
plot_smear_XA(Tpa_LG_F_TPM, "Tpa_LG_F_TPM", "Tpa_LG_meanTPM", TTT_LG_Tpa_sex_bias_TPM$table, 0.05, "Tpa LG", "TPM", LG_max_FC, "hard", -1),
plot_smear_XA(Tps_LG_F_TPM, "Tps_LG_F_TPM", "Tps_LG_meanTPM", TTT_LG_Tps_sex_bias_TPM$table, 0.05, "Tps LG", "TPM", LG_max_FC, "hard", -1), ncol = 1)

pdf("SB_RT_TPM_smear_XA_hard.pdf",  width = 9, height = 20)
SB_RT_TPM_smear_XA_hard
dev.off()

pdf("SB_HD_TPM_smear_XA_hard.pdf",  width = 9, height = 20)
SB_HD_TPM_smear_XA_hard
dev.off()

pdf("SB_LG_TPM_smear_XA_hard.pdf",  width = 9, height = 20)
SB_LG_TPM_smear_XA_hard
dev.off()



### wFC


RT_max_FC <- ceiling(max(abs(c(TTT_RT_Tbi_sex_bias_FPKM$table$logFC, TTT_RT_Tce_sex_bias_FPKM$table$logFC, TTT_RT_Tcm_sex_bias_FPKM$table$logFC,TTT_RT_Tpa_sex_bias_FPKM$table$logFC,TTT_RT_Tps_sex_bias_FPKM$table$logFC))))
SB_RT_FPKM_smear_wFC <- plot_grid(
plot_smear(Tbi_RT_F_FPKM, "Tbi_RT_F_FPKM", "Tbi_RT_meanFPKM", TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, "Tbi RT", "FPKM", RT_max_FC, 1),
plot_smear(Tce_RT_F_FPKM, "Tce_RT_F_FPKM", "Tce_RT_meanFPKM", TTT_RT_Tce_sex_bias_FPKM$table, 0.05, "Tce RT", "FPKM", RT_max_FC, 1),
plot_smear(Tcm_RT_F_FPKM, "Tcm_RT_F_FPKM", "Tcm_RT_meanFPKM", TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, "Tcm RT", "FPKM", RT_max_FC, 1),
plot_smear(Tpa_RT_F_FPKM, "Tpa_RT_F_FPKM", "Tpa_RT_meanFPKM", TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, "Tpa RT", "FPKM", RT_max_FC, 1),
plot_smear(Tps_RT_F_FPKM, "Tps_RT_F_FPKM", "Tps_RT_meanFPKM", TTT_RT_Tps_sex_bias_FPKM$table, 0.05, "Tps RT", "FPKM", RT_max_FC, 1), ncol = 2)

HD_max_FC <- ceiling(max(abs(c(TTT_HD_Tbi_sex_bias_FPKM$table$logFC, TTT_HD_Tce_sex_bias_FPKM$table$logFC, TTT_HD_Tcm_sex_bias_FPKM$table$logFC,TTT_HD_Tpa_sex_bias_FPKM$table$logFC,TTT_HD_Tps_sex_bias_FPKM$table$logFC))))
SB_HD_FPKM_smear_wFC <- plot_grid(
plot_smear(Tbi_HD_F_FPKM, "Tbi_HD_F_FPKM", "Tbi_HD_meanFPKM", TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, "Tbi HD", "FPKM", HD_max_FC, 1),
plot_smear(Tce_HD_F_FPKM, "Tce_HD_F_FPKM", "Tce_HD_meanFPKM", TTT_HD_Tce_sex_bias_FPKM$table, 0.05, "Tce HD", "FPKM", HD_max_FC, 1),
plot_smear(Tcm_HD_F_FPKM, "Tcm_HD_F_FPKM", "Tcm_HD_meanFPKM", TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, "Tcm HD", "FPKM", HD_max_FC, 1),
plot_smear(Tpa_HD_F_FPKM, "Tpa_HD_F_FPKM", "Tpa_HD_meanFPKM", TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, "Tpa HD", "FPKM", HD_max_FC, 1),
plot_smear(Tps_HD_F_FPKM, "Tps_HD_F_FPKM", "Tps_HD_meanFPKM", TTT_HD_Tps_sex_bias_FPKM$table, 0.05, "Tps HD", "FPKM", HD_max_FC, 1), ncol = 2)

LG_max_FC <- ceiling(max(abs(c(TTT_LG_Tbi_sex_bias_FPKM$table$logFC, TTT_LG_Tce_sex_bias_FPKM$table$logFC, TTT_LG_Tcm_sex_bias_FPKM$table$logFC,TTT_LG_Tpa_sex_bias_FPKM$table$logFC,TTT_LG_Tps_sex_bias_FPKM$table$logFC))))
SB_LG_FPKM_smear_wFC <- plot_grid(
plot_smear(Tbi_LG_F_FPKM, "Tbi_LG_F_FPKM", "Tbi_LG_meanFPKM", TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, "Tbi LG", "FPKM", LG_max_FC, 1),
plot_smear(Tce_LG_F_FPKM, "Tce_LG_F_FPKM", "Tce_LG_meanFPKM", TTT_LG_Tce_sex_bias_FPKM$table, 0.05, "Tce LG", "FPKM", LG_max_FC, 1),
plot_smear(Tcm_LG_F_FPKM, "Tcm_LG_F_FPKM", "Tcm_LG_meanFPKM", TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, "Tcm LG", "FPKM", LG_max_FC, 1),
plot_smear(Tpa_LG_F_FPKM, "Tpa_LG_F_FPKM", "Tpa_LG_meanFPKM", TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, "Tpa LG", "FPKM", LG_max_FC, 1),
plot_smear(Tps_LG_F_FPKM, "Tps_LG_F_FPKM", "Tps_LG_meanFPKM", TTT_LG_Tps_sex_bias_FPKM$table, 0.05, "Tps LG", "FPKM", LG_max_FC, 1), ncol = 2)

pdf("SB_RT_FPKM_smear_wFC.pdf",  width = 10, height = 15)
SB_RT_FPKM_smear_wFC
dev.off()
getwd() ## where has my plot gone....

pdf("SB_HD_FPKM_smear_wFC.pdf",  width = 10, height = 15)
SB_HD_FPKM_smear_wFC
dev.off()
getwd() ## where has my plot gone....

pdf("SB_LG_FPKM_smear_wFC.pdf",  width = 10, height = 15)
SB_LG_FPKM_smear_wFC
dev.off()
getwd() ## where has my plot gone....



RT_max_FC <- ceiling(max(abs(c(TTT_RT_Tbi_sex_bias_TPM$table$logFC, TTT_RT_Tce_sex_bias_TPM$table$logFC, TTT_RT_Tcm_sex_bias_TPM$table$logFC,TTT_RT_Tpa_sex_bias_TPM$table$logFC,TTT_RT_Tps_sex_bias_TPM$table$logFC))))
SB_RT_TPM_smear_wFC <- plot_grid(
plot_smear(Tbi_RT_F_TPM, "Tbi_RT_F_TPM", "Tbi_RT_meanTPM", TTT_RT_Tbi_sex_bias_TPM$table, 0.05, "Tbi RT", "TPM", RT_max_FC, 1),
plot_smear(Tce_RT_F_TPM, "Tce_RT_F_TPM", "Tce_RT_meanTPM", TTT_RT_Tce_sex_bias_TPM$table, 0.05, "Tce RT", "TPM", RT_max_FC, 1),
plot_smear(Tcm_RT_F_TPM, "Tcm_RT_F_TPM", "Tcm_RT_meanTPM", TTT_RT_Tcm_sex_bias_TPM$table, 0.05, "Tcm RT", "TPM", RT_max_FC, 1),
plot_smear(Tpa_RT_F_TPM, "Tpa_RT_F_TPM", "Tpa_RT_meanTPM", TTT_RT_Tpa_sex_bias_TPM$table, 0.05, "Tpa RT", "TPM", RT_max_FC, 1),
plot_smear(Tps_RT_F_TPM, "Tps_RT_F_TPM", "Tps_RT_meanTPM", TTT_RT_Tps_sex_bias_TPM$table, 0.05, "Tps RT", "TPM", RT_max_FC, 1), ncol = 2)

HD_max_FC <- ceiling(max(abs(c(TTT_HD_Tbi_sex_bias_TPM$table$logFC, TTT_HD_Tce_sex_bias_TPM$table$logFC, TTT_HD_Tcm_sex_bias_TPM$table$logFC,TTT_HD_Tpa_sex_bias_TPM$table$logFC,TTT_HD_Tps_sex_bias_TPM$table$logFC))))
SB_HD_TPM_smear_wFC <- plot_grid(
plot_smear(Tbi_HD_F_TPM, "Tbi_HD_F_TPM", "Tbi_HD_meanTPM", TTT_HD_Tbi_sex_bias_TPM$table, 0.05, "Tbi HD", "TPM", HD_max_FC, 1),
plot_smear(Tce_HD_F_TPM, "Tce_HD_F_TPM", "Tce_HD_meanTPM", TTT_HD_Tce_sex_bias_TPM$table, 0.05, "Tce HD", "TPM", HD_max_FC, 1),
plot_smear(Tcm_HD_F_TPM, "Tcm_HD_F_TPM", "Tcm_HD_meanTPM", TTT_HD_Tcm_sex_bias_TPM$table, 0.05, "Tcm HD", "TPM", HD_max_FC, 1),
plot_smear(Tpa_HD_F_TPM, "Tpa_HD_F_TPM", "Tpa_HD_meanTPM", TTT_HD_Tpa_sex_bias_TPM$table, 0.05, "Tpa HD", "TPM", HD_max_FC, 1),
plot_smear(Tps_HD_F_TPM, "Tps_HD_F_TPM", "Tps_HD_meanTPM", TTT_HD_Tps_sex_bias_TPM$table, 0.05, "Tps HD", "TPM", HD_max_FC, 1), ncol = 2)

LG_max_FC <- ceiling(max(abs(c(TTT_LG_Tbi_sex_bias_TPM$table$logFC, TTT_LG_Tce_sex_bias_TPM$table$logFC, TTT_LG_Tcm_sex_bias_TPM$table$logFC,TTT_LG_Tpa_sex_bias_TPM$table$logFC,TTT_LG_Tps_sex_bias_TPM$table$logFC))))
SB_LG_TPM_smear_wFC <- plot_grid(
plot_smear(Tbi_LG_F_TPM, "Tbi_LG_F_TPM", "Tbi_LG_meanTPM", TTT_LG_Tbi_sex_bias_TPM$table, 0.05, "Tbi LG", "TPM", LG_max_FC, 1),
plot_smear(Tce_LG_F_TPM, "Tce_LG_F_TPM", "Tce_LG_meanTPM", TTT_LG_Tce_sex_bias_TPM$table, 0.05, "Tce LG", "TPM", LG_max_FC, 1),
plot_smear(Tcm_LG_F_TPM, "Tcm_LG_F_TPM", "Tcm_LG_meanTPM", TTT_LG_Tcm_sex_bias_TPM$table, 0.05, "Tcm LG", "TPM", LG_max_FC, 1),
plot_smear(Tpa_LG_F_TPM, "Tpa_LG_F_TPM", "Tpa_LG_meanTPM", TTT_LG_Tpa_sex_bias_TPM$table, 0.05, "Tpa LG", "TPM", LG_max_FC, 1),
plot_smear(Tps_LG_F_TPM, "Tps_LG_F_TPM", "Tps_LG_meanTPM", TTT_LG_Tps_sex_bias_TPM$table, 0.05, "Tps LG", "TPM", LG_max_FC, 1), ncol = 2)

pdf("SB_RT_TPM_smear_wFC.pdf",  width = 10, height = 15)
SB_RT_TPM_smear_wFC
dev.off()
getwd() ## where has my plot gone....

pdf("SB_HD_TPM_smear_wFC.pdf",  width = 10, height = 15)
SB_HD_TPM_smear_wFC
dev.off()
getwd() ## where has my plot gone....

pdf("SB_LG_TPM_smear_wFC.pdf",  width = 10, height = 15)
SB_LG_TPM_smear_wFC
dev.off()
getwd() ## where has my plot gone....

				


### XA

SB_RT_FPKM_smear_XA_soft_wFC <- plot_grid(
plot_smear_XA(Tbi_RT_F_FPKM, "Tbi_RT_F_FPKM", "Tbi_RT_meanFPKM", TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, "Tbi RT", "FPKM", RT_max_FC, "soft", 1),
plot_smear_XA(Tce_RT_F_FPKM, "Tce_RT_F_FPKM", "Tce_RT_meanFPKM", TTT_RT_Tce_sex_bias_FPKM$table, 0.05, "Tce RT", "FPKM", RT_max_FC, "soft", 1),
plot_smear_XA(Tcm_RT_F_FPKM, "Tcm_RT_F_FPKM", "Tcm_RT_meanFPKM", TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, "Tcm RT", "FPKM", RT_max_FC, "soft", 1),
plot_smear_XA(Tpa_RT_F_FPKM, "Tpa_RT_F_FPKM", "Tpa_RT_meanFPKM", TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, "Tpa RT", "FPKM", RT_max_FC, "soft", 1),
plot_smear_XA(Tps_RT_F_FPKM, "Tps_RT_F_FPKM", "Tps_RT_meanFPKM", TTT_RT_Tps_sex_bias_FPKM$table, 0.05, "Tps RT", "FPKM", RT_max_FC, "soft", 1), ncol = 1)

SB_HD_FPKM_smear_XA_soft_wFC <- plot_grid(
plot_smear_XA(Tbi_HD_F_FPKM, "Tbi_HD_F_FPKM", "Tbi_HD_meanFPKM", TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, "Tbi HD", "FPKM", HD_max_FC, "soft", 1),
plot_smear_XA(Tce_HD_F_FPKM, "Tce_HD_F_FPKM", "Tce_HD_meanFPKM", TTT_HD_Tce_sex_bias_FPKM$table, 0.05, "Tce HD", "FPKM", HD_max_FC, "soft", 1),
plot_smear_XA(Tcm_HD_F_FPKM, "Tcm_HD_F_FPKM", "Tcm_HD_meanFPKM", TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, "Tcm HD", "FPKM", HD_max_FC, "soft", 1),
plot_smear_XA(Tpa_HD_F_FPKM, "Tpa_HD_F_FPKM", "Tpa_HD_meanFPKM", TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, "Tpa HD", "FPKM", HD_max_FC, "soft", 1),
plot_smear_XA(Tps_HD_F_FPKM, "Tps_HD_F_FPKM", "Tps_HD_meanFPKM", TTT_HD_Tps_sex_bias_FPKM$table, 0.05, "Tps HD", "FPKM", HD_max_FC, "soft", 1), ncol = 1)

SB_LG_FPKM_smear_XA_soft_wFC <- plot_grid(
plot_smear_XA(Tbi_LG_F_FPKM, "Tbi_LG_F_FPKM", "Tbi_LG_meanFPKM", TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, "Tbi LG", "FPKM", LG_max_FC, "soft", 1),
plot_smear_XA(Tce_LG_F_FPKM, "Tce_LG_F_FPKM", "Tce_LG_meanFPKM", TTT_LG_Tce_sex_bias_FPKM$table, 0.05, "Tce LG", "FPKM", LG_max_FC, "soft", 1),
plot_smear_XA(Tcm_LG_F_FPKM, "Tcm_LG_F_FPKM", "Tcm_LG_meanFPKM", TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, "Tcm LG", "FPKM", LG_max_FC, "soft", 1),
plot_smear_XA(Tpa_LG_F_FPKM, "Tpa_LG_F_FPKM", "Tpa_LG_meanFPKM", TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, "Tpa LG", "FPKM", LG_max_FC, "soft", 1),
plot_smear_XA(Tps_LG_F_FPKM, "Tps_LG_F_FPKM", "Tps_LG_meanFPKM", TTT_LG_Tps_sex_bias_FPKM$table, 0.05, "Tps LG", "FPKM", LG_max_FC, "soft", 1), ncol = 1)

pdf("SB_RT_FPKM_smear_XA_soft_wFC.pdf",  width = 9, height = 20)
SB_RT_FPKM_smear_XA_soft_wFC
dev.off()
getwd() ## where has my plot gone....

pdf("SB_HD_FPKM_smear_XA_soft_wFC.pdf",  width = 9, height = 20)
SB_HD_FPKM_smear_XA_soft_wFC
dev.off()
getwd() ## where has my plot gone....

pdf("SB_LG_FPKM_smear_XA_soft_wFC.pdf",  width = 9, height = 20)
SB_LG_FPKM_smear_XA_soft_wFC
dev.off()
getwd() ## where has my plot gone....


SB_RT_FPKM_smear_XA_hard_wFC <- plot_grid(
plot_smear_XA(Tbi_RT_F_FPKM, "Tbi_RT_F_FPKM", "Tbi_RT_meanFPKM", TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, "Tbi RT", "FPKM", RT_max_FC, "hard", 1),
plot_smear_XA(Tce_RT_F_FPKM, "Tce_RT_F_FPKM", "Tce_RT_meanFPKM", TTT_RT_Tce_sex_bias_FPKM$table, 0.05, "Tce RT", "FPKM", RT_max_FC, "hard", 1),
plot_smear_XA(Tcm_RT_F_FPKM, "Tcm_RT_F_FPKM", "Tcm_RT_meanFPKM", TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, "Tcm RT", "FPKM", RT_max_FC, "hard", 1),
plot_smear_XA(Tpa_RT_F_FPKM, "Tpa_RT_F_FPKM", "Tpa_RT_meanFPKM", TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, "Tpa RT", "FPKM", RT_max_FC, "hard", 1),
plot_smear_XA(Tps_RT_F_FPKM, "Tps_RT_F_FPKM", "Tps_RT_meanFPKM", TTT_RT_Tps_sex_bias_FPKM$table, 0.05, "Tps RT", "FPKM", RT_max_FC, "hard", 1), ncol = 1)

SB_HD_FPKM_smear_XA_hard_wFC <- plot_grid(
plot_smear_XA(Tbi_HD_F_FPKM, "Tbi_HD_F_FPKM", "Tbi_HD_meanFPKM", TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, "Tbi HD", "FPKM", HD_max_FC, "hard", 1),
plot_smear_XA(Tce_HD_F_FPKM, "Tce_HD_F_FPKM", "Tce_HD_meanFPKM", TTT_HD_Tce_sex_bias_FPKM$table, 0.05, "Tce HD", "FPKM", HD_max_FC, "hard", 1),
plot_smear_XA(Tcm_HD_F_FPKM, "Tcm_HD_F_FPKM", "Tcm_HD_meanFPKM", TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, "Tcm HD", "FPKM", HD_max_FC, "hard", 1),
plot_smear_XA(Tpa_HD_F_FPKM, "Tpa_HD_F_FPKM", "Tpa_HD_meanFPKM", TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, "Tpa HD", "FPKM", HD_max_FC, "hard", 1),
plot_smear_XA(Tps_HD_F_FPKM, "Tps_HD_F_FPKM", "Tps_HD_meanFPKM", TTT_HD_Tps_sex_bias_FPKM$table, 0.05, "Tps HD", "FPKM", HD_max_FC, "hard", 1), ncol = 1)

SB_LG_FPKM_smear_XA_hard_wFC <- plot_grid(
plot_smear_XA(Tbi_LG_F_FPKM, "Tbi_LG_F_FPKM", "Tbi_LG_meanFPKM", TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, "Tbi LG", "FPKM", LG_max_FC, "hard", 1),
plot_smear_XA(Tce_LG_F_FPKM, "Tce_LG_F_FPKM", "Tce_LG_meanFPKM", TTT_LG_Tce_sex_bias_FPKM$table, 0.05, "Tce LG", "FPKM", LG_max_FC, "hard", 1),
plot_smear_XA(Tcm_LG_F_FPKM, "Tcm_LG_F_FPKM", "Tcm_LG_meanFPKM", TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, "Tcm LG", "FPKM", LG_max_FC, "hard", 1),
plot_smear_XA(Tpa_LG_F_FPKM, "Tpa_LG_F_FPKM", "Tpa_LG_meanFPKM", TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, "Tpa LG", "FPKM", LG_max_FC, "hard", 1),
plot_smear_XA(Tps_LG_F_FPKM, "Tps_LG_F_FPKM", "Tps_LG_meanFPKM", TTT_LG_Tps_sex_bias_FPKM$table, 0.05, "Tps LG", "FPKM", LG_max_FC, "hard", 1), ncol = 1)

pdf("SB_RT_FPKM_smear_XA_hard_wFC.pdf",  width = 9, height = 20)
SB_RT_FPKM_smear_XA_hard_wFC
dev.off()

pdf("SB_HD_FPKM_smear_XA_hard_wFC.pdf",  width = 9, height = 20)
SB_HD_FPKM_smear_XA_hard_wFC
dev.off()

pdf("SB_LG_FPKM_smear_XA_hard_wFC.pdf",  width = 9, height = 20)
SB_LG_FPKM_smear_XA_hard_wFC
dev.off()


SB_RT_TPM_smear_XA_soft_wFC <- plot_grid(
plot_smear_XA(Tbi_RT_F_TPM, "Tbi_RT_F_TPM", "Tbi_RT_meanTPM", TTT_RT_Tbi_sex_bias_TPM$table, 0.05, "Tbi RT", "TPM", RT_max_FC, "soft", 1),
plot_smear_XA(Tce_RT_F_TPM, "Tce_RT_F_TPM", "Tce_RT_meanTPM", TTT_RT_Tce_sex_bias_TPM$table, 0.05, "Tce RT", "TPM", RT_max_FC, "soft", 1),
plot_smear_XA(Tcm_RT_F_TPM, "Tcm_RT_F_TPM", "Tcm_RT_meanTPM", TTT_RT_Tcm_sex_bias_TPM$table, 0.05, "Tcm RT", "TPM", RT_max_FC, "soft", 1),
plot_smear_XA(Tpa_RT_F_TPM, "Tpa_RT_F_TPM", "Tpa_RT_meanTPM", TTT_RT_Tpa_sex_bias_TPM$table, 0.05, "Tpa RT", "TPM", RT_max_FC, "soft", 1),
plot_smear_XA(Tps_RT_F_TPM, "Tps_RT_F_TPM", "Tps_RT_meanTPM", TTT_RT_Tps_sex_bias_TPM$table, 0.05, "Tps RT", "TPM", RT_max_FC, "soft", 1), ncol = 1)

SB_HD_TPM_smear_XA_soft_wFC <- plot_grid(
plot_smear_XA(Tbi_HD_F_TPM, "Tbi_HD_F_TPM", "Tbi_HD_meanTPM", TTT_HD_Tbi_sex_bias_TPM$table, 0.05, "Tbi HD", "TPM", HD_max_FC, "soft", 1),
plot_smear_XA(Tce_HD_F_TPM, "Tce_HD_F_TPM", "Tce_HD_meanTPM", TTT_HD_Tce_sex_bias_TPM$table, 0.05, "Tce HD", "TPM", HD_max_FC, "soft", 1),
plot_smear_XA(Tcm_HD_F_TPM, "Tcm_HD_F_TPM", "Tcm_HD_meanTPM", TTT_HD_Tcm_sex_bias_TPM$table, 0.05, "Tcm HD", "TPM", HD_max_FC, "soft", 1),
plot_smear_XA(Tpa_HD_F_TPM, "Tpa_HD_F_TPM", "Tpa_HD_meanTPM", TTT_HD_Tpa_sex_bias_TPM$table, 0.05, "Tpa HD", "TPM", HD_max_FC, "soft", 1),
plot_smear_XA(Tps_HD_F_TPM, "Tps_HD_F_TPM", "Tps_HD_meanTPM", TTT_HD_Tps_sex_bias_TPM$table, 0.05, "Tps HD", "TPM", HD_max_FC, "soft", 1), ncol = 1)

SB_LG_TPM_smear_XA_soft_wFC <- plot_grid(
plot_smear_XA(Tbi_LG_F_TPM, "Tbi_LG_F_TPM", "Tbi_LG_meanTPM", TTT_LG_Tbi_sex_bias_TPM$table, 0.05, "Tbi LG", "TPM", LG_max_FC, "soft", 1),
plot_smear_XA(Tce_LG_F_TPM, "Tce_LG_F_TPM", "Tce_LG_meanTPM", TTT_LG_Tce_sex_bias_TPM$table, 0.05, "Tce LG", "TPM", LG_max_FC, "soft", 1),
plot_smear_XA(Tcm_LG_F_TPM, "Tcm_LG_F_TPM", "Tcm_LG_meanTPM", TTT_LG_Tcm_sex_bias_TPM$table, 0.05, "Tcm LG", "TPM", LG_max_FC, "soft", 1),
plot_smear_XA(Tpa_LG_F_TPM, "Tpa_LG_F_TPM", "Tpa_LG_meanTPM", TTT_LG_Tpa_sex_bias_TPM$table, 0.05, "Tpa LG", "TPM", LG_max_FC, "soft", 1),
plot_smear_XA(Tps_LG_F_TPM, "Tps_LG_F_TPM", "Tps_LG_meanTPM", TTT_LG_Tps_sex_bias_TPM$table, 0.05, "Tps LG", "TPM", LG_max_FC, "soft", 1), ncol = 1)

pdf("SB_RT_TPM_smear_XA_soft_wFC.pdf",  width = 9, height = 20)
SB_RT_TPM_smear_XA_soft_wFC
dev.off()

pdf("SB_HD_TPM_smear_XA_soft_wFC.pdf",  width = 9, height = 20)
SB_HD_TPM_smear_XA_soft_wFC
dev.off()

pdf("SB_LG_TPM_smear_XA_soft_wFC.pdf",  width = 9, height = 20)
SB_LG_TPM_smear_XA_soft_wFC
dev.off()



SB_RT_TPM_smear_XA_hard_wFC <- plot_grid(
plot_smear_XA(Tbi_RT_F_TPM, "Tbi_RT_F_TPM", "Tbi_RT_meanTPM", TTT_RT_Tbi_sex_bias_TPM$table, 0.05, "Tbi RT", "TPM", RT_max_FC, "hard", 1),
plot_smear_XA(Tce_RT_F_TPM, "Tce_RT_F_TPM", "Tce_RT_meanTPM", TTT_RT_Tce_sex_bias_TPM$table, 0.05, "Tce RT", "TPM", RT_max_FC, "hard", 1),
plot_smear_XA(Tcm_RT_F_TPM, "Tcm_RT_F_TPM", "Tcm_RT_meanTPM", TTT_RT_Tcm_sex_bias_TPM$table, 0.05, "Tcm RT", "TPM", RT_max_FC, "hard", 1),
plot_smear_XA(Tpa_RT_F_TPM, "Tpa_RT_F_TPM", "Tpa_RT_meanTPM", TTT_RT_Tpa_sex_bias_TPM$table, 0.05, "Tpa RT", "TPM", RT_max_FC, "hard", 1),
plot_smear_XA(Tps_RT_F_TPM, "Tps_RT_F_TPM", "Tps_RT_meanTPM", TTT_RT_Tps_sex_bias_TPM$table, 0.05, "Tps RT", "TPM", RT_max_FC, "hard", 1), ncol = 1)

SB_HD_TPM_smear_XA_hard_wFC <- plot_grid(
plot_smear_XA(Tbi_HD_F_TPM, "Tbi_HD_F_TPM", "Tbi_HD_meanTPM", TTT_HD_Tbi_sex_bias_TPM$table, 0.05, "Tbi HD", "TPM", HD_max_FC, "hard", 1),
plot_smear_XA(Tce_HD_F_TPM, "Tce_HD_F_TPM", "Tce_HD_meanTPM", TTT_HD_Tce_sex_bias_TPM$table, 0.05, "Tce HD", "TPM", HD_max_FC, "hard", 1),
plot_smear_XA(Tcm_HD_F_TPM, "Tcm_HD_F_TPM", "Tcm_HD_meanTPM", TTT_HD_Tcm_sex_bias_TPM$table, 0.05, "Tcm HD", "TPM", HD_max_FC, "hard", 1),
plot_smear_XA(Tpa_HD_F_TPM, "Tpa_HD_F_TPM", "Tpa_HD_meanTPM", TTT_HD_Tpa_sex_bias_TPM$table, 0.05, "Tpa HD", "TPM", HD_max_FC, "hard", 1),
plot_smear_XA(Tps_HD_F_TPM, "Tps_HD_F_TPM", "Tps_HD_meanTPM", TTT_HD_Tps_sex_bias_TPM$table, 0.05, "Tps HD", "TPM", HD_max_FC, "hard", 1), ncol = 1)

SB_LG_TPM_smear_XA_hard_wFC <- plot_grid(
plot_smear_XA(Tbi_LG_F_TPM, "Tbi_LG_F_TPM", "Tbi_LG_meanTPM", TTT_LG_Tbi_sex_bias_TPM$table, 0.05, "Tbi LG", "TPM", LG_max_FC, "hard", 1),
plot_smear_XA(Tce_LG_F_TPM, "Tce_LG_F_TPM", "Tce_LG_meanTPM", TTT_LG_Tce_sex_bias_TPM$table, 0.05, "Tce LG", "TPM", LG_max_FC, "hard", 1),
plot_smear_XA(Tcm_LG_F_TPM, "Tcm_LG_F_TPM", "Tcm_LG_meanTPM", TTT_LG_Tcm_sex_bias_TPM$table, 0.05, "Tcm LG", "TPM", LG_max_FC, "hard", 1),
plot_smear_XA(Tpa_LG_F_TPM, "Tpa_LG_F_TPM", "Tpa_LG_meanTPM", TTT_LG_Tpa_sex_bias_TPM$table, 0.05, "Tpa LG", "TPM", LG_max_FC, "hard", 1),
plot_smear_XA(Tps_LG_F_TPM, "Tps_LG_F_TPM", "Tps_LG_meanTPM", TTT_LG_Tps_sex_bias_TPM$table, 0.05, "Tps LG", "TPM", LG_max_FC, "hard", 1), ncol = 1)

pdf("SB_RT_TPM_smear_XA_hard_wFC.pdf",  width = 9, height = 20)
SB_RT_TPM_smear_XA_hard_wFC
dev.off()

pdf("SB_HD_TPM_smear_XA_hard_wFC.pdf",  width = 9, height = 20)
SB_HD_TPM_smear_XA_hard_wFC
dev.off()

pdf("SB_LG_TPM_smear_XA_hard_wFC.pdf",  width = 9, height = 20)
SB_LG_TPM_smear_XA_hard_wFC
dev.off()





########################################################################################################
##### are MB or FB genes enriched on the X?


MB_FB_on_X <- function(TTT_df, FDR, logFC_cut, full_info_table, chr_class){
	
	TTT_df_sig = subset(TTT_df, TTT_df$FDR < 0.05)
	TTT_df_sig_MB = subset(TTT_df_sig, TTT_df_sig$logFC < 0)
	TTT_df_sig_FB = subset(TTT_df_sig, TTT_df_sig$logFC > 0)
	TTT_df_sig_MB_FC = subset(TTT_df_sig_MB, TTT_df_sig_MB$logFC <= -1 * logFC_cut )	
	TTT_df_sig_FB_FC = subset(TTT_df_sig_FB, TTT_df_sig_FB$logFC >= logFC_cut)		


	FB_genes = as.character(TTT_df_sig_FB_FC$genes)
	MB_genes = as.character(TTT_df_sig_MB_FC$genes)
	SB_genes = c(FB_genes, MB_genes)
	all_exp_genes = as.character(TTT_df$genes)
	
	print("N female-biased")
	print(length(FB_genes))
	print("N male-biased")
	print(length(MB_genes))	

	## only consider genes that wee expressed
	exp_info_table <- subset(full_info_table, full_info_table$gene_id %in% all_exp_genes)

	FB_genes_df = subset(exp_info_table,   exp_info_table$gene_id %in% FB_genes)
	MB_genes_df = subset(exp_info_table,   exp_info_table$gene_id %in% MB_genes)
	UB_genes_df = subset(exp_info_table, !(exp_info_table$gene_id %in% SB_genes))
	
	
	#N_FB_X = length(subset(FB_genes_df, FB_genes_df$chr_hard == "X")[,1])
	N_FB_X = length(subset(FB_genes_df, eval(parse(text=paste("FB_genes_df",'$',chr_class,sep=''))) == "X")[,1])
	N_FB_A = length(subset(FB_genes_df, eval(parse(text=paste("FB_genes_df",'$',chr_class,sep=''))) == "A")[,1])
	N_MB_X = length(subset(MB_genes_df, eval(parse(text=paste("MB_genes_df",'$',chr_class,sep=''))) == "X")[,1])
	N_MB_A = length(subset(MB_genes_df, eval(parse(text=paste("MB_genes_df",'$',chr_class,sep=''))) == "A")[,1])
	N_All_exp_X = length(subset(exp_info_table, eval(parse(text=paste("exp_info_table",'$',chr_class,sep=''))) == "X")[,1])
	N_All_exp_A = length(subset(exp_info_table, eval(parse(text=paste("exp_info_table",'$',chr_class,sep=''))) == "A")[,1])

	FB_FT_mat <- matrix(c(N_FB_X,(N_All_exp_X - N_FB_X),N_FB_A,(N_All_exp_A - N_FB_A)), nrow = 2)
	MB_FT_mat <- matrix(c(N_MB_X,(N_All_exp_X - N_MB_X),N_MB_A,(N_All_exp_A - N_MB_A)), nrow = 2)
	
	FB_fisher <- fisher.test(FB_FT_mat, alternative="two.sided") ### enrichment (over or under-rep)
	MB_fisher <- fisher.test(MB_FT_mat, alternative="two.sided") ### enrichment (over or under-rep)

	x_name <- deparse(substitute(TTT_df))
	x_name2 <- strsplit(x_name, "_")[[1]]
	sp <- x_name2[3]
	tiss <- x_name2[2]

	out_df <- as.data.frame(cbind(sp, tiss, chr_class, paste("SB_FDR", FDR, "logFC", logFC_cut, sep = "_"), N_FB_X, N_FB_A, N_MB_X, N_MB_A, N_All_exp_X, N_All_exp_A, FB_fisher$estimate,  MB_fisher$estimate,  FB_fisher$p.value,  MB_fisher$p.value))
	colnames(out_df ) <- c("sp", "tiss", "chr_class", "SB_class", "N_FB_X", "N_FB_A", "N_MB_X", "N_MB_A", "N_All_exp_X", "N_All_exp_A", "FB_fisher_OR",  "MB_fisher_OR",  "FB_fisher_p",  "MB_fisher_p")
	
	return(out_df)
}


FT_MB_FB_X_A_hard_FPKM_df <- as.data.frame(rbind(
MB_FB_on_X(TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, 0, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_RT_Tce_sex_bias_FPKM$table, 0.05, 0, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, 0, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, 0, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_RT_Tps_sex_bias_FPKM$table, 0.05, 0, dat1_Tps, "chr_hard"),
MB_FB_on_X(TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, 0, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_HD_Tce_sex_bias_FPKM$table, 0.05, 0, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, 0, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, 0, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_HD_Tps_sex_bias_FPKM$table, 0.05, 0, dat1_Tps, "chr_hard"),
MB_FB_on_X(TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, 0, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_LG_Tce_sex_bias_FPKM$table, 0.05, 0, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, 0, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, 0, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_LG_Tps_sex_bias_FPKM$table, 0.05, 0, dat1_Tps, "chr_hard")))

FT_MB_FB_X_A_soft_FPKM_df <- as.data.frame(rbind(
MB_FB_on_X(TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, 0, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_RT_Tce_sex_bias_FPKM$table, 0.05, 0, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, 0, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, 0, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_RT_Tps_sex_bias_FPKM$table, 0.05, 0, dat1_Tps, "chr_soft"),
MB_FB_on_X(TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, 0, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_HD_Tce_sex_bias_FPKM$table, 0.05, 0, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, 0, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, 0, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_HD_Tps_sex_bias_FPKM$table, 0.05, 0, dat1_Tps, "chr_soft"),
MB_FB_on_X(TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, 0, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_LG_Tce_sex_bias_FPKM$table, 0.05, 0, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, 0, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, 0, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_LG_Tps_sex_bias_FPKM$table, 0.05, 0, dat1_Tps, "chr_soft")))

FT_MB_FB_X_A_hard_TPM_df <- as.data.frame(rbind(
MB_FB_on_X(TTT_RT_Tbi_sex_bias_TPM$table, 0.05, 0, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_RT_Tce_sex_bias_TPM$table, 0.05, 0, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_RT_Tcm_sex_bias_TPM$table, 0.05, 0, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_RT_Tpa_sex_bias_TPM$table, 0.05, 0, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_RT_Tps_sex_bias_TPM$table, 0.05, 0, dat1_Tps, "chr_hard"),
MB_FB_on_X(TTT_HD_Tbi_sex_bias_TPM$table, 0.05, 0, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_HD_Tce_sex_bias_TPM$table, 0.05, 0, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_HD_Tcm_sex_bias_TPM$table, 0.05, 0, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_HD_Tpa_sex_bias_TPM$table, 0.05, 0, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_HD_Tps_sex_bias_TPM$table, 0.05, 0, dat1_Tps, "chr_hard"),
MB_FB_on_X(TTT_LG_Tbi_sex_bias_TPM$table, 0.05, 0, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_LG_Tce_sex_bias_TPM$table, 0.05, 0, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_LG_Tcm_sex_bias_TPM$table, 0.05, 0, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_LG_Tpa_sex_bias_TPM$table, 0.05, 0, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_LG_Tps_sex_bias_TPM$table, 0.05, 0, dat1_Tps, "chr_hard")))

FT_MB_FB_X_A_soft_TPM_df <- as.data.frame(rbind(
MB_FB_on_X(TTT_RT_Tbi_sex_bias_TPM$table, 0.05, 0, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_RT_Tce_sex_bias_TPM$table, 0.05, 0, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_RT_Tcm_sex_bias_TPM$table, 0.05, 0, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_RT_Tpa_sex_bias_TPM$table, 0.05, 0, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_RT_Tps_sex_bias_TPM$table, 0.05, 0, dat1_Tps, "chr_soft"),
MB_FB_on_X(TTT_HD_Tbi_sex_bias_TPM$table, 0.05, 0, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_HD_Tce_sex_bias_TPM$table, 0.05, 0, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_HD_Tcm_sex_bias_TPM$table, 0.05, 0, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_HD_Tpa_sex_bias_TPM$table, 0.05, 0, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_HD_Tps_sex_bias_TPM$table, 0.05, 0, dat1_Tps, "chr_soft"),
MB_FB_on_X(TTT_LG_Tbi_sex_bias_TPM$table, 0.05, 0, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_LG_Tce_sex_bias_TPM$table, 0.05, 0, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_LG_Tcm_sex_bias_TPM$table, 0.05, 0, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_LG_Tpa_sex_bias_TPM$table, 0.05, 0, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_LG_Tps_sex_bias_TPM$table, 0.05, 0, dat1_Tps, "chr_soft")))

FT_MB_FB_X_A_hard_FPKM_wFC_df <- as.data.frame(rbind(
MB_FB_on_X(TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, 1, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_RT_Tce_sex_bias_FPKM$table, 0.05, 1, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, 1, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, 1, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_RT_Tps_sex_bias_FPKM$table, 0.05, 1, dat1_Tps, "chr_hard"),
MB_FB_on_X(TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, 1, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_HD_Tce_sex_bias_FPKM$table, 0.05, 1, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, 1, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, 1, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_HD_Tps_sex_bias_FPKM$table, 0.05, 1, dat1_Tps, "chr_hard"),
MB_FB_on_X(TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, 1, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_LG_Tce_sex_bias_FPKM$table, 0.05, 1, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, 1, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, 1, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_LG_Tps_sex_bias_FPKM$table, 0.05, 1, dat1_Tps, "chr_hard")))

FT_MB_FB_X_A_soft_FPKM_wFC_df <- as.data.frame(rbind(
MB_FB_on_X(TTT_RT_Tbi_sex_bias_FPKM$table, 0.05, 1, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_RT_Tce_sex_bias_FPKM$table, 0.05, 1, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_RT_Tcm_sex_bias_FPKM$table, 0.05, 1, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_RT_Tpa_sex_bias_FPKM$table, 0.05, 1, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_RT_Tps_sex_bias_FPKM$table, 0.05, 1, dat1_Tps, "chr_soft"),
MB_FB_on_X(TTT_HD_Tbi_sex_bias_FPKM$table, 0.05, 1, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_HD_Tce_sex_bias_FPKM$table, 0.05, 1, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_HD_Tcm_sex_bias_FPKM$table, 0.05, 1, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_HD_Tpa_sex_bias_FPKM$table, 0.05, 1, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_HD_Tps_sex_bias_FPKM$table, 0.05, 1, dat1_Tps, "chr_soft"),
MB_FB_on_X(TTT_LG_Tbi_sex_bias_FPKM$table, 0.05, 1, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_LG_Tce_sex_bias_FPKM$table, 0.05, 1, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_LG_Tcm_sex_bias_FPKM$table, 0.05, 1, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_LG_Tpa_sex_bias_FPKM$table, 0.05, 1, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_LG_Tps_sex_bias_FPKM$table, 0.05, 1, dat1_Tps, "chr_soft")))

FT_MB_FB_X_A_hard_TPM_wFC_df <- as.data.frame(rbind(
MB_FB_on_X(TTT_RT_Tbi_sex_bias_TPM$table, 0.05, 1, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_RT_Tce_sex_bias_TPM$table, 0.05, 1, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_RT_Tcm_sex_bias_TPM$table, 0.05, 1, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_RT_Tpa_sex_bias_TPM$table, 0.05, 1, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_RT_Tps_sex_bias_TPM$table, 0.05, 1, dat1_Tps, "chr_hard"),
MB_FB_on_X(TTT_HD_Tbi_sex_bias_TPM$table, 0.05, 1, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_HD_Tce_sex_bias_TPM$table, 0.05, 1, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_HD_Tcm_sex_bias_TPM$table, 0.05, 1, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_HD_Tpa_sex_bias_TPM$table, 0.05, 1, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_HD_Tps_sex_bias_TPM$table, 0.05, 1, dat1_Tps, "chr_hard"),
MB_FB_on_X(TTT_LG_Tbi_sex_bias_TPM$table, 0.05, 1, dat1_Tbi, "chr_hard"),
MB_FB_on_X(TTT_LG_Tce_sex_bias_TPM$table, 0.05, 1, dat1_Tce, "chr_hard"),
MB_FB_on_X(TTT_LG_Tcm_sex_bias_TPM$table, 0.05, 1, dat1_Tcm, "chr_hard"),
MB_FB_on_X(TTT_LG_Tpa_sex_bias_TPM$table, 0.05, 1, dat1_Tpa, "chr_hard"),
MB_FB_on_X(TTT_LG_Tps_sex_bias_TPM$table, 0.05, 1, dat1_Tps, "chr_hard")))

FT_MB_FB_X_A_soft_TPM_wFC_df <- as.data.frame(rbind(
MB_FB_on_X(TTT_RT_Tbi_sex_bias_TPM$table, 0.05, 1, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_RT_Tce_sex_bias_TPM$table, 0.05, 1, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_RT_Tcm_sex_bias_TPM$table, 0.05, 1, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_RT_Tpa_sex_bias_TPM$table, 0.05, 1, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_RT_Tps_sex_bias_TPM$table, 0.05, 1, dat1_Tps, "chr_soft"),
MB_FB_on_X(TTT_HD_Tbi_sex_bias_TPM$table, 0.05, 1, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_HD_Tce_sex_bias_TPM$table, 0.05, 1, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_HD_Tcm_sex_bias_TPM$table, 0.05, 1, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_HD_Tpa_sex_bias_TPM$table, 0.05, 1, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_HD_Tps_sex_bias_TPM$table, 0.05, 1, dat1_Tps, "chr_soft"),
MB_FB_on_X(TTT_LG_Tbi_sex_bias_TPM$table, 0.05, 1, dat1_Tbi, "chr_soft"),
MB_FB_on_X(TTT_LG_Tce_sex_bias_TPM$table, 0.05, 1, dat1_Tce, "chr_soft"),
MB_FB_on_X(TTT_LG_Tcm_sex_bias_TPM$table, 0.05, 1, dat1_Tcm, "chr_soft"),
MB_FB_on_X(TTT_LG_Tpa_sex_bias_TPM$table, 0.05, 1, dat1_Tpa, "chr_soft"),
MB_FB_on_X(TTT_LG_Tps_sex_bias_TPM$table, 0.05, 1, dat1_Tps, "chr_soft")))


####### FDR correct

FT_MB_FB_X_A_hard_FPKM_df$FB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_hard_FPKM_df$FB_fisher_p)), method = "BH")
FT_MB_FB_X_A_hard_FPKM_df$MB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_hard_FPKM_df$MB_fisher_p)), method = "BH")
FT_MB_FB_X_A_soft_FPKM_df$FB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_soft_FPKM_df$FB_fisher_p)), method = "BH")
FT_MB_FB_X_A_soft_FPKM_df$MB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_soft_FPKM_df$MB_fisher_p)), method = "BH")
FT_MB_FB_X_A_hard_TPM_df$FB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_hard_TPM_df$FB_fisher_p)), method = "BH")
FT_MB_FB_X_A_hard_TPM_df$MB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_hard_TPM_df$MB_fisher_p)), method = "BH")
FT_MB_FB_X_A_soft_TPM_df$FB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_soft_TPM_df$FB_fisher_p)), method = "BH")
FT_MB_FB_X_A_soft_TPM_df$MB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_soft_TPM_df$MB_fisher_p)), method = "BH")

FT_MB_FB_X_A_hard_FPKM_wFC_df$FB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_hard_FPKM_wFC_df$FB_fisher_p)), method = "BH")
FT_MB_FB_X_A_hard_FPKM_wFC_df$MB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_hard_FPKM_wFC_df$MB_fisher_p)), method = "BH")
FT_MB_FB_X_A_soft_FPKM_wFC_df$FB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_soft_FPKM_wFC_df$FB_fisher_p)), method = "BH")
FT_MB_FB_X_A_soft_FPKM_wFC_df$MB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_soft_FPKM_wFC_df$MB_fisher_p)), method = "BH")
FT_MB_FB_X_A_hard_TPM_wFC_df$FB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_hard_TPM_wFC_df$FB_fisher_p)), method = "BH")
FT_MB_FB_X_A_hard_TPM_wFC_df$MB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_hard_TPM_wFC_df$MB_fisher_p)), method = "BH")
FT_MB_FB_X_A_soft_TPM_wFC_df$FB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_soft_TPM_wFC_df$FB_fisher_p)), method = "BH")
FT_MB_FB_X_A_soft_TPM_wFC_df$MB_fisher_FDR <- p.adjust(as.numeric(as.character(FT_MB_FB_X_A_soft_TPM_wFC_df$MB_fisher_p)), method = "BH")

write.csv(FT_MB_FB_X_A_hard_FPKM_df, "FT_MB_FB_X_A_hard_FPKM_df.csv", row.names = F)
write.csv(FT_MB_FB_X_A_soft_FPKM_df, "FT_MB_FB_X_A_soft_FPKM_df.csv", row.names = F)
write.csv(FT_MB_FB_X_A_hard_FPKM_wFC_df, "FT_MB_FB_X_A_hard_FPKM_wFC_df.csv", row.names = F)
write.csv(FT_MB_FB_X_A_soft_FPKM_wFC_df, "FT_MB_FB_X_A_soft_FPKM_wFC_df.csv", row.names = F)
write.csv(FT_MB_FB_X_A_hard_TPM_df, "FT_MB_FB_X_A_hard_TPM_df.csv", row.names = F)
write.csv(FT_MB_FB_X_A_soft_TPM_df, "FT_MB_FB_X_A_soft_TPM_df.csv", row.names = F)
write.csv(FT_MB_FB_X_A_hard_TPM_wFC_df, "FT_MB_FB_X_A_hard_TPM_wFC_df.csv", row.names = F)
write.csv(FT_MB_FB_X_A_soft_TPM_wFC_df, "FT_MB_FB_X_A_soft_TPM_wFC_df.csv", row.names = F)




##################################################################################################
### plot SB X A


prop_FB_plot <- function(df,ymax, tit_txt){
	ggplot(df, aes(factor(sp), prop_FB, fill = chr_type)) + 
	geom_bar(stat="identity", position = "dodge", colour="black") + 
	theme_bw() +
	scale_fill_manual(values=c("grey", "firebrick2")) +
	xlab ("Species") + 
	ylab ("Proportion of genes")  +
	scale_y_continuous(expand = c(0,0), labels = scales::number_format(accuracy = 0.01)) + 
	coord_cartesian(ylim=c(-0,ymax)) + 
	ggtitle(tit_txt)
}

prop_MB_plot <- function(df,ymax, tit_txt){
	ggplot(df, aes(factor(sp), prop_MB, fill = chr_type)) + 
	geom_bar(stat="identity", position = "dodge", colour="black") + 
	theme_bw() +
	scale_fill_manual(values=c("white", "royalblue2")) +
	xlab ("Species") + 
	ylab ("Proportion of genes")  +
	scale_y_continuous(expand = c(0,0), labels = scales::number_format(accuracy = 0.01)) + 
	coord_cartesian(ylim=c(-0,ymax)) + 
	ggtitle(tit_txt)
}


plot_SB_XA <-  function(FT_df, run_type){
	
	FT_df$N_FB_X       <-  as.numeric(as.character(FT_df$N_FB_X))
	FT_df$N_FB_A       <-  as.numeric(as.character(FT_df$N_FB_A))
	FT_df$N_MB_X       <-  as.numeric(as.character(FT_df$N_MB_X))
	FT_df$N_MB_A       <-  as.numeric(as.character(FT_df$N_MB_A))
	FT_df$N_All_exp_X      <- as.numeric(as.character(FT_df$N_All_exp_X))
	FT_df$N_All_exp_A      <- as.numeric(as.character(FT_df$N_All_exp_A))
	FT_df$FB_fisher_OR <- as.numeric(as.character(FT_df$FB_fisher_OR))
	FT_df$MB_fisher_OR <- as.numeric(as.character(FT_df$MB_fisher_OR))

	FT_df$tiss_ord <- ordered(FT_df$tiss, levels = c("HD", "LG", "RT"))

	prop_FB_X_A_df <- 
		as.data.frame(
		cbind(
		c(FT_df$N_FB_X, FT_df$N_FB_A),
		c(FT_df$N_All_exp_X, FT_df$N_All_exp_A),
		c(rep("X", length(FT_df[,1])),rep("A", length(FT_df[,1]))),
		c(as.character(FT_df$tiss),as.character(FT_df$tiss)),
		c(as.character(FT_df$sp),as.character(FT_df$sp))
		))



	colnames(prop_FB_X_A_df) <- c("N_FB", "N_all", "chr_type", "tiss", "sp")
	prop_FB_X_A_df$N_FB <- as.numeric(as.character(prop_FB_X_A_df$N_FB))
	prop_FB_X_A_df$N_all <- as.numeric(as.character(prop_FB_X_A_df$N_all))

	prop_FB_X_A_df$prop_FB <- prop_FB_X_A_df$N_FB / prop_FB_X_A_df$N_all 

	prop_FB_X_A_df_RT <- subset(prop_FB_X_A_df,prop_FB_X_A_df$tiss == "RT")
	prop_FB_X_A_df_HD <- subset(prop_FB_X_A_df,prop_FB_X_A_df$tiss == "HD")
	prop_FB_X_A_df_LG <- subset(prop_FB_X_A_df,prop_FB_X_A_df$tiss == "LG")

	prop_MB_X_A_df <- 
		as.data.frame(
		cbind(
		c(FT_df$N_MB_X, FT_df$N_MB_A),
		c(FT_df$N_All_exp_X, FT_df$N_All_exp_A),
		c(rep("X", length(FT_df[,1])),rep("A", length(FT_df[,1]))),
		c(as.character(FT_df$tiss),as.character(FT_df$tiss)),
		c(as.character(FT_df$sp),as.character(FT_df$sp))
		))

	colnames(prop_MB_X_A_df) <- c("N_MB", "N_all", "chr_type", "tiss", "sp")
	prop_MB_X_A_df$N_MB <- as.numeric(as.character(prop_MB_X_A_df$N_MB))
	prop_MB_X_A_df$N_all <- as.numeric(as.character(prop_MB_X_A_df$N_all))

	prop_MB_X_A_df$prop_MB <- prop_MB_X_A_df$N_MB / prop_MB_X_A_df$N_all 

	prop_MB_X_A_df_RT <- subset(prop_MB_X_A_df,prop_MB_X_A_df$tiss == "RT")
	prop_MB_X_A_df_HD <- subset(prop_MB_X_A_df,prop_MB_X_A_df$tiss == "HD")
	prop_MB_X_A_df_LG <- subset(prop_MB_X_A_df,prop_MB_X_A_df$tiss == "LG")

	#max(prop_FB_X_A_df_WB$prop_FB, prop_MB_X_A_df_WB$prop_MB) * 1.02
	
	P_FB_RT <- prop_FB_plot(prop_FB_X_A_df_RT, max(prop_FB_X_A_df_RT$prop_FB, prop_MB_X_A_df_RT$prop_MB) * 1.05, paste("RT | FB", run_type, sep = " | "))
	P_FB_HD <- prop_FB_plot(prop_FB_X_A_df_HD, max(prop_FB_X_A_df_HD$prop_FB, prop_MB_X_A_df_HD$prop_MB) * 1.05, paste("HD | FB", run_type, sep = " | "))
	P_FB_LG <- prop_FB_plot(prop_FB_X_A_df_LG, max(prop_FB_X_A_df_LG$prop_FB, prop_MB_X_A_df_LG$prop_MB) * 1.05, paste("LG | FB", run_type, sep = " | "))

	P_MB_RT <- prop_MB_plot(prop_MB_X_A_df_RT, max(prop_FB_X_A_df_RT$prop_FB, prop_MB_X_A_df_RT$prop_MB) * 1.05, paste("RT | MB", run_type, sep = " | "))
	P_MB_HD <- prop_MB_plot(prop_MB_X_A_df_HD, max(prop_FB_X_A_df_HD$prop_FB, prop_MB_X_A_df_HD$prop_MB) * 1.05, paste("HD | MB", run_type, sep = " | "))
	P_MB_LG <- prop_MB_plot(prop_MB_X_A_df_LG, max(prop_FB_X_A_df_LG$prop_FB, prop_MB_X_A_df_LG$prop_MB) * 1.05, paste("LG | MB", run_type, sep = " | "))
	
	P_out_prop <- plot_grid(P_FB_HD, P_MB_HD, P_FB_LG, P_MB_LG, P_FB_RT, P_MB_RT, ncol = 2)
	
	
	
	### plot OR
	
	OR_MB_P1 <- ggplot(FT_df, aes(tiss_ord, MB_fisher_OR, fill = sp)) + 
		geom_bar(stat="identity", position = "dodge", colour="black") + 
		theme_bw() +
		scale_fill_brewer(palette="Blues")+
		xlab ("Tissue") + 
		ylab ("Odds Ratio")  +
		geom_hline(yintercept=1) +
		scale_y_continuous(expand = c(0,0)) + 
		coord_cartesian(ylim=c(-0, max(FT_df$MB_fisher_OR, FT_df$FB_fisher_OR) * 1.02)) + ggtitle(run_type)
	
	OR_FB_P1 <-  ggplot(FT_df, aes(tiss_ord, FB_fisher_OR, fill = sp)) + 
		geom_bar(stat="identity", position = "dodge", colour="black") + 
		theme_bw() +
		scale_fill_brewer(palette="Reds")+
		xlab ("Tissue") + 
		ylab ("Odds Ratio")  +
		geom_hline(yintercept=1) +
		scale_y_continuous(expand = c(0,0)) + 
		coord_cartesian(ylim=c(-0, max(FT_df$MB_fisher_OR, FT_df$FB_fisher_OR) * 1.02)) + ggtitle(run_type)
	
	P_out_OR <- plot_grid(OR_FB_P1, OR_MB_P1)
	
	outlist <- list("P_out_prop" = P_out_prop, "P_out_OR" = P_out_OR)
	return(outlist)
	
}


## with fold-change cutoff

pdf("OR_SB_FT_MB_FB_X_A_hard_FPKM_wFC.pdf", width = 6, height = 4)
plot_SB_XA(FT_MB_FB_X_A_hard_FPKM_wFC_df, "hard_FPKM_wFC")$P_out_OR
dev.off()
getwd() ## where has my plot gone....

pdf("OR_SB_FT_MB_FB_X_A_soft_FPKM_wFC.pdf", width = 6, height = 4)
plot_SB_XA(FT_MB_FB_X_A_soft_FPKM_wFC_df, "soft_FPKM_wFC")$P_out_OR
dev.off()
getwd() ## where has my plot gone....

pdf("prop_SB_FT_MB_FB_X_A_hard_FPKM_wFC.pdf", width = 8, height = 9)
plot_SB_XA(FT_MB_FB_X_A_hard_FPKM_wFC_df, "hard_FPKM_wFC")$P_out_prop
dev.off()
getwd() ## where has my plot gone....

pdf("prop_SB_FT_MB_FB_X_A_soft_FPKM_wFC.pdf", width = 8, height = 9)
plot_SB_XA(FT_MB_FB_X_A_soft_FPKM_wFC_df, "soft_FPKM_wFC")$P_out_prop
dev.off()
getwd() ## where has my plot gone....

pdf("OR_SB_FT_MB_FB_X_A_hard_TPM_wFC.pdf", width = 6, height = 4)
plot_SB_XA(FT_MB_FB_X_A_hard_TPM_wFC_df, "hard_TPM_wFC")$P_out_OR
dev.off()
getwd() ## where has my plot gone....

pdf("OR_SB_FT_MB_FB_X_A_soft_TPM_wFC.pdf", width = 6, height = 4)
plot_SB_XA(FT_MB_FB_X_A_soft_TPM_wFC_df, "soft_TPM_wFC")$P_out_OR
dev.off()
getwd() ## where has my plot gone....

pdf("prop_SB_FT_MB_FB_X_A_hard_TPM_wFC.pdf", width = 8, height = 9)
plot_SB_XA(FT_MB_FB_X_A_hard_TPM_wFC_df, "hard_TPM_wFC")$P_out_prop
dev.off()
getwd() ## where has my plot gone....

pdf("prop_SB_FT_MB_FB_X_A_soft_TPM_wFC.pdf", width = 8, height = 9)
plot_SB_XA(FT_MB_FB_X_A_soft_TPM_wFC_df, "soft_TPM_wFC")$P_out_prop
dev.off()
getwd() ## where has my plot gone....

## no fold-change cutoff

pdf("OR_SB_FT_MB_FB_X_A_hard_FPKM.pdf", width = 6, height = 4)
plot_SB_XA(FT_MB_FB_X_A_hard_FPKM_df, "hard_FPKM")$P_out_OR
dev.off()
getwd() ## where has my plot gone....

pdf("OR_SB_FT_MB_FB_X_A_soft_FPKM.pdf", width = 6, height = 4)
plot_SB_XA(FT_MB_FB_X_A_soft_FPKM_df, "soft_FPKM")$P_out_OR
dev.off()
getwd() ## where has my plot gone....

pdf("prop_SB_FT_MB_FB_X_A_hard_FPKM.pdf", width = 8, height = 9)
plot_SB_XA(FT_MB_FB_X_A_hard_FPKM_df, "hard_FPKM")$P_out_prop
dev.off()
getwd() ## where has my plot gone....

pdf("prop_SB_FT_MB_FB_X_A_soft_FPKM.pdf", width = 8, height = 9)
plot_SB_XA(FT_MB_FB_X_A_soft_FPKM_df, "soft_FPKM")$P_out_prop
dev.off()
getwd() ## where has my plot gone....

pdf("OR_SB_FT_MB_FB_X_A_hard_TPM.pdf", width = 6, height = 4)
plot_SB_XA(FT_MB_FB_X_A_hard_TPM_df, "hard_TPM")$P_out_OR
dev.off()
getwd() ## where has my plot gone....

pdf("OR_SB_FT_MB_FB_X_A_soft_TPM.pdf", width = 6, height = 4)
plot_SB_XA(FT_MB_FB_X_A_soft_TPM_df, "soft_TPM")$P_out_OR
dev.off()
getwd() ## where has my plot gone....

pdf("prop_SB_FT_MB_FB_X_A_hard_TPM.pdf", width = 8, height = 9)
plot_SB_XA(FT_MB_FB_X_A_hard_TPM_df, "hard_TPM")$P_out_prop
dev.off()
getwd() ## where has my plot gone....

pdf("prop_SB_FT_MB_FB_X_A_soft_TPM.pdf", width = 8, height = 9)
plot_SB_XA(FT_MB_FB_X_A_soft_TPM_df, "soft_TPM")$P_out_prop
dev.off()
getwd() ## where has my plot gone....





#######################################################################################################################################################################
##### plot log2MF along chromosome 


plot_MF_along_chr = function(RT_df, HD_df, LG_df, lg_want, sp){
	
	# add Tiss = 
	RT_df$tiss = rep("RT", length(RT_df[,1]))
	HD_df$tiss = rep("HD", length(HD_df[,1]))	
	LG_df$tiss = rep("LG", length(LG_df[,1]))
	
	### join 
	
	join_df = as.data.frame(cbind(
		c(RT_df$tiss, HD_df$tiss, LG_df$tiss),
		c(as.character(RT_df$lg_pos), as.character(HD_df$lg_pos ), as.character(LG_df$lg_pos )),
		c(RT_df$gene_rel_midpoint, HD_df$gene_rel_midpoint, LG_df$gene_rel_midpoint),
		c(RT_df[,(length(RT_df) - 2)],   HD_df[,(length(HD_df) - 2)],  LG_df[,(length(LG_df) - 2)]),
		c(RT_df[,(length(RT_df) - 1)],   HD_df[,(length(HD_df) - 1)],  LG_df[,(length(LG_df) - 1)])
	))
	
	colnames(join_df) <- c("tiss", "lg", "gene_rel_midpoint", "F_AvFPKM", "M_AvFPKM")
	
	
	join_df$tiss = as.character(join_df$tiss)	
	#join_df$lg   = as.character(join_df$lg)	
	join_df$gene_rel_midpoint = as.numeric(as.character(join_df$gene_rel_midpoint))
	join_df$M_AvFPKM = as.numeric(as.character(join_df$M_AvFPKM))
	join_df$F_AvFPKM = as.numeric(as.character(join_df$F_AvFPKM))
	join_df$AvFPKM_log2MF <- log2(join_df$M_AvFPKM /join_df$F_AvFPKM)

	
	## sub to linkage group
	
	print(str(join_df))
	join_df_lg = subset(join_df, join_df$lg == lg_want)	
	
	print(head(join_df_lg))		
	P1 = ggplot(join_df_lg, aes(gene_rel_midpoint)) + 
 		 geom_point(aes(y = AvFPKM_log2MF, colour = tiss)) + geom_line(aes(y = AvFPKM_log2MF, colour = tiss))
 		 
 		 
 	## per tissue
 	
 	join_df_lg_RT = subset(join_df_lg, join_df_lg$tiss == "RT")
 	join_df_lg_HD = subset(join_df_lg, join_df_lg$tiss == "HD")
 	join_df_lg_LG = subset(join_df_lg, join_df_lg$tiss == "LG")	
 	
	P1_RT = ggplot(join_df_lg_RT, aes(gene_rel_midpoint)) + 
			theme_classic() +
 		 	geom_point(aes(y = AvFPKM_log2MF, colour = tiss)) + geom_line(aes(y = AvFPKM_log2MF, colour = tiss)) +
 		 	coord_cartesian(ylim=c(-5,5)) +
 		    geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
 		    scale_colour_manual(values=c("slateblue1")) + xlab("Position") + ylab("log2(M/F exp)")

	P1_HD = ggplot(join_df_lg_HD, aes(gene_rel_midpoint)) + 
			theme_classic() +
 		 	geom_point(aes(y = AvFPKM_log2MF, colour = tiss)) + geom_line(aes(y = AvFPKM_log2MF, colour = tiss)) +
 		 	coord_cartesian(ylim=c(-5,5)) +
 		    geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
 		    scale_colour_manual(values=c("forestgreen")) + xlab("Position") + ylab("log2(M/F exp)")


	P1_LG = ggplot(join_df_lg_LG, aes(gene_rel_midpoint)) + 
			theme_classic() +
 		 	geom_point(aes(y = AvFPKM_log2MF, colour = tiss)) + geom_line(aes(y = AvFPKM_log2MF, colour = tiss)) +
 		 	coord_cartesian(ylim=c(-5,5)) +
 		    geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
 		    scale_colour_manual(values=c("deepskyblue1")) + xlab("Position") + ylab("log2(M/F exp)")

	P1_together = plot_grid(P1_HD,P1_LG,P1_RT, ncol = 1, nrow = 3)

	title <- ggdraw() + draw_label(paste(sp, "|", lg_want), fontface='bold')
	P1_together2 <-	plot_grid(title, P1_together, ncol=1, rel_heights=c(0.1, 2)) 



 	return(P1_together2) 
}


### for X can plot the longest (in terms on N genes)
head(dat1_Tbi)
head(dat1_Tbi$lg)

sort(table(subset(dat1_Tbi, dat1_Tbi$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]
sort(table(subset(dat1_Tce, dat1_Tce$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]
sort(table(subset(dat1_Tcm, dat1_Tcm$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]
sort(table(subset(dat1_Tpa, dat1_Tpa$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]
sort(table(subset(dat1_Tps, dat1_Tps$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]

Tbi_percent_gene_top_3 <- sum(sort(table(subset(dat1_Tbi, dat1_Tbi$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]) / sum(sort(table(subset(dat1_Tbi, dat1_Tbi$lg =="lgX")$lg_pos),decreasing=TRUE)) * 100
Tce_percent_gene_top_3 <- sum(sort(table(subset(dat1_Tce, dat1_Tce$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]) / sum(sort(table(subset(dat1_Tce, dat1_Tce$lg =="lgX")$lg_pos),decreasing=TRUE)) * 100
Tcm_percent_gene_top_3 <- sum(sort(table(subset(dat1_Tcm, dat1_Tcm$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]) / sum(sort(table(subset(dat1_Tcm, dat1_Tcm$lg =="lgX")$lg_pos),decreasing=TRUE)) * 100
Tpa_percent_gene_top_3 <- sum(sort(table(subset(dat1_Tpa, dat1_Tpa$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]) / sum(sort(table(subset(dat1_Tpa, dat1_Tpa$lg =="lgX")$lg_pos),decreasing=TRUE)) * 100
Tps_percent_gene_top_3 <- sum(sort(table(subset(dat1_Tps, dat1_Tps$lg =="lgX")$lg_pos),decreasing=TRUE)[1:3]) / sum(sort(table(subset(dat1_Tps, dat1_Tps$lg =="lgX")$lg_pos),decreasing=TRUE)) * 100


Tbi_lgX1_MF <- plot_MF_along_chr(Tbi_RT_F_FPKM, Tbi_HD_F_FPKM, Tbi_LG_F_FPKM, "lgX_ordNA_scaf1319", "Tbi")
Tce_lgX1_MF <- plot_MF_along_chr(Tce_RT_F_FPKM, Tce_HD_F_FPKM, Tce_LG_F_FPKM, "lgX_ordNA_scaf1319", "Tce")
Tcm_lgX1_MF <- plot_MF_along_chr(Tcm_RT_F_FPKM, Tcm_HD_F_FPKM, Tcm_LG_F_FPKM, "lgX_ordNA_scaf1319", "Tcm")
Tpa_lgX1_MF <- plot_MF_along_chr(Tpa_RT_F_FPKM, Tpa_HD_F_FPKM, Tpa_LG_F_FPKM, "lgX_ordNA_scaf1319", "Tpa")
Tps_lgX1_MF <- plot_MF_along_chr(Tps_RT_F_FPKM, Tps_HD_F_FPKM, Tps_LG_F_FPKM, "lgX_ordNA_scaf1319", "Tps")

Tbi_lgX2_MF <- plot_MF_along_chr(Tbi_RT_F_FPKM, Tbi_HD_F_FPKM, Tbi_LG_F_FPKM, "lgX_ordNA_scaf2537", "Tbi")
Tce_lgX2_MF <- plot_MF_along_chr(Tce_RT_F_FPKM, Tce_HD_F_FPKM, Tce_LG_F_FPKM, "lgX_ordNA_scaf2537", "Tce")
Tcm_lgX2_MF <- plot_MF_along_chr(Tcm_RT_F_FPKM, Tcm_HD_F_FPKM, Tcm_LG_F_FPKM, "lgX_ordNA_scaf2537", "Tcm")
Tpa_lgX2_MF <- plot_MF_along_chr(Tpa_RT_F_FPKM, Tpa_HD_F_FPKM, Tpa_LG_F_FPKM, "lgX_ordNA_scaf2537", "Tpa")
Tps_lgX2_MF <- plot_MF_along_chr(Tps_RT_F_FPKM, Tps_HD_F_FPKM, Tps_LG_F_FPKM, "lgX_ordNA_scaf2537", "Tps")

Tbi_lgX3_MF <- plot_MF_along_chr(Tbi_RT_F_FPKM, Tbi_HD_F_FPKM, Tbi_LG_F_FPKM, "lgX_ordNA_scaf1899", "Tbi")
Tce_lgX3_MF <- plot_MF_along_chr(Tce_RT_F_FPKM, Tce_HD_F_FPKM, Tce_LG_F_FPKM, "lgX_ordNA_scaf1899", "Tce")
Tcm_lgX3_MF <- plot_MF_along_chr(Tcm_RT_F_FPKM, Tcm_HD_F_FPKM, Tcm_LG_F_FPKM, "lgX_ordNA_scaf699",  "Tcm")
Tpa_lgX3_MF <- plot_MF_along_chr(Tpa_RT_F_FPKM, Tpa_HD_F_FPKM, Tpa_LG_F_FPKM, "lgX_ordNA_scaf699",  "Tpa")
Tps_lgX3_MF <- plot_MF_along_chr(Tps_RT_F_FPKM, Tps_HD_F_FPKM, Tps_LG_F_FPKM, "lgX_ordNA_scaf699",  "Tps")


Tbi_percent_gene_top_3_title <- ggdraw() + draw_label(paste("% lgX genes =", Tbi_percent_gene_top_3), fontface='bold')
Tce_percent_gene_top_3_title <- ggdraw() + draw_label(paste("% lgX genes =", Tce_percent_gene_top_3), fontface='bold')
Tcm_percent_gene_top_3_title <- ggdraw() + draw_label(paste("% lgX genes =", Tcm_percent_gene_top_3), fontface='bold')
Tpa_percent_gene_top_3_title <- ggdraw() + draw_label(paste("% lgX genes =", Tpa_percent_gene_top_3), fontface='bold')
Tps_percent_gene_top_3_title <- ggdraw() + draw_label(paste("% lgX genes =", Tps_percent_gene_top_3), fontface='bold')
	
Tbi_lgX1to3_MF   <- plot_grid(Tbi_lgX1_MF, Tbi_lgX2_MF, Tbi_lgX3_MF, ncol=3) 
Tbi_lgX1to3_MF_t <-	plot_grid(Tbi_percent_gene_top_3_title, Tbi_lgX1to3_MF, ncol=1, rel_heights=c(0.1, 2)) 
Tce_lgX1to3_MF   <- plot_grid(Tce_lgX1_MF, Tce_lgX2_MF, Tce_lgX3_MF, ncol=3) 
Tce_lgX1to3_MF_t <-	plot_grid(Tce_percent_gene_top_3_title, Tce_lgX1to3_MF, ncol=1, rel_heights=c(0.1, 2)) 
Tcm_lgX1to3_MF   <- plot_grid(Tcm_lgX1_MF, Tcm_lgX2_MF, Tcm_lgX3_MF, ncol=3) 
Tcm_lgX1to3_MF_t <-	plot_grid(Tcm_percent_gene_top_3_title, Tcm_lgX1to3_MF, ncol=1, rel_heights=c(0.1, 2)) 
Tpa_lgX1to3_MF   <- plot_grid(Tpa_lgX1_MF, Tpa_lgX2_MF, Tpa_lgX3_MF, ncol=3) 
Tpa_lgX1to3_MF_t <-	plot_grid(Tpa_percent_gene_top_3_title, Tpa_lgX1to3_MF, ncol=1, rel_heights=c(0.1, 2)) 
Tps_lgX1to3_MF   <- plot_grid(Tps_lgX1_MF, Tps_lgX2_MF, Tps_lgX3_MF, ncol=3) 
Tps_lgX1to3_MF_t <-	plot_grid(Tps_percent_gene_top_3_title, Tps_lgX1to3_MF, ncol=1, rel_heights=c(0.1, 2)) 

Tbi_lg1_MF <- plot_MF_along_chr(Tbi_RT_F_FPKM, Tbi_HD_F_FPKM, Tbi_LG_F_FPKM, "lg1", "Tbi")
Tce_lg1_MF <- plot_MF_along_chr(Tce_RT_F_FPKM, Tce_HD_F_FPKM, Tce_LG_F_FPKM, "lg1", "Tce")
Tcm_lg1_MF <- plot_MF_along_chr(Tcm_RT_F_FPKM, Tcm_HD_F_FPKM, Tcm_LG_F_FPKM, "lg1", "Tcm")
Tpa_lg1_MF <- plot_MF_along_chr(Tpa_RT_F_FPKM, Tpa_HD_F_FPKM, Tpa_LG_F_FPKM, "lg1", "Tpa")
Tps_lg1_MF <- plot_MF_along_chr(Tps_RT_F_FPKM, Tps_HD_F_FPKM, Tps_LG_F_FPKM, "lg1", "Tps")



## Tbi lgX
pdf("Tbi_lgX_MF.pdf", width = 12, height = 7)
Tbi_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tbi_lgX_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tbi_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....


## Tce lgX
pdf("Tce_lgX_MF.pdf", width = 12, height = 7)
Tce_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tce_lgX_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tce_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....

## Tcm lgX
pdf("Tcm_lgX_MF.pdf", width = 12, height = 7)
Tcm_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tcm_lgX_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tcm_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....


## Tpa lgX
pdf("Tpa_lgX_MF.pdf", width = 12, height = 7)
Tpa_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tpa_lgX_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tpa_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....


## Tps lgX
pdf("Tps_lgX_MF.pdf", width = 12, height = 7)
Tps_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tps_lgX_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tps_lgX1to3_MF_t
dev.off()
getwd() ## where has my plot gone....


## Tbi lg1
pdf("Tbi_lg1_MF.pdf", width = 12, height = 7)
Tbi_lg1_MF 
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tbi_lg1_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tbi_lg1_MF 
dev.off()
getwd() ## where has my plot gone....


## Tce lg1
pdf("Tce_lg1_MF.pdf", width = 12, height = 7)
Tce_lg1_MF 
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tce_lg1_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tce_lg1_MF 
dev.off()
getwd() ## where has my plot gone....

## Tcm lg1
pdf("Tcm_lg1_MF.pdf", width = 12, height = 7)
Tcm_lg1_MF 
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tcm_lg1_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tcm_lg1_MF 
dev.off()
getwd() ## where has my plot gone....


## Tpa lg1
pdf("Tpa_lg1_MF.pdf", width = 12, height = 7)
Tpa_lg1_MF 
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tpa_lg1_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tpa_lg1_MF 
dev.off()
getwd() ## where has my plot gone....


## Tps lg1
pdf("Tps_lg1_MF.pdf", width = 12, height = 7)
Tps_lg1_MF 
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tps_lg1_MF.png", width = 12, height = 7, units = "in", bg = "white", res = 300)
Tps_lg1_MF 
dev.off()
getwd() ## where has my plot gone....



#######################################################################################################################################################################
##### Ortholog heatmap























































