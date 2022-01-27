# selection_on_the_X.R

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

dat1_Tbi_raw$lg <- ifelse(dat1_Tbi_raw$lg == "", NA, dat1_Tbi_raw$lg)
dat1_Tce_raw$lg <- ifelse(dat1_Tce_raw$lg == "", NA, dat1_Tce_raw$lg)
dat1_Tcm_raw$lg <- ifelse(dat1_Tcm_raw$lg == "", NA, dat1_Tcm_raw$lg)
dat1_Tpa_raw$lg <- ifelse(dat1_Tpa_raw$lg == "", NA, dat1_Tpa_raw$lg)
dat1_Tps_raw$lg <- ifelse(dat1_Tps_raw$lg == "", NA, dat1_Tps_raw$lg)

## from Expression_analyses.R
## sex-bias
Tbi_RT <- read.table ("../output/Exp_out/TTT_RT_Tbi_sex_bias_FPKM.csv", header = T, sep = ',')
Tce_RT <- read.table ("../output/Exp_out/TTT_RT_Tce_sex_bias_FPKM.csv", header = T, sep = ',')
Tcm_RT <- read.table ("../output/Exp_out/TTT_RT_Tcm_sex_bias_FPKM.csv", header = T, sep = ',')
Tpa_RT <- read.table ("../output/Exp_out/TTT_RT_Tpa_sex_bias_FPKM.csv", header = T, sep = ',')
Tps_RT <- read.table ("../output/Exp_out/TTT_RT_Tps_sex_bias_FPKM.csv", header = T, sep = ',')

Tbi_HD <- read.table ("../output/Exp_out/TTT_HD_Tbi_sex_bias_FPKM.csv", header = T, sep = ',')
Tce_HD <- read.table ("../output/Exp_out/TTT_HD_Tce_sex_bias_FPKM.csv", header = T, sep = ',')
Tcm_HD <- read.table ("../output/Exp_out/TTT_HD_Tcm_sex_bias_FPKM.csv", header = T, sep = ',')
Tpa_HD <- read.table ("../output/Exp_out/TTT_HD_Tpa_sex_bias_FPKM.csv", header = T, sep = ',')
Tps_HD <- read.table ("../output/Exp_out/TTT_HD_Tps_sex_bias_FPKM.csv", header = T, sep = ',')

Tbi_LG <- read.table ("../output/Exp_out/TTT_LG_Tbi_sex_bias_FPKM.csv", header = T, sep = ',')
Tce_LG <- read.table ("../output/Exp_out/TTT_LG_Tce_sex_bias_FPKM.csv", header = T, sep = ',')
Tcm_LG <- read.table ("../output/Exp_out/TTT_LG_Tcm_sex_bias_FPKM.csv", header = T, sep = ',')
Tpa_LG <- read.table ("../output/Exp_out/TTT_LG_Tpa_sex_bias_FPKM.csv", header = T, sep = ',')
Tps_LG <- read.table ("../output/Exp_out/TTT_LG_Tps_sex_bias_FPKM.csv", header = T, sep = ',')

## average exp
Tbi_RT_FKPM <- read.table ("../output/Exp_out/Tbi_RT_F_FPKM.csv", header = T, sep = ',')
Tce_RT_FKPM <- read.table ("../output/Exp_out/Tce_RT_F_FPKM.csv", header = T, sep = ',')
Tcm_RT_FKPM <- read.table ("../output/Exp_out/Tcm_RT_F_FPKM.csv", header = T, sep = ',')
Tpa_RT_FKPM <- read.table ("../output/Exp_out/Tpa_RT_F_FPKM.csv", header = T, sep = ',')
Tps_RT_FKPM <- read.table ("../output/Exp_out/Tps_RT_F_FPKM.csv", header = T, sep = ',')

Tbi_HD_FKPM <- read.table ("../output/Exp_out/Tbi_HD_F_FPKM.csv", header = T, sep = ',')
Tce_HD_FKPM <- read.table ("../output/Exp_out/Tce_HD_F_FPKM.csv", header = T, sep = ',')
Tcm_HD_FKPM <- read.table ("../output/Exp_out/Tcm_HD_F_FPKM.csv", header = T, sep = ',')
Tpa_HD_FKPM <- read.table ("../output/Exp_out/Tpa_HD_F_FPKM.csv", header = T, sep = ',')
Tps_HD_FKPM <- read.table ("../output/Exp_out/Tps_HD_F_FPKM.csv", header = T, sep = ',')

Tbi_LG_FKPM <- read.table ("../output/Exp_out/Tbi_LG_F_FPKM.csv", header = T, sep = ',')
Tce_LG_FKPM <- read.table ("../output/Exp_out/Tce_LG_F_FPKM.csv", header = T, sep = ',')
Tcm_LG_FKPM <- read.table ("../output/Exp_out/Tcm_LG_F_FPKM.csv", header = T, sep = ',')
Tpa_LG_FKPM <- read.table ("../output/Exp_out/Tpa_LG_F_FPKM.csv", header = T, sep = ',')
Tps_LG_FKPM <- read.table ("../output/Exp_out/Tps_LG_F_FPKM.csv", header = T, sep = ',')


#### GC and CUB from GC_CUB.R
GC_all  <- read.table ("../output/sel_out/GC.csv", header = T, sep = ',')
ENC_all <- read.table ("../output/sel_out/ENC.csv", header = T, sep = ',')

head(GC_all)
head(ENC_all)


dat1_orth_raw <- read.table ("TbiTceTcmTpaTps_orth_1000_chr_info_and_counts.csv", header = T, sep = ',')

#### change wd to output folder
dir.create("../output/sel_out")
setwd("../output/sel_out")


#########################################################################################################################
## gene info into dicts
## Note hash tables in R studio give an error:
# Error: protect(): protection stack overflow
# Error: no more error handlers available (recursive errors?); invoking 'abort' restart
# this is a bug https://github.com/rstudio/rstudio/issues/5546


all_raw_dfs <- c("dat1_Tbi_raw", "dat1_Tce_raw", "dat1_Tcm_raw", "dat1_Tpa_raw", "dat1_Tps_raw")

levels(as.factor(dat1_Tbi_raw$lg))



Orth_dict <- hash()
hard_chr_dict <- hash()
soft_chr_dict <- hash()
all_orths = c()
gene_to_orth_dict <- hash()
gene_to_lg_dict <- hash()

for(d in all_raw_dfs){
	sp <- strsplit(d, "_")[[1]][2]
	test_df <- eval(parse(text=paste(d,sep='')))

	for(i in seq(1:length(test_df[,1]))){
		gene_n <- test_df$gene_id[i]
		HOG_n  <- test_df$HOG[i]
		hard_n <- test_df$chr_hard[i]
		soft_n <- test_df$chr_soft[i]
		lg_n   <- test_df$lg[i]
		if(! is.na(HOG_n)){
			Orth_dict[[paste(sp, HOG_n, sep = "___")]] <- gene_n
			all_orths = c(all_orths, HOG_n)
			gene_to_orth_dict[[gene_n]] <- HOG_n
		}
	  
		if(! is.na(hard_n)){
		  
			hard_chr_dict[[gene_n]] <- hard_n
		}

		if(! is.na(soft_n)){
			soft_chr_dict[[gene_n]] <- soft_n
		}
		if(! is.na(lg_n)){
		  gene_to_lg_dict[[gene_n]] <- lg_n
		}
	}
}

all_orths <- unique(all_orths)

soft_chr_dict[["TBI_07206"]]
Orth_dict[["Tpa___HOG_17469"]]
length(Orth_dict)
length(hard_chr_dict)
length(soft_chr_dict)


######################################################################################################
### SB gene data

all_RT_dfs <- c("Tbi_RT", "Tce_RT", "Tcm_RT", "Tpa_RT", "Tps_RT")

RT_SB_FDR_dict <- hash()
RT_SB_logFC_dict <- hash()

for(d in all_RT_dfs){
  test_df <- eval(parse(text=paste(d,sep='')))
  print(head(test_df))
  
  for(i in seq(1:length(test_df[,1]))){
    gene_n   <- test_df$genes[i]
    logFC_n  <- test_df$logFC[i]
    FDR_n    <- test_df$FDR[i]
    
    if(! is.na(FDR_n)){
      RT_SB_FDR_dict[[gene_n]] <- FDR_n
    }
    
    if(! is.na(logFC_n)){
      RT_SB_logFC_dict[[gene_n]] <- logFC_n
    }
  }
}


all_HD_dfs <- c("Tbi_HD", "Tce_HD", "Tcm_HD", "Tpa_HD", "Tps_HD")

HD_SB_FDR_dict <- hash()
HD_SB_logFC_dict <- hash()

for(d in all_HD_dfs){
  test_df <- eval(parse(text=paste(d,sep='')))
  
  for(i in seq(1:length(test_df[,1]))){
    gene_n   <- test_df$genes[i]
    logFC_n  <- test_df$logFC[i]
    FDR_n    <- test_df$FDR[i]
    
    if(! is.na(FDR_n)){
      HD_SB_FDR_dict[[gene_n]] <- FDR_n
    }
    
    if(! is.na(logFC_n)){
      HD_SB_logFC_dict[[gene_n]] <- logFC_n
    }
  }
}


all_LG_dfs <- c("Tbi_LG", "Tce_LG", "Tcm_LG", "Tpa_LG", "Tps_LG")

LG_SB_FDR_dict <- hash()
LG_SB_logFC_dict <- hash()

for(d in all_LG_dfs){
  test_df <- eval(parse(text=paste(d,sep='')))
  
  for(i in seq(1:length(test_df[,1]))){
    gene_n   <- test_df$genes[i]
    logFC_n  <- test_df$logFC[i]
    FDR_n    <- test_df$FDR[i]
    
    if(! is.na(FDR_n)){
      LG_SB_FDR_dict[[gene_n]] <- FDR_n
    }
    
    if(! is.na(logFC_n)){
      LG_SB_logFC_dict[[gene_n]] <- logFC_n
    }
  }
}

RT_SB_FDR_dict[["TPS_13859"]]
RT_SB_logFC_dict[["TPA_06181"]]


##### exp lev

all_RT_FKPM_dict_dfs <- c("Tbi_RT_FKPM ", "Tce_RT_FKPM ", "Tcm_RT_FKPM ", "Tpa_RT_FKPM ", "Tps_RT_FKPM ")

RT_SF_FKPM_dict <- hash()
RT_SM_FKPM_dict <- hash()

for(d in all_RT_FKPM_dict_dfs ){
  test_df <- eval(parse(text=paste(d,sep='')))
  colnames(test_df) <- gsub("T.._", "", colnames(test_df))
  print(head(test_df))
  
  for(i in seq(1:length(test_df[,1]))){
    gene_n   <- test_df$gene_id[i]
    SF_n     <- test_df$SF_RT_meanFPKM[i]
    SM_n     <- test_df$SM_RT_meanFPKM[i]
    
    if(! is.na(SF_n)){
      RT_SF_FKPM_dict[[gene_n]] <- SF_n
    }
    
    if(! is.na(SM_n)){
      RT_SM_FKPM_dict[[gene_n]] <- SM_n
    }
    
  }
}



all_HD_FKPM_dict_dfs <- c("Tbi_HD_FKPM ", "Tce_HD_FKPM ", "Tcm_HD_FKPM ", "Tpa_HD_FKPM ", "Tps_HD_FKPM ")

HD_SF_FKPM_dict <- hash()
HD_SM_FKPM_dict <- hash()

for(d in all_HD_FKPM_dict_dfs ){
  test_df <- eval(parse(text=paste(d,sep='')))
  colnames(test_df) <- gsub("T.._", "", colnames(test_df))
  print(head(test_df))
  
  for(i in seq(1:length(test_df[,1]))){
    gene_n   <- test_df$gene_id[i]
    SF_n     <- test_df$SF_HD_meanFPKM[i]
    SM_n     <- test_df$SM_HD_meanFPKM[i]
    
    if(! is.na(SF_n)){
      HD_SF_FKPM_dict[[gene_n]] <- SF_n
    }
    
    if(! is.na(SM_n)){
      HD_SM_FKPM_dict[[gene_n]] <- SM_n
    }
    
  }
}

all_LG_FKPM_dict_dfs <- c("Tbi_LG_FKPM ", "Tce_LG_FKPM ", "Tcm_LG_FKPM ", "Tpa_LG_FKPM ", "Tps_LG_FKPM ")

LG_SF_FKPM_dict <- hash()
LG_SM_FKPM_dict <- hash()

for(d in all_LG_FKPM_dict_dfs ){
  test_df <- eval(parse(text=paste(d,sep='')))
  colnames(test_df) <- gsub("T.._", "", colnames(test_df))
  print(head(test_df))
  
  for(i in seq(1:length(test_df[,1]))){
    gene_n   <- test_df$gene_id[i]
    SF_n     <- test_df$SF_LG_meanFPKM[i]
    SM_n     <- test_df$SM_LG_meanFPKM[i]
    
    if(! is.na(SF_n)){
      LG_SF_FKPM_dict[[gene_n]] <- SF_n
    }
    
    if(! is.na(SM_n)){
      LG_SM_FKPM_dict[[gene_n]] <- SM_n
    }
    
  }
}



RT_SF_FKPM_dict[["TPS_13859"]]
RT_SM_FKPM_dict[["TPS_13859"]]
RT_SF_FKPM_dict[["TPS_00008"]]
RT_SM_FKPM_dict[["TPS_00008"]]

HD_SM_FKPM_dict[["TPS_00008"]]
LG_SM_FKPM_dict[["TPS_00008"]]



#### GC data
GC_dict <- hash()

head(GC_all )

for(i in seq(1:length(GC_all[,1]))){
    gene_n   <- GC_all$g_name[i]
    GC_n     <- GC_all$GC_all[i]
    
    if(! is.na(GC_n)){
      GC_dict[[gene_n]] <- GC_n
    }
}


GC_dict[["Tbi_HOG_746"]]

#### CUB (ENC) data

ENC_dict <- hash()
for(i in seq(1:length(ENC_all[,1]))){
  gene_n   <- ENC_all$Gene[i]
  ENC_n   <- ENC_all$ENC[i]
  
  if(! is.na(ENC_n)){
    ENC_dict[[gene_n]] <- ENC_n
  }
}

ENC_dict[["Tbi_HOG_746"]]








######################################################################################################
### positive selection data

pos_sel_dat <- read.table("../../selection/timema_543_branches_with-ncat-codon-rate_sites_with_h0_wgenename.tsv", sep = "\t", header = T)
pos_sel_dat$gene_name<- as.character(pos_sel_dat$gene_name )

head(pos_sel_dat)

### add chr and SB to positive selection results  

positive_sel_by_chr <- function(pos_sel_df){
  
  print(length(pos_sel_df[,1]))
  chr_soft_class <- c()
  chr_hard_class <- c()
  RT_logFC       <- c()
  RT_FDR         <- c()
  HD_logFC       <- c()
  HD_FDR         <- c()
  LG_logFC       <- c()
  LG_FDR         <- c()
  RT_SF_FKPM     <- c()
  RT_SM_FKPM     <- c()  
  HD_SF_FKPM     <- c()
  HD_SM_FKPM     <- c()  
  LG_SF_FKPM     <- c()
  LG_SM_FKPM     <- c()  
  GC             <- c()   
  ENC            <- c()
  lg             <- c()

  for(i in seq(1, length(pos_sel_df[,1]))){
    gene_n <- pos_sel_df$gene_name[i]
    hog_name <- paste(pos_sel_df$branch_name[i], pos_sel_df$gene[i], sep = "_")
    chr_class_s <- soft_chr_dict[[gene_n]]
    if(length(chr_class_s) == 0){chr_class_s = NA}   
    chr_soft_class <- c(chr_soft_class, chr_class_s)
    
    chr_class_h <- hard_chr_dict[[gene_n]]
    if(length(chr_class_h) == 0){chr_class_h = NA}   
    chr_hard_class <- c(chr_hard_class, chr_class_h)

    RT_logFC_n <- RT_SB_logFC_dict[[gene_n]]
    if(length(RT_logFC_n) == 0){RT_logFC_n = NA}   
    RT_logFC <- c(RT_logFC, RT_logFC_n)
    
    RT_FDR_n <- RT_SB_FDR_dict[[gene_n]]
    if(length(RT_FDR_n) == 0){RT_FDR_n = NA}   
    RT_FDR <- c(RT_FDR, RT_FDR_n)
    
    HD_logFC_n <- HD_SB_logFC_dict[[gene_n]]
    if(length(HD_logFC_n) == 0){HD_logFC_n = NA}   
    HD_logFC <- c(HD_logFC, HD_logFC_n)

    HD_FDR_n <- HD_SB_FDR_dict[[gene_n]]
    if(length(HD_FDR_n) == 0){HD_FDR_n = NA}   
    HD_FDR <- c(HD_FDR, HD_FDR_n)
    
    LG_logFC_n <- LG_SB_logFC_dict[[gene_n]]
    if(length(LG_logFC_n) == 0){LG_logFC_n = NA}   
    LG_logFC <- c(LG_logFC, LG_logFC_n)
    
    LG_FDR_n <- LG_SB_FDR_dict[[gene_n]]
    if(length(LG_FDR_n) == 0){LG_FDR_n = NA}   
    LG_FDR <- c(LG_FDR, LG_FDR_n)
    
    RT_SF_FKPM_n <- RT_SF_FKPM_dict[[gene_n]]
    if(length(RT_SF_FKPM_n) == 0){RT_SF_FKPM_n = NA}   
    RT_SF_FKPM <- c(RT_SF_FKPM, RT_SF_FKPM_n)
    
    RT_SM_FKPM_n <- RT_SM_FKPM_dict[[gene_n]]
    if(length(RT_SM_FKPM_n) == 0){RT_SM_FKPM_n = NA}   
    RT_SM_FKPM <- c(RT_SM_FKPM, RT_SM_FKPM_n)

    HD_SF_FKPM_n <- HD_SF_FKPM_dict[[gene_n]]
    if(length(HD_SF_FKPM_n) == 0){HD_SF_FKPM_n = NA}   
    HD_SF_FKPM <- c(HD_SF_FKPM, HD_SF_FKPM_n)
    
    HD_SM_FKPM_n <- HD_SM_FKPM_dict[[gene_n]]
    if(length(HD_SM_FKPM_n) == 0){HD_SM_FKPM_n = NA}   
    HD_SM_FKPM <- c(HD_SM_FKPM, HD_SM_FKPM_n)    
    
    LG_SF_FKPM_n <- LG_SF_FKPM_dict[[gene_n]]
    if(length(LG_SF_FKPM_n) == 0){LG_SF_FKPM_n = NA}   
    LG_SF_FKPM <- c(LG_SF_FKPM, LG_SF_FKPM_n)
    
    LG_SM_FKPM_n <- LG_SM_FKPM_dict[[gene_n]]
    if(length(LG_SM_FKPM_n) == 0){LG_SM_FKPM_n = NA}   
    LG_SM_FKPM <- c(LG_SM_FKPM, LG_SM_FKPM_n)
    
    GC_n <- GC_dict[[hog_name]]
    if(length(GC_n) == 0){GC_n = NA}   
    GC <- c(GC, GC_n)   
    
    ENC_n <- ENC_dict[[hog_name]]
    if(length(ENC_n) == 0){ENC_n = NA}   
    ENC <- c(ENC, ENC_n)    
    
    lg_n <- gene_to_lg_dict[[gene_n]]
    if(length(lg_n) == 0){lg_n = NA}   
    lg <- c(lg, lg_n) 
  }

  pos_sel_df$chr_soft_class <- chr_soft_class
  pos_sel_df$chr_hard_class <- chr_hard_class
  pos_sel_df$RT_logFC <- RT_logFC
  pos_sel_df$RT_FDR   <- RT_FDR 
  pos_sel_df$HD_logFC <- HD_logFC
  pos_sel_df$HD_FDR   <- HD_FDR 
  pos_sel_df$LG_logFC <- LG_logFC
  pos_sel_df$LG_FDR   <- LG_FDR 
  pos_sel_df$RT_SF_FKPM <- RT_SF_FKPM  
  pos_sel_df$RT_SM_FKPM <- RT_SM_FKPM  
  pos_sel_df$HD_SF_FKPM <- HD_SF_FKPM  
  pos_sel_df$HD_SM_FKPM <- HD_SM_FKPM  
  pos_sel_df$LG_SF_FKPM <- LG_SF_FKPM  
  pos_sel_df$LG_SM_FKPM <- LG_SM_FKPM  
  pos_sel_df$GC         <- GC
  pos_sel_df$ENC        <- ENC
  pos_sel_df$lg         <- lg 
  pos_sel_df <- subset(pos_sel_df, !is.na(pos_sel_df$chr_soft_class))
  
  print(length(pos_sel_df[,1]))  
  return(pos_sel_df)
}

pos_sel_dat_2 <- positive_sel_by_chr(pos_sel_dat)
head(pos_sel_dat_2)

##### are genes showing +ve sel overrep on the X?

test_possel_overrep_X <- function(df_all, pos_q_threh, chr_type){
  df_all$qvalue <- as.numeric(df_all$qvalue)
  print(length(df_all[,1]))
  df_all <- subset(df_all, eval(parse(text=paste('df_all', '$chr_', chr_type, '_class', sep = ''))) != "NA")
  print(length(df_all[,1]))
  
  df_X 	<- subset(df_all, eval(parse(text=paste('df_all', '$chr_', chr_type, '_class', sep = ''))) == "X")
  print(length(df_X[,1]))
  
  X_genes_sel      = length(subset(df_X, df_X$qvalue  <  pos_q_threh)[,1])
  X_genes_no_sel   = length(subset(df_X, df_X$qvalue  >= pos_q_threh)[,1])
  
  All_genes_sel    = length(subset(df_all,   df_all$qvalue  <  pos_q_threh)[,1])
  All_genes_no_sel = length(subset(df_all,   df_all$qvalue  >= pos_q_threh)[,1])
  
  FT_mat <- matrix(c(X_genes_sel,(All_genes_sel - X_genes_sel), X_genes_no_sel,(All_genes_no_sel - X_genes_no_sel)), nrow = 2)
  
  print(fisher.test(FT_mat, alternative="two.sided")) 
  
  ## plot props by sp
  
  get_prop <- function(sp){
    X_genes_sel_sp      = length(subset(df_X, df_X$qvalue  <  pos_q_threh & df_X$branch_name == sp)[,1])
    X_genes_no_sel_sp   = length(subset(df_X, df_X$qvalue  >= pos_q_threh & df_X$branch_name == sp)[,1])
    All_genes_sel_sp    = length(subset(df_all,   df_all$qvalue  <  pos_q_threh & df_all$branch_name == sp)[,1])
    All_genes_no_sel_sp = length(subset(df_all,   df_all$qvalue  >= pos_q_threh & df_all$branch_name == sp)[,1])  
    sp_X_prop_sel <- X_genes_sel_sp / (X_genes_sel_sp + X_genes_no_sel_sp) 
    sp_A_prop_sel <- (All_genes_sel_sp - X_genes_sel_sp) / ((All_genes_no_sel_sp - X_genes_no_sel_sp) + (All_genes_sel_sp - X_genes_sel_sp))
  
    return(as.data.frame(rbind(c(sp, "X", sp_X_prop_sel), c(sp, "A", sp_A_prop_sel))))
    
    }
  
  prop_sel_df <- as.data.frame(rbind(get_prop("Tbi"), get_prop("Tce"),get_prop("Tcm"),get_prop("Tpa"),get_prop("Tps")))
  colnames(prop_sel_df) <- c("sp", "chr", "prop_possel")
  prop_sel_df$prop_possel <- as.numeric(prop_sel_df$prop_possel)
  
  print(str(prop_sel_df))

  max_y = max(prop_sel_df$prop_possel * 1.05)
  
  if(chr_type == "soft"){
    P1b <- ggplot(prop_sel_df, aes(x = factor(sp), y = prop_possel, fill = chr)) + 
      geom_col(width = 0.5, colour="black", position=position_dodge(width=0.6)) + 
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      
      scale_fill_manual(values = c("darkgrey", "darkorange2")) + 
      xlab ("Species pair") + 
      ylab ("Prop of positive selected genes")  + 
      ggtitle(paste(chr_type, ", qval thresh = ",pos_q_threh))  + ylim(0,max_y )
  }

  if(chr_type == "hard"){
    P1b <- ggplot(prop_sel_df, aes(x = factor(sp), y = prop_possel, fill = chr)) + 
      geom_col(width = 0.5, colour="black", position=position_dodge(width=0.6)) + 
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      
      scale_fill_manual(values = c("darkgrey", "red3")) + 
      xlab ("Species pair") + 
      ylab ("Prop  of positive selected genes")  + 
      ggtitle(paste(chr_type, ", qval thresh = ",pos_q_threh))  + ylim(0,max_y )
  }    
  return(P1b)
}


test_possel_overrep_X(pos_sel_dat_2, 0.05, "soft")


pdf(paste("prop_possel_soft" ,".pdf", sep = ""), width = 6, height = 8)
test_possel_overrep_X(pos_sel_dat_2, 0.05, "soft")
dev.off()
getwd() ## where has my plot gone....

pdf(paste("prop_possel_hard" ,".pdf", sep = ""), width = 6, height = 8)
test_possel_overrep_X(pos_sel_dat_2, 0.05, "hard")
dev.off()
getwd() ## where has my plot gone....


##########################################################################################################
### dNdS

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

### get dN/dS for all non-selected sites (p0 + p1) 

## soft
pos_sel_dat_2$omega0_1 <- (pos_sel_dat_2$omega0 * pos_sel_dat_2$p0 + 1 * pos_sel_dat_2$p1) / (pos_sel_dat_2$p0 + pos_sel_dat_2$p1)
pos_sel_dat_2$sp_soft  <- paste(pos_sel_dat_2$branch_name, pos_sel_dat_2$chr_soft_class, sep = "_") 

## get means
omega0_1_means_soft  <- summarySE(pos_sel_dat_2, measurevar="omega0_1", groupvars=c("sp_soft"))
omega0_1_means_soft$upper <- omega0_1_means_soft$omega0_1 +  omega0_1_means_soft$se
omega0_1_means_soft$lower <- omega0_1_means_soft$omega0_1 -  omega0_1_means_soft$se
omega0_1_means_soft$chr_soft <- str_split_fixed(omega0_1_means_soft$sp_soft, "_", 2)[,2]

## plot
P_mean_omega0_1_soft <- ggplot() + 
  geom_errorbar(data=omega0_1_means_soft, mapping=aes(x=sp_soft, ymin=upper, ymax=lower,  color= chr_soft), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=omega0_1_means_soft, mapping=aes(x = sp_soft, y=omega0_1,  fill=chr_soft), size=4, shape=21) + 
  scale_color_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE") +
  scale_fill_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE")
  
pdf(paste("P_mean_omega0_1_soft" ,".pdf", sep = ""), width = 6, height = 8)
P_mean_omega0_1_soft
dev.off()
getwd() ## where has my plot gone....


## test wilcox

wilcox_omega0_1_soft <- as.data.frame(cbind(
c("Tbi", "Tce", "Tcm", "Tpa", "Tps"),
c(
wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_A")$omega0_1, paired = FALSE)$p.value,
wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_A")$omega0_1, paired = FALSE)$p.value,
wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_A")$omega0_1, paired = FALSE)$p.value,
wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_A")$omega0_1, paired = FALSE)$p.value,
wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_A")$omega0_1, paired = FALSE)$p.value
)))
colnames(wilcox_omega0_1_soft) <- c("sp", "wilcox_p")
wilcox_omega0_1_soft$wilcox_p <- as.numeric(wilcox_omega0_1_soft$wilcox_p)
wilcox_omega0_1_soft$wilcox_FDR <- p.adjust(wilcox_omega0_1_soft$wilcox_p, method = "BH")
write.csv(wilcox_omega0_1_soft, "wilcox_omega0_1_soft.csv", row.names = FALSE)

## hard
pos_sel_dat_2$omega0_1 <- (pos_sel_dat_2$omega0 * pos_sel_dat_2$p0 + 1 * pos_sel_dat_2$p1) / (pos_sel_dat_2$p0 + pos_sel_dat_2$p1)
pos_sel_dat_2$sp_hard  <- paste(pos_sel_dat_2$branch_name, pos_sel_dat_2$chr_hard_class, sep = "_") 

## get means
omega0_1_means_hard  <- summarySE(pos_sel_dat_2, measurevar="omega0_1", groupvars=c("sp_hard"))
omega0_1_means_hard$upper <- omega0_1_means_hard$omega0_1 +  omega0_1_means_hard$se
omega0_1_means_hard$lower <- omega0_1_means_hard$omega0_1 -  omega0_1_means_hard$se
omega0_1_means_hard$chr_hard <- str_split_fixed(omega0_1_means_hard$sp_hard, "_", 2)[,2]

## plot
P_mean_omega0_1_hard <- ggplot() + 
  geom_errorbar(data=omega0_1_means_hard, mapping=aes(x=sp_hard, ymin=upper, ymax=lower,  color= chr_hard), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=omega0_1_means_hard, mapping=aes(x = sp_hard, y=omega0_1,  fill=chr_hard), size=4, shape=21) + 
  scale_color_manual(values=c("darkgrey", "red3")) + ggtitle("hard chr, mean +/- SE") +
  scale_fill_manual(values=c("darkgrey", "red3")) + ggtitle("hard chr, mean +/- SE")

pdf(paste("P_mean_omega0_1_hard" ,".pdf", sep = ""), width = 6, height = 8)
P_mean_omega0_1_hard
dev.off()
getwd() ## where has my plot gone....

## test wilcox

wilcox_omega0_1_hard <- as.data.frame(cbind(
  c("Tbi", "Tce", "Tcm", "Tpa", "Tps"),
  c(
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tbi_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tbi_A")$omega0_1, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tce_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tce_A")$omega0_1, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tcm_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tcm_A")$omega0_1, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tpa_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tpa_A")$omega0_1, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tps_X")$omega0_1, subset(pos_sel_dat_2, pos_sel_dat_2$sp_hard == "Tps_A")$omega0_1, paired = FALSE)$p.value
  )))
colnames(wilcox_omega0_1_hard) <- c("sp", "wilcox_p")
wilcox_omega0_1_hard$wilcox_p <- as.numeric(wilcox_omega0_1_hard$wilcox_p)
wilcox_omega0_1_hard$wilcox_FDR <- p.adjust(wilcox_omega0_1_hard$wilcox_p, method = "BH")
write.csv(wilcox_omega0_1_hard, "wilcox_omega0_1_hard.csv", row.names = FALSE)

write.csv(pos_sel_dat_2, "pos_sel_dat_2_ttt.csv")





##################################################################################################################
###


### get mean expersion

## make missing expression 0
GE_all <- as.data.frame(cbind(pos_sel_dat_2$RT_SF_FKPM, pos_sel_dat_2$RT_SM_FKPM, pos_sel_dat_2$HD_SF_FKPM, pos_sel_dat_2$HD_SM_FKPM, pos_sel_dat_2$LG_SF_FKPM, pos_sel_dat_2$LG_SM_FKPM ))
GE_all[is.na(GE_all )] <- 0
pos_sel_dat_2$mean_FKPM <- rowMeans(GE_all)
head(pos_sel_dat_2)



pos_sel_dat_2_A <- subset(pos_sel_dat_2, pos_sel_dat_2$chr_soft_class == "A")
pos_sel_dat_2_X <- subset(pos_sel_dat_2, pos_sel_dat_2$chr_soft_class == "X")

cor.test(pos_sel_dat_2_A$omega0_1, pos_sel_dat_2_A$GC, method = "spearman") ### neg corr - rho -0.2580059  p-value < 2.2e-16
cor.test(pos_sel_dat_2_X$omega0_1, pos_sel_dat_2_X$GC, method = "spearman") ### neg corr - but much weaker rho -0.04380179 p-value = 0.05062

pos_sel_dat_2_Tbi_X <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_X")
pos_sel_dat_2_Tce_X <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_X")
pos_sel_dat_2_Tcm_X <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_X")
pos_sel_dat_2_Tpa_X <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_X")
pos_sel_dat_2_Tps_X <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_X")
pos_sel_dat_2_Tbi_A <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_A")
pos_sel_dat_2_Tce_A <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_A")
pos_sel_dat_2_Tcm_A <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_A")
pos_sel_dat_2_Tpa_A <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_A")
pos_sel_dat_2_Tps_A <- subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_A")

cor.test(pos_sel_dat_2_Tbi_X$omega0_1, pos_sel_dat_2_Tbi_X$GC, method = "spearman")
cor.test(pos_sel_dat_2_Tce_X$omega0_1, pos_sel_dat_2_Tce_X$GC, method = "spearman")
cor.test(pos_sel_dat_2_Tcm_X$omega0_1, pos_sel_dat_2_Tcm_X$GC, method = "spearman")
cor.test(pos_sel_dat_2_Tpa_X$omega0_1, pos_sel_dat_2_Tpa_X$GC, method = "spearman")
cor.test(pos_sel_dat_2_Tps_X$omega0_1, pos_sel_dat_2_Tps_X$GC, method = "spearman")
cor.test(pos_sel_dat_2_Tbi_A$omega0_1, pos_sel_dat_2_Tbi_A$GC, method = "spearman")
cor.test(pos_sel_dat_2_Tce_A$omega0_1, pos_sel_dat_2_Tce_A$GC, method = "spearman")
cor.test(pos_sel_dat_2_Tcm_A$omega0_1, pos_sel_dat_2_Tcm_A$GC, method = "spearman")
cor.test(pos_sel_dat_2_Tpa_A$omega0_1, pos_sel_dat_2_Tpa_A$GC, method = "spearman")
cor.test(pos_sel_dat_2_Tps_A$omega0_1, pos_sel_dat_2_Tps_A$GC, method = "spearman")




head(pos_sel_dat_2)

#### perm anova

m1_omega0_1 <- glm(pos_sel_dat_2$omega0_1 ~ pos_sel_dat_2$mean_FKPM + pos_sel_dat_2$GC * pos_sel_dat_2$chr_soft_class)
m1_omega0_1_out <- drop1(m1_omega0_1,~.,test="F") 
m1_omega0_1_out

pos_sel_dat_2_Tbi <- subset(pos_sel_dat_2, pos_sel_dat_2$branch_name == "Tbi")
pos_sel_dat_2_Tce <- subset(pos_sel_dat_2, pos_sel_dat_2$branch_name == "Tce")
pos_sel_dat_2_Tcm <- subset(pos_sel_dat_2, pos_sel_dat_2$branch_name == "Tcm")
pos_sel_dat_2_Tpa <- subset(pos_sel_dat_2, pos_sel_dat_2$branch_name == "Tpa")
pos_sel_dat_2_Tps <- subset(pos_sel_dat_2, pos_sel_dat_2$branch_name == "Tps")

m1_omega0_1_Tbi <- glm(pos_sel_dat_2_Tbi$omega0_1 ~ pos_sel_dat_2_Tbi$mean_FKPM + pos_sel_dat_2_Tbi$GC * pos_sel_dat_2_Tbi$chr_soft_class)
m1_omega0_1_Tbi_out <- drop1(m1_omega0_1_Tbi,~.,test="F") 
m1_omega0_1_Tbi_out

m1_omega0_1_Tce <- glm(pos_sel_dat_2_Tce$omega0_1 ~ pos_sel_dat_2_Tce$mean_FKPM + pos_sel_dat_2_Tce$GC * pos_sel_dat_2_Tce$chr_soft_class)
m1_omega0_1_Tce_out <- drop1(m1_omega0_1_Tce,~.,test="F") 
m1_omega0_1_Tce_out

m1_omega0_1_Tcm <- glm(pos_sel_dat_2_Tcm$omega0_1 ~ pos_sel_dat_2_Tcm$mean_FKPM + pos_sel_dat_2_Tcm$GC * pos_sel_dat_2_Tcm$chr_soft_class)
m1_omega0_1_Tcm_out <- drop1(m1_omega0_1_Tcm,~.,test="F") 
m1_omega0_1_Tcm_out

m1_omega0_1_Tpa <- glm(pos_sel_dat_2_Tpa$omega0_1 ~ pos_sel_dat_2_Tpa$mean_FKPM + pos_sel_dat_2_Tpa$GC * pos_sel_dat_2_Tpa$chr_soft_class)
m1_omega0_1_Tpa_out <- drop1(m1_omega0_1_Tpa,~.,test="F") 
m1_omega0_1_Tpa_out

m1_omega0_1_Tps <- glm(pos_sel_dat_2_Tps$omega0_1 ~ pos_sel_dat_2_Tps$mean_FKPM + pos_sel_dat_2_Tps$GC * pos_sel_dat_2_Tps$chr_soft_class)
m1_omega0_1_Tps_out <- drop1(m1_omega0_1_Tps,~.,test="F") 
m1_omega0_1_Tps_out

### overall test
## permutes all vals
perm_all_glm <- function(df){
  
 as.character(df$chr_soft_class)
  
  dd <- as.data.frame(cbind(
    sample(as.character(df$omega0_1)),
    as.character(df$chr_soft_class),
    as.character(df$GC),
    as.character(df$mean_FKPM)    
  ))
  
  colnames(dd) <- c("omega0_1 ", "chr_soft_class", "GC", "mean_FKPM")
  
  dd$omega0_1          <- as.numeric(as.character(dd$omega0_1))
  dd$chr_soft_class    <- as.factor(as.character(dd$chr_soft_class))
  dd$GC                <- as.numeric(as.character(dd$GC))
  dd$mean_FKPM         <- as.numeric(as.character(dd$mean_FKPM))

  model_1 <- glm(dd$omega0_1 ~ dd$mean_FKPM + dd$GC * dd$chr_soft_class)
  out <- drop1(model_1,~.,test="F") 

  F_mean_FKPM = out$F[2]  
  F_GC = out$F[3]
  F_chr_soft_class = out$F[4]
  F_int = out$F[5]
  
  P_mean_FKPM = out$P[2]  
  P_GC = out$P[3]
  P_chr_soft_class = out$P[4]
  P_int = out$P[5]
  
  out_vals <- c(F_mean_FKPM, F_GC, F_chr_soft_class, F_int, P_mean_FKPM, P_GC, P_chr_soft_class, P_int)
  return(out_vals)

}

N_perm = 20000

Tbi_all_perm_out <- data.frame()
Tce_all_perm_out <- data.frame()
Tcm_all_perm_out <- data.frame()
Tpa_all_perm_out <- data.frame()
Tps_all_perm_out <- data.frame()

for(i in seq(1:N_perm)){
  run_all_Tbi <- perm_all_glm(pos_sel_dat_2_Tbi)	
  Tbi_all_perm_out <- rbind(Tbi_all_perm_out, run_all_Tbi)	

  run_all_Tce <- perm_all_glm(pos_sel_dat_2_Tce)	
  Tce_all_perm_out <- rbind(Tce_all_perm_out, run_all_Tce)
  
  run_all_Tcm <- perm_all_glm(pos_sel_dat_2_Tcm)	
  Tcm_all_perm_out <- rbind(Tcm_all_perm_out, run_all_Tcm)
  
  run_all_Tpa <- perm_all_glm(pos_sel_dat_2_Tpa)	
  Tpa_all_perm_out <- rbind(Tpa_all_perm_out, run_all_Tpa)	
  
  run_all_Tps <- perm_all_glm(pos_sel_dat_2_Tps)	
  Tps_all_perm_out <- rbind(Tps_all_perm_out, run_all_Tps)	
  
}

colnames(Tbi_all_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")
colnames(Tce_all_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")
colnames(Tcm_all_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")
colnames(Tpa_all_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")
colnames(Tps_all_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")


#### output

write.csv(Tbi_all_perm_out, file=paste("Tbi_all_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(Tce_all_perm_out, file=paste("Tce_all_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(Tcm_all_perm_out, file=paste("Tcm_all_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(Tpa_all_perm_out, file=paste("Tpa_all_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(Tps_all_perm_out, file=paste("Tps_all_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)


#################################################################################################################
##### just flip chr type


perm_chr_glm <- function(df){
  as.character(df$chr_soft_class)
  
  dd <- as.data.frame(cbind(
    as.character(df$omega0_1),
    sample(as.character(df$chr_soft_class)),
    as.character(df$GC),
    as.character(df$mean_FKPM)    
  ))
  
  colnames(dd) <- c("omega0_1 ", "chr_soft_class", "GC", "mean_FKPM")
  
  dd$omega0_1          <- as.numeric(as.character(dd$omega0_1))
  dd$chr_soft_class    <- as.factor(as.character(dd$chr_soft_class))
  dd$GC                <- as.numeric(as.character(dd$GC))
  dd$mean_FKPM         <- as.numeric(as.character(dd$mean_FKPM))
  
  model_1 <- glm(dd$omega0_1 ~ dd$mean_FKPM + dd$GC * dd$chr_soft_class)
  out <- drop1(model_1,~.,test="F") 
  
  #print(head(dd))
  
  F_mean_FKPM = out$F[2]  
  F_GC = out$F[3]
  F_chr_soft_class = out$F[4]
  F_int = out$F[5]
  
  P_mean_FKPM = out$P[2]  
  P_GC = out$P[3]
  P_chr_soft_class = out$P[4]
  P_int = out$P[5]
  
  out_vals <- c(F_mean_FKPM, F_GC, F_chr_soft_class, F_int, P_mean_FKPM, P_GC, P_chr_soft_class, P_int)
  return(out_vals)
  
}


Tbi_chr_perm_out <- data.frame()
Tce_chr_perm_out <- data.frame()
Tcm_chr_perm_out <- data.frame()
Tpa_chr_perm_out <- data.frame()
Tps_chr_perm_out <- data.frame()

for(i in seq(1:N_perm)){
  run_chr_Tbi <- perm_chr_glm(pos_sel_dat_2_Tbi)	
  Tbi_chr_perm_out <- rbind(Tbi_chr_perm_out, run_chr_Tbi)	
  
  run_chr_Tce <- perm_chr_glm(pos_sel_dat_2_Tce)	
  Tce_chr_perm_out <- rbind(Tce_chr_perm_out, run_chr_Tce)
  
  run_chr_Tcm <- perm_chr_glm(pos_sel_dat_2_Tcm)	
  Tcm_chr_perm_out <- rbind(Tcm_chr_perm_out, run_chr_Tcm)
  
  run_chr_Tpa <- perm_chr_glm(pos_sel_dat_2_Tpa)	
  Tpa_chr_perm_out <- rbind(Tpa_chr_perm_out, run_chr_Tpa)	
  
  run_chr_Tps <- perm_chr_glm(pos_sel_dat_2_Tps)	
  Tps_chr_perm_out <- rbind(Tps_chr_perm_out, run_chr_Tps)	
  
}

colnames(Tbi_chr_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")
colnames(Tce_chr_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")
colnames(Tcm_chr_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")
colnames(Tpa_chr_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")
colnames(Tps_chr_perm_out) <- c("F_mean_FKPM", "F_GC", "F_chr_soft_class", "F_int", "P_mean_FKPM", "P_GC", "P_chr_soft_class", "P_int")

write.csv(Tbi_chr_perm_out, file=paste("Tbi_chr_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(Tce_chr_perm_out, file=paste("Tce_chr_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(Tcm_chr_perm_out, file=paste("Tcm_chr_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(Tpa_chr_perm_out, file=paste("Tpa_chr_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)
write.csv(Tps_chr_perm_out, file=paste("Tps_chr_perm_out_nperm=", N_perm, ".csv", sep = ""), row.names=FALSE)



##############################################################################################################
## get p vals

get_P_val <- function(perm_vector, orig_TS){
  v1 <- ifelse(perm_vector > orig_TS , perm_vector, NA)
  v2 <- v1[!is.na(v1)]
  N_over = length(v2)
  P = N_over / length(perm_vector)
 # print(N_over)
  return(P)	
}

Tbi_F_mean_FKPM      = m1_omega0_1_Tbi_out$F[2]  
Tbi_F_GC             = m1_omega0_1_Tbi_out$F[3]
Tbi_F_chr_soft_class = m1_omega0_1_Tbi_out$F[4]
Tbi_F_int            = m1_omega0_1_Tbi_out$F[5]

get_P_val(Tbi_all_perm_out$F_mean_FKPM, Tbi_F_mean_FKPM )
get_P_val(Tbi_all_perm_out$F_GC, Tbi_F_GC )
get_P_val(Tbi_all_perm_out$F_chr_soft_class, Tbi_F_chr_soft_class)
get_P_val(Tbi_all_perm_out$F_int, Tbi_F_int)


Tce_F_mean_FKPM      = m1_omega0_1_Tce_out$F[2]  
Tce_F_GC             = m1_omega0_1_Tce_out$F[3]
Tce_F_chr_soft_class = m1_omega0_1_Tce_out$F[4]
Tce_F_int            = m1_omega0_1_Tce_out$F[5]

get_P_val(Tce_all_perm_out$F_mean_FKPM, Tce_F_mean_FKPM )
get_P_val(Tce_all_perm_out$F_GC, Tce_F_GC )
get_P_val(Tce_all_perm_out$F_chr_soft_class, Tce_F_chr_soft_class)
get_P_val(Tce_all_perm_out$F_int, Tce_F_int)


Tcm_F_mean_FKPM      = m1_omega0_1_Tcm_out$F[2]  
Tcm_F_GC             = m1_omega0_1_Tcm_out$F[3]
Tcm_F_chr_soft_class = m1_omega0_1_Tcm_out$F[4]
Tcm_F_int            = m1_omega0_1_Tcm_out$F[5]

get_P_val(Tcm_all_perm_out$F_mean_FKPM, Tcm_F_mean_FKPM )
get_P_val(Tcm_all_perm_out$F_GC, Tcm_F_GC )
get_P_val(Tcm_all_perm_out$F_chr_soft_class, Tcm_F_chr_soft_class)
get_P_val(Tcm_all_perm_out$F_int, Tcm_F_int)

Tpa_F_mean_FKPM      = m1_omega0_1_Tpa_out$F[2]  
Tpa_F_GC             = m1_omega0_1_Tpa_out$F[3]
Tpa_F_chr_soft_class = m1_omega0_1_Tpa_out$F[4]
Tpa_F_int            = m1_omega0_1_Tpa_out$F[5]

get_P_val(Tpa_all_perm_out$F_mean_FKPM, Tpa_F_mean_FKPM )
get_P_val(Tpa_all_perm_out$F_GC, Tpa_F_GC )
get_P_val(Tpa_all_perm_out$F_chr_soft_class, Tpa_F_chr_soft_class)
get_P_val(Tpa_all_perm_out$F_int, Tpa_F_int)

Tps_F_mean_FKPM      = m1_omega0_1_Tps_out$F[2]  
Tps_F_GC             = m1_omega0_1_Tps_out$F[3]
Tps_F_chr_soft_class = m1_omega0_1_Tps_out$F[4]
Tps_F_int            = m1_omega0_1_Tps_out$F[5]

get_P_val(Tps_all_perm_out$F_mean_FKPM, Tps_F_mean_FKPM )
get_P_val(Tps_all_perm_out$F_GC, Tps_F_GC )
get_P_val(Tps_all_perm_out$F_chr_soft_class, Tps_F_chr_soft_class)
get_P_val(Tps_all_perm_out$F_int, Tps_F_int)



get_P_val(Tbi_all_perm_out$F_chr_soft_class, Tbi_F_chr_soft_class)
get_P_val(Tce_all_perm_out$F_chr_soft_class, Tce_F_chr_soft_class)
get_P_val(Tcm_all_perm_out$F_chr_soft_class, Tcm_F_chr_soft_class)
get_P_val(Tpa_all_perm_out$F_chr_soft_class, Tpa_F_chr_soft_class)
get_P_val(Tps_all_perm_out$F_chr_soft_class, Tps_F_chr_soft_class)

### very similar results
get_P_val(Tbi_chr_perm_out$F_chr_soft_class, Tbi_F_chr_soft_class)
get_P_val(Tce_chr_perm_out$F_chr_soft_class, Tce_F_chr_soft_class)
get_P_val(Tcm_chr_perm_out$F_chr_soft_class, Tcm_F_chr_soft_class)
get_P_val(Tpa_chr_perm_out$F_chr_soft_class, Tpa_F_chr_soft_class)
get_P_val(Tps_chr_perm_out$F_chr_soft_class, Tps_F_chr_soft_class)




#################################################################################
## LG

pos_sel_dat_2$sp_lg  <- paste(pos_sel_dat_2$branch_name, pos_sel_dat_2$lg, sep = "_") 
pos_sel_dat_2$sp_lg  <- ordered(pos_sel_dat_2$sp_lg, levels=c(
  "Tbi_lg1", "Tbi_lg2", "Tbi_lg3", "Tbi_lg4", "Tbi_lg5", "Tbi_lg6", "Tbi_lg7", "Tbi_lg8", "Tbi_lg9", "Tbi_lg10", "Tbi_lg11", "Tbi_lg12", "Tbi_lgX", 
  "Tce_lg1", "Tce_lg2", "Tce_lg3", "Tce_lg4", "Tce_lg5", "Tce_lg6", "Tce_lg7", "Tce_lg8", "Tce_lg9", "Tce_lg10", "Tce_lg11", "Tce_lg12", "Tce_lgX",        
  "Tcm_lg1", "Tcm_lg2", "Tcm_lg3", "Tcm_lg4", "Tcm_lg5", "Tcm_lg6", "Tcm_lg7", "Tcm_lg8", "Tcm_lg9", "Tcm_lg10", "Tcm_lg11", "Tcm_lg12", "Tcm_lgX", 
  "Tpa_lg1", "Tpa_lg2", "Tpa_lg3", "Tpa_lg4", "Tpa_lg5", "Tpa_lg6", "Tpa_lg7", "Tpa_lg8", "Tpa_lg9", "Tpa_lg10", "Tpa_lg11", "Tpa_lg12", "Tpa_lgX", 
  "Tps_lg1", "Tps_lg2", "Tps_lg3", "Tps_lg4", "Tps_lg5", "Tps_lg6", "Tps_lg7", "Tps_lg8", "Tps_lg9", "Tps_lg10", "Tps_lg11", "Tps_lg12", "Tps_lgX" ))

## get means
omega0_1_means_lg  <- summarySE(subset(pos_sel_dat_2, ! is.na(pos_sel_dat_2$lg)), measurevar="omega0_1", groupvars=c("sp_lg"))
omega0_1_means_lg$upper <- omega0_1_means_lg$omega0_1 +  omega0_1_means_lg$se
omega0_1_means_lg$lower <- omega0_1_means_lg$omega0_1 -  omega0_1_means_lg$se
omega0_1_means_lg$lg <- str_split_fixed(omega0_1_means_lg$sp_lg, "_", 2)[,2]

### exclude lgs with <X genes
min_genes = 50
omega0_1_means_lg <- subset(omega0_1_means_lg, omega0_1_means_lg$N > min_genes )

## plot
P_mean_omega0_1_lg <- ggplot() + 
  geom_errorbar(data=omega0_1_means_lg, mapping=aes(x=sp_lg, ymin=upper, ymax=lower,  color= lg), width=0.2, size=1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point(data=omega0_1_means_lg, mapping=aes(x = sp_lg, y=omega0_1,  fill=lg), size=4, shape=21) + 
  scale_color_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) +
  scale_fill_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) + ggtitle(paste("lg, mean +/- SE, min gene N = ", min_genes))

pdf(paste("P_mean_omega0_1_lg_mingeneN", min_genes ,".pdf", sep = ""), width = 12, height = 5)
P_mean_omega0_1_lg
dev.off()
getwd() ## where has my plot gone....



length(subset(pos_sel_dat_2, pos_sel_dat_2$sp_lg == "Tbi_lg12")[,1])
length(subset(pos_sel_dat_2, pos_sel_dat_2$sp_lg == "Tpa_lgX")[,1])

######################################################################################################################################################
### what about dN, dS, branch lens?

head(pos_sel_dat_2)


##################### Branch lengths first (this is all sites )
## get means
branch_length_means_soft  <- summarySE(pos_sel_dat_2, measurevar="branch_length", groupvars=c("sp_soft"))
branch_length_means_soft$upper <- branch_length_means_soft$branch_length +  branch_length_means_soft$se
branch_length_means_soft$lower <- branch_length_means_soft$branch_length -  branch_length_means_soft$se
branch_length_means_soft$chr_soft <- str_split_fixed(branch_length_means_soft$sp_soft, "_", 2)[,2]

## plot
P_mean_branch_length_soft <- ggplot() + 
  geom_errorbar(data=branch_length_means_soft, mapping=aes(x=sp_soft, ymin=upper, ymax=lower,  color= chr_soft), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=branch_length_means_soft, mapping=aes(x = sp_soft, y=branch_length,  fill=chr_soft), size=4, shape=21) + 
  scale_color_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE") +
  scale_fill_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE")

pdf(paste("P_mean_branch_length_soft" ,".pdf", sep = ""), width = 6, height = 8)
P_mean_branch_length_soft
dev.off()
getwd() ## where has my plot gone....


wilcox_branch_length_soft <- as.data.frame(cbind(
  c("Tbi", "Tce", "Tcm", "Tpa", "Tps"),
  c(
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_X")$branch_length, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_A")$branch_length, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_X")$branch_length, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_A")$branch_length, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_X")$branch_length, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_A")$branch_length, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_X")$branch_length, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_A")$branch_length, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_X")$branch_length, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_A")$branch_length, paired = FALSE)$p.value
  )))
colnames(wilcox_branch_length_soft) <- c("sp", "wilcox_p")
wilcox_branch_length_soft$wilcox_p <- as.numeric(wilcox_branch_length_soft$wilcox_p)
wilcox_branch_length_soft$wilcox_FDR <- p.adjust(wilcox_branch_length_soft$wilcox_p, method = "BH")
write.csv(wilcox_branch_length_soft, "wilcox_branch_length_soft.csv", row.names = FALSE)



############## dN and dS sep (excludes +ve sites as want to decompose my dN/dS to look which of these drives it)

get_dn_and_ds <- function(branch_len, dnds){
  options(scipen = 999)
  dnds = noquote(format(dnds , digits=10, nsmall=10))
  int_bit <- strsplit(as.character(dnds), "\\.")[[1]][1]
  dec_bit <- strsplit(as.character(dnds), "\\.")[[1]][2]
  denom = 10 ^ nchar(dec_bit)
  numer = as.numeric(int_bit) *  denom + as.numeric(dec_bit)
  # 
  # print(dnds)
  # print( dec_bit )
  # 
  # print( denom)
  # print( numer)  
  
  dN = (branch_len / (numer + denom)) * numer 
  dS = (branch_len / (numer + denom)) * denom  
  # print(dN)
  # print(dS) 
  output = list("dN" = dN, "dS" = dS)
  return(output)
}

dN_all <- c()
dS_all <- c()
for(i in 1:length(pos_sel_dat_2[,1])){
  branch_val <- pos_sel_dat_2$branch_length[i]
  dnds_val   <- pos_sel_dat_2$omega0_1[i]
  dn_val <- get_dn_and_ds (branch_val, dnds_val)$dN
  ds_val <- get_dn_and_ds (branch_val, dnds_val)$dS
  dN_all <- c(dN_all, dn_val)
  dS_all <- c(dS_all, ds_val)
}

pos_sel_dat_2$dN <- dN_all
pos_sel_dat_2$dS <- dS_all

# pos_sel_dat_2$dNdS_check <- pos_sel_dat_2$dN / pos_sel_dat_2$dS 
# plot (pos_sel_dat_2$dNdS_check, pos_sel_dat_2$omega0_1 )
head(pos_sel_dat_2, n = 30)

#### plot dN
## get means
dN_means_soft  <- summarySE(pos_sel_dat_2, measurevar="dN", groupvars=c("sp_soft"))
dN_means_soft$upper <- dN_means_soft$dN +  dN_means_soft$se
dN_means_soft$lower <- dN_means_soft$dN -  dN_means_soft$se
dN_means_soft$chr_soft <- str_split_fixed(dN_means_soft$sp_soft, "_", 2)[,2]

## plot
P_mean_dN_soft <- ggplot() + 
  geom_errorbar(data=dN_means_soft, mapping=aes(x=sp_soft, ymin=upper, ymax=lower,  color= chr_soft), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=dN_means_soft, mapping=aes(x = sp_soft, y=dN,  fill=chr_soft), size=4, shape=21) + 
  scale_color_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE") +
  scale_fill_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE")

pdf(paste("P_mean_dN_soft" ,".pdf", sep = ""), width = 6, height = 8)
P_mean_dN_soft
dev.off()
getwd() ## where has my plot gone....

wilcox_dN_soft <- as.data.frame(cbind(
  c("Tbi", "Tce", "Tcm", "Tpa", "Tps"),
  c(
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_X")$dN, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_A")$dN, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_X")$dN, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_A")$dN, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_X")$dN, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_A")$dN, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_X")$dN, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_A")$dN, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_X")$dN, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_A")$dN, paired = FALSE)$p.value
  )))
colnames(wilcox_dN_soft) <- c("sp", "wilcox_p")
wilcox_dN_soft$wilcox_p <- as.numeric(wilcox_dN_soft$wilcox_p)
wilcox_dN_soft$wilcox_FDR <- p.adjust(wilcox_dN_soft$wilcox_p, method = "BH")
write.csv(wilcox_dN_soft, "wilcox_dN_soft.csv", row.names = FALSE)

#################################################################################
## LG

## get means
branch_length_means_lg  <- summarySE(subset(pos_sel_dat_2, ! is.na(pos_sel_dat_2$lg)), measurevar="branch_length", groupvars=c("sp_lg"))
branch_length_means_lg$upper <- branch_length_means_lg$branch_length +  branch_length_means_lg$se
branch_length_means_lg$lower <- branch_length_means_lg$branch_length -  branch_length_means_lg$se
branch_length_means_lg$lg <- str_split_fixed(branch_length_means_lg$sp_lg, "_", 2)[,2]

### exclude lgs with <X genes
min_genes = 50
branch_length_means_lg <- subset(branch_length_means_lg, branch_length_means_lg$N > min_genes )

## plot
P_mean_branch_length_lg <- ggplot() + 
  geom_errorbar(data=branch_length_means_lg, mapping=aes(x=sp_lg, ymin=upper, ymax=lower,  color= lg), width=0.2, size=1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point(data=branch_length_means_lg, mapping=aes(x = sp_lg, y=branch_length,  fill=lg), size=4, shape=21) + 
  scale_color_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) +
  scale_fill_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) + ggtitle(paste("lg, mean +/- SE, min gene N = ", min_genes))

pdf(paste("P_mean_branch_length_lg_mingeneN", min_genes ,".pdf", sep = ""), width = 12, height = 5)
P_mean_branch_length_lg
dev.off()
getwd() ## where has my plot gone....











#### plot ds
## get means
dS_means_soft  <- summarySE(pos_sel_dat_2, measurevar="dS", groupvars=c("sp_soft") , na.rm = T)
dS_means_soft$upper <- dS_means_soft$dS +  dS_means_soft$se
dS_means_soft$lower <- dS_means_soft$dS -  dS_means_soft$se
dS_means_soft$chr_soft <- str_split_fixed(dS_means_soft$sp_soft, "_", 2)[,2]

## plot
P_mean_dS_soft <- ggplot() + 
  geom_errorbar(data=dS_means_soft, mapping=aes(x=sp_soft, ymin=upper, ymax=lower,  color= chr_soft), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=dS_means_soft, mapping=aes(x = sp_soft, y=dS,  fill=chr_soft), size=4, shape=21) + 
  scale_color_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE") +
  scale_fill_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE")

pdf(paste("P_mean_dS_soft" ,".pdf", sep = ""), width = 6, height = 8)
P_mean_dS_soft
dev.off()
getwd() ## where has my plot gone....

wilcox_dS_soft <- as.data.frame(cbind(
  c("Tbi", "Tce", "Tcm", "Tpa", "Tps"),
  c(
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_X")$dS, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tbi_A")$dS, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_X")$dS, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tce_A")$dS, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_X")$dS, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tcm_A")$dS, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_X")$dS, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tpa_A")$dS, paired = FALSE)$p.value,
    wilcox.test(subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_X")$dS, subset(pos_sel_dat_2, pos_sel_dat_2$sp_soft == "Tps_A")$dS, paired = FALSE)$p.value
  )))
colnames(wilcox_dS_soft) <- c("sp", "wilcox_p")
wilcox_dS_soft$wilcox_p <- as.numeric(wilcox_dS_soft$wilcox_p)
wilcox_dS_soft$wilcox_FDR <- p.adjust(wilcox_dS_soft$wilcox_p, method = "BH")
write.csv(wilcox_dS_soft, "wilcox_dS_soft.csv", row.names = FALSE)



#### ratio of dNA to dNX

subset(dN_means_soft, dN_means_soft$sp_soft == "Tbi_A")$dN / subset(dN_means_soft, dN_means_soft$sp_soft == "Tbi_X")$dN
subset(dN_means_soft, dN_means_soft$sp_soft == "Tce_A")$dN / subset(dN_means_soft, dN_means_soft$sp_soft == "Tce_X")$dN
subset(dN_means_soft, dN_means_soft$sp_soft == "Tcm_A")$dN / subset(dN_means_soft, dN_means_soft$sp_soft == "Tcm_X")$dN
subset(dN_means_soft, dN_means_soft$sp_soft == "Tpa_A")$dN / subset(dN_means_soft, dN_means_soft$sp_soft == "Tpa_X")$dN
subset(dN_means_soft, dN_means_soft$sp_soft == "Tps_A")$dN / subset(dN_means_soft, dN_means_soft$sp_soft == "Tps_X")$dN

#### ratio of dSA to dSX

subset(dS_means_soft, dS_means_soft$sp_soft == "Tbi_A")$dS / subset(dS_means_soft, dS_means_soft$sp_soft == "Tbi_X")$dS
subset(dS_means_soft, dS_means_soft$sp_soft == "Tce_A")$dS / subset(dS_means_soft, dS_means_soft$sp_soft == "Tce_X")$dS
subset(dS_means_soft, dS_means_soft$sp_soft == "Tcm_A")$dS / subset(dS_means_soft, dS_means_soft$sp_soft == "Tcm_X")$dS
subset(dS_means_soft, dS_means_soft$sp_soft == "Tpa_A")$dS / subset(dS_means_soft, dS_means_soft$sp_soft == "Tpa_X")$dS
subset(dS_means_soft, dS_means_soft$sp_soft == "Tps_A")$dS / subset(dS_means_soft, dS_means_soft$sp_soft == "Tps_X")$dS





####################################################################################################
#### GC

## get means
GC_means_soft  <- summarySE(pos_sel_dat_2, measurevar="GC", groupvars=c("sp_soft"))
GC_means_soft$upper <- GC_means_soft$GC +  GC_means_soft$se
GC_means_soft$lower <- GC_means_soft$GC -  GC_means_soft$se
GC_means_soft$chr_soft <- str_split_fixed(GC_means_soft$sp_soft, "_", 2)[,2]

## plot
P_mean_GC_soft <- ggplot() + 
  geom_errorbar(data=GC_means_soft, mapping=aes(x=sp_soft, ymin=upper, ymax=lower,  color= chr_soft), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=GC_means_soft, mapping=aes(x = sp_soft, y=GC,  fill=chr_soft), size=4, shape=21) + 
  scale_color_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE") +
  scale_fill_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE")

pdf(paste("P_mean_GC_soft" ,".pdf", sep = ""), width = 6, height = 8)
P_mean_GC_soft
dev.off()
getwd() ## where has my plot gone....



#################################################################################
## LG

## get means
GC_means_lg  <- summarySE(subset(pos_sel_dat_2, ! is.na(pos_sel_dat_2$lg)), measurevar="GC", groupvars=c("sp_lg"))
GC_means_lg$upper <- GC_means_lg$GC +  GC_means_lg$se
GC_means_lg$lower <- GC_means_lg$GC -  GC_means_lg$se
GC_means_lg$lg <- str_split_fixed(GC_means_lg$sp_lg, "_", 2)[,2]

### exclude lgs with <X genes
min_genes = 50
GC_means_lg <- subset(GC_means_lg, GC_means_lg$N > min_genes )

## plot
P_mean_GC_lg <- ggplot() + 
  geom_errorbar(data=GC_means_lg, mapping=aes(x=sp_lg, ymin=upper, ymax=lower,  color= lg), width=0.2, size=1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point(data=GC_means_lg, mapping=aes(x = sp_lg, y=GC,  fill=lg), size=4, shape=21) + 
  scale_color_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) +
  scale_fill_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) + ggtitle(paste("lg, mean +/- SE, min gene N = ", min_genes))

pdf(paste("P_mean_GC_lg_mingeneN", min_genes ,".pdf", sep = ""), width = 12, height = 5)
P_mean_GC_lg
dev.off()
getwd() ## where has my plot gone....



#############################################
### ENC

## get means
ENC_means_soft  <- summarySE(pos_sel_dat_2, measurevar="ENC", groupvars=c("sp_soft"))
ENC_means_soft$upper <- ENC_means_soft$ENC +  ENC_means_soft$se
ENC_means_soft$lower <- ENC_means_soft$ENC -  ENC_means_soft$se
ENC_means_soft$chr_soft <- str_split_fixed(ENC_means_soft$sp_soft, "_", 2)[,2]

## plot
P_mean_ENC_soft <- ggplot() + 
  geom_errorbar(data=ENC_means_soft, mapping=aes(x=sp_soft, ymin=upper, ymax=lower,  color= chr_soft), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=ENC_means_soft, mapping=aes(x = sp_soft, y=ENC,  fill=chr_soft), size=4, shape=21) + 
  scale_color_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE") +
  scale_fill_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft chr, mean +/- SE")

pdf(paste("P_mean_ENC_soft" ,".pdf", sep = ""), width = 6, height = 8)
P_mean_ENC_soft
dev.off()
getwd() ## where has my plot gone....



## LG

## get means
ENC_means_lg  <- summarySE(subset(pos_sel_dat_2, ! is.na(pos_sel_dat_2$lg)), measurevar="ENC", groupvars=c("sp_lg"))
ENC_means_lg$upper <- ENC_means_lg$ENC +  ENC_means_lg$se
ENC_means_lg$lower <- ENC_means_lg$ENC -  ENC_means_lg$se
ENC_means_lg$lg <- str_split_fixed(ENC_means_lg$sp_lg, "_", 2)[,2]

### exclude lgs with <X genes
min_genes = 50
ENC_means_lg <- subset(ENC_means_lg, ENC_means_lg$N > min_genes )

## plot
P_mean_ENC_lg <- ggplot() + 
  geom_errorbar(data=ENC_means_lg, mapping=aes(x=sp_lg, ymin=upper, ymax=lower,  color= lg), width=0.2, size=1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point(data=ENC_means_lg, mapping=aes(x = sp_lg, y=ENC,  fill=lg), size=4, shape=21) + 
  scale_color_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) +
  scale_fill_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) + ggtitle(paste("lg, mean +/- SE, min gene N = ", min_genes))

pdf(paste("P_mean_ENC_lg_mingeneN", min_genes ,".pdf", sep = ""), width = 12, height = 5)
P_mean_ENC_lg
dev.off()
getwd() ## where has my plot gone....





head(pos_sel_dat_2)
library(car)
m1 <- glm(pos_sel_dat_2$ENC ~ pos_sel_dat_2$GC + pos_sel_dat_2$RT_SF_FKPM + pos_sel_dat_2$chr_soft_class + pos_sel_dat_2$branch_name)
summary(m1)


m1 <- glm(pos_sel_dat_2$omega0_1 ~  pos_sel_dat_2$GC + pos_sel_dat_2$RT_SF_FKPM + pos_sel_dat_2$chr_soft_class + pos_sel_dat_2$branch_name)
summary(m1)

Anova(m1, type = 3)

m1 <- glm(pos_sel_dat_2$omega0_1 ~  pos_sel_dat_2$GC +pos_sel_dat_2$ENC + pos_sel_dat_2$RT_SF_FKPM + pos_sel_dat_2$chr_soft_class + pos_sel_dat_2$branch_name)
summary(m1)
Anova(m1, type = 3)

m1 <- glm(pos_sel_dat_2$dN~  pos_sel_dat_2$GC +pos_sel_dat_2$ENC + pos_sel_dat_2$RT_SF_FKPM + pos_sel_dat_2$chr_soft_class + pos_sel_dat_2$branch_name)
summary(m1)
Anova(m1, type = 3)

m1 <- glm(pos_sel_dat_2$dS~  pos_sel_dat_2$GC +pos_sel_dat_2$ENC + pos_sel_dat_2$RT_SF_FKPM + pos_sel_dat_2$chr_soft_class + pos_sel_dat_2$branch_name)
summary(m1)
Anova(m1, type = 3)

### par cor

### data


 

library(ppcor)
pcor.test(X1, X2, X3, method = c("spearman")) ## partial corr of prop_off_sire and inv_rel_order accounting for spd_size  (corr of X1 and X2 accounting for X3)

### this only works for 3 vars - what if I want more?
## do this: 

q1 = as.data.frame(cbind(X1,X2,X3,X4,X5))
colnames(q1) <- c("X1", "X2", "X3", "X4", "X5")
pcor(q1,method = c("spearman"))

## can plot with added variable plots

library(car)
o8 <- glm(X1 ~ X2 + X3)
avPlots(o8)




################################################################################################################
### SB


## SB +ve vals = malebiased

pos_sel_dat_2$RT_SB   <- ifelse(pos_sel_dat_2$RT_FDR < 0.05 & pos_sel_dat_2$RT_logFC > 0 ,   "MB", 
                         ifelse(pos_sel_dat_2$RT_FDR < 0.05 & pos_sel_dat_2$RT_logFC < 0 ,   "FB",  "UB"))
pos_sel_dat_2$HD_SB   <- ifelse(pos_sel_dat_2$HD_FDR < 0.05 & pos_sel_dat_2$HD_logFC > 0 ,   "MB", 
                         ifelse(pos_sel_dat_2$HD_FDR < 0.05 & pos_sel_dat_2$HD_logFC < 0 ,   "FB",  "UB"))
pos_sel_dat_2$LG_SB   <- ifelse(pos_sel_dat_2$LG_FDR < 0.05 & pos_sel_dat_2$LG_logFC > 0 ,   "MB", 
                         ifelse(pos_sel_dat_2$LG_FDR < 0.05 & pos_sel_dat_2$LG_logFC < 0 ,   "FB",  "UB"))
head(pos_sel_dat_2, n = 100)

## soft HD
pos_sel_dat_2$sp_soft_HDSB  <- paste(pos_sel_dat_2$branch_name, pos_sel_dat_2$chr_soft_class, pos_sel_dat_2$HD_SB , sep = "_") 

## get means
omega0_1_means_soft_HDSB  <- summarySE(pos_sel_dat_2, measurevar="omega0_1", groupvars=c("sp_soft_HDSB"))
omega0_1_means_soft_HDSB$upper <- omega0_1_means_soft_HDSB$omega0_1 +  omega0_1_means_soft_HDSB$se
omega0_1_means_soft_HDSB$lower <- omega0_1_means_soft_HDSB$omega0_1 -  omega0_1_means_soft_HDSB$se
omega0_1_means_soft_HDSB$chr_soft <- str_split_fixed(omega0_1_means_soft_HDSB$sp_soft, "_", 2)[,2]
omega0_1_means_soft_HDSB$SB       <- str_split_fixed(omega0_1_means_soft_HDSB$sp_soft, "_", 3)[,3]
omega0_1_means_soft_HDSB_2        <- subset(omega0_1_means_soft_HDSB, omega0_1_means_soft_HDSB$SB != "NA") 


omega0_1_means_soft_HDSB$chr_soft <- ordered(omega0_1_means_soft_HDSB$chr_soft, levels = c("A_UB", "X_UB", "A_FB", "X_FB", "A_MB", "X_MB"))

omega0_1_means_soft_HDSB_2$sp_soft_HDSB <- ordered(omega0_1_means_soft_HDSB_2$sp_soft_HDSB, levels = c(
  "Tbi_A_UB", "Tbi_X_UB", "Tbi_A_FB", "Tbi_X_FB", "Tbi_A_MB", "Tbi_X_MB",
  "Tce_A_UB", "Tce_X_UB", "Tce_A_FB", "Tce_X_FB", "Tce_A_MB", "Tce_X_MB",
  "Tcm_A_UB", "Tcm_X_UB", "Tcm_A_FB", "Tcm_X_FB", "Tcm_A_MB", "Tcm_X_MB",
  "Tpa_A_UB", "Tpa_X_UB", "Tpa_A_FB", "Tpa_X_FB", "Tpa_A_MB", "Tpa_X_MB",
  "Tps_A_UB", "Tps_X_UB", "Tps_A_FB", "Tps_X_FB", "Tps_A_MB", "Tps_X_MB"))

## plot
P_mean_omega0_1_soft_HDSB <- ggplot() + 
  geom_errorbar(data=omega0_1_means_soft_HDSB_2, mapping=aes(x=sp_soft_HDSB, ymin=upper, ymax=lower,  color= chr_soft), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=omega0_1_means_soft_HDSB_2  , mapping=aes(x = sp_soft_HDSB, y=omega0_1,  fill=chr_soft), size=4, shape=21) + 
  scale_color_manual(values=c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + 
  scale_fill_manual(values=c("pink", "lightblue", "darkgrey","red", "darkblue", "black"))+ ggtitle("HD soft chr, mean +/- SE")
P_mean_omega0_1_soft_HDSB <- P_mean_omega0_1_soft_HDSB + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## soft LG
pos_sel_dat_2$sp_soft_LGSB  <- paste(pos_sel_dat_2$branch_name, pos_sel_dat_2$chr_soft_class, pos_sel_dat_2$LG_SB , sep = "_") 

## get means
omega0_1_means_soft_LGSB  <- summarySE(pos_sel_dat_2, measurevar="omega0_1", groupvars=c("sp_soft_LGSB"))
omega0_1_means_soft_LGSB$upper <- omega0_1_means_soft_LGSB$omega0_1 +  omega0_1_means_soft_LGSB$se
omega0_1_means_soft_LGSB$lower <- omega0_1_means_soft_LGSB$omega0_1 -  omega0_1_means_soft_LGSB$se
omega0_1_means_soft_LGSB$chr_soft <- str_split_fixed(omega0_1_means_soft_LGSB$sp_soft, "_", 2)[,2]
omega0_1_means_soft_LGSB$SB       <- str_split_fixed(omega0_1_means_soft_LGSB$sp_soft, "_", 3)[,3]
omega0_1_means_soft_LGSB_2        <- subset(omega0_1_means_soft_LGSB, omega0_1_means_soft_LGSB$SB != "NA") 


omega0_1_means_soft_LGSB$chr_soft <- ordered(omega0_1_means_soft_LGSB$chr_soft, levels = c("A_UB", "X_UB", "A_FB", "X_FB", "A_MB", "X_MB"))

omega0_1_means_soft_LGSB_2$sp_soft_LGSB <- ordered(omega0_1_means_soft_LGSB_2$sp_soft_LGSB, levels = c(
  "Tbi_A_UB", "Tbi_X_UB", "Tbi_A_FB", "Tbi_X_FB", "Tbi_A_MB", "Tbi_X_MB",
  "Tce_A_UB", "Tce_X_UB", "Tce_A_FB", "Tce_X_FB", "Tce_A_MB", "Tce_X_MB",
  "Tcm_A_UB", "Tcm_X_UB", "Tcm_A_FB", "Tcm_X_FB", "Tcm_A_MB", "Tcm_X_MB",
  "Tpa_A_UB", "Tpa_X_UB", "Tpa_A_FB", "Tpa_X_FB", "Tpa_A_MB", "Tpa_X_MB",
  "Tps_A_UB", "Tps_X_UB", "Tps_A_FB", "Tps_X_FB", "Tps_A_MB", "Tps_X_MB"))

## plot
P_mean_omega0_1_soft_LGSB <- ggplot() + 
  geom_errorbar(data=omega0_1_means_soft_LGSB_2, mapping=aes(x=sp_soft_LGSB, ymin=upper, ymax=lower,  color= chr_soft), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=omega0_1_means_soft_LGSB_2  , mapping=aes(x = sp_soft_LGSB, y=omega0_1,  fill=chr_soft), size=4, shape=21) + 
  scale_color_manual(values=c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + 
  scale_fill_manual(values=c("pink", "lightblue", "darkgrey","red", "darkblue", "black"))+ ggtitle("LG soft chr, mean +/- SE")
P_mean_omega0_1_soft_LGSB <- P_mean_omega0_1_soft_LGSB + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## soft
pos_sel_dat_2$sp_soft_RTSB  <- paste(pos_sel_dat_2$branch_name, pos_sel_dat_2$chr_soft_class, pos_sel_dat_2$RT_SB , sep = "_") 

## get means
omega0_1_means_soft_RTSB  <- summarySE(pos_sel_dat_2, measurevar="omega0_1", groupvars=c("sp_soft_RTSB"))
omega0_1_means_soft_RTSB$upper <- omega0_1_means_soft_RTSB$omega0_1 +  omega0_1_means_soft_RTSB$se
omega0_1_means_soft_RTSB$lower <- omega0_1_means_soft_RTSB$omega0_1 -  omega0_1_means_soft_RTSB$se
omega0_1_means_soft_RTSB$chr_soft <- str_split_fixed(omega0_1_means_soft_RTSB$sp_soft, "_", 2)[,2]
omega0_1_means_soft_RTSB$SB       <- str_split_fixed(omega0_1_means_soft_RTSB$sp_soft, "_", 3)[,3]
omega0_1_means_soft_RTSB_2        <- subset(omega0_1_means_soft_RTSB, omega0_1_means_soft_RTSB$SB != "NA") 


omega0_1_means_soft_RTSB$chr_soft <- ordered(omega0_1_means_soft_RTSB$chr_soft, levels = c("A_UB", "X_UB", "A_FB", "X_FB", "A_MB", "X_MB"))

omega0_1_means_soft_RTSB_2$sp_soft_RTSB <- ordered(omega0_1_means_soft_RTSB_2$sp_soft_RTSB, levels = c(
  "Tbi_A_UB", "Tbi_X_UB", "Tbi_A_FB", "Tbi_X_FB", "Tbi_A_MB", "Tbi_X_MB",
  "Tce_A_UB", "Tce_X_UB", "Tce_A_FB", "Tce_X_FB", "Tce_A_MB", "Tce_X_MB",
  "Tcm_A_UB", "Tcm_X_UB", "Tcm_A_FB", "Tcm_X_FB", "Tcm_A_MB", "Tcm_X_MB",
  "Tpa_A_UB", "Tpa_X_UB", "Tpa_A_FB", "Tpa_X_FB", "Tpa_A_MB", "Tpa_X_MB",
  "Tps_A_UB", "Tps_X_UB", "Tps_A_FB", "Tps_X_FB", "Tps_A_MB", "Tps_X_MB"))

## plot
P_mean_omega0_1_soft_RTSB <- ggplot() + 
  geom_errorbar(data=omega0_1_means_soft_RTSB_2, mapping=aes(x=sp_soft_RTSB, ymin=upper, ymax=lower,  color= chr_soft), width=0.2, size=1) + 
  theme_bw() +
  geom_point(data=omega0_1_means_soft_RTSB_2  , mapping=aes(x = sp_soft_RTSB, y=omega0_1,  fill=chr_soft), size=4, shape=21) + 
  scale_color_manual(values=c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + 
  scale_fill_manual(values=c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + ggtitle("RT soft chr, mean +/- SE")
P_mean_omega0_1_soft_RTSB <- P_mean_omega0_1_soft_RTSB + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf(paste("P_mean_omega0_1_soft_HDSB" ,".pdf", sep = ""), width = 12, height = 8)
P_mean_omega0_1_soft_HDSB 
dev.off()
getwd() ## where has my plot gone....

pdf(paste("P_mean_omega0_1_soft_LGSB" ,".pdf", sep = ""), width = 12, height = 8)
P_mean_omega0_1_soft_LGSB 
dev.off()
getwd() ## where has my plot gone....

pdf(paste("P_mean_omega0_1_soft_RTSB" ,".pdf", sep = ""), width = 12, height = 8)
P_mean_omega0_1_soft_RTSB 
dev.off()
getwd() ## where has my plot gone....



###### +ve sel by SB

head(pos_sel_dat_2)

possel_X_SB <- function(df_all, pos_q_threh){
  
  df_all$soft_SB_HD <- paste(df_all$branch_name, df_all$chr_soft_class, df_all$HD_SB, "HD", sep = "_")
  df_all$soft_SB_LG <- paste(df_all$branch_name, df_all$chr_soft_class, df_all$LG_SB, "LG", sep = "_")
  df_all$soft_SB_RT <- paste(df_all$branch_name, df_all$chr_soft_class, df_all$RT_SB, "RT", sep = "_")

  ## N total genes in each cat
  N_HD <- as.data.frame(table(df_all$soft_SB_HD))
  N_LG <- as.data.frame(table(df_all$soft_SB_LG))
  N_RT <- as.data.frame(table(df_all$soft_SB_RT))
  
  N_HDLGRT          <- rbind(N_HD,  N_LG,  N_RT)
  N_HDLGRT$sp       <- str_split_fixed(as.character(N_HDLGRT$Var1), "_", 4)[,1]
  N_HDLGRT$chr_soft <- str_split_fixed(as.character(N_HDLGRT$Var1), "_", 4)[,2]
  N_HDLGRT$SB       <- str_split_fixed(as.character(N_HDLGRT$Var1), "_", 4)[,3]
  N_HDLGRT$tiss     <- str_split_fixed(as.character(N_HDLGRT$Var1), "_", 4)[,4]
  
  N_HDLGRT_2 <- subset(N_HDLGRT, N_HDLGRT$SB != "NA")
  
  ## N pos sel genes in each cat
  
  df_pos <- subset(df_all, df_all$qvalue < pos_q_threh)
  
  N_pos_HD <- as.data.frame(table(df_pos$soft_SB_HD))
  N_pos_LG <- as.data.frame(table(df_pos$soft_SB_LG))
  N_pos_RT <- as.data.frame(table(df_pos$soft_SB_RT))
  
  N_pos_HDLGRT          <- rbind(N_pos_HD,  N_pos_LG,  N_pos_RT)
  N_pos_HDLGRT$sp       <- str_split_fixed(as.character(N_pos_HDLGRT$Var1), "_", 4)[,1]
  N_pos_HDLGRT$chr_soft <- str_split_fixed(as.character(N_pos_HDLGRT$Var1), "_", 4)[,2]
  N_pos_HDLGRT$SB       <- str_split_fixed(as.character(N_pos_HDLGRT$Var1), "_", 4)[,3]
  N_pos_HDLGRT$tiss     <- str_split_fixed(as.character(N_pos_HDLGRT$Var1), "_", 4)[,4]
  
  N_pos_HDLGRT_2 <- subset(N_pos_HDLGRT, N_pos_HDLGRT$SB != "NA")
  
  out_df <- merge(N_HDLGRT_2, N_pos_HDLGRT_2, by.x = 1, by.y = 1, all.x = TRUE)
  out_df$Freq_B_0 <- ifelse(is.na(out_df$Freq.y), 0, out_df$Freq.y )
  
  out_df <- out_df[,c(1:6,12)] 
  
  colnames(out_df) <- c("Group", "N_tot", "sp", "chr_soft", "SB", "tiss", "N_pos")
  out_df$prop_pos <- out_df$N_pos / out_df$N_tot
  out_df$chr_SB  <- paste(out_df$chr_soft, out_df$SB, sep = "_")
  out_df$sp_chr_SB <- paste(out_df$sp, out_df$chr_soft, out_df$SB, sep = "_")  
  return(out_df  )
  
}


q_thesh = 0.05

N_pos_sel_SB <- possel_X_SB(pos_sel_dat_2, q_thesh) 
N_pos_sel_SB$sp_chr_SB <- ordered(N_pos_sel_SB$sp_chr_SB, levels = c(
  "Tbi_A_UB", "Tbi_X_UB", "Tbi_A_FB", "Tbi_X_FB", "Tbi_A_MB", "Tbi_X_MB",
  "Tce_A_UB", "Tce_X_UB", "Tce_A_FB", "Tce_X_FB", "Tce_A_MB", "Tce_X_MB",
  "Tcm_A_UB", "Tcm_X_UB", "Tcm_A_FB", "Tcm_X_FB", "Tcm_A_MB", "Tcm_X_MB",
  "Tpa_A_UB", "Tpa_X_UB", "Tpa_A_FB", "Tpa_X_FB", "Tpa_A_MB", "Tpa_X_MB",
  "Tps_A_UB", "Tps_X_UB", "Tps_A_FB", "Tps_X_FB", "Tps_A_MB", "Tps_X_MB"))

N_pos_sel_SB_HD <- subset(N_pos_sel_SB, N_pos_sel_SB$tiss == "HD")

prop_pos_HD <- ggplot(N_pos_sel_SB_HD, aes(x = factor(sp_chr_SB), y = prop_pos, fill = chr_SB)) + 
  geom_col(width = 0.5, colour="black", position=position_dodge(width=0.6)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + 
  xlab ("Group") + 
  ylab ("Prop positive selected genes")  + 
  ggtitle(paste("soft, HD, qval thresh = ",q_thesh))
prop_pos_HD <- prop_pos_HD  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

N_pos_HD <- ggplot(N_pos_sel_SB_HD, aes(x = factor(sp_chr_SB), y = N_pos, fill = chr_SB)) + 
  geom_col(width = 0.5, colour="black", position=position_dodge(width=0.6)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + 
  xlab ("Group") + 
  ylab ("N positive selected genes")  + 
  ggtitle(paste("soft, HD, qval thresh = ",q_thesh))
N_pos_HD  <- N_pos_HD + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


N_pos_sel_SB_LG <- subset(N_pos_sel_SB, N_pos_sel_SB$tiss == "LG")

prop_pos_LG <- ggplot(N_pos_sel_SB_LG, aes(x = factor(sp_chr_SB), y = prop_pos, fill = chr_SB)) + 
  geom_col(width = 0.5, colour="black", position=position_dodge(width=0.6)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + 
  xlab ("Group") + 
  ylab ("Prop positive selected genes")  + 
  ggtitle(paste("soft, LG, qval thresh = ",q_thesh))
prop_pos_LG <- prop_pos_LG  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

N_pos_LG <- ggplot(N_pos_sel_SB_LG, aes(x = factor(sp_chr_SB), y = N_pos, fill = chr_SB)) + 
  geom_col(width = 0.5, colour="black", position=position_dodge(width=0.6)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + 
  xlab ("Group") + 
  ylab ("N positive selected genes")  + 
  ggtitle(paste("soft, LG, qval thresh = ",q_thesh))
N_pos_LG  <- N_pos_LG + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


N_pos_sel_SB_RT <- subset(N_pos_sel_SB, N_pos_sel_SB$tiss == "RT")

prop_pos_RT <- ggplot(N_pos_sel_SB_RT, aes(x = factor(sp_chr_SB), y = prop_pos, fill = chr_SB)) + 
  geom_col(width = 0.5, colour="black", position=position_dodge(width=0.6)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + 
  xlab ("Group") + 
  ylab ("Prop positive selected genes")  + 
  ggtitle(paste("soft, RT, qval thresh = ",q_thesh))
prop_pos_RT <- prop_pos_RT  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

N_pos_RT <- ggplot(N_pos_sel_SB_RT, aes(x = factor(sp_chr_SB), y = N_pos, fill = chr_SB)) + 
  geom_col(width = 0.5, colour="black", position=position_dodge(width=0.6)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values = c("pink", "lightblue", "darkgrey","red", "darkblue", "black")) + 
  xlab ("Group") + 
  ylab ("N positive selected genes")  + 
  ggtitle(paste("soft, RT, qval thresh = ",q_thesh))
N_pos_RT  <- N_pos_RT + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



pdf(paste("P_mean_omega0_1_soft_HDLGRTSB" ,".pdf", sep = ""), width = 8, height = 16)
plot_grid(P_mean_omega0_1_soft_HDSB, P_mean_omega0_1_soft_LGSB ,P_mean_omega0_1_soft_RTSB, ncol = 1 )
dev.off()
getwd() ## where has my plot gone....


pdf(paste("prop_pos_HDLGRT" ,".pdf", sep = ""), width = 8, height = 16)
plot_grid(prop_pos_HD, prop_pos_LG ,prop_pos_RT, ncol = 1 )
dev.off()
getwd() ## where has my plot gone....

pdf(paste("N_pos_HDLGRT" ,".pdf", sep = ""), width = 8, height = 16)
plot_grid(N_pos_HD, N_pos_LG ,N_pos_RT, ncol = 1 )
dev.off()
getwd() ## where has my plot gone....






















print (sessionInfo())
# 
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12 hash_2.2.6.1    edgeR_3.32.1    limma_3.46.0    cowplot_1.1.1   Rtsne_0.15      plyr_1.8.6      modeest_2.4.0   stringr_1.4.0   ggplot2_3.3.5  
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.6          RColorBrewer_1.1-2  pillar_1.4.7        compiler_4.0.3      tools_4.0.3         timeSeries_3062.100 digest_0.6.27       rpart_4.1-15       
# [9] lattice_0.20-41     lifecycle_0.2.0     tibble_3.0.6        gtable_0.3.0        stable_1.1.4        clue_0.3-58         pkgconfig_2.0.3     rlang_0.4.10       
# [17] rmutil_1.1.5        statip_0.2.3        withr_2.4.1         dplyr_1.0.3         cluster_2.1.0       generics_0.1.0      vctrs_0.3.6         locfit_1.5-9.4     
# [25] grid_4.0.3          tidyselect_1.1.0    glue_1.4.2          R6_2.5.0            spatial_7.3-13      farver_2.0.3        purrr_0.3.4         magrittr_2.0.1     
# [33] fBasics_3042.89.1   scales_1.1.1        ellipsis_0.3.1      stabledist_0.7-1    timeDate_3043.102   colorspace_2.0-0    labeling_0.4.2      stringi_1.5.3      
# [41] munsell_0.5.0       crayon_1.4.0  

writeLines(capture.output(sessionInfo()), "sel_analy_sess_info.txt")

