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


## from Expression_analyses.R
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

Orth_dict <- hash()
hard_chr_dict <- hash()
soft_chr_dict <- hash()
all_orths = c()
gene_to_orth_dict <- hash()

for(d in all_raw_dfs){
	sp <- strsplit(d, "_")[[1]][2]
	test_df <- eval(parse(text=paste(d,sep='')))

	for(i in seq(1:length(test_df[,1]))){
		gene_n <- test_df$gene_id[i]
		HOG_n  <- test_df$HOG[i]
		hard_n <- test_df$chr_hard[i]
		soft_n <- test_df$chr_soft[i]
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
	}
}

all_orths <- unique(all_orths)

soft_chr_dict[["TBI_07206"]]

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
  
  test <- c()
  for(i in seq(1, length(pos_sel_df[,1]))){
    gene_n <- pos_sel_df$gene_name[i]
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
  }

  pos_sel_df$chr_soft_class <- chr_soft_class
  pos_sel_df$chr_hard_class <- chr_hard_class
  pos_sel_df$RT_logFC <- RT_logFC
  pos_sel_df$RT_FDR   <- RT_FDR 
  pos_sel_df$HD_logFC <- HD_logFC
  pos_sel_df$HD_FDR   <- HD_FDR 
  pos_sel_df$LG_logFC <- LG_logFC
  pos_sel_df$LG_FDR   <- LG_FDR 
  pos_sel_df <- subset(pos_sel_df, !is.na(pos_sel_df$chr_soft_class))
  
  print(length(pos_sel_df[,1]))  
  return(pos_sel_df)
}

pos_sel_dat_2 <- positive_sel_by_chr(pos_sel_dat)




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


test_possel_overrep_X(pos_sel_dat_2, 0.05, "hard")


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

