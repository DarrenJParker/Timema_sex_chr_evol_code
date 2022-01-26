# GC_genomic.R

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
Tbi_dat <- read.csv(file="data/output/MF_cov/Tbi_MFcov_filt_1000_wGC.csv")
Tce_dat <- read.csv(file="data/output/MF_cov/Tce_MFcov_filt_1000_wGC.csv")
Tcm_dat <- read.csv(file="data/output/MF_cov/Tcm_MFcov_filt_1000_wGC.csv")
Tpa_dat <- read.csv(file="data/output/MF_cov/Tpa_MFcov_filt_1000_wGC.csv")
Tps_dat <- read.csv(file="data/output/MF_cov/Tps_MFcov_filt_1000_wGC.csv")





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



plot_GC <- function(df){
  
  sp <- gsub("_dat", "", deparse(substitute(df)))
  print(sp)
  ###################################################################################
  ### soft - means
  GC_means_soft       <- summarySE(df , measurevar="GC", groupvars=c("class_soft"))
  GC_means_soft$upper <- GC_means_soft$GC +  GC_means_soft$se
  GC_means_soft$lower <- GC_means_soft$GC -  GC_means_soft$se
  
  ## plot
  P_mean_GC_soft <- ggplot() + 
    geom_errorbar(data=GC_means_soft, mapping=aes(x=class_soft, ymin=upper, ymax=lower,  color= class_soft), width=0.2, size=1) + 
    theme_bw() +
    geom_point(data=GC_means_soft, mapping=aes(x = class_soft, y=GC,  fill=class_soft), size=4, shape=21) + 
    scale_color_manual(values=c("A" = "darkgrey", "X" = "darkorange2")) +
    scale_fill_manual(values=c("A" = "darkgrey", "X" = "darkorange2"))  + ggtitle(paste(sp, ",soft chr, mean +/- SE", sep = ""))
  
  
  ###################################################################################
  ### LG - means
  df_LG <- subset(df, ! is.na(as.character(df$LG)))
  df_LG$LG_ord <- ordered(df_LG$LG, levels=c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8", "lg9", "lg10", "lg11", "lg12", "lgX" ))
  
  GC_means_LG       <- summarySE(df_LG, measurevar="GC", groupvars=c("LG_ord"))
  GC_means_LG$upper <- GC_means_LG$GC +  GC_means_LG$se
  GC_means_LG$lower <- GC_means_LG$GC -  GC_means_LG$se
  
  ## plot
  P_mean_GC_LG <- ggplot() + 
    geom_errorbar(data=GC_means_LG, mapping=aes(x=LG_ord, ymin=upper, ymax=lower,  color= LG_ord), width=0.2, size=1) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_point(data=GC_means_LG, mapping=aes(x = LG_ord, y=GC,  fill=LG_ord), size=4, shape=21) + 
    scale_color_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) +
    scale_fill_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) + ggtitle(paste(sp, ", LG, mean +/- SE", sep = ""))
  
  
  ###################################################################################
  ### soft wt med

  wt_med_A <- weighted.median(subset(df, df$class_soft == "A")$GC, subset(df, df$class_soft == "A")$length)  
  wt_med_X <- weighted.median(subset(df, df$class_soft == "X")$GC, subset(df, df$class_soft == "X")$length)
  # print(wt_med_A )
  # print(wt_med_X )
  
  wt_med_soft_df <- as.data.frame(cbind(c(wt_med_A, wt_med_X), c("A", "X")))
  colnames(wt_med_soft_df) <- c("wtmed_GC", "chr_soft")
  wt_med_soft_df$wtmed_GC <- as.numeric(wt_med_soft_df$wtmed_GC )
  P_wt_med_GC_soft <- ggplot(wt_med_soft_df, aes(chr_soft, wtmed_GC, fill = chr_soft)) +
    geom_bar(position="dodge",stat="identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values=c("A" = "darkgrey", "X" = "darkorange2")) + ggtitle(paste(sp, ", soft, wt med", sep = ""))
  

  ###################################################################################
  ### LG wt med

  wt_med_lg1  <- weighted.median(subset(df, df$LG == "lg1" )$GC, subset(df, df$LG == "lg1" )$length)  
  wt_med_lg2  <- weighted.median(subset(df, df$LG == "lg2" )$GC, subset(df, df$LG == "lg2" )$length)  
  wt_med_lg3  <- weighted.median(subset(df, df$LG == "lg3" )$GC, subset(df, df$LG == "lg3" )$length)  
  wt_med_lg4  <- weighted.median(subset(df, df$LG == "lg4" )$GC, subset(df, df$LG == "lg4" )$length)  
  wt_med_lg5  <- weighted.median(subset(df, df$LG == "lg5" )$GC, subset(df, df$LG == "lg5" )$length)  
  wt_med_lg6  <- weighted.median(subset(df, df$LG == "lg6" )$GC, subset(df, df$LG == "lg6" )$length)  
  wt_med_lg7  <- weighted.median(subset(df, df$LG == "lg7" )$GC, subset(df, df$LG == "lg7" )$length)  
  wt_med_lg8  <- weighted.median(subset(df, df$LG == "lg8" )$GC, subset(df, df$LG == "lg8" )$length)  
  wt_med_lg9  <- weighted.median(subset(df, df$LG == "lg9" )$GC, subset(df, df$LG == "lg9" )$length)  
  wt_med_lg10 <- weighted.median(subset(df, df$LG == "lg10")$GC, subset(df, df$LG == "lg10")$length)  
  wt_med_lg11 <- weighted.median(subset(df, df$LG == "lg11")$GC, subset(df, df$LG == "lg11")$length)  
  wt_med_lg12 <- weighted.median(subset(df, df$LG == "lg12")$GC, subset(df, df$LG == "lg12")$length)  
  wt_med_lgX  <- weighted.median(subset(df, df$LG == "lgX" )$GC, subset(df, df$LG == "lgX" )$length)  
  
  wt_med_LG_df <- as.data.frame(cbind(c(wt_med_lg1, wt_med_lg2, wt_med_lg3, wt_med_lg4, wt_med_lg5, wt_med_lg6, wt_med_lg7, wt_med_lg8, wt_med_lg9, wt_med_lg10, wt_med_lg11, wt_med_lg12, wt_med_lgX), 
                  c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8", "lg9", "lg10", "lg11", "lg12", "lgX")))
  colnames(wt_med_LG_df) <- c("wtmed_GC", "LG")
  wt_med_LG_df$LG <- ordered( wt_med_LG_df$LG, levels=c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8", "lg9", "lg10", "lg11", "lg12", "lgX" ))
  
  wt_med_LG_df$wtmed_GC <- as.numeric(wt_med_LG_df$wtmed_GC )  
  #$print(  wt_med_LG_df)
  
  P_wt_med_GC_LG <- ggplot(wt_med_LG_df, aes(LG, wtmed_GC, fill = LG)) +
    geom_bar(position="dodge",stat="identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values=c("A" = "darkgrey", "lgX" = "darkorange2")) + ggtitle(paste(sp, ", LG, wt med", sep = "")) #+ coord_cartesian(ylim=c(0.3,0.37))
  

  out_list <- list("P_mean_GC_soft" = P_mean_GC_soft, "P_mean_GC_LG" = P_mean_GC_LG, "P_wt_med_GC_soft" = P_wt_med_GC_soft, "P_wt_med_GC_LG" = P_wt_med_GC_LG )
  return(out_list)
  
}


setwd("data/output/sel_out")

pdf("P_mean_GC_g_LG.pdf", width = 	8, height = 10)
plot_grid(
  plot_GC(Tbi_dat)$P_mean_GC_LG,
  plot_GC(Tce_dat)$P_mean_GC_LG,
  plot_GC(Tcm_dat)$P_mean_GC_LG,
  plot_GC(Tpa_dat)$P_mean_GC_LG,
  plot_GC(Tps_dat)$P_mean_GC_LG, ncol = 2
)
dev.off()
getwd() ## where has my plot gone....?



pdf("P_mean_GC_g_soft.pdf", width = 	12, height = 5)
plot_grid(
  plot_GC(Tbi_dat)$P_mean_GC_soft,
  plot_GC(Tce_dat)$P_mean_GC_soft,
  plot_GC(Tcm_dat)$P_mean_GC_soft,
  plot_GC(Tpa_dat)$P_mean_GC_soft,
  plot_GC(Tps_dat)$P_mean_GC_soft, ncol = 5
)
dev.off()
getwd() ## where has my plot gone....?



pdf("P_wt_med_GC_g_LG.pdf", width = 	8, height = 10)
plot_grid(
  plot_GC(Tbi_dat)$P_wt_med_GC_LG,
  plot_GC(Tce_dat)$P_wt_med_GC_LG,
  plot_GC(Tcm_dat)$P_wt_med_GC_LG,
  plot_GC(Tpa_dat)$P_wt_med_GC_LG,
  plot_GC(Tps_dat)$P_wt_med_GC_LG, ncol = 2
)
dev.off()
getwd() ## where has my plot gone....?




pdf("P_wt_med_GC_g_soft.pdf", width = 	12, height = 5)
plot_grid(
  plot_GC(Tbi_dat)$P_wt_med_GC_soft,
  plot_GC(Tce_dat)$P_wt_med_GC_soft,
  plot_GC(Tcm_dat)$P_wt_med_GC_soft,
  plot_GC(Tpa_dat)$P_wt_med_GC_soft,
  plot_GC(Tps_dat)$P_wt_med_GC_soft, ncol = 5
)

dev.off()
getwd() ## where has my plot gone....?




