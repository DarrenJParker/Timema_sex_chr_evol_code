## Male_and_female_coverage.R

### libs

library(ggplot2)
library(stringr)
library(modeest)
library(plyr)
library(cowplot)
library(grid)
library(matrixStats)
library(spatstat)
library(RColorBrewer)

########################################################################################################################################################################
### DATA

# file ext

infileext = "_pairedcov_minlen=1000_contig_cov_KF2_wLGinfopos_withangsDnucl.csv"

# raw
dat_Tbi <- read.table (paste("data/coverage_and_heterozygosity/Tbi", infileext, sep = ""), header = T, sep = ',')
dat_Tce <- read.table (paste("data/coverage_and_heterozygosity/Tce", infileext, sep = ""), header = T, sep = ',')
dat_Tcm <- read.table (paste("data/coverage_and_heterozygosity/Tcm", infileext, sep = ""), header = T, sep = ',')
dat_Tpa <- read.table (paste("data/coverage_and_heterozygosity/Tpa", infileext, sep = ""), header = T, sep = ',')
dat_Tps <- read.table (paste("data/coverage_and_heterozygosity/Tps", infileext, sep = ""), header = T, sep = ',')

### tidy header
colnames(dat_Tbi) <- gsub("_to_T.._v8_pe_BWA_mapqfilt_30aDR_contig_cov.txt", "", colnames(dat_Tbi))
colnames(dat_Tce) <- gsub("_to_T.._v8_pe_BWA_mapqfilt_30aDR_contig_cov.txt", "", colnames(dat_Tce))
colnames(dat_Tcm) <- gsub("_to_T.._v8_pe_BWA_mapqfilt_30aDR_contig_cov.txt", "", colnames(dat_Tcm))
colnames(dat_Tpa) <- gsub("_to_T.._v8_pe_BWA_mapqfilt_30aDR_contig_cov.txt", "", colnames(dat_Tpa))
colnames(dat_Tps) <- gsub("_to_T.._v8_pe_BWA_mapqfilt_30aDR_contig_cov.txt", "", colnames(dat_Tps))

head(dat_Tps)

#### set min contig length
cutoff_len = 1000 ### set to desired min contig length 

dat1_Tbi <- subset(dat_Tbi, dat_Tbi$length >= cutoff_len)
dat1_Tce <- subset(dat_Tce, dat_Tce$length >= cutoff_len)
dat1_Tcm <- subset(dat_Tcm, dat_Tcm$length >= cutoff_len)
dat1_Tpa <- subset(dat_Tpa, dat_Tpa$length >= cutoff_len)
dat1_Tps <- subset(dat_Tps, dat_Tps$length >= cutoff_len)

c(length(dat_Tbi[,1]), length(dat1_Tbi[,1]))
c(length(dat_Tce[,1]), length(dat1_Tce[,1]))
c(length(dat_Tcm[,1]), length(dat1_Tcm[,1]))
c(length(dat_Tpa[,1]), length(dat1_Tpa[,1]))
c(length(dat_Tps[,1]), length(dat1_Tps[,1]))


## get LG
dat1_Tbi$LG <- gsub("_.*", "", as.character(dat1_Tbi$lg_biggest_block))
dat1_Tce$LG <- gsub("_.*", "", as.character(dat1_Tce$lg_biggest_block))
dat1_Tcm$LG <- gsub("_.*", "", as.character(dat1_Tcm$lg_biggest_block))
dat1_Tpa$LG <- gsub("_.*", "", as.character(dat1_Tpa$lg_biggest_block))
dat1_Tps$LG <- gsub("_.*", "", as.character(dat1_Tps$lg_biggest_block))

head(dat1_Tbi)

### output

dir.create("data/output")
dir.create("data/output/MF_cov_nucl")
setwd("data/output/MF_cov_nucl")

########################################################################################################################################################################
### Functions

peakfinder <- function(d){
  dh <- hist(d,plot=FALSE, breaks=1000)
  ins <- dh[["counts"]]
  nbins <- length(ins)
  ss <- which(rank(ins)%in%seq(from=nbins,to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}


peakfinder_lessprecise <- function(d){
  dh <- hist(d,plot=FALSE, breaks=500)
  ins <- dh[["counts"]]
  nbins <- length(ins)
  ss <- which(rank(ins)%in%seq(from=nbins,to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}



### male to female cov

MF_cov_sum <- function(df, sp){
	
	female_covs <- df[,grep(paste("^",sp,"_F", sep = ""),colnames(df))]
	male_covs   <- df[,grep(paste("^",sp,"_M", sep = ""),colnames(df))]
	df$Female_cov_sum <- rowSums(female_covs)
	
	print(head(male_covs))
	
	print(length(male_covs))
	if (length(male_covs) > 100){
    	df$Male_cov_sum   <- male_covs
	} else {
   		df$Male_cov_sum   <- rowSums(male_covs)
	}	
	
	##### filter contigs with 0 cov in females or males 
	
	df_filt  <- subset(df , df$Male_cov_sum   > 0 & df$Female_cov_sum > 0)
	
	print(length(df[,1]))
	print(length(df_filt[,1]))

	##### contigs with 0 cov in females or males 
	
	MF_0   <- subset(df , df$Male_cov_sum   == 0 | df$Female_cov_sum == 0)
	print(length(	MF_0[,1]))
	
	## norm by modal cov
	
	male_mode_cov   = mlv(df_filt$Male_cov_sum, method = "shorth")
	female_mode_cov = mlv(df_filt$Female_cov_sum, method = "shorth")
	df_filt$Female_cov_sum_norm_mode = df_filt$Female_cov_sum / female_mode_cov
	df_filt$Male_cov_sum_norm_mode   = df_filt$Male_cov_sum /   male_mode_cov

	#### calc M to F ratio
	### Scaffolds were considered to be X candidates if they had Log2(M/F coverage) within the range [A_coord -1.1, A_coord -0.9];

	df_filt$M_F_mode <- log2(df_filt$Male_cov_sum_norm_mode / df_filt$Female_cov_sum_norm_mode)

	#### adjusted peak grab
	df_filt_cut <- subset(df_filt, df_filt$M_F_mode < -0.5) 

	Sex_chr_peak <- peakfinder(df_filt_cut$M_F_mode)
	Auto_peak    <- peakfinder(df_filt$M_F_mode)

	if(length(Sex_chr_peak) == 0){
		Sex_chr_peak <- peakfinder_lessprecise(df_filt_cut$M_F_mode)	
	}
	
	print(Sex_chr_peak)


	df_filt$class_soft <- ifelse(df_filt$M_F_mode < Auto_peak - 0.5, "X", "A")
	df_filt$class_hard <- ifelse(df_filt$M_F_mode > Sex_chr_peak - 0.1 & df_filt$M_F_mode < Sex_chr_peak + 0.1, "X", "A")
	df_filt$to_fill <- paste(df_filt$class_soft, df_filt$class_hard, sep = "")
	
	print(head(df_filt))
	
	p2_adj <- ggplot(df_filt, aes(x=M_F_mode, fill = to_fill)) +
    	theme_bw() +
    	geom_histogram(binwidth=0.01, show.legend = FALSE) +
    	coord_cartesian(xlim=c(-3,3))  +
    	geom_hline(yintercept = 0) +
    	xlab("log2(Male cov / Female cov)") + 
    ylab("Number of contigs") +
    	ggtitle(sp) + 
    	geom_vline(xintercept = Sex_chr_peak, linetype="dashed", color = "black", size=0.2) +
    	geom_vline(xintercept = Auto_peak, linetype="dotted", color = "black", size=0.2) + 
    	scale_fill_manual(values = c("darkgrey", "darkorange2", "red3")) 

	out_list = list("p2_adj" = p2_adj, "Sex_chr_peak" = Sex_chr_peak,  "Auto_peak" = Auto_peak,  "df_filt" = df_filt)
	return(out_list)	
}


Tbi_out <- MF_cov_sum(dat1_Tbi, "Tbi")
#Tbi_out$p2_adj
Tce_out <- MF_cov_sum(dat1_Tce, "Tce")
Tcm_out <- MF_cov_sum(dat1_Tcm, "Tcm")

#### I know that "Tps_F_ReSeq_Ps08", "Tps_F_ReSeq_Ps12" have very bad coverage for some reason. Here I exclude.
Tps_drop <- c("Tps_F_ReSeq_Ps08", "Tps_F_ReSeq_Ps12")
dat1_Tps_c <- dat1_Tps[ , !(names(dat1_Tps) %in% Tps_drop)]
Tps_out <- MF_cov_sum(dat1_Tps_c, "Tps")

#### I know that "Tpa_F_H56" has very bad coverage for some reason. Here I exclude
Tpa_drop <- c("Tpa_F_H56")
dat1_Tpa_c <- dat1_Tpa[ , !(names(dat1_Tpa) %in% Tpa_drop)]
Tpa_out <- MF_cov_sum(dat1_Tpa_c, "Tpa")



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### plot Nucl

Tbi_df_filt <- Tbi_out$df_filt
Tce_df_filt <- Tce_out$df_filt
Tcm_df_filt <- Tcm_out$df_filt
Tpa_df_filt <- Tpa_out$df_filt
Tps_df_filt <- Tps_out$df_filt


Nucl_plot <- function(df){
	
	df_filt <- subset(df, df$nSites >= 1000) ## filt by N covered sites
	
	## calc pi
	
	df_filt$nucl_div <- df_filt$tP / df_filt$nSites 

	X_soft <- subset(df_filt, df_filt$class_soft == 'X')
	A_soft <- subset(df_filt, df_filt$class_soft == 'A')
	wt_med_X_nucl_soft <- weighted.median(X_soft$nucl_div, X_soft$nSites)
	wt_med_A_nucl_soft <- weighted.median(A_soft$nucl_div, A_soft$nSites)
	XA_ratio_soft <- wt_med_X_nucl_soft / wt_med_A_nucl_soft
	
	X_hard <- subset(df_filt, df_filt$class_hard == 'X')
	A_hard <- subset(df_filt, df_filt$class_hard == 'A')
	wt_med_X_nucl_hard <- weighted.median(X_hard$nucl_div, X_hard$nSites)
	wt_med_A_nucl_hard <- weighted.median(A_hard$nucl_div, A_hard$nSites)
	XA_ratio_hard <- wt_med_X_nucl_hard / wt_med_A_nucl_hard
	
	out_list <- list("wt_med_X_nucl_soft" = wt_med_X_nucl_soft, "wt_med_A_nucl_soft" = wt_med_A_nucl_soft, "XA_ratio_soft" = XA_ratio_soft, "wt_med_X_nucl_hard" = wt_med_X_nucl_hard, "wt_med_A_nucl_hard" = wt_med_A_nucl_hard, "XA_ratio_hard" = XA_ratio_hard)
	
	return(out_list)
}


### ratio X A pi

soft_nucl_XA_df <- as.data.frame(
rbind(
c("Tbi", Nucl_plot(Tbi_df_filt)$XA_ratio_soft),
c("Tce", Nucl_plot(Tce_df_filt)$XA_ratio_soft),
c("Tcm", Nucl_plot(Tcm_df_filt)$XA_ratio_soft),
c("Tpa", Nucl_plot(Tpa_df_filt)$XA_ratio_soft),
c("Tps", Nucl_plot(Tps_df_filt)$XA_ratio_soft)
))

colnames(soft_nucl_XA_df) <- c("sp", "X_A_nucl")
soft_nucl_XA_df$X_A_nucl  <- as.numeric(as.character(soft_nucl_XA_df$X_A_nucl))


XA_nucl_ratio_soft <- ggplot(soft_nucl_XA_df , aes(sp, X_A_nucl, fill = sp)) + 
		geom_bar(position="dodge",stat="identity", colour="black") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		ylim(c(0,1)) + geom_hline(yintercept= 0.75,  linetype="dashed") + 
		xlab ("Species") + 
		ylab ("X:A Nucl") + scale_fill_brewer(palette = "Set3") + 
		scale_y_continuous(expand = c(0,0), limits = c(0,1)) + ggtitle("soft")

pdf(paste("XA_Pnucl_ratio_soft_class_", cutoff_len, ".pdf", sep = ""), width = 	3, height = 10)
XA_nucl_ratio_soft 
dev.off()
getwd() ## where has my plot gone....?


### ratio X A pi

hard_nucl_XA_df <- as.data.frame(
rbind(
c("Tbi", Nucl_plot(Tbi_df_filt)$XA_ratio_hard),
c("Tce", Nucl_plot(Tce_df_filt)$XA_ratio_hard),
c("Tcm", Nucl_plot(Tcm_df_filt)$XA_ratio_hard),
c("Tpa", Nucl_plot(Tpa_df_filt)$XA_ratio_hard),
c("Tps", Nucl_plot(Tps_df_filt)$XA_ratio_hard)
))

colnames(hard_nucl_XA_df) <- c("sp", "X_A_nucl")
hard_nucl_XA_df$X_A_nucl  <- as.numeric(as.character(hard_nucl_XA_df$X_A_nucl))


XA_nucl_ratio_hard <- ggplot(hard_nucl_XA_df , aes(sp, X_A_nucl, fill = sp)) + 
		geom_bar(position="dodge",stat="identity", colour="black") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		ylim(c(0,1)) + geom_hline(yintercept= 0.75,  linetype="dashed") + 
		xlab ("Species") + 
		ylab ("X:A Nucl") + scale_fill_brewer(palette = "Set3") + 
		scale_y_continuous(expand = c(0,0), limits = c(0,1)) + ggtitle("hard")

pdf(paste("XA_Pnucl_ratio_hard_class_", cutoff_len, ".pdf", sep = ""), width = 	3, height = 10)
XA_nucl_ratio_hard 
dev.off()
getwd() ## where has my plot gone....?






##############################################################################################################################
### actual values

nucl_XA_soft <- as.data.frame(cbind(
c("Tbi", "Tbi", "Tce", "Tce", "Tcm", "Tcm", "Tpa", "Tpa", "Tps", "Tps"),
c("A", "X","A", "X","A", "X","A", "X","A", "X"),
c(
Nucl_plot(Tbi_df_filt)$wt_med_A_nucl_soft,
Nucl_plot(Tbi_df_filt)$wt_med_X_nucl_soft,
Nucl_plot(Tce_df_filt)$wt_med_A_nucl_soft,
Nucl_plot(Tce_df_filt)$wt_med_X_nucl_soft,
Nucl_plot(Tcm_df_filt)$wt_med_A_nucl_soft,
Nucl_plot(Tcm_df_filt)$wt_med_X_nucl_soft,
Nucl_plot(Tpa_df_filt)$wt_med_A_nucl_soft,
Nucl_plot(Tpa_df_filt)$wt_med_X_nucl_soft,
Nucl_plot(Tps_df_filt)$wt_med_A_nucl_soft,
Nucl_plot(Tps_df_filt)$wt_med_X_nucl_soft)))


colnames(nucl_XA_soft) <- c("sp", "chr", "pi")

nucl_XA_soft$pi <- as.numeric(nucl_XA_soft$pi)
str(nucl_XA_soft)


nucl_XA_soft_p1 <- ggplot(nucl_XA_soft, aes(sp, pi, fill = chr)) + 
		geom_bar(position="dodge",stat="identity") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		xlab ("Sample") + 
		ylab ("Nucl div") +
		scale_fill_manual(values=c("darkgrey", "darkorange2")) + ggtitle("soft") 


nucl_XA_hard <- as.data.frame(cbind(
c("Tbi", "Tbi", "Tce", "Tce", "Tcm", "Tcm", "Tpa", "Tpa", "Tps", "Tps"),
c("A", "X","A", "X","A", "X","A", "X","A", "X"),
c(
Nucl_plot(Tbi_df_filt)$wt_med_A_nucl_hard,
Nucl_plot(Tbi_df_filt)$wt_med_X_nucl_hard,
Nucl_plot(Tce_df_filt)$wt_med_A_nucl_hard,
Nucl_plot(Tce_df_filt)$wt_med_X_nucl_hard,
Nucl_plot(Tcm_df_filt)$wt_med_A_nucl_hard,
Nucl_plot(Tcm_df_filt)$wt_med_X_nucl_hard,
Nucl_plot(Tpa_df_filt)$wt_med_A_nucl_hard,
Nucl_plot(Tpa_df_filt)$wt_med_X_nucl_hard,
Nucl_plot(Tps_df_filt)$wt_med_A_nucl_hard,
Nucl_plot(Tps_df_filt)$wt_med_X_nucl_hard)))


colnames(nucl_XA_hard) <- c("sp", "chr", "pi")

nucl_XA_hard$pi <- as.numeric(nucl_XA_hard$pi)
str(nucl_XA_hard)


nucl_XA_hard_p1 <- ggplot(nucl_XA_hard, aes(sp, pi, fill = chr)) + 
		geom_bar(position="dodge",stat="identity") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		xlab ("Sample") + 
		ylab ("Nucl div") +
		scale_fill_manual(values=c("darkgrey", "red3")) + ggtitle("hard") 



pdf(paste("nucl_XA_soft_p1_", cutoff_len, ".pdf", sep = ""), width = 	4, height = 5)
nucl_XA_soft_p1
dev.off()
getwd() ## where has my plot gone....?



pdf(paste("nucl_XA_hard_p1_", cutoff_len, ".pdf", sep = ""), width = 	4, height = 5)
nucl_XA_hard_p1
dev.off()
getwd() ## where has my plot gone....?



#### as table with Ne

nucl_XA_soft_Ne <- as.data.frame(rbind(
c(Nucl_plot(Tbi_df_filt)$wt_med_X_nucl_soft, Nucl_plot(Tbi_df_filt)$wt_med_A_nucl_soft, Nucl_plot(Tbi_df_filt)$XA_ratio_soft),
c(Nucl_plot(Tce_df_filt)$wt_med_X_nucl_soft, Nucl_plot(Tce_df_filt)$wt_med_A_nucl_soft, Nucl_plot(Tce_df_filt)$XA_ratio_soft),
c(Nucl_plot(Tcm_df_filt)$wt_med_X_nucl_soft, Nucl_plot(Tcm_df_filt)$wt_med_A_nucl_soft, Nucl_plot(Tcm_df_filt)$XA_ratio_soft),
c(Nucl_plot(Tpa_df_filt)$wt_med_X_nucl_soft, Nucl_plot(Tpa_df_filt)$wt_med_A_nucl_soft, Nucl_plot(Tpa_df_filt)$XA_ratio_soft),
c(Nucl_plot(Tps_df_filt)$wt_med_X_nucl_soft, Nucl_plot(Tps_df_filt)$wt_med_A_nucl_soft, Nucl_plot(Tps_df_filt)$XA_ratio_soft)))

colnames(nucl_XA_soft_Ne) <- c("wt_med_X_nucl_soft", "wt_med_A_nucl_soft", "XA_ratio_soft")
nucl_XA_soft_Ne$sp <- c("Tbi", "Tce", "Tcm", "Tpa", "Tps")

mutation_rate_hel = 2.9E-09	## Heliconius Keightley PD, Pinharanda A, Ness RW, Simpson F, Dasmahapatra KK, Mallet J, et al. (2015). Estimation of the spontaneous mutation rate in Heliconius melpomene. Mol Biol Evol 32: 239–243.
mutation_rate_dro = 2.8E-09 ## Drosophila Keightley PD, Ness RW, Halligan DL, Haddrill PR (2014). Estimation of the spontaneous mutation rate per nucleotide site in a Drosophila melanogaster full-sib family. Genetics 196: 313–320.
## pi = 4Neu
nucl_XA_soft_Ne$NeX_hel <- nucl_XA_soft_Ne$wt_med_X_nucl_soft / 4 / mutation_rate_hel
nucl_XA_soft_Ne$NeA_hel <- nucl_XA_soft_Ne$wt_med_A_nucl_soft / 4 / mutation_rate_hel
nucl_XA_soft_Ne$NeX_dro <- nucl_XA_soft_Ne$wt_med_X_nucl_soft / 4 / mutation_rate_dro
nucl_XA_soft_Ne$NeA_dro <- nucl_XA_soft_Ne$wt_med_A_nucl_soft / 4 / mutation_rate_dro


write.csv(nucl_XA_soft_Ne, "nucl_XA_soft_Ne.csv")

##########################################################################################################
### by LG  



pi_per_LG <- function(df){
  
  df_filt <- subset(df, df$nSites >= 1000) ## filt by N covered sites
  
  ## calc pi
  
  df_filt$nucl_div <- df_filt$tP / df_filt$nSites 

  df_filt_lg1   <- subset(df_filt, df_filt$LG == "lg1")	
  df_filt_lg2   <- subset(df_filt, df_filt$LG == "lg2")
  df_filt_lg3   <- subset(df_filt, df_filt$LG == "lg3")
  df_filt_lg4   <- subset(df_filt, df_filt$LG == "lg4")
  df_filt_lg5   <- subset(df_filt, df_filt$LG == "lg5")
  df_filt_lg6   <- subset(df_filt, df_filt$LG == "lg6")	
  df_filt_lg7   <- subset(df_filt, df_filt$LG == "lg7")
  df_filt_lg8   <- subset(df_filt, df_filt$LG == "lg8")
  df_filt_lg9   <- subset(df_filt, df_filt$LG == "lg9")
  df_filt_lg10  <- subset(df_filt, df_filt$LG == "lg10")
  df_filt_lg11  <- subset(df_filt, df_filt$LG == "lg11")
  df_filt_lg12  <- subset(df_filt, df_filt$LG == "lg12")
  df_filt_lgX   <- subset(df_filt, df_filt$LG == "lgX")
  
  
  wt_med_lg1_nucl  <- weighted.median(df_filt_lg1$nucl_div,  df_filt_lg1$nSites)
  wt_med_lg2_nucl  <- weighted.median(df_filt_lg2$nucl_div,  df_filt_lg2$nSites)
  wt_med_lg3_nucl  <- weighted.median(df_filt_lg3$nucl_div,  df_filt_lg3$nSites)
  wt_med_lg4_nucl  <- weighted.median(df_filt_lg4$nucl_div,  df_filt_lg4$nSites)
  wt_med_lg5_nucl  <- weighted.median(df_filt_lg5$nucl_div,  df_filt_lg5$nSites)
  wt_med_lg6_nucl  <- weighted.median(df_filt_lg6$nucl_div,  df_filt_lg6$nSites)
  wt_med_lg7_nucl  <- weighted.median(df_filt_lg7$nucl_div,  df_filt_lg7$nSites)
  wt_med_lg8_nucl  <- weighted.median(df_filt_lg8$nucl_div,  df_filt_lg8$nSites)
  wt_med_lg9_nucl  <- weighted.median(df_filt_lg9$nucl_div,  df_filt_lg9$nSites)
  wt_med_lg10_nucl <- weighted.median(df_filt_lg10$nucl_div, df_filt_lg10$nSites)
  wt_med_lg11_nucl <- weighted.median(df_filt_lg11$nucl_div, df_filt_lg11$nSites)
  wt_med_lg12_nucl <- weighted.median(df_filt_lg12$nucl_div, df_filt_lg12$nSites)
  wt_med_lgX_nucl  <- weighted.median(df_filt_lgX$nucl_div,  df_filt_lgX$nSites)
  
  out_list <- list(
  "wt_med_lg1_nucl"  = wt_med_lg1_nucl,
  "wt_med_lg2_nucl"  = wt_med_lg2_nucl,
  "wt_med_lg3_nucl"  = wt_med_lg3_nucl,
  "wt_med_lg4_nucl"  = wt_med_lg4_nucl,
  "wt_med_lg5_nucl"  = wt_med_lg5_nucl,
  "wt_med_lg6_nucl"  = wt_med_lg6_nucl,
  "wt_med_lg7_nucl"  = wt_med_lg7_nucl,
  "wt_med_lg8_nucl"  = wt_med_lg8_nucl,
  "wt_med_lg9_nucl"  = wt_med_lg9_nucl,
  "wt_med_lg10_nucl" = wt_med_lg10_nucl,
  "wt_med_lg11_nucl" = wt_med_lg11_nucl,
  "wt_med_lg12_nucl" = wt_med_lg12_nucl,
  "wt_med_lgX_nucl"  = wt_med_lgX_nucl)
  
  return(out_list)
}

plot_pi_LG <- function(pi_LG_vals, tit_text){
  
  LG <- c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8", "lg9", "lg10", "lg11", "lg12", "lgX")
  
  wt_med_pi <- c(
    pi_LG_vals$wt_med_lg1_nucl,
    pi_LG_vals$wt_med_lg2_nucl,
    pi_LG_vals$wt_med_lg3_nucl,
    pi_LG_vals$wt_med_lg4_nucl,
    pi_LG_vals$wt_med_lg5_nucl,
    pi_LG_vals$wt_med_lg6_nucl,
    pi_LG_vals$wt_med_lg7_nucl,
    pi_LG_vals$wt_med_lg8_nucl,
    pi_LG_vals$wt_med_lg9_nucl,
    pi_LG_vals$wt_med_lg10_nucl,
    pi_LG_vals$wt_med_lg11_nucl,
    pi_LG_vals$wt_med_lg12_nucl,
    pi_LG_vals$wt_med_lgX_nucl)
  
  df1 <- as.data.frame(cbind(LG,  wt_med_pi))
  df1$wt_med_pi <- as.numeric(df1$wt_med_pi)
  df1$LG_ord <- ordered(df1$LG, levels=c("lg1","lg2","lg3","lg4","lg5","lg6","lg7","lg8","lg9","lg10","lg11","lg12","lgX")) 
  
   p1 <- ggplot(df1, aes(LG_ord, wt_med_pi, fill = LG_ord)) + 
    geom_bar(position="dodge",stat="identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab ("Linkage group") + 
    ylab ("nucleotide diversity") + 
    scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "yellow2","#666666","lightblue","royalblue2", "darkorchid", "red3")) + ggtitle(tit_text) 
  return(p1)

}  

Tbi_pi_LG <- pi_per_LG(Tbi_df_filt)
Tce_pi_LG <- pi_per_LG(Tce_df_filt)
Tcm_pi_LG <- pi_per_LG(Tcm_df_filt)
Tpa_pi_LG <- pi_per_LG(Tpa_df_filt)
Tps_pi_LG <- pi_per_LG(Tps_df_filt)


pdf(paste("Pi_LG_", cutoff_len, ".pdf", sep = ""), width = 	9, height = 15)
plot_grid(
  plot_pi_LG(Tbi_pi_LG, "Tbi"),
  plot_pi_LG(Tce_pi_LG, "Tce"),
  plot_pi_LG(Tcm_pi_LG, "Tcm"),
  plot_pi_LG(Tpa_pi_LG, "Tpa"),
  plot_pi_LG(Tbi_pi_LG, "Tps"),
  ncol = 2)
dev.off()
getwd() ## where has my plot gone....?











########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), paste("Male_and_female_coverage.R_sessionInfo_len_", cutoff_len, ".txt"), sep = "")



































