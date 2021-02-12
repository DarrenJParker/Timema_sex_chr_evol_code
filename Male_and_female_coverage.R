## Male_and_female_coverage.R

### libs

library(ggplot2)
library(stringr)
library(modeest)
library(plyr)
library(cowplot)
library(grid)



########################################################################################################################################################################
### DATA

# file ext

infileext = "_pairedcov_minlen=1000_contig_cov_KF2_wLGinfopos_withangsD.csv"

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

### output

dir.create("data/output")
dir.create("data/output/MF_cov")
setwd("data/output/MF_cov")

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




##########################################################################################################################
##### cov of each sample

cov_plot <- function(df, samp, min_cov, max_cov) {
	
	median_cov <- median(eval(parse(text=paste(df,"$",samp, sep = ""))))
	modal_cov  <- mlv(eval(parse(text=paste(df,"$",samp, sep = ""))), method = "hsm")
	q99 <- quantile(eval(parse(text=paste(df,"$",samp, sep = ""))), .99)

	p1 <- ggplot(eval(parse(text=df)), aes(x=eval(parse(text=samp)))) +
   		theme_bw() +
    	geom_histogram(color="darkgrey", fill="darkgrey", binwidth=0.2,alpha = 1) +
    xlim(c(min_cov, max_cov)) + 
    	geom_hline(yintercept = 0) +
    	xlab("Coverage") + 
    ylab("N contigs") + 
    	ggtitle(paste(samp, " | len >= ", cutoff_len, sep = "")) +
    	geom_vline(xintercept = median_cov, color = "blue", size=1) + 
    	geom_vline(xintercept = q99, linetype="dotted", color = "black", size=1)+    	
     geom_vline(xintercept = modal_cov, linetype="dashed", color = "red", size=1.5)    	

	out_df <- as.data.frame(rbind(c(samp, round(median_cov, digits=2), round(modal_cov, digits=2), round(q99, digits=2))))
	colnames(out_df) <- c("sample", "med_cov", "mod_cov", "q99")

	output = list("p1" = p1, "out_df" = out_df)
	return(output)


}


cov_plot_LG_plus_ove <- function(df, samp) {

	df1 <- eval(parse(text=df))
	df1$LG <- as.character(df1$LG)

	p1 <- ggplot(df1, aes(x=eval(parse(text=samp)))) +
		theme_bw() +
		geom_line(data=subset(df1,LG == "lg1"), color="#1B9E77", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg2"), color="#D95F02", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg3"), color="#7570B3", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg4"), color="#E7298A", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg5"), color="#66A61E", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg6"), color="#E6AB02", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg7"), color="#A6761D", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg8"), color="yellow2", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg9"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg10"), color="lightblue", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg11"), color="royalblue2", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg12"), color="darkorchid", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lgX"), color="red3", size = 1, stat="density") +
		geom_line(data=subset(df1,LG != "NA"), color="black", size = 1, stat="density") +	
		ggtitle(paste(samp, sep = "")) +
		xlim(c(0,40)) + 
		xlab("Coverage")
	
	fake_dat <- as.data.frame(cbind(
		seq(1,14),
		c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8", "lg9", "lg10", "lg11", "lg12", "lgX", "all_lg")
		))

	p222 <- ggplot(fake_dat, aes(V2, V1)) + 
		theme_bw() +
		geom_point(aes(colour = factor(V2)), shape = 95, size = 15) +
		scale_colour_manual(values=c(
		"lg1" ="#1B9E77",
		"lg2" ="#D95F02",
		"lg3" ="#7570B3",
		"lg4" ="#E7298A", 
		"lg5" ="#66A61E",
		"lg6" ="#E6AB02", 
		"lg7" ="#A6761D", 
		"lg8" ="yellow2",
		"lg9" ="#666666",
		"lg10" ="lightblue",
		"lg11" ="royalblue2",
		"lg12" ="darkorchid",
		"lgX"  = "red3",
		"all_lg" = "black"
		)) + theme(legend.title=element_blank())

	legend_out <- cowplot::get_legend(p222)
	
	out_list = list("p1" = p1, "legend" = legend_out)
	return(out_list )	
}


cov_plot_LG <- function(df, samp) {

	df1 <- eval(parse(text=df))
	df1$LG <- as.character(df1$LG)

	p1 <- ggplot(df1, aes(x=eval(parse(text=samp)))) +
		theme_bw() +
		geom_line(data=subset(df1,LG == "lg1"), color="#1B9E77", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg2"), color="#D95F02", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg3"), color="#7570B3", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg4"), color="#E7298A", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg5"), color="#66A61E", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg6"), color="#E6AB02", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg7"), color="#A6761D", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg8"), color="yellow2", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg9"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg10"), color="lightblue", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg11"), color="royalblue2", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg12"), color="darkorchid", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lgX"), color="red3", size = 1, stat="density") +
		#geom_line(data=subset(df1,LG != "NA"), color="black", size = 1, stat="density") +	
		ggtitle(paste(samp, sep = "")) +
		xlim(c(0,40)) + 
		#geom_hline(yintercept=0, colour="black", size=1) +
		xlab("Coverage")
		
	fake_dat <- as.data.frame(cbind(
		seq(1,13),
		c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8", "lg9", "lg10", "lg11", "lg12", "lgX")
		))

	p222 <- ggplot(fake_dat, aes(V2, V1)) + 
		theme_bw() +
		geom_point(aes(colour = factor(V2)), shape = 95, size = 15) +
		scale_colour_manual(values=c(
		"lg1" ="#1B9E77",
		"lg2" ="#D95F02",
		"lg3" ="#7570B3",
		"lg4" ="#E7298A", 
		"lg5" ="#66A61E",
		"lg6" ="#E6AB02", 
		"lg7" ="#A6761D", 
		"lg8" ="yellow2",
		"lg9" ="#666666",
		"lg10" ="lightblue",
		"lg11" ="royalblue2",
		"lg12" ="darkorchid",
		"lgX"  = "red3"
		)) + theme(legend.title=element_blank())

	legend_out <- cowplot::get_legend(p222)
	
	out_list = list("p1" = p1, "legend" = legend_out)
	return(out_list )	

}

cov_plot_LG_XA <- function(df, samp) {

	df1 <- eval(parse(text=df))
	df1$LG <- as.character(df1$LG)

	p1 <- ggplot(df1, aes(x=eval(parse(text=samp)))) +
		theme_bw() +

		geom_line(data=subset(df1,LG == "lg1"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg2"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg3"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg4"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg5"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg6"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg7"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg8"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg9"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg10"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg11"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lg12"), color="#666666", size = 1, stat="density") +
		geom_line(data=subset(df1,LG == "lgX"), color="red3", size = 1, stat="density") +
		#geom_line(data=subset(df1,LG != "NA"), color="black", size = 1, stat="density") +	
		ggtitle(paste(samp, sep = "")) +
		xlim(c(0,40)) + 
		xlab("Coverage")

	fake_dat <- as.data.frame(cbind(
		seq(1,2),
		c("lgX", "A")
		))

	p222 <- ggplot(fake_dat, aes(V2, V1)) + 
		theme_bw() +
		geom_point(aes(colour = factor(V2)), shape = 95, size = 15) +
		scale_colour_manual(values=c(
		"A" ="#666666",
		"lgX"  = "red3"
		)) + theme(legend.title=element_blank())

	legend_out <- cowplot::get_legend(p222)
	
	out_list = list("p1" = p1, "legend" = legend_out)
	return(out_list )	


}

 
#################################################################################################################
#################################################################################################################
### Plot contig coverage

### Tbi

min_cov_1 = -1
max_cov_1 = 100


Tbi_a1 <- plot_grid(
cov_plot("dat1_Tbi", "Tbi_M_13_Tbi", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tbi", "Tbi_M_14_Tbi", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tbi", "Tbi_M_15_Tbi", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tbi", "Tbi_M_16_Tbi", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tbi", "Tbi_F_CC86B", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tbi", "Tbi_F_CC86C", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tbi", "Tbi_F_CC87B", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tbi", "Tbi_F_CC87C", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tbi", "Tbi_F_CC88B", min_cov_1, max_cov_1)$p1,
ncol = 2, nrow = 5)

pdf(paste("Tbi_cov_min=", min_cov_1, "_max=", max_cov_1, ".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tbi_a1, ncol = 1)
dev.off()
getwd() ## where has my plot gone....?

 
### Tce

Tce_a1 <- plot_grid(
cov_plot("dat1_Tce", "Tce_M_05_HM15", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce", "Tce_M_06_HM16", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce", "Tce_M_07_HM33", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce", "Tce_M_08_HM61", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce", "Tce_F_CC22B", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce", "Tce_F_CC22C", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce", "Tce_F_CC24B", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce", "Tce_F_CC24C", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce", "Tce_F_CC25B", min_cov_1, max_cov_1)$p1,
ncol = 2, nrow = 5)

pdf(paste("Tce_cov_min=", min_cov_1, "_max=", max_cov_1, ".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tce_a1, ncol = 1)
dev.off()
getwd() ## where has my plot gone....?

### Tcm

Tcm_a1 <- plot_grid(
cov_plot("dat1_Tcm", "Tcm_M_01_HM148", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tcm", "Tcm_M_02_HM149", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tcm", "Tcm_M_03_HM150", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tcm", "Tcm_M_04_HM151", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tcm", "Tcm_F_HM217", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tcm", "Tcm_F_HM218", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tcm", "Tcm_F_HM219", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tcm", "Tcm_F_HM220", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tcm", "Tcm_F_HM221", min_cov_1, max_cov_1)$p1,
ncol = 2, nrow = 5)

pdf(paste("Tcm_cov_min=", min_cov_1, "_max=", max_cov_1, ".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tcm_a1, ncol = 1)
dev.off()
getwd() ## where has my plot gone....?

### Tpa


Tpa_a1 <- plot_grid(
cov_plot("dat1_Tpa", "Tpa_M_09_Tpa", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tpa", "Tpa_M_10_Tpa", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tpa", "Tpa_M_11_Tpa", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tpa", "Tpa_M_12_Tpa", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tpa", "Tpa_F_H54", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tpa", "Tpa_F_H56", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tpa", "Tpa_F_PA_CD", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tpa", "Tpa_F_PA_E", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tpa", "Tpa_F_Pa_AB", min_cov_1, max_cov_1)$p1,
ncol = 2, nrow = 5)

pdf(paste("Tpa_cov_min=", min_cov_1, "_max=", max_cov_1, ".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tpa_a1, ncol = 1)
dev.off()
getwd() ## where has my plot gone....?



### Tps


Tps_a1 <- plot_grid(
cov_plot("dat1_Tps", "Tps_M_17_HM99", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tps", "Tps_M_18_HM100", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tps", "Tps_M_19_HM101", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tps", "Tps_M_20_15255", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tps", "Tps_F_ReSeq_Ps08", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tps", "Tps_F_ReSeq_Ps12", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tps", "Tps_F_ReSeq_Ps14", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tps", "Tps_F_ReSeq_Ps16", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tps", "Tps_F_ReSeq_Ps18", min_cov_1, max_cov_1)$p1,
ncol = 2, nrow = 5)

pdf(paste("Tps_cov_min=", min_cov_1, "_max=", max_cov_1, ".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tps_a1, ncol = 1)
dev.off()
getwd() ## where has my plot gone....?


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
#Tpa_out$p2_adj


######################################################################################################################
##### export tables and plots

### sex chromosome adjusted peak
sex_chr_peaks_df <- as.data.frame(
cbind(c("Tbi","Tce","Tcm","Tpa","Tps"),
c(
Tbi_out$Sex_chr_peak,
Tce_out$Sex_chr_peak,
Tcm_out$Sex_chr_peak,
Tpa_out$Sex_chr_peak,
Tps_out$Sex_chr_peak
)))

colnames(sex_chr_peaks_df) <- c("Sp", "sex_chr_peaks")
write.csv(sex_chr_peaks_df, file=paste("sex_chr_peaks_", cutoff_len,  ".csv", sep = ""), row.names=FALSE)

###### full (filtered) table

head(Tbi_out$df_filt)

write.csv(Tbi_out$df_filt, file=paste("Tbi_MFcov_filt_", cutoff_len,  ".csv", sep = ""), row.names=FALSE)
write.csv(Tce_out$df_filt, file=paste("Tce_MFcov_filt_", cutoff_len,  ".csv", sep = ""), row.names=FALSE)
write.csv(Tcm_out$df_filt, file=paste("Tcm_MFcov_filt_", cutoff_len,  ".csv", sep = ""), row.names=FALSE)
write.csv(Tpa_out$df_filt, file=paste("Tpa_MFcov_filt_", cutoff_len,  ".csv", sep = ""), row.names=FALSE)
write.csv(Tps_out$df_filt, file=paste("Tps_MFcov_filt_", cutoff_len,  ".csv", sep = ""), row.names=FALSE)

####### plots

png(filename = paste("Allsp_MFcov_filt_", cutoff_len, "_wAdjpeak2.png", sep = ""), width = 6, height = 8, units = "in", bg = "white", res = 300)
plot_grid(
Tbi_out$p2_adj,
Tce_out$p2_adj,
Tcm_out$p2_adj,
Tpa_out$p2_adj,
Tps_out$p2_adj,
ncol = 2, nrow = 3)
dev.off()
getwd() ## where has my plot gone....

pdf(paste("Allsp_MFcov_filt_", cutoff_len, "_wAdjpeak2.pdf", sep = ""), width = 6, height = 8, )
plot_grid(
Tbi_out$p2_adj,
Tce_out$p2_adj,
Tcm_out$p2_adj,
Tpa_out$p2_adj,
Tps_out$p2_adj,
ncol = 2, nrow = 3)
dev.off()
getwd() ## where has my plot gone....

########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), paste("Male_and_female_coverage.R_sessionInfo_len_", cutoff_len, ".txt"), sep = "")



































