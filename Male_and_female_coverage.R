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





#######################################################################################################################################################################
########## Coverage and linkage groups


## get contigs assingned to a single linkage group
dat1_Tbi_inLG <- subset(dat1_Tbi, dat1_Tbi$LG != "NA" & dat1_Tbi$multi_linkage_groups == "NO")
dat1_Tce_inLG <- subset(dat1_Tce, dat1_Tce$LG != "NA" & dat1_Tce$multi_linkage_groups == "NO")
dat1_Tcm_inLG <- subset(dat1_Tcm, dat1_Tcm$LG != "NA" & dat1_Tcm$multi_linkage_groups == "NO")
dat1_Tpa_inLG <- subset(dat1_Tpa, dat1_Tpa$LG != "NA" & dat1_Tpa$multi_linkage_groups == "NO")
dat1_Tps_inLG <- subset(dat1_Tps, dat1_Tps$LG != "NA" & dat1_Tps$multi_linkage_groups == "NO")


### Plot by linkage group

### Tbi

Tbi_a1 <- plot_grid(
cov_plot_LG("dat1_Tbi", "Tbi_M_13_Tbi")$p1,
cov_plot_LG("dat1_Tbi", "Tbi_M_14_Tbi")$p1,
cov_plot_LG("dat1_Tbi", "Tbi_M_15_Tbi")$p1,
cov_plot_LG("dat1_Tbi", "Tbi_M_16_Tbi")$p1,
cov_plot_LG("dat1_Tbi", "Tbi_F_CC86B")$p1,
cov_plot_LG("dat1_Tbi", "Tbi_F_CC86C")$p1,
cov_plot_LG("dat1_Tbi", "Tbi_F_CC87B")$p1,
cov_plot_LG("dat1_Tbi", "Tbi_F_CC87C")$p1,
cov_plot_LG("dat1_Tbi", "Tbi_F_CC88B")$p1,
ncol = 2, nrow = 5)

Tbi_a2 <- plot_grid(
cov_plot_LG("dat1_Tbi", "Tbi_M_13_Tbi")$legend, 
ncol = 1)

pdf(paste("Tbi_cov_LG",".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tbi_a1, Tbi_a2, ncol = 2, rel_widths = c(4,1))
dev.off()
getwd() ## where has my plot gone....?

 
### Tce

Tce_a1 <- plot_grid(
cov_plot_LG("dat1_Tce", "Tce_M_05_HM15")$p1,
cov_plot_LG("dat1_Tce", "Tce_M_06_HM16")$p1,
cov_plot_LG("dat1_Tce", "Tce_M_07_HM33")$p1,
cov_plot_LG("dat1_Tce", "Tce_M_08_HM61")$p1,
cov_plot_LG("dat1_Tce", "Tce_F_CC22B")$p1,
cov_plot_LG("dat1_Tce", "Tce_F_CC22C")$p1,
cov_plot_LG("dat1_Tce", "Tce_F_CC24B")$p1,
cov_plot_LG("dat1_Tce", "Tce_F_CC24C")$p1,
cov_plot_LG("dat1_Tce", "Tce_F_CC25B")$p1,
ncol = 2, nrow = 5)

Tce_a2 <- plot_grid(
cov_plot_LG("dat1_Tce", "Tce_F_CC22B")$legend, 
ncol = 1)

pdf(paste("Tce_cov_LG", ".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tce_a1, Tce_a2, ncol = 2, rel_widths = c(4,1))
dev.off()
getwd() ## where has my plot gone....?

### Tcm

Tcm_a1 <- plot_grid(
cov_plot_LG("dat1_Tcm", "Tcm_M_01_HM148")$p1,
cov_plot_LG("dat1_Tcm", "Tcm_M_02_HM149")$p1,
cov_plot_LG("dat1_Tcm", "Tcm_M_03_HM150")$p1,
cov_plot_LG("dat1_Tcm", "Tcm_M_04_HM151")$p1,
cov_plot_LG("dat1_Tcm", "Tcm_F_HM217")$p1,
cov_plot_LG("dat1_Tcm", "Tcm_F_HM218")$p1,
cov_plot_LG("dat1_Tcm", "Tcm_F_HM219")$p1,
cov_plot_LG("dat1_Tcm", "Tcm_F_HM220")$p1,
cov_plot_LG("dat1_Tcm", "Tcm_F_HM221")$p1,
ncol = 2, nrow = 5)

Tcm_a2 <- plot_grid(
cov_plot_LG("dat1_Tcm", "Tcm_F_HM221")$legend, 
ncol = 1)

pdf(paste("Tcm_cov_LG", ".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tcm_a1, Tcm_a2, ncol = 2, rel_widths = c(4,1))
dev.off()
getwd() ## where has my plot gone....?

### Tpa


Tpa_a1 <- plot_grid(
cov_plot_LG("dat1_Tpa", "Tpa_M_09_Tpa")$p1,
cov_plot_LG("dat1_Tpa", "Tpa_M_10_Tpa")$p1,
cov_plot_LG("dat1_Tpa", "Tpa_M_11_Tpa")$p1,
cov_plot_LG("dat1_Tpa", "Tpa_M_12_Tpa")$p1,
cov_plot_LG("dat1_Tpa", "Tpa_F_H54")$p1,
#cov_plot_LG("dat1_Tpa", "Tpa_F_H56")$p1,
cov_plot_LG("dat1_Tpa", "Tpa_F_PA_CD")$p1,
cov_plot_LG("dat1_Tpa", "Tpa_F_PA_E")$p1,
cov_plot_LG("dat1_Tpa", "Tpa_F_Pa_AB")$p1,
ncol = 2, nrow = 5)

Tpa_a2 <- plot_grid(
cov_plot_LG("dat1_Tpa", "Tpa_F_Pa_AB")$legend, 
ncol = 1)

pdf(paste("Tpa_cov_LG",".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tpa_a1, Tpa_a2, ncol = 2, rel_widths = c(4,1))
dev.off()
getwd() ## where has my plot gone....?



### Tps


Tps_a1 <- plot_grid(
cov_plot_LG("dat1_Tps", "Tps_M_17_HM99")$p1,
cov_plot_LG("dat1_Tps", "Tps_M_18_HM100")$p1,
cov_plot_LG("dat1_Tps", "Tps_M_19_HM101")$p1,
cov_plot_LG("dat1_Tps", "Tps_M_20_15255")$p1,
#cov_plot_LG("dat1_Tps", "Tps_F_ReSeq_Ps08")$p1,
#cov_plot_LG("dat1_Tps", "Tps_F_ReSeq_Ps12")$p1,
cov_plot_LG("dat1_Tps", "Tps_F_ReSeq_Ps14")$p1,
cov_plot_LG("dat1_Tps", "Tps_F_ReSeq_Ps16")$p1,
cov_plot_LG("dat1_Tps", "Tps_F_ReSeq_Ps18")$p1,
ncol = 2, nrow = 5)

Tps_a2 <- plot_grid(
cov_plot_LG("dat1_Tps", "Tps_M_20_15255")$legend, 
ncol = 1)

pdf(paste("Tps_cov_LG", ".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tps_a1, Tps_a2, ncol = 2, rel_widths = c(4,1))
dev.off()
getwd() ## where has my plot gone....?













#######################################################################################################################################################################
########## heterozy on the X

calc_prop_het <- function(df, samp_names) {
	df_out <- df
	out_head_1 <- colnames(df_out )
	out_head_2 <- c()
	for(s in samp_names){
		out_head_2 <- c(out_head_2, paste("Phet_", s, sep = ""),  paste("TcovS_", s, sep = ""))
		print(s)
		Phet      <- eval(parse(text=paste('df_out$hetero_mlest_',s, sep=''))) / (eval(parse(text=paste('df_out$homo_mlest_',s, sep=''))) + eval(parse(text=paste('df_out$hetero_mlest_',s, sep=''))))
		TcovSites <- eval(parse(text=paste('df_out$homo_mlest_',s, sep=''))) + eval(parse(text=paste('df_out$hetero_mlest_',s, sep='')))
		df_out <- cbind(df_out, Phet, TcovSites)
	}
	
	colnames(df_out) <- c(out_head_1, out_head_2)
	print(out_head_2)
	return(df_out)
}


Tbi_df_filt <- Tbi_out$df_filt
Tce_df_filt <- Tce_out$df_filt
Tcm_df_filt <- Tcm_out$df_filt
Tpa_df_filt <- Tpa_out$df_filt
Tps_df_filt <- Tps_out$df_filt

Tbi_samp_names <- c("Tbi_F_CC86B", "Tbi_F_CC86C", "Tbi_F_CC87B", "Tbi_F_CC87C", "Tbi_F_CC88B", "Tbi_M_13_Tbi", "Tbi_M_14_Tbi", "Tbi_M_15_Tbi", "Tbi_M_16_Tbi")
Tce_samp_names <- c("Tce_F_CC22B", "Tce_F_CC22C", "Tce_F_CC24B", "Tce_F_CC24C", "Tce_F_CC25B", "Tce_M_05_HM15", "Tce_M_06_HM16", "Tce_M_07_HM33", "Tce_M_08_HM61")
Tcm_samp_names <- c("Tcm_F_HM217", "Tcm_F_HM218", "Tcm_F_HM219", "Tcm_F_HM220", "Tcm_F_HM221", "Tcm_M_01_HM148", "Tcm_M_02_HM149", "Tcm_M_03_HM150", "Tcm_M_04_HM151")
Tpa_samp_names <- c("Tpa_F_H54", "Tpa_F_PA_CD", "Tpa_F_PA_E", "Tpa_F_Pa_AB", "Tpa_M_09_Tpa", "Tpa_M_10_Tpa", "Tpa_M_11_Tpa", "Tpa_M_12_Tpa")
Tps_samp_names <- c("Tps_F_ReSeq_Ps14", "Tps_F_ReSeq_Ps16", "Tps_F_ReSeq_Ps18", "Tps_M_17_HM99", "Tps_M_18_HM100", "Tps_M_19_HM101", "Tps_M_20_15255")

Tbi_F_samp_names <- c("Tbi_F_CC86B", "Tbi_F_CC86C", "Tbi_F_CC87B", "Tbi_F_CC87C", "Tbi_F_CC88B")
Tce_F_samp_names <- c("Tce_F_CC22B", "Tce_F_CC22C", "Tce_F_CC24B", "Tce_F_CC24C", "Tce_F_CC25B")
Tcm_F_samp_names <- c("Tcm_F_HM217", "Tcm_F_HM218", "Tcm_F_HM219", "Tcm_F_HM220", "Tcm_F_HM221")
Tpa_F_samp_names <- c("Tpa_F_H54", "Tpa_F_PA_CD", "Tpa_F_PA_E", "Tpa_F_Pa_AB")
Tps_F_samp_names <- c("Tps_F_ReSeq_Ps14", "Tps_F_ReSeq_Ps16", "Tps_F_ReSeq_Ps18")

Tbi_df_filt_het <- calc_prop_het(Tbi_df_filt, Tbi_samp_names) 
Tce_df_filt_het <- calc_prop_het(Tce_df_filt, Tce_samp_names) 
Tcm_df_filt_het <- calc_prop_het(Tcm_df_filt, Tcm_samp_names) 
Tpa_df_filt_het <- calc_prop_het(Tpa_df_filt, Tpa_samp_names) 
Tps_df_filt_het <- calc_prop_het(Tps_df_filt, Tps_samp_names) 

head(Tbi_df_filt_het )

X_A_Phet <-  function(df, samp_names){
	
	### hard class	
	
	df_hard_X  <- subset(df, df$class_hard == "X")
	df_hard_A  <- subset(df, df$class_hard == "A")		
	
	out_df_hard <- c()
	out_df_XA_ratio_hard <- c()
	for(s in samp_names){
		print(s)
		wt_med_X <- weighted.median(eval(parse(text=paste('df_hard_X$Phet_',s, sep=''))), eval(parse(text=paste('df_hard_X$TcovS_',s, sep=''))))
		wt_med_A <- weighted.median(eval(parse(text=paste('df_hard_A$Phet_',s, sep=''))), eval(parse(text=paste('df_hard_A$TcovS_',s, sep=''))))
	
		out_line <- c(s, "X", wt_med_X)
		out_df_hard <- rbind(out_df_hard, out_line)
		out_line <- c(s, "A", wt_med_A)
		out_df_hard <- rbind(out_df_hard, out_line)		
		
		X_A_ratio <- wt_med_X / wt_med_A
		out_line <- c(s, "XA_ratio_hard", X_A_ratio)
		out_df_XA_ratio_hard <- 	rbind(out_df_XA_ratio_hard, out_line)		
	}	

	out_df_hard <- as.data.frame(out_df_hard)
	colnames(out_df_hard) <- c("samp", "class", "wt_med_Phet")
	out_df_hard$wt_med_Phet <- as.numeric(as.character(out_df_hard$wt_med_Phet))

	out_df_XA_ratio_hard <- as.data.frame(out_df_XA_ratio_hard)
	colnames(out_df_XA_ratio_hard) <- c("samp", "class", "wt_med_Phet_XA_ratio")
	out_df_XA_ratio_hard$wt_med_Phet_XA_ratio <- as.numeric(as.character(out_df_XA_ratio_hard$wt_med_Phet_XA_ratio))	
	
	### soft class	
	
	df_soft_X  <- subset(df, df$class_soft == "X")
	df_soft_A  <- subset(df, df$class_soft == "A")		
	
	out_df_soft <- c()
	out_df_XA_ratio_soft <- c()
	for(s in samp_names){
		print(s)
		wt_med_X <- weighted.median(eval(parse(text=paste('df_soft_X$Phet_',s, sep=''))), eval(parse(text=paste('df_soft_X$TcovS_',s, sep=''))))
		wt_med_A <- weighted.median(eval(parse(text=paste('df_soft_A$Phet_',s, sep=''))), eval(parse(text=paste('df_soft_A$TcovS_',s, sep=''))))
	
		out_line <- c(s, "X", wt_med_X)
		out_df_soft <- rbind(out_df_soft, out_line)
		out_line <- c(s, "A", wt_med_A)
		out_df_soft <- rbind(out_df_soft, out_line)	
		
		X_A_ratio <- wt_med_X / wt_med_A
		out_line <- c(s, "XA_ratio_soft", X_A_ratio)
		out_df_XA_ratio_soft <- 	rbind(out_df_XA_ratio_soft, out_line)						
	}	

	out_df_soft <- as.data.frame(out_df_soft)
	colnames(out_df_soft) <- c("samp", "class", "wt_med_Phet")
	out_df_soft$wt_med_Phet <- as.numeric(as.character(out_df_soft$wt_med_Phet))

	out_df_XA_ratio_soft <- as.data.frame(out_df_XA_ratio_soft)
	colnames(out_df_XA_ratio_soft) <- c("samp", "class", "wt_med_Phet_XA_ratio")
	out_df_XA_ratio_soft$wt_med_Phet_XA_ratio <- as.numeric(as.character(out_df_XA_ratio_soft$wt_med_Phet_XA_ratio))	
	

	### split in 3 cats
	
	df_XX  <- subset(df, df$to_fill == "XX")
	df_XA  <- subset(df, df$to_fill == "XA")
	df_AA  <- subset(df, df$to_fill == "AA")		
	
	print(c(length(df_XX[,1]), length(df_XA[,1]), length(df_AA[,1])))
	
	out_df3cat <- c()
	for(s in samp_names){
		print(s)
		wt_med_XX <- weighted.median(eval(parse(text=paste('df_XX$Phet_',s, sep=''))), eval(parse(text=paste('df_XX$TcovS_',s, sep=''))))
		wt_med_XA <- weighted.median(eval(parse(text=paste('df_XA$Phet_',s, sep=''))), eval(parse(text=paste('df_XA$TcovS_',s, sep=''))))		
		wt_med_AA <- weighted.median(eval(parse(text=paste('df_AA$Phet_',s, sep=''))), eval(parse(text=paste('df_AA$TcovS_',s, sep=''))))
	
		out_line <- c(s, "XX", wt_med_XX)
		out_df3cat <- rbind(out_df3cat, out_line)
		out_line <- c(s, "XA", wt_med_XA)
		out_df3cat <- rbind(out_df3cat, out_line)		
		out_line <- c(s, "AA", wt_med_AA)
		out_df3cat <- rbind(out_df3cat, out_line)					
	}	

	out_df3cat <- as.data.frame(out_df3cat)
	colnames(out_df3cat) <- c("samp", "class", "wt_med_Phet")
	out_df3cat$wt_med_Phet <- as.numeric(as.character(out_df3cat$wt_med_Phet))

	
	out_list = list("out_df_soft" = out_df_soft, "out_df_hard" = out_df_hard, "out_df3cat" = out_df3cat, "out_df_XA_ratio_hard" = 	out_df_XA_ratio_hard,  "out_df_XA_ratio_soft" = 	out_df_XA_ratio_soft)
	return(out_list)	
	
}

head(Tbi_df_filt_het)

Tbi_Phet <- X_A_Phet(Tbi_df_filt_het, Tbi_samp_names)
Tce_Phet <- X_A_Phet(Tce_df_filt_het, Tce_samp_names)
Tcm_Phet <- X_A_Phet(Tcm_df_filt_het, Tcm_samp_names)
Tpa_Phet <- X_A_Phet(Tpa_df_filt_het, Tpa_samp_names)
Tps_Phet <- X_A_Phet(Tps_df_filt_het, Tps_samp_names)

#### plot


plot_phet_1 <- function(df,tit_text){
	p1 <- ggplot(df, aes(samp, wt_med_Phet, fill = class)) + 
		geom_bar(position="dodge",stat="identity") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		xlab ("Sample") + 
		ylab ("Heterozygosity") +
		scale_fill_manual(values=c("darkgrey", "red3")) + ggtitle(tit_text) 
	return(p1)
}	

plot_phet_2 <- function(df,tit_text){
  max_y <- max(df$wt_med_Phet)
  
  df$sex <- str_split_fixed(df$samp, "_", 3)[,2]
  df$class_2 <- paste(df$class, df$sex , sep = "")
  print(df)
  
	p1 <- ggplot(df, aes(samp, wt_med_Phet, fill = class)) + 
		geom_bar(position="dodge",stat="identity") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		xlab ("Sample") + 
		ylab ("Heterozygosity") +  ylim(c(0, max_y * 1.2)) +
	  scale_fill_manual(values=c("darkgrey", "darkorange2")) + ggtitle(tit_text)  
		#scale_fill_manual(values=c("darkgrey",  "black", "darkorange2", "yellow")) + ggtitle(tit_text)  
	return(p1)
}	

plot_phet_2(Tbi_Phet$out_df_soft, "Tbi")

plot_phet_3 <- function(df,tit_text){
	p1 <- ggplot(df, aes(samp, wt_med_Phet, fill = class)) + 
		geom_bar(position="dodge",stat="identity") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		xlab ("Sample") + 
		ylab ("Heterozygosity") +
		scale_fill_manual(values=c("darkgrey", "darkorange2", "red3")) + ggtitle(tit_text) 
	return(p1)
}	


pdf(paste("Phet_hard_class_", cutoff_len, ".pdf", sep = ""), width = 9, height = 12)
plot_grid(
plot_phet_1(Tbi_Phet$out_df_hard, "Tbi"),
plot_phet_1(Tce_Phet$out_df_hard, "Tce"),
plot_phet_1(Tcm_Phet$out_df_hard, "Tcm"),
plot_phet_1(Tpa_Phet$out_df_hard, "Tpa"),
plot_phet_1(Tps_Phet$out_df_hard, "Tps"), ncol = 2 
)
dev.off()
getwd() ## where has my plot gone....?


pdf(paste("Phet_soft_class_", cutoff_len, ".pdf", sep = ""), width = 	9, height = 12)
plot_grid(
plot_phet_2(Tbi_Phet$out_df_soft, "Tbi"),
plot_phet_2(Tce_Phet$out_df_soft, "Tce"),
plot_phet_2(Tcm_Phet$out_df_soft, "Tcm"),
plot_phet_2(Tpa_Phet$out_df_soft, "Tpa"),
plot_phet_2(Tps_Phet$out_df_soft, "Tps"), ncol = 2 
)
dev.off()
getwd() ## where has my plot gone....?



# plot_phet_3(Tpa_Phet$out_df3cat, "Tpa")


###################################################################################################
### XA ratio

XA_ratio_all <- as.data.frame(rbind(
Tbi_Phet$out_df_XA_ratio_hard,
Tce_Phet$out_df_XA_ratio_hard,
Tcm_Phet$out_df_XA_ratio_hard,
Tpa_Phet$out_df_XA_ratio_hard,
Tps_Phet$out_df_XA_ratio_hard,
Tbi_Phet$out_df_XA_ratio_soft,
Tce_Phet$out_df_XA_ratio_soft,
Tcm_Phet$out_df_XA_ratio_soft,
Tpa_Phet$out_df_XA_ratio_soft,
Tps_Phet$out_df_XA_ratio_soft))


str(XA_ratio_all)



XA_ratio_all$sex <- str_split_fixed(XA_ratio_all$samp, "_", 3)[,2]
XA_ratio_all$sp <- str_split_fixed(XA_ratio_all$samp, "_", 3)[,1]
XA_ratio_F_hard <- subset(XA_ratio_all, XA_ratio_all$sex == "F" & XA_ratio_all$class == "XA_ratio_hard")
XA_ratio_F_soft <- subset(XA_ratio_all, XA_ratio_all$sex == "F" & XA_ratio_all$class == "XA_ratio_soft")


XA_het_ratio_hard <- ggplot(XA_ratio_F_hard, aes(samp, wt_med_Phet_XA_ratio, fill = sp)) + 
		geom_bar(position="dodge",stat="identity", colour="black") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		ylim(c(0,1)) + geom_hline(yintercept= 0.75,  linetype="dashed") + 
		
		xlab ("Sample") + 
		ylab ("X:A Heterozygosity") + scale_fill_brewer(palette = "Set3") + 
		scale_y_continuous(expand = c(0,0), limits = c(0,1)) 
		

XA_het_ratio_soft <- ggplot(XA_ratio_F_soft, aes(samp, wt_med_Phet_XA_ratio, fill = sp)) + 
		geom_bar(position="dodge",stat="identity", colour="black") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		ylim(c(0,1)) + geom_hline(yintercept= 0.75,  linetype="dashed") + 
		
		xlab ("Sample") + 
		ylab ("X:A Heterozygosity") + scale_fill_brewer(palette = "Set3") + 
		scale_y_continuous(expand = c(0,0), limits = c(0,1)) 


pdf(paste("XA_Phet_ratio_hard_class_", cutoff_len, ".pdf", sep = ""), width = 	7, height = 10)
XA_het_ratio_hard 
dev.off()
getwd() ## where has my plot gone....?


pdf(paste("XA_Phet_ratio_soft_class_", cutoff_len, ".pdf", sep = ""), width = 	7, height = 10)
XA_het_ratio_soft 
dev.off()
getwd() ## where has my plot gone....?

min(XA_ratio_F_hard$wt_med_Phet_XA_ratio)
max(XA_ratio_F_hard$wt_med_Phet_XA_ratio)

##########################################################################################################
### by LG  


head(Tbi_df_filt_het)

LG_Phet <-  function(df, samp_names){

	df_filt <- subset(df, df$multi_linkage_groups == "NO")
	print(head(df_filt))		
	### LG 	

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
	
	print(head(df_filt_lg6 ))
	
		
	out_df_LG <- c()
	for(s in samp_names){
		print(s)
		wt_med_lg1  <- weighted.median(eval(parse(text=paste('df_filt_lg1$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lg1$TcovS_',s, sep=''))))
		wt_med_lg2  <- weighted.median(eval(parse(text=paste('df_filt_lg2$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lg2$TcovS_',s, sep=''))))
		wt_med_lg3  <- weighted.median(eval(parse(text=paste('df_filt_lg3$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lg3$TcovS_',s, sep=''))))
		wt_med_lg4  <- weighted.median(eval(parse(text=paste('df_filt_lg4$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lg4$TcovS_',s, sep=''))))
		wt_med_lg5  <- weighted.median(eval(parse(text=paste('df_filt_lg5$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lg5$TcovS_',s, sep=''))))
		wt_med_lg6  <- weighted.median(eval(parse(text=paste('df_filt_lg6$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lg6$TcovS_',s, sep=''))))
		wt_med_lg7  <- weighted.median(eval(parse(text=paste('df_filt_lg7$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lg7$TcovS_',s, sep=''))))
		wt_med_lg8  <- weighted.median(eval(parse(text=paste('df_filt_lg8$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lg8$TcovS_',s, sep=''))))
		wt_med_lg9  <- weighted.median(eval(parse(text=paste('df_filt_lg9$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lg9$TcovS_',s, sep=''))))
		wt_med_lg10 <- weighted.median(eval(parse(text=paste('df_filt_lg10$Phet_',s, sep=''))), eval(parse(text=paste('df_filt_lg10$TcovS_',s, sep=''))))
		wt_med_lg11 <- weighted.median(eval(parse(text=paste('df_filt_lg11$Phet_',s, sep=''))), eval(parse(text=paste('df_filt_lg11$TcovS_',s, sep=''))))
		wt_med_lg12 <- weighted.median(eval(parse(text=paste('df_filt_lg12$Phet_',s, sep=''))), eval(parse(text=paste('df_filt_lg12$TcovS_',s, sep=''))))
		wt_med_lgX  <- weighted.median(eval(parse(text=paste('df_filt_lgX$Phet_',s, sep=''))),  eval(parse(text=paste('df_filt_lgX$TcovS_',s, sep=''))))
		
		out_line <- c(s, "lg1", wt_med_lg1)
		out_df_LG <- rbind(out_df_LG, out_line)	
		out_line <- c(s, "lg2", wt_med_lg2)
		out_df_LG <- rbind(out_df_LG, out_line)
		out_line <- c(s, "lg3", wt_med_lg3)
		out_df_LG <- rbind(out_df_LG, out_line)				
		out_line <- c(s, "lg4", wt_med_lg4)
		out_df_LG <- rbind(out_df_LG, out_line)	
		out_line <- c(s, "lg5", wt_med_lg5)
		out_df_LG <- rbind(out_df_LG, out_line)
		out_line <- c(s, "lg6", wt_med_lg6)
		out_df_LG <- rbind(out_df_LG, out_line)	
		out_line <- c(s, "lg7", wt_med_lg7)
		out_df_LG <- rbind(out_df_LG, out_line)	
		out_line <- c(s, "lg8", wt_med_lg8)
		out_df_LG <- rbind(out_df_LG, out_line)
		out_line <- c(s, "lg9", wt_med_lg9)
		out_df_LG <- rbind(out_df_LG, out_line)	
		out_line <- c(s, "lg10", wt_med_lg10)
		out_df_LG <- rbind(out_df_LG, out_line)	
		out_line <- c(s, "lg11", wt_med_lg11)
		out_df_LG <- rbind(out_df_LG, out_line)
		out_line <- c(s, "lg12", wt_med_lg12)
		out_df_LG <- rbind(out_df_LG, out_line)
		out_line <- c(s, "lgX", wt_med_lgX)
		out_df_LG <- rbind(out_df_LG, out_line)			
	}	

	out_df_LG <- as.data.frame(out_df_LG)
	colnames(out_df_LG) <- c("samp", "LG", "wt_med_Phet")
	out_df_LG$wt_med_Phet <- as.numeric(as.character(out_df_LG$wt_med_Phet))
	out_df_LG$LG_ord <- ordered(out_df_LG$LG, levels=c("lg1","lg2","lg3","lg4","lg5","lg6","lg7","lg8","lg9","lg10","lg11","lg12","lgX")) 

	#out_list = list("out_df_soft" = out_df_soft, "out_df_hard" = out_df_hard, "out_df3cat" = out_df3cat)
	return(out_df_LG)	
	
}

Tbi_LG_Phet <- LG_Phet(Tbi_df_filt_het, Tbi_samp_names)
Tce_LG_Phet <- LG_Phet(Tce_df_filt_het, Tce_samp_names)
Tcm_LG_Phet <- LG_Phet(Tcm_df_filt_het, Tcm_samp_names)
Tpa_LG_Phet <- LG_Phet(Tpa_df_filt_het, Tpa_samp_names)
Tps_LG_Phet <- LG_Phet(Tps_df_filt_het, Tps_samp_names)

plot_phet_LG_1 <- function(df,tit_text){
	p1 <- ggplot(df, aes(samp, wt_med_Phet, fill = LG_ord)) + 
		geom_bar(position="dodge",stat="identity") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90)) +
		xlab ("Sample") + 
		ylab ("Heterozygosity") + 
		scale_fill_manual(values=c("#1B9E77", "#D95F02","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "yellow2","#666666","lightblue","royalblue2", "darkorchid", "red3")) + ggtitle(tit_text) 
	return(p1)
}	



pdf(paste("Phet_LG_", cutoff_len, ".pdf", sep = ""), width = 	9, height = 20)
plot_grid(
plot_phet_LG_1(Tbi_LG_Phet, "Tbi"),
plot_phet_LG_1(Tce_LG_Phet, "Tce"),
plot_phet_LG_1(Tcm_LG_Phet, "Tcm"),
plot_phet_LG_1(Tpa_LG_Phet, "Tpa"),
plot_phet_LG_1(Tps_LG_Phet, "Tps"),
ncol = 1)
dev.off()
getwd() ## where has my plot gone....?


## F only

Tbi_F_LG_Phet <- LG_Phet(Tbi_df_filt_het, Tbi_F_samp_names)
Tce_F_LG_Phet <- LG_Phet(Tce_df_filt_het, Tce_F_samp_names)
Tcm_F_LG_Phet <- LG_Phet(Tcm_df_filt_het, Tcm_F_samp_names)
Tpa_F_LG_Phet <- LG_Phet(Tpa_df_filt_het, Tpa_F_samp_names)
Tps_F_LG_Phet <- LG_Phet(Tps_df_filt_het, Tps_F_samp_names)

pdf(paste("Phet_LG_F_", cutoff_len, ".pdf", sep = ""), width = 	8, height = 20)
plot_grid(
  plot_phet_LG_1(Tbi_F_LG_Phet, "Tbi"),
  plot_phet_LG_1(Tce_F_LG_Phet, "Tce"),
  plot_phet_LG_1(Tcm_F_LG_Phet, "Tcm"),
  plot_phet_LG_1(Tpa_F_LG_Phet, "Tpa"),
  plot_phet_LG_1(Tps_F_LG_Phet, "Tps"),
  ncol = 1)
dev.off()
getwd() ## where has my plot gone....?




############################################# ############################################# ############################################# 
############################################# hist XA

### Hard to see much here
hist_min_len = 20000



Tbi_df_filt_het_l <- as.data.frame(cbind(
c(
Tbi_df_filt_het$Phet_Tbi_F_CC86B, 
Tbi_df_filt_het$Phet_Tbi_F_CC86C,
Tbi_df_filt_het$Phet_Tbi_F_CC87B,
Tbi_df_filt_het$Phet_Tbi_F_CC87C,
Tbi_df_filt_het$Phet_Tbi_F_CC88B,
Tbi_df_filt_het$Phet_Tbi_M_13_Tbi, 
Tbi_df_filt_het$Phet_Tbi_M_14_Tbi,
Tbi_df_filt_het$Phet_Tbi_M_15_Tbi,
Tbi_df_filt_het$Phet_Tbi_M_16_Tbi),
c(
rep("Tbi_F_CC86B", length(Tbi_df_filt_het[,1])),
rep("Tbi_F_CC86C", length(Tbi_df_filt_het[,1])),
rep("Tbi_F_CC87B", length(Tbi_df_filt_het[,1])),
rep("Tbi_F_CC87C", length(Tbi_df_filt_het[,1])),
rep("Tbi_F_CC88B", length(Tbi_df_filt_het[,1])),
rep("Tbi_M_13_Tbi", length(Tbi_df_filt_het[,1])),
rep("Tbi_M_14_Tbi", length(Tbi_df_filt_het[,1])),
rep("Tbi_M_15_Tbi", length(Tbi_df_filt_het[,1])),
rep("Tbi_M_16_Tbi", length(Tbi_df_filt_het[,1]))),
c(
Tbi_df_filt_het$length, 
Tbi_df_filt_het$length, 
Tbi_df_filt_het$length, 
Tbi_df_filt_het$length, 
Tbi_df_filt_het$length, 
Tbi_df_filt_het$length, 
Tbi_df_filt_het$length, 
Tbi_df_filt_het$length, 
Tbi_df_filt_het$length),
c(
Tbi_df_filt_het$class_hard, 
Tbi_df_filt_het$class_hard, 
Tbi_df_filt_het$class_hard, 
Tbi_df_filt_het$class_hard, 
Tbi_df_filt_het$class_hard, 
Tbi_df_filt_het$class_hard, 
Tbi_df_filt_het$class_hard, 
Tbi_df_filt_het$class_hard, 
Tbi_df_filt_het$class_hard)))

colnames(Tbi_df_filt_het_l) <- c("Phet", "samp", "length", "class_hard")
Tbi_df_filt_het_l$Phet <- as.numeric(Tbi_df_filt_het_l$Phet)
Tbi_df_filt_het_l$length <- as.numeric(Tbi_df_filt_het_l$length)



Tbi_df_filt_het_l_s <- subset(Tbi_df_filt_het_l, Tbi_df_filt_het_l$length >= hist_min_len)
head(Tbi_df_filt_het_l_s )

length(Tbi_df_filt_het_l_s[,1])

# ggplot(Tbi_df_filt_het_l_s, aes(log(Phet))) + 
	# theme_bw() +
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC86B"), color="darkgrey", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC86B"), color="red3", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC86C"), color="darkgrey", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC86C"), color="red3", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC87B"), color="darkgrey", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC87B"), color="red3", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC87C"), color="darkgrey", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC87C"), color="red3", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC88B"), color="darkgrey", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC88B"), color="red3", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_M_13_Tbi"), color="darkgrey", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_M_13_Tbi"), color="red3", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_M_14_Tbi"), color="darkgrey", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_M_14_Tbi"), color="red3", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_M_15_Tbi"), color="darkgrey", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_M_15_Tbi"), color="red3", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_M_16_Tbi"), color="darkgrey", size = 1, stat="density") + 
	# geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_M_16_Tbi"), color="red3", size = 1, stat="density") 



Tbi_Phet_hist_F <- ggplot(Tbi_df_filt_het_l_s, aes(log(Phet + 1))) + 
	theme_bw() +
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC86B"), color="darkgrey", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC86B"), color="red3", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC86C"), color="darkgrey", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC86C"), color="red3", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC87B"), color="darkgrey", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC87B"), color="red3", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC87C"), color="darkgrey", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC87C"), color="red3", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_F_CC88B"), color="darkgrey", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_F_CC88B"), color="red3", size = 1, stat="density") + xlim(c(-.0001, 0.01)) + ggtitle(paste("Tbi Females min len = ",hist_min_len))


Tbi_Phet_hist_M <- ggplot(Tbi_df_filt_het_l_s, aes(log10(Phet + 1))) + 
	theme_bw() +
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_M_13_Tbi"), color="darkgrey", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_M_13_Tbi"), color="red3", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_M_14_Tbi"), color="darkgrey", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_M_14_Tbi"), color="red3", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_M_15_Tbi"), color="darkgrey", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_M_15_Tbi"), color="red3", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "A" & samp == "Tbi_M_16_Tbi"), color="darkgrey", size = 1, stat="density") + 
	geom_line(data=subset(Tbi_df_filt_het_l_s,class_hard  == "X" & samp == "Tbi_M_16_Tbi"), color="red3", size = 1, stat="density") + xlim(c(-.0001, 0.01)) + ggtitle(paste("Tbi Males, min len = ",hist_min_len))


plot_grid(Tbi_Phet_hist_F, Tbi_Phet_hist_M, ncol = 1)

pdf(paste("Tbi_Phet_hist", hist_min_len, ".pdf", sep = ""), width = 8, height = 7)
plot_grid(Tbi_Phet_hist_F, Tbi_Phet_hist_M, ncol = 1)
dev.off()
getwd() ## where has my plot gone....?



########################################################################################################################################################################
####### output session info
print (sessionInfo())

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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-2  spatstat_1.64-1     rpart_4.1-15        nlme_3.1-151        spatstat.data_1.7-0 matrixStats_0.58.0  cowplot_1.1.1      
# [8] plyr_1.8.6          modeest_2.4.0       stringr_1.4.0       ggplot2_3.3.5      
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.6           pillar_1.4.7         compiler_4.0.3       tools_4.0.3          timeSeries_3062.100  digest_0.6.27       
# [7] goftest_1.2-2        lattice_0.20-41      lifecycle_0.2.0      tibble_3.0.6         gtable_0.3.0         stable_1.1.4        
# [13] clue_0.3-58          mgcv_1.8-33          pkgconfig_2.0.3      rlang_0.4.10         Matrix_1.3-2         rmutil_1.1.5        
# [19] statip_0.2.3         withr_2.4.1          dplyr_1.0.3          cluster_2.1.0        generics_0.1.0       vctrs_0.3.6         
# [25] spatstat.utils_2.0-0 tidyselect_1.1.0     glue_1.4.2           R6_2.5.0             spatial_7.3-13       polyclip_1.10-0     
# [31] farver_2.0.3         tensor_1.5           deldir_0.2-9         purrr_0.3.4          magrittr_2.0.1       splines_4.0.3       
# [37] fBasics_3042.89.1    scales_1.1.1         ellipsis_0.3.1       abind_1.4-5          stabledist_0.7-1     timeDate_3043.102   
# [43] colorspace_2.0-0     labeling_0.4.2       stringi_1.5.3        munsell_0.5.0        crayon_1.4.0     

writeLines(capture.output(sessionInfo()), paste("Male_and_female_coverage.R_sessionInfo_len_", cutoff_len, ".txt"), sep = "")



































