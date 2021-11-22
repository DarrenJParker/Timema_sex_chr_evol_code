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

dat1_Tce_Nosil_raw <- read.table ("data/Nosil_mapping/Tce_Nosil_30aDR_cov.csv", header = T, sep = ',')
length(dat1_Tce_Nosil_raw[,1])

### filter 

dat1_Tce_Nosil <- dat1_Tce_Nosil_raw
dat1_Tce_Nosil_nolgNA <- subset(dat1_Tce_Nosil_raw, dat1_Tce_Nosil_raw$lg != "lgNA")
length(dat1_Tce_Nosil_nolgNA[,1])

subset(dat1_Tce_Nosil_raw, dat1_Tce_Nosil_raw$lg == "lg13")

# keep <- rowSums(dat1_Tce_Nosil_raw[,3:11] > 5) >= 6 ## remove contigs with 0 coverage in 2 or more libs
# dat1_Tce_Nosil <- dat1_Tce_Nosil_raw[keep,]
# length(dat1_Tce_Nosil[,1])

#dat1_Tce_Nosil <- subset(dat1_Tce_Nosil, dat1_Tce_Nosil$lg != "lgNA")

### output

dir.create("data/Nosil_mapping/plots")
setwd("data/Nosil_mapping/plots")

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
    	ggtitle(paste(samp, " | len >= ", sep = "")) +
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
	df1$LG <- as.character(df1$lg)

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
		geom_line(data=subset(df1,LG == "lg13"), color="red3", size = 1, stat="density") +
		geom_line(data=subset(df1,LG != "lgNA"), color="black", size = 1, stat="density") +	
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
		"lg13"  = "red3",
		"lgNA"  = "black"
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

### Tce

Tce_a1 <- plot_grid(
cov_plot("dat1_Tce_Nosil", "Tce_M_05_HM15", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce_Nosil", "Tce_M_06_HM16", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce_Nosil", "Tce_M_07_HM33", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce_Nosil", "Tce_M_08_HM61", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce_Nosil", "Tce_F_CC22B", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce_Nosil", "Tce_F_CC22C", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce_Nosil", "Tce_F_CC24B", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce_Nosil", "Tce_F_CC24C", min_cov_1, max_cov_1)$p1,
cov_plot("dat1_Tce_Nosil", "Tce_F_CC25B", min_cov_1, max_cov_1)$p1,
ncol = 2, nrow = 5)

# pdf(paste("Tce_cov_min=", min_cov_1, "_max=", max_cov_1, ".pdf", sep = ""), width = 	10, height = 15)
# plot_grid(Tce_a1, ncol = 1)
# dev.off()
# getwd() ## where has my plot gone....?



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



MF_cov_sum_nolgNA <- function(df, sp){
	
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

	out_list = list("df_filt" = df_filt)
	return(out_list)	
}


Tce_Nosil_out <- MF_cov_sum(dat1_Tce_Nosil, "Tce")
Tce_Nosil_out_nolgNA <- MF_cov_sum_nolgNA(dat1_Tce_Nosil_nolgNA, "Tce")

##########################################################################################################################################################################
### PLOT

head(Tce_Nosil_out$df_filt )

tapply(Tce_Nosil_out$df_filt$M_F_mode, Tce_Nosil_out$df_filt$lg,  quantile, probs=0.95)
tapply(Tce_Nosil_out$df_filt$M_F_mode, Tce_Nosil_out$df_filt$lg,  quantile, probs=0.5)

Tce_Nosil_out_nolgNA$df_filt$lg_ord <- ordered(Tce_Nosil_out_nolgNA$df_filt$lg, levels = c("lg1", "lg2", "lg3", "lg4", "lg5", "lg6", "lg7", "lg8", "lg9", "lg10", "lg11", "lg12", "lg13"))

Nosil_MFcov_nolgNA <- ggplot(Tce_Nosil_out_nolgNA$df_filt, aes(lg_ord, M_F_mode)) + 
	theme_bw() +
	geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), aes(fill = factor(lg_ord)), stroke = 1, dotsize=2 , binwidth = 0.01) + 
	coord_cartesian(ylim=c(-1.5,1.5)) +
	ylab ("log2(Male cov / Female cov)") +
	xlab ("Linkage group") + 
	scale_fill_manual(values=c("grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "grey30", "red2")) + theme(legend.position = "none")

pdf("Nosil_MFcov_nolgNA.pdf", width = 	5, height = 7.5)
Nosil_MFcov_nolgNA 
dev.off()
getwd() ## where has my plot gone....?

write.csv(
tapply(Tce_Nosil_out_nolgNA$df_filt$M_F_mode, Tce_Nosil_out_nolgNA$df_filt$lg, median), "Tce_Nosil_out_nolgNA_meds.csv"
)


########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), paste("Male_and_female_coverage.R_sessionInfo_len_", cutoff_len, ".txt"), sep = "")
































