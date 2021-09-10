## libs

library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library("SuperExactTest")


### read in data

setwd("data/counts")

Tbi_dat_1000 <- read.csv("Tbi_1000_chr_info_and_counts.csv")
Tce_dat_1000 <- read.csv("Tce_1000_chr_info_and_counts.csv")
Tcm_dat_1000 <- read.csv("Tcm_1000_chr_info_and_counts.csv")
Tpa_dat_1000 <- read.csv("Tpa_1000_chr_info_and_counts.csv")
Tps_dat_1000 <- read.csv("Tps_1000_chr_info_and_counts.csv")

Tbi_dat_5000 <- read.csv("Tbi_5000_chr_info_and_counts.csv")
Tce_dat_5000 <- read.csv("Tce_5000_chr_info_and_counts.csv")
Tcm_dat_5000 <- read.csv("Tcm_5000_chr_info_and_counts.csv")
Tpa_dat_5000 <- read.csv("Tpa_5000_chr_info_and_counts.csv")
Tps_dat_5000 <- read.csv("Tps_5000_chr_info_and_counts.csv")


### functions

X_hogs <- function(df){
	df_HOGs <- subset(df, df$HOG != "NA")
	N_HOGs <- length(df_HOGs[,1])
	X_HOG_hard <- subset(df_HOGs, df_HOGs$chr_hard == "X")$HOG
	X_HOG_soft <- subset(df_HOGs, df_HOGs$chr_soft == "X")$HOG
	
	print("N X orth (soft):")
	print(length(X_HOG_soft))
	
	print("N X orth (hard):")
	print(length(X_HOG_hard))
	
	out_list = list("X_HOG_hard" = X_HOG_hard, "X_HOG_soft" = X_HOG_soft, "N_HOGs" = N_HOGs)
	return(out_list)
}

five_sp_DE_venn <- function(Tbi,Tce,Tcm,Tpa,Tps,title){
	venny.plot <- venn.diagram(
	list("Tbi" = Tbi, "Tce" = Tce, "Tcm" = Tcm, "Tpa" = Tpa, "Tps" = Tps ), filename = NULL,
                            fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                            margin = 0.6, cat.dist = 0.23, main = title, main.pos = c(0.5,0.8), main.cex = 2, main.fontface = "bold", cat.cex = 2)
	return(venny.plot)
}


#### output

dir.create("../output/X_orths")
setwd("../output/X_orths")

###### make a venn

X_HOG_venn_hard_1000 <- five_sp_DE_venn(X_hogs(Tbi_dat_1000)$X_HOG_hard, X_hogs(Tce_dat_1000)$X_HOG_hard, X_hogs(Tcm_dat_1000)$X_HOG_hard, X_hogs(Tpa_dat_1000)$X_HOG_hard, X_hogs(Tps_dat_1000)$X_HOG_hard , "X-linked genes")
X_HOG_venn_soft_1000 <- five_sp_DE_venn(X_hogs(Tbi_dat_1000)$X_HOG_soft, X_hogs(Tce_dat_1000)$X_HOG_soft, X_hogs(Tcm_dat_1000)$X_HOG_soft, X_hogs(Tpa_dat_1000)$X_HOG_soft, X_hogs(Tps_dat_1000)$X_HOG_soft , "X-linked genes")
X_HOG_venn_hard_5000 <- five_sp_DE_venn(X_hogs(Tbi_dat_5000)$X_HOG_hard, X_hogs(Tce_dat_5000)$X_HOG_hard, X_hogs(Tcm_dat_5000)$X_HOG_hard, X_hogs(Tpa_dat_5000)$X_HOG_hard, X_hogs(Tps_dat_5000)$X_HOG_hard , "X-linked genes")
X_HOG_venn_soft_5000 <- five_sp_DE_venn(X_hogs(Tbi_dat_5000)$X_HOG_soft, X_hogs(Tce_dat_5000)$X_HOG_soft, X_hogs(Tcm_dat_5000)$X_HOG_soft, X_hogs(Tpa_dat_5000)$X_HOG_soft, X_hogs(Tps_dat_5000)$X_HOG_soft , "X-linked genes")


# grid.arrange(gTree(children=X_HOG_venn_hard_5000),ncol = 1)
# grid.arrange(gTree(children=X_HOG_venn_soft_5000),ncol = 1)

ggsave(file="X_HOG_venn_hard_1000.pdf", X_HOG_venn_hard_1000)
ggsave(file="X_HOG_venn_hard_5000.pdf", X_HOG_venn_hard_5000)
ggsave(file="X_HOG_venn_soft_1000.pdf", X_HOG_venn_soft_1000)
ggsave(file="X_HOG_venn_soft_5000.pdf", X_HOG_venn_soft_5000)


### Super exact test 

SET_wFDR <- function(list_of_lists, total_N ){
	res=supertest(list_of_lists, n= total_N)
	res_out <- as.data.frame(cbind(
		summary(res)$Table$Intersections,
		summary(res)$Table$Degree ,
		summary(res)$Table$Observed.Overlap,
		summary(res)$Table$Expected.Overlap ,
		summary(res)$Table$FE,
		summary(res)$Table$P.value))

	colnames(res_out) <- c("Intersections","Degree" , "Observed.Overlap","Expected.Overlap" ,"FE","P.value")
	res_out$P.value <- as.numeric(res_out$P.value)
	res_out$FDR     <- p.adjust(res_out$P.value)

	return(res_out)	
}


total_N_HOGs_1000 <- X_hogs(Tbi_dat_1000)$N_HOGs
l_dat_hard_1000 <- list("Tbi" = X_hogs(Tbi_dat_1000)$X_HOG_hard, "Tce" = X_hogs(Tce_dat_1000)$X_HOG_hard, "Tcm" = X_hogs(Tcm_dat_1000)$X_HOG_hard, "Tpa" = X_hogs(Tpa_dat_1000)$X_HOG_hard, "Tps" = X_hogs(Tps_dat_1000)$X_HOG_hard )
l_dat_soft_1000 <- list("Tbi" = X_hogs(Tbi_dat_1000)$X_HOG_soft, "Tce" = X_hogs(Tce_dat_1000)$X_HOG_soft, "Tcm" = X_hogs(Tcm_dat_1000)$X_HOG_soft, "Tpa" = X_hogs(Tpa_dat_1000)$X_HOG_soft, "Tps" = X_hogs(Tps_dat_1000)$X_HOG_soft )


write.csv(SET_wFDR(l_dat_hard_1000, total_N_HOGs_1000), file="SET_orths_X_summary_table_hard_1000.csv", row.names=FALSE)
write.csv(SET_wFDR(l_dat_hard_1000, total_N_HOGs_1000), file="SET_orths_X_summary_table_hard_5000.csv", row.names=FALSE)
write.csv(SET_wFDR(l_dat_soft_1000, total_N_HOGs_1000), file="SET_orths_X_summary_table_soft_1000.csv", row.names=FALSE)
write.csv(SET_wFDR(l_dat_soft_1000, total_N_HOGs_1000), file="SET_orths_X_summary_table_soft_5000.csv", row.names=FALSE)



### X in all species


Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

X_soft_HOG_list <-  Intersect(l_dat_soft_1000)
write.table(X_soft_HOG_list, file="X_soft_HOG_list.tsv", row.names = F, col.names = F, sep = "\t")
