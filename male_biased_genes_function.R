# male_biased_genes_function.R

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
library(topGO)
library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library("SuperExactTest")
require(dplyr)


print (sessionInfo())

### read data

setwd("data/Exp_out")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### data

Tbi_SB_RT_FPKM_TTT <- read.table ("TTT_RT_Tbi_sex_bias_FPKM.csv", header = T, sep = ',')
Tce_SB_RT_FPKM_TTT <- read.table ("TTT_RT_Tce_sex_bias_FPKM.csv", header = T, sep = ',')
Tcm_SB_RT_FPKM_TTT <- read.table ("TTT_RT_Tcm_sex_bias_FPKM.csv", header = T, sep = ',')
Tpa_SB_RT_FPKM_TTT <- read.table ("TTT_RT_Tpa_sex_bias_FPKM.csv", header = T, sep = ',')
Tps_SB_RT_FPKM_TTT <- read.table ("TTT_RT_Tps_sex_bias_FPKM.csv", header = T, sep = ',')

#### change wd to output folder
dir.create("../male_biased_gene_function")
setwd("../male_biased_gene_function")

#### load annotation

geneID2GO_Tbi_Droso <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Droso/Tbi_b3v08.max_arth_b2g_droso_b2g.gff_droso_ontology_term_fortopgo.txt")
geneID2GO_Tce_Droso <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Droso/Tce_b3v08.max_arth_b2g_droso_b2g.gff_droso_ontology_term_fortopgo.txt")
geneID2GO_Tcm_Droso <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Droso/Tcm_b3v08.max_arth_b2g_droso_b2g.gff_droso_ontology_term_fortopgo.txt")
geneID2GO_Tpa_Droso <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Droso/Tpa_b3v08.max_arth_b2g_droso_b2g.gff_droso_ontology_term_fortopgo.txt")
geneID2GO_Tps_Droso <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Droso/Tps_b3v08.max_arth_b2g_droso_b2g.gff_droso_ontology_term_fortopgo.txt")

geneID2GO_Tbi_Arth <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Arth/Tbi_b3v08.max_arth_b2g_droso_b2g.gff_arth_ontology_term_fortopgo.txt")
geneID2GO_Tce_Arth <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Arth/Tce_b3v08.max_arth_b2g_droso_b2g.gff_arth_ontology_term_fortopgo.txt")
geneID2GO_Tcm_Arth <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Arth/Tcm_b3v08.max_arth_b2g_droso_b2g.gff_arth_ontology_term_fortopgo.txt")
geneID2GO_Tpa_Arth <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Arth/Tpa_b3v08.max_arth_b2g_droso_b2g.gff_arth_ontology_term_fortopgo.txt")
geneID2GO_Tps_Arth <- readMappings(file = "/Users/dparker/Documents/University/Lausanne/Timema_genomes/GO_terms/Arth/Tps_b3v08.max_arth_b2g_droso_b2g.gff_arth_ontology_term_fortopgo.txt")







#########################################################################################################################
## get X-linked MB genes

make_named_numeric_vector_MB_X <- function(TTT){
  # class MB genes on the X as sig (0) and rest as non-sig (1)
  TTT$MB_sig   <- ifelse(TTT$soft_chr_class == "X" & TTT$FDR < 0.05 & TTT$logFC > 0, 0, 1)
  print(length(subset(TTT, TTT$MB_sig  == 0)[,1]))
  TTT_2        <- as.data.frame(cbind(TTT$genes,TTT$MB_sig))
  TTT_2$V2 <- as.numeric(as.character(TTT_2$V2))
  full_list <- as.list(TTT_2)
  full_list_GL <- full_list$V2
  names(full_list_GL) <- full_list$V1  
  return(full_list_GL)
  }

Tbi_GL <- make_named_numeric_vector_MB_X(Tbi_SB_RT_FPKM_TTT) 
Tce_GL <- make_named_numeric_vector_MB_X(Tce_SB_RT_FPKM_TTT) 
Tcm_GL <- make_named_numeric_vector_MB_X(Tcm_SB_RT_FPKM_TTT) 
Tpa_GL <- make_named_numeric_vector_MB_X(Tpa_SB_RT_FPKM_TTT) 
Tps_GL <- make_named_numeric_vector_MB_X(Tps_SB_RT_FPKM_TTT) 


make_named_numeric_vector_MB_Xonly <- function(TTT){
  # class MB genes on the X as sig (0) and rest as non-sig (1)
  TTT <- subset(TTT, TTT$soft_chr_class  == "X")
  TTT$MB_sig   <- ifelse(TTT$soft_chr_class == "X" & TTT$FDR < 0.05 & TTT$logFC > 0, 0, 1)
  print(length(subset(TTT, TTT$MB_sig  == 0)[,1]))
  TTT_2        <- as.data.frame(cbind(TTT$genes,TTT$MB_sig))
  TTT_2$V2 <- as.numeric(as.character(TTT_2$V2))
  full_list <- as.list(TTT_2)
  full_list_GL <- full_list$V2
  names(full_list_GL) <- full_list$V1  
  return(full_list_GL)
}

Tbi_GL_Xonly  <- make_named_numeric_vector_MB_Xonly(Tbi_SB_RT_FPKM_TTT) 
Tce_GL_Xonly  <- make_named_numeric_vector_MB_Xonly(Tce_SB_RT_FPKM_TTT) 
Tcm_GL_Xonly  <- make_named_numeric_vector_MB_Xonly(Tcm_SB_RT_FPKM_TTT) 
Tpa_GL_Xonly  <- make_named_numeric_vector_MB_Xonly(Tpa_SB_RT_FPKM_TTT) 
Tps_GL_Xonly  <- make_named_numeric_vector_MB_Xonly(Tps_SB_RT_FPKM_TTT)




run_enrichment <- function(genelist, ref, sig_for_GO){
  
  ### make rule for classing sig / non-sig 
  
  topDiffGenes <- function(allScore) {return(allScore < sig_for_GO)}
  
  #### make GOdata object
  #### setting node size as 10 so at least 10 genes must be annot per GO terms 
  #### do enrichment test
  
  GODATA_BP = new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 10)
  
  ### get N GOs used
  
  GO_term_use_BP_list = GODATA_BP@graph@nodes
  N_GO_term_use_BP = length(GODATA_BP@graph@nodes)
  resultFisher <- runTest(GODATA_BP, algorithm = "weight01", statistic = "fisher")
  
  ### combined tables
  allRes1_BP <- GenTable(GODATA_BP, Fisher_w01 = resultFisher, ranksOf = "Fisher_w01", topNodes = length(GODATA_BP@graph@nodes), numChar = 200)
  
  sig_fisher_BP_GO     = subset(allRes1_BP, allRes1_BP$Fisher_w01 < sig_for_GO)$GO.ID
  
  ## return everything!
  out_list = list("N_GO_term_use_BP" = N_GO_term_use_BP, 
                  "GO_term_use_BP_list" = GO_term_use_BP_list, 
                  "allRes1_BP" = allRes1_BP, 
                  "sig_fisher_BP_GO" = sig_fisher_BP_GO,
                  "GODATA_BP" = GODATA_BP) 
  return(out_list)
  
}

#### run the enrichment stuff (0.05)

Tbi_enrich  <- run_enrichment(Tbi_GL, geneID2GO_Tbi_Droso, 0.05)
Tce_enrich  <- run_enrichment(Tce_GL, geneID2GO_Tce_Droso, 0.05)
Tcm_enrich  <- run_enrichment(Tcm_GL, geneID2GO_Tcm_Droso, 0.05)
Tpa_enrich  <- run_enrichment(Tpa_GL, geneID2GO_Tpa_Droso, 0.05)
Tps_enrich  <- run_enrichment(Tps_GL, geneID2GO_Tps_Droso, 0.05)
head(Tbi_enrich$allRes1_BP, n = 20)
head(Tce_enrich$allRes1_BP, n = 20)
head(Tcm_enrich$allRes1_BP, n = 20)
head(Tpa_enrich$allRes1_BP, n = 20)
head(Tps_enrich$allRes1_BP, n = 20)


Tbi_enrich  <- run_enrichment(Tbi_GL, geneID2GO_Tbi_Arth, 0.05)
Tce_enrich  <- run_enrichment(Tce_GL, geneID2GO_Tce_Arth, 0.05)
Tcm_enrich  <- run_enrichment(Tcm_GL, geneID2GO_Tcm_Arth, 0.05)
Tpa_enrich  <- run_enrichment(Tpa_GL, geneID2GO_Tpa_Arth, 0.05)
Tps_enrich  <- run_enrichment(Tps_GL, geneID2GO_Tps_Arth, 0.05)
head(Tbi_enrich$allRes1_BP, n = 20)
head(Tce_enrich$allRes1_BP, n = 20)
head(Tcm_enrich$allRes1_BP, n = 20)
head(Tpa_enrich$allRes1_BP, n = 20)
head(Tps_enrich$allRes1_BP, n = 20)



### X_only

Tbi_enrich_Xonly  <- run_enrichment(Tbi_GL_Xonly, geneID2GO_Tbi_Droso, 0.05)
Tce_enrich_Xonly  <- run_enrichment(Tce_GL_Xonly, geneID2GO_Tce_Droso, 0.05)
Tcm_enrich_Xonly  <- run_enrichment(Tcm_GL_Xonly, geneID2GO_Tcm_Droso, 0.05)
Tpa_enrich_Xonly  <- run_enrichment(Tpa_GL_Xonly, geneID2GO_Tpa_Droso, 0.05)
Tps_enrich_Xonly  <- run_enrichment(Tps_GL_Xonly, geneID2GO_Tps_Droso, 0.05)
head(Tbi_enrich_Xonly$allRes1_BP, n = 20)
head(Tce_enrich_Xonly$allRes1_BP, n = 20)
head(Tcm_enrich_Xonly$allRes1_BP, n = 20)
head(Tpa_enrich_Xonly$allRes1_BP, n = 20)
head(Tps_enrich_Xonly$allRes1_BP, n = 20)



