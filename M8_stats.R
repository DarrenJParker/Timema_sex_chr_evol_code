## data
dat_M8 <- read.csv("data/selection/selectome_M8_wchr.csv")
head(dat_M8)
str(dat_M8)

### get pvals
## final D is 2*LR from Godon
## lookup pval with 1.d.f and divide by 2 as an approximation to take into account the test statistic does not have an asymptotic χ2 distribution. 
## Instead, the correct null distribution is a 1:1 mixture of point mass 0 and χ2 (Chernoff 1954; see also Self and Liang 1987, case 5)". https://academic.oup.com/mbe/article/28/3/1217/993378

dat_M8$p <- pchisq(dat_M8$finalD, df=1, lower.tail=F) / 2
dat_M8$qvalue <- p.adjust(dat_M8$p, method = "BH")

dat_M8$chr_class <- ifelse(dat_M8$soft_chr == "AAAA", "A",
                           ifelse(dat_M8$soft_chr == "AAAAA", "A",
                                  ifelse(dat_M8$soft_chr == "XXXX", "X",
                                         ifelse(dat_M8$soft_chr == "XXXXX", "X", "ERROR"))))

levels(as.factor(dat_M8$chr_class))
# [1] "A" "X"

##### are genes showing +ve sel overrep on the X?


dat_M8_X <- subset(dat_M8, dat_M8$chr_class == "X")
dat_M8_A <- subset(dat_M8, dat_M8$chr_class == "A")
pos_q_threh = 0.05
X_genes_sel      = length(subset(dat_M8_X, dat_M8_X$qvalue  <  pos_q_threh)[,1])
X_genes_no_sel   = length(subset(dat_M8_X, dat_M8_X$qvalue  >= pos_q_threh)[,1])
All_genes_sel    = length(subset(dat_M8,   dat_M8$qvalue  <  pos_q_threh)[,1])
All_genes_no_sel = length(subset(dat_M8,   dat_M8$qvalue  >= pos_q_threh)[,1])
  
FT_mat <- matrix(c(X_genes_sel,(All_genes_sel - X_genes_sel), X_genes_no_sel,(All_genes_no_sel - X_genes_no_sel)), nrow = 2)
  
fisher.test(FT_mat, alternative="two.sided")


# Fisher's Exact Test for Count Data
# 
# data:  FT_mat
# p-value = 0.1077
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5338199 10.3063915
# sample estimates:
# odds ratio 
#    2.89281 




