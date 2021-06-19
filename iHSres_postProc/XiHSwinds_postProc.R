library(dplyr)
library(ggplot2)
library(ggdendro)

# Set parent result folder as the working directory 
setwd("~/Documents/kaaj/selectSign/selSweeps/selScan/resultsUS/iHS/set1/")

XiHSwind_1cM_SAS = read.csv("XiHSwind_1cM_plusB_full_SAS.csv", header = TRUE)


XiHSwind_1cM_SAS_filtrd = XiHSwind_1cM_SAS %>% filter(nSNPs > 20)
ggplot(data=XiHSwind_1cM_SAS_filtrd, aes(x=propTopMax)) + geom_histogram(binwidth = 0.001, boundary = 0, fill = "blue", colour = "blue") + facet_grid(pOp ~ .) +
  geom_vline(xintercept=0.1, colour = "Red", size = .2) + theme_bw() + xlab("Proportion of extreme iHS SNPs in the windows") + 
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=10), axis.title = element_text(size = 12))

XiHSwind_1cM_SAS_sel = XiHSwind_1cM_SAS %>% filter(propTopMax > 0.1 & nSNPs > 20) 
XiHSwind_1cM_SAS_summary = XiHSwind_1cM_SAS %>% filter(propTop1 > 0.1 & nSNPs > 20) %>% group_by(pOp) %>% summarise(Num_windows = n())
XiHSwind_1cM_SAS %>% filter(propTop1 > 0.2 & nSNPs > 20) %>% group_by(pOp) %>% summarise(Num_windows = n())

ovlapInputDf = XiHSwind_1cM_SAS_sel %>% select(startWin,endWin,chromosme,pOp)
source("~/Documents/kaaj/selectSign/selSweeps/scripts/iHSres_postProc/overlaps_genomRegions.R")
ovelap_1cm_SAS = findOverlaps(ovlapInputDf)
rownames(ovelap_1cm_SAS) = population

ovlp_1cM_dist = as.dist(t(1 - ovelap_1cm_SAS/colSums(ovelap_1cm_SAS)))
ovlp_1cM_hc = hclust(ovlp_1cM_dist, method = "average")

ggdendrogram(ovlp_1cM_hc, rotate = TRUE, theme_dendro = FALSE)

XiHSwind_halfcM_SAS = read.csv("XiHSwind_halfcM_plusB_full_SAS.csv", header = TRUE)


XiHSwind_halfcM_SAS_filtrd = XiHSwind_halfcM_SAS %>% filter(nSNPs > 20)
ggplot(data=XiHSwind_halfcM_SAS_filtrd, aes(x=propTopMax)) + geom_histogram(binwidth = 0.001, boundary = 0, fill = "blue", colour = "blue") + facet_grid(pOp ~ .) +
  geom_vline(xintercept=0.1, colour = "Red", size = .2) + theme_bw() + xlab("Proportion of extreme iHS SNPs in the windows") + 
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=10), axis.title = element_text(size = 12))

XiHSwind_halfcM_SAS_sel = XiHSwind_halfcM_SAS %>% filter(propTopMax > 0.1 & nSNPs > 20)
XiHSwind_halfcM_SAS_sumry = XiHSwind_halfcM_SAS_sel %>% group_by(pOp) %>% summarise(windNum = n())

XiHSwind_halfcM_SAS_selB = XiHSwind_halfcM_SAS %>% filter(propTop1 > 0.1 & nSNPs > 20)
XiHSwind_halfcM_SAS %>% filter(propTopMax > 0.1 & nSNPs > 20) %>% group_by(pOp) %>% summarise(Num_windows = n())

ovlapInputDf_halfcM = XiHSwind_halfcM_SAS_sel %>% select(startWin,endWin,chromosme,pOp)
ovelap_halfcM_SAS = findOverlaps(ovlapInputDf_halfcM)

rownames(ovelap_halfcM_SAS) = population

ovlp_halfcM_dist = as.dist(t(1 - ovelap_halfcM_SAS/colSums(ovelap_halfcM_SAS)))
ovlp_halfcM_hc = hclust(ovlp_halfcM_dist, method = "average")
ggdendrogram(ovlp_halfcM_hc, rotate = TRUE, theme_dendro = FALSE)

nonovlp_halfcM_dist = as.dist((XiHSwind_halfcM_SAS_sumry$windNum - ovelap_halfcM_SAS)/XiHSwind_halfcM_SAS_sumry$windNum)
nonovlp_halfcM_hc = hclust(nonovlp_halfcM_dist, method = "average")

ggdendrogram(nonovlp_halfcM_hc, rotate = TRUE, theme_dendro = FALSE)
  