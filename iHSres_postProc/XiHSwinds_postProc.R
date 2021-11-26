library(dplyr)
library(ggplot2)
library(ggdendro)


# Store path to functions
# pathDesk=c("~/Documents/kaaj/selectSign/selSweeps/scripts/iHSres_postProc/")
pathMac=c("~/OneDrive/Documents/kaaj/selectSign/scripts/selSweeps/iHSres_postProc/")

# Call the function files
# source(paste(pathMac,"propHiIHS_wind.R",sep=""))
source(paste(pathMac,"overlaps_genomRegions.R", sep = ""))


# Set parent result folder as the working directory 
# setwd("~/Documents/kaaj/selectSign/selSweeps/selScan/resultsUS/iHS/set1/")
setwd("~/OneDrive/Documents/kaaj/selectSign/iHS/dataFiles/XiHSwinds/")

# Load the file containing xtreme iHS windows of 26 SAS populations 
XiHSwind_1cM_SAS = read.csv("XiHSwind_1cM_plusB_full_SAS.csv", header = TRUE)



# Plot histograms for proportions of xtreme iHS SNPs separately for the populations
XiHSwind_1cM_SAS %>% filter(nSNPs > 20) %>% 
  ggplot(aes(x=propTopMax)) + 
  geom_histogram(binwidth = 0.001, boundary = 0, fill = "blue", colour = "blue") + 
  facet_grid(pOp ~ .) + geom_vline(xintercept=0.1, colour = "Red", size = .2) + 
  theme_bw() + xlab("Proportion of extreme iHS SNPs in the windows") + 
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=8), axis.title = element_text(size = 12))

# Plot histograms for proportions of xtreme iHS SNPs separately for the populations
XiHSwind_1cM_SAS %>% filter(nSNPs > 20) %>% mutate(winSz = endWin - startWin) %>%
  ggplot(aes(x=winSz/1000000)) + geom_histogram(binwidth = 0.25)  + theme_bw() + 
  xlab("Size of windows (Mb)") + scale_x_continuous(breaks = seq(0,15,1))
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=12), axis.title = element_text(size = 15))

XiHSwind_1cM_SAS %>% mutate(winSz = endWin - startWin) %>% 
  group_by(pOp) %>% summarise(median(winSz))  


# Filter the main df for windows with > 10% extreme iHS SNPs and > 20 SNps
XiHSwind_1cM_SAS_sel = XiHSwind_1cM_SAS %>% filter(propTopMax > 0.1 & nSNPs > 20) 
# Estimate number of xtreme iHS windows per population
XiHSwind_1cM_SAS %>% filter(propTopMax > 0.1 & nSNPs > 20) %>% group_by(pOp) %>% summarise(Num_windows = n())
# Estimate number of xtreme iHS windows (20% xtreme iHS SNPs) per population
XiHSwind_1cM_SAS %>% filter(propTop1 > 0.2 & nSNPs > 20) %>% group_by(pOp) %>% summarise(Num_windows = n())


XiHSwind_1cM_SAS_sel %>% mutate(winSz = endWin - startWin) %>%
  ggplot(aes(x=winSz/1000000)) + geom_histogram(binwidth = 0.25, boundary = 0)  + 
  theme_bw() + xlab("Size of windows (Mb)") + scale_x_continuous(breaks = seq(0,4,.5)) +
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=12), axis.title = element_text(size = 15))

XiHSwind_1cM_SAS_sel %>% mutate(winSz = endWin - startWin) %>% 
  group_by(pOp) %>% summarise(median(winSz))

XiHSwind_1cM_SAS_sel %>% mutate(winSz = endWin - startWin) %>%
  ggplot(aes(x=pOp, y=winSz/1000000)) + geom_boxplot()   + 
  theme_bw() + xlab("Size of windows (Mb)") + scale_y_continuous(breaks = seq(0,4,.5)) +
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=12), axis.title = element_text(size = 15))

# Input df for overlap estimation: only 4 columns from the filtered df
ovlapInputDf = XiHSwind_1cM_SAS_sel %>% select(startWin,endWin,chromosme,pOp)
# change pop IDs into factors
ovlapInputDf$pOp=as.factor(ovlapInputDf$pOp)
# Run the overlap estimation program
ovelap_1cm_SAS = findOverlaps(ovlapInputDf)

# Store the population IDs
# populations = unname(levels(ovlapInputDf$pOp))
# Name th rows of overlap matrix
# rownames(ovelap_1cm_SAS) = populations

# Storing the diagonals = number of xtreme iHS windows for that population
ovelap_diag = diag(as.matrix(ovelap_1cm_SAS))

# Calculate distance from the overlap window data
ovlp_1cM_dist = as.dist(t(1 - ovelap_1cm_SAS/ovelap_diag))
ovlp_1cM_dist = dist(ovelap_1cm_SAS/ovelap_diag)
# Run hierarchical clustering of the overlap distance matrix
ovlp_1cM_hc = hclust(ovlp_1cM_dist, method = "average")
# Plot a dendogram 
ggdendrogram(ovlp_1cM_hc, rotate = FALSE, theme_dendro = FALSE)

metDat = read.csv("pop_metadata.csv", header=TRUE)
colnames(metDat)[7]="windows"
metDat$Demography = as.factor(metDat$Demography)
metDat$Lifestyle = as.factor(metDat$Lifestyle)
levels(metDat$Demography)[which(! levels(metDat$Demography) %in% "Bottleneck")] = "noBottleneck" 


ggplot(data=metDat, aes(x=Samples, y=windows, colour = Demography)) + geom_point(size=5) + theme_bw() + 
  xlab("Number of samples per population") + ylab("Number of extreme iHS windows") +
  theme(strip.background = element_rect(fill="white", color="black"),
        axis.text = element_text(size=12), axis.title = element_text(size = 15),
        legend.text = element_text(size = 12))

cor.test(metDat$Samples,metDat$windows, method = c("spearman"))

ggplot(data=metDat, aes(x=Demography, y=windows)) + geom_boxplot() + theme_bw() + 
  xlab("Demographic history") + ylab("Number of extreme iHS windows") +
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=12), axis.title = element_text(size = 15))

actual_diff = diff(tapply(metDat$windows, metDat$Demography, median))


## Permutation test function
permTest = function(var,treatmnt,iter,statistic){
  
  # Initializing empty dataframes and lists
  pVal=0
  effSz=c()
  
  for (i in 1:iter){
    effSz[i]=diff(by(var, sample(treatmnt, length(treatmnt), FALSE), statistic))
  }
  
  actDiff = diff(tapply(var, treatmnt, statistic))
  pVal=sum(abs(effSz) >= abs(actDiff))/(iter)
  return(list(actDiff,effSz,pVal))
}

res_Dem = permTest(metDat$windows,metDat$Demography,10000,statistic = median)

hist(res_Dem[[2]], breaks=50, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=res_Dem[[1]], lwd=3, col="red")



ggplot(data=metDat, aes(x=Lifestyle, y=windows)) + geom_boxplot() + theme_bw() + 
  xlab("Lifestyle") + ylab("Number of extreme iHS windows") +
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=12), axis.title = element_text(size = 15))



annotWind = data.frame(read.table("./annotXiHSwinds/annotXiHSwindsSel_1cM_SAS.bed", header=FALSE, skip =0))

colnames(annotWind)=c("chrom","start","end","pop","annotations","fChrom","fStart","fEnd")

annotWind$pop = as.factor(annotWind$pop)
annotWind$annotations = factor(annotWind$annotations, 
                                  levels=c("exon","intron","intergenic"))

annotWindFin = annotWind %>% mutate(bpLength = (fEnd - fStart))
annotWind_summary = annotWindFin %>% group_by(pop,annotations) %>% 
  summarise(bpSize=sum(bpLength)) %>% mutate(relFreq = bpSize/sum(bpSize))

ggplot(data = annotWind_summary, aes(x = pop, y = relFreq, 
                                     fill = annotations, order = desc(annotations))) + 
  geom_bar(stat = "identity") + theme_bw() + xlab("SAS populations") + 
  ylab("Relative proportions") + theme(strip.background = element_rect(fill="white", color="black"),
        axis.text = element_text(size=12), axis.title = element_text(size = 15),
        legend.text = element_text(size=12))

ggplot(data = annotWind_summary, aes(x = pop, y = bpSize/1000000, 
                                     fill = annotations, order = desc(annotations))) + 
  geom_bar(stat = "identity") + theme_bw() + xlab("SAS populations") + 
  ylab("Size (Mbp)") + theme(strip.background = element_rect(fill="white", color="black"),
                             axis.text = element_text(size=12), axis.title = element_text(size = 15),
                             legend.text = element_text(size=12))



annotWindFull = data.frame(read.table("./annotXiHSwinds/annotXiHSwindsF_1cM_SAS.bed", header=FALSE, skip =0))

colnames(annotWindFull)=c("chrom","start","end","pop","annotations","fChrom","fStart","fEnd")

annotWindFull$pop = as.factor(annotWindFull$pop)
annotWindFull$annotations = factor(annotWindFull$annotations, 
                               levels=c("exon","intron","intergenic"))

annotWindFull = annotWindFull %>% mutate(bpLength = (fEnd - fStart))
annotWindFull_summary = annotWindFull %>% group_by(pop,annotations) %>% 
  summarise(bpSize=sum(bpLength)) %>% mutate(relFreq = bpSize/sum(bpSize))

ggplot(data = annotWindFull_summary, aes(x = pop, y = relFreq, 
                                     fill = annotations, order = desc(annotations))) + 
  geom_bar(stat = "identity") + theme_bw() + xlab("SAS populations") + 
  ylab("Relative proportions") + theme(strip.background = element_rect(fill="white", color="black"),
                                       axis.text = element_text(size=12), axis.title = element_text(size = 15),
                                       legend.text = element_text(size=12))
---------------------------------------------------------------------------------------------------

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
  