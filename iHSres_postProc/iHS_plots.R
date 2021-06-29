## Packages
library(dplyr)
library(ggplot2)

# Paths on the popgen_ws1@desktop
#source("~/Documents/kaaj/selectSign/selSweeps/scripts/iHSres_postProc/propHiIHS_wind.R")
#source("~/Documents/kaaj/selectSign/selSweeps/scripts/iHSres_postProc/propXiHS_clusters.R")

# Paths on the GI cluster
source("/ndata/home/genomeindia/gi-advance/gi3/dipto/aNalysis/selectSign/selSweeps/iHSres_postProc/propHiIHS_wind.R")
source("/ndata/home/genomeindia/gi-advance/gi3/dipto/aNalysis/selectSign/selSweeps/iHSres_postProc/propXiHS_clusters.R")


# Set parent result folder as the working directory 
setwd("/ndata/home/genomeindia/gi-advance/gi3/dipto/aNalysis/selectSign/selscan/resultS/iHS/setFinal/")

# Load the text file containing the selected population IDs 
popID = readLines("selectSA.pops")

# Initiate lists, dataframes etc
ls_iHS_top1 = list()
propHiIHSWind = list()
hiIHSwind_full_ls = list()

# Plots top 1% iHS scores of a chromosome of all the populations in a single file
for (chr in 1:22) {
  
  print(paste("Started chromosome", chr))
  # initiate counter used of list assigments
  iter = 1
  for (pID in popID){
    
    print(paste("of", pID, "populations"))
    # Load the nomralised iHS scores 
    iHS_chr = data.frame(read.table(paste(pID,"/","GA100k_chr",chr,"_",pID,".ihs.out.100bins.norm.final",sep=""), header = FALSE, skip = 0))
    colnames(iHS_chr) <- c("locID", "chrPos", "freq_d","iHH_d","iHH_a","iHS_unst","iHS","qScore","cmPos")
    
    # Create an extra-column containing |iHS| scores
    iHS_chr = iHS_chr %>% mutate(iHS_mod = Mod(iHS))
    # Calculate the 99% quantile value of the |iHS|
    #val_top1<-unname(quantile(iHS_chr$iHS_mod,.99))
    # Select all the |iHS| scores that are higher than the 99% quantile value (top 1%)
    #iHS_top1 = iHS_chr %>% filter(iHS_mod > val_top1) %>% mutate(pOp = pID)
    # Save the object into a list 
    #ls_iHS_top1[[iter]] = iHS_top1
    
    propHiIHSWind_temp = findXiHS_SNPclusters(iHS_chr,1,0.2)
  
    propHiIHSWind[[iter]] =  propHiIHSWind_temp  %>% mutate(pOp = pID)
    
    iter = iter + 1
  }
  # Name the dataframe that would store the top 1% |iHS| scores of a given chromosome for all the populations
  #dfNam = paste("iHS_top1_chr", chr, sep="")
  dfNam0 = paste("propHiiHSwind_chr", chr, sep="")
  # Convert the list object to a dataframe 
  #assign(dfNam, do.call(rbind.data.frame, ls_iHS_top1))
  assign(dfNam0, do.call(rbind.data.frame, propHiIHSWind))
  # Create a temporary dataframe to use for plotting
  #iHS_top1 = get(dfNam)
  
  # Main dataframe conraining the xtreme iHS windows
  propHiiHSwind = get(dfNam0)
  
  
  
  # Estimates the window containing max proportion of xIHS SNPs and adds those columns 
  XiHSwind_full = propHiiHSwind %>% rowwise() %>% 
    mutate(propTopMax = max(c(propTop1,propTop1_pB,propTop1_nB)), 
           startWin = case_when(
             which.max(c(propTop1,propTop1_pB,propTop1_nB)) == 3 ~ fSNP_nB,
             which.max(c(propTop1,propTop1_pB,propTop1_nB)) == 1 | 2 ~ fSNP),
           endWin = case_when(
             which.max(c(propTop1,propTop1_pB,propTop1_nB)) == 2 ~ lSNP_pB,
             which.max(c(propTop1,propTop1_pB,propTop1_nB)) == 1 | 3 ~ lSNP))
  
  
  #XiIHSwind_temp = XiHSwind_full %>% filter(propTop1 > 0.1 & nSNPs > 20) %>% 
  #  select(propTop1,fSNP,lSNP,pOp)
  
  
  # Filters the windows that contains at least 10% hiIHS SNPs and 20 SNPs
  #XiIHSwind_fin = XiHSwind_full %>% filter(propTopMax > 0.1 & nSNPs > 20) %>% 
  #  select(propTopMax,startWin,endWin,pOp)
  
  
  
  # Filters top 1% of the windows based on the proportion of hiIHS SNPs
  # hiIHSwind_top1 = propHiiHSwind %>% group_by(pOp) %>% filter(propTop1 > quantile(propTop1,.99) & nSNPs > 20)
  
  # collating the hiIHSwind data for all chromosomes 
  hiIHSwind_full_ls[[chr]] = XiHSwind_full %>% mutate(chromosme = chr)
  
  
  # create the ggplot Object with chromosome positions as X and top 1% |iHS| as Y
  # plotObj <- ggplot(data=tempDat, aes(x = chrPos/1000000, y = iHS_mod, colour=pOp)) 
  #ggplot() +  geom_rect(data = XiIHSwind_fin, 
  #                      aes(xmin = startWin/1000000, xmax = endWin/1000000, 
  #                          ymin = -Inf, ymax = Inf, fill = propTopMax), alpha = 0.33) +
  #  geom_point(data=iHS_top1, aes(x = chrPos/1000000, y = iHS_mod, colour=pOp), size=.1, show.legend = FALSE) + 
  #  facet_grid(pOp ~ .) + theme_bw() + xlab("Position (Mb)") + ylab("|iHS|") +
  #  theme(strip.background = element_rect(fill="white", color="black"),
  #        legend.title = element_blank()) +
  #  scale_fill_gradient(low = "#132B4385", high = "#56B1F785")
  
  # Generate the scatter plot
  # plotObj + geom_point(size=.1) + facet_grid(pOp ~ .) + theme_bw() +  
  #   xlab("Position (Mb)") + ylab("|iHS|")
  
  # Create the output file name
  #outFname = paste("XiHS_1plusBcM_10p_chr", chr, ".pdf", sep = "")
  # Save the plot for a given chromosome
  #ggsave(outFname, path =  "figures/XiHS_1cM_plusB_10p_plots/")
  
  #ggplot() +  geom_rect(data = XiIHSwind_temp, 
  #                      aes(xmin = fSNP/1000000, xmax = lSNP/1000000, 
  #                          ymin = -Inf, ymax = Inf, fill = propTop1), alpha = 0.33) +
  #  geom_point(data=iHS_top1, aes(x = chrPos/1000000, y = iHS_mod, colour=pOp), size=.1, show.legend = FALSE) + 
  #  facet_grid(pOp ~ .) + theme_bw() + xlab("Position (Mb)") + ylab("|iHS|") +
  #  theme(strip.background = element_rect(fill="white", color="black"),
  #        legend.title = element_blank()) +
  #  scale_fill_gradient(low = "#132B4385", high = "#56B1F785")
  
  # Create the output file name
  #outF1name = paste("XiHS_1cM_10p_chr", chr, ".pdf", sep = "")
  # Save the plot for a given chromosome
  #ggsave(outF1name, path =  "figures/XiHS_1cM_10p_plots/")
  
  print(paste("Finished all populations for chromosome", chr))
}

# Collated df conatining all the xtreme IHS windows on all chromosomes of all the populations
XiHSwind_full_SAS = do.call(rbind.data.frame, hiIHSwind_full_ls)
write.csv(XiHSwind_full_SAS,"/ndata/home/genomeindia/gi-advance/gi3/dipto/aNalysis/selectSign/selscan/resultS/iHS/XiHSwinds/XiHSwind_1cM_plusB_full_SAS.csv")



#XiHSwind_SAS_summary = XiHSwind_full_SAS %>% filter(propTop1 > 0.1 & nSNPs > 20) %>% group_by(pOp,chromosme) %>% summarise(Num_windows = n()) 

#hist(hiIHSwind_full$propTop1, xlab = "Proportion of extreme IHS SNPs", xlim = c(0,.8), ylim = c(0,200), main = "Distribution of extreme IHS SNPs proportions")


