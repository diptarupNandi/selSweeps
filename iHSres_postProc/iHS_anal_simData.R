## Packages
library(dplyr)
library(ggplot2)

source("~/Documents/kaaj/selectSign/selSweeps/scripts/iHSres_postProc/propHiIHS_wind.R")
source("~/Documents/kaaj/selectSign/selSweeps/scripts/iHSres_postProc/propXiHS_clusters.R")


# Set parent result folder as the working directory 
setwd("~/Documents/kaaj/selectSign/selSweeps/selScan/simulations/results/sweep_drift_40Mb/")

# Load the text file containing the selected population IDs 
simParams = readLines("simParams.txt")


# Initiate lists, dataframes etc
ls_iHS_top1 = list()
iHS_top1_chr = list()
propHiIHSWind = list()
hiIHSwind_top1_full_ls = list()
propHiIHSwind_full_ls = list()

# Plots top 1% iHS scores of a chromosome of all the populations in a single file
for (chr in 1:2) {
  
  print(paste("Started chromosome", chr))
  # initiate counter used of list assigments
  iter = 1
  for (sPar in simParams){
    
    print(paste("for", sPar))
    
    # Load the nomralised iHS scores 
    iHS_chr = data.frame(read.table(paste("sweep_chr",chr,"_",sPar,"_100ind_40Mb.biAll.ihs.out.100bins.norm.final",sep=""), header = FALSE, skip = 0))
    colnames(iHS_chr) <- c("locID", "chrPos", "freq_d","iHH_d","iHH_a","iHS_unst","iHS","qScore","cmPos")
    
    # Create an extra-column containing |iHS| scores
    iHS_chr = iHS_chr %>% mutate(iHS_mod = Mod(iHS))
    # Calculate the 99% quantile value of the |iHS|
    val_top1<-unname(quantile(iHS_chr$iHS_mod,.99))
    # Select all the |iHS| scores that are higher than the 99% quantile value (top 1%)
    iHS_top1 = iHS_chr %>% filter(iHS_mod > val_top1) %>% mutate(simNe = sPar)
    # Save the object into a list 
    ls_iHS_top1[[iter]] = iHS_top1
    
    propHiIHSWind_temp = findXiHS_SNPclusters(iHS_chr,1,0.1)
    propHiIHSWind[[iter]] =  propHiIHSWind_temp  %>% mutate(simNe = sPar)
    
    iter = iter + 1
  }
  # Name the dataframe that would store the top 1% |iHS| scores of a given chromosome for all the populations
  dfNam = paste("iHS_top1_chr", chr, sep="")
  dfNam0 = paste("propHiiHSwind_chr", chr, sep="")
  # Convert the list object to a dataframe 
  assign(dfNam, do.call(rbind.data.frame, ls_iHS_top1))
  assign(dfNam0, do.call(rbind.data.frame, propHiIHSWind))
  # Create dataframes 
  iHS_top1 = get(dfNam)
  propHiiHSwind = get(dfNam0)
  
  # Filters the windows that contains at least 10% hiIHS SNPs and 20 SNPs
  # hiIHSwind_10p = propHiiHSwind %>% filter (propTop1 > 0.1 & nSNPs > 10)
  
  # Filters top 1% of the windows based on the proportion of hiIHS SNPs
  # hiIHSwind_top1 = propHiiHSwind %>% group_by(simNe) %>% filter(propTop1 > quantile(propTop1,.99) & nSNPs > 10)
  
  
  # collating the top1% iHS values for all chromosomes 
  iHS_top1_chr[[chr]] = iHS_top1 %>% mutate(chromosome = chr)
  
  # collating the top1% hiIHSwind data for all chromosomes 
  # hiIHSwind_top1_full_ls[[chr]] = hiIHSwind_top1 %>% mutate(chromosme = chr)
  
  # collating the hiIHS window data for all chromosomes 
  propHiIHSwind_full_ls[[chr]] = propHiiHSwind %>% mutate(chromosome = chr)
  
  # create the ggplot Object with chromosome positions as X and top 1% |iHS| as Y
  # ggplot() + geom_point(data=iHS_top1, aes(x = chrPos/1000000, y = iHS_mod, colour=simNe), size=.1) + facet_grid(simNe ~ .) +
    # theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

  
  
  # Create the output file name
  # outFname = paste("xIHS_1cM_drift_20Mb_chr", chr, ".jpeg", sep = "")
  # Save the plot for a given chromosome
  # ggsave(outFname, path =  "/home/popgen_ws1/Documents/kaaj/selectSign/selSweeps/selScan/resultsUS/iHS/set1/simPlots/driftOnly_20Mb/")
  
  print(paste("Finished all populations for chromosome", chr))
}

# Collated df conatining all the top1% hiIHS windows on all chromosomes 
# hiIHSwind_top1_full = do.call(rbind.data.frame, hiIHSwind_top1_full_ls)

# Collated df conatining all the hiIHS windows on all chromosomes 
iHS_top1_full = do.call(rbind.data.frame, iHS_top1_chr)

ggplot() + geom_point(data=iHS_top1_full, aes(x = chrPos/1000000, y = iHS_mod, colour=simNe), size=.5) + facet_grid(simNe ~ chromosome) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|") + 
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=10), axis.title = element_text(size = 12))


# Collated df conatining all the hiIHS windows on all chromosomes 
propHiIHSwind_full = do.call(rbind.data.frame, propHiIHSwind_full_ls)


ggplot(data=propHiIHSwind_full, aes(x=propTop1_pB)) + geom_histogram(binwidth = 0.015) + facet_grid(simNe ~ chromosome) +
  theme_bw() + xlab("Proportion of extreme iHS SNPs in the windows") + 
  theme(strip.background = element_rect(fill="white", color="black"), legend.position = "none",
        axis.text = element_text(size=10), axis.title = element_text(size = 12))


ggplot(data=propHiIHSwind_full, aes(x=propTop1_nB)) + geom_boxplot() + facet_grid(simNe ~ chromosome) +
  theme_bw() + xlab("Proportion of extreme IHS SNPs in the windows")


  
