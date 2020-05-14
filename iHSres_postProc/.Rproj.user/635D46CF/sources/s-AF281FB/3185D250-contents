## Packages
library(dplyr)
library(ggplot2)

source("~/Documents/kaaj/selectSign/selScan/scripts/iHS/propHiIHS_wind.R")


# Set parent result folder as the working directory 
setwd("~/Documents/kaaj/selectSign/selScan/resultsUS/iHS/set1/")

# Load the text file containing the selected population IDs 
popID = readLines("popSelect.txt")

# Initiate lists, dataframes etc
ls_iHS_top1 = list()
propHiIHSWind = list()

# Plots top 1% iHS scores of a chromosome of all the populations in a single file
for (chr in 1:22) {
  
  print(paste("Started chromosome", chr))
  # initiate counter used of list assigments
  iter = 1
  for (pID in popID){
    
    print(paste("for", pID, "populations"))
    # Load the nomralised iHS scores 
    iHS_chr = data.frame(read.table(paste(pID,"/","GA100k_chr",chr,"_",pID,".ihs.out.100bins.norm.final",sep=""), header = FALSE, skip = 0))
    colnames(iHS_chr) <- c("locID", "chrPos", "freq_d","iHH_d","iHH_a","iHS_unst","iHS","qScore","cmPos")
    
    # Create an extra-column containing |iHS| scores
    iHS_chr = iHS_chr %>% mutate(iHS_mod = Mod(iHS))
    # Calculate the 99% quantile value of the |iHS|
    val_top1<-unname(quantile(iHS_chr$iHS_mod,.99))
    # Select all the |iHS| scores that are higher than the 99% quantile value (top 1%)
    iHS_top1 = iHS_chr %>% filter(iHS_mod > val_top1) %>% mutate(pOp = pID)
    # Save the object into a list 
    ls_iHS_top1[[iter]] = iHS_top1
    
    propHiIHSWind_temp = find_HiiHS_SNPclusters(iHS_chr,.5)
    propHiIHSWind[[iter]] =  propHiIHSWind_temp  %>% mutate(pOp = pID)
    
    iter = iter + 1
  }
  # Name the dataframe that would store the top 1% |iHS| scores of a given chromosome for all the populations
  dfNam = paste("iHS_top1_chr", chr, sep="")
  dfNam0 = paste("propHiiHSwind_chr", chr, sep="")
  # Convert the list object to a dataframe 
  assign(dfNam, do.call(rbind.data.frame, ls_iHS_top1))
  assign(dfNam0, do.call(rbind.data.frame, propHiIHSWind))
  # Create a temporary dataframe to use for plotting
  iHS_top1 = get(dfNam)
  propHiiHSwind = get(dfNam0)
  
  hiIHSwind = propHiiHSwind %>% filter (propTop1 > 0.1 & nSNPs > 20)
  
  # create the ggplot Object with chromosome positions as X and top 1% |iHS| as Y
  # plotObj <- ggplot(data=tempDat, aes(x = chrPos/1000000, y = iHS_mod, colour=pOp)) 
  ggplot() +  geom_rect(data = hiIHSwind, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf, fill = pOp), alpha = 0.4) +
    geom_point(data=iHS_top1, aes(x = chrPos/1000000, y = iHS_mod, colour=pOp), size=.1) + facet_grid(pOp ~ .) +
    theme_bw() + xlab("Position (Mb)") + ylab("|iHS|") 
  
  # Generate the scatter plot
  # plotObj + geom_point(size=.1) + facet_grid(pOp ~ .) + theme_bw() +  
  #   xlab("Position (Mb)") + ylab("|iHS|")
  
  # Create the output file name
  outFname = paste("iHS_top1_plot_chr", chr, ".pdf", sep = "")
  # Save the plot for a given chromosome
  ggsave(outFname, path =  "iHS_plots/")
  
  print(paste("Finished all populations for chromosome", chr))
}
