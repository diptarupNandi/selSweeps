## Function to quantify gnomic windows that contain the highest proportion of snps with the top 1% of |iHS| score.
# Inputs modified selscan normalized output file (.norm) with an extra column containing the genetic distance (in cm) and 
# Size of the window in centimorgans. 


findXiHS_SNPclusters<-function(iHSdf,winSz,buf){
  # Set colnames to iHSdf 
  colnames(iHSdf) <- c("locID", "chrPos", "freq_d","iHH_d","iHH_a","iHS_unst","iHS","qScore","cmPos","iHS_mod")
  # Range of genetic distances (in Centimorgans)
  rang_cmPos = range(iHSdf$cmPos)
  # Start position of the scan = nearest multiple of the window size lower than the minimum cmPos
  startScan=winSz*floor(rang_cmPos[1]/winSz)
  # End position of the scan = nearest multiple of the window size higher than the maximum cmPos
  endScan=winSz*ceiling(rang_cmPos[2]/winSz)
  # Calculate the 99% quantile value of the |iHS|
  val_top1<-unname(quantile(iHSdf$iHS_mod,.99))
  
  # Initialize matrix, vectors
  gNomWind_summaryDf_cm1 = c()
  iter=1
  for (i in seq(startScan+winSz,endScan,winSz)){
    # Filter the snps that lie within scanning window
    iHS_SNPs = iHSdf %>% filter(cmPos >= (i - winSz) & cmPos < (i))
    
    # Filter the snps that lie within scanning window + buffer
    iHS_SNPs_pB = iHSdf %>% filter(cmPos >= (i - winSz) & cmPos < (i + buf))
    
    # Filter the snps that lie within scanning window - buffer
    iHS_SNPs_nB = iHSdf %>% filter(cmPos >= (i - winSz - buf) & cmPos < (i))
    
    if (nrow(iHS_SNPs) == 0){
      # If empty then move to the next iteration
      next
    } else {
      
      WindNo = iter
      # Start and end positions on scanning window in centimorgans
      startP = i - winSz
      endP = i
      # First and the last SNPs to lie within the window
      fSNP = iHS_SNPs$chrPos[1]
      lSNP = tail(iHS_SNPs$chrPos,1)
      # Total number of SNPs within the scanning window
      nSNPs = nrow(iHS_SNPs)
      # Proportion of SNPs in the window that have top 1% iHS score
      propTop1 = length(iHS_SNPs$iHS_mod[iHS_SNPs$iHS_mod > val_top1]) / nrow(iHS_SNPs)
      
      
      # Last SNP to lie within the window + buffer
      lSNP_pB = tail(iHS_SNPs_pB$chrPos,1)
      # Total number of SNPs within the scanning window + buffer
      nSNPs_pB = nrow(iHS_SNPs_pB)
      
      # Proportion of SNPs in the window+buffer that have top 1% iHS score
      propTop1pB = length(iHS_SNPs_pB$iHS_mod[iHS_SNPs_pB$iHS_mod > val_top1]) / nrow(iHS_SNPs_pB)
      
      # First SNP to lie within the window - buffer
      fSNP_nB = iHS_SNPs_nB$chrPos[1]
    
      # Total number of SNPs within the scanning window + buffer
      nSNPs_nB = nrow(iHS_SNPs_nB)
      
      # Proportion of SNPs in the window+buffer that have top 1% iHS score
      propTop1nB = length(iHS_SNPs_nB$iHS_mod[iHS_SNPs_nB$iHS_mod > val_top1]) / nrow(iHS_SNPs_nB)
      
      # Put all the parameters in a matrix
      gNomWind_summaryDf_cm1 = rbind(gNomWind_summaryDf_cm1,c(WindNo,startP,endP,fSNP,lSNP,lSNP_pB,fSNP_nB,nSNPs,nSNPs_pB,nSNPs_nB,propTop1,propTop1pB,propTop1nB))
      
      iter=iter+1
    }
  }
  
  # Convert the matrix to a dataframe
  propHiIHSwind_summary = data.frame(gNomWind_summaryDf_cm1)
  # Assign column names to the 
  colnames(propHiIHSwind_summary) = c("WindNo","startP","endP","fSNP","lSNP","lSNP_pB", "fSNP_nB","nSNPs","nSNPs_pB","nSNPs_nB","propTop1","propTop1_pB","propTop1_nB")
  propHiIHSwind_summary = propHiIHSwind_summary %>% arrange(desc(propTop1))
  return(propHiIHSwind_summary)
  
}