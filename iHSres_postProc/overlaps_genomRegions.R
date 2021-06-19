findOverlaps <- function(posDf){
  colnames(posDf) = c("startPos", "endPos", "chrom", "pop")
  population = unname(levels(posDf$pop))
  chromosome = unique(posDf$chrom)
  
  # Initiate lists, dataframes etc
  
  
  ovlapDf = data.frame(matrix(ncol = length(population), nrow = length(population)))
  colnames(ovlapDf) = population
  colnames(ovlapDf) = population
  
  for (p in population){
    iter = 1
    tempLs = list()
    
    for (chr in chromosome){
      objDF = posDf %>% filter(pop %in% p, chrom == chr)
      targDf = posDf %>% filter(! pop %in% p, chrom == chr)
      if (nrow(objDF) == 0 | nrow(targDf) == 0){
        next
      } else {
        for (ind in 1:nrow(objDF)){
          startP = objDF$startPos[ind]
          endP = objDF$endPos[ind]
          tempCoord = targDf %>% filter((startPos >= startP & startPos <= endP) | (endPos >= startP & endPos <= endP) | (startPos <= startP & endPos >= endP))
          if (nrow(tempCoord) == 0){
            # If empty then move to the next iteration
            next
          } else {
            tempLs[[iter]] = tempCoord 
            iter = iter + 1
          }
        }
      }
        
      }
      tempOvDf = do.call(rbind.data.frame, tempLs)
      ovSummary = tempOvDf %>% group_by(pop) %>% summarise(ovlpNum = n())
      
      if (nrow(ovSummary) < length(population)){
        missPop = population[! population %in% ovSummary$pop]
        tempPop = data.frame(pop=missPop, ovlpNum=rep(0,length(missPop)))
        ovSummary_temp = rbind(ovSummary, tempPop)
        ovSummary_temp$pop = as.character(ovSummary_temp$pop)
        ovSummaryF = ovSummary_temp %>% arrange(pop)
        ovlapDf[p] = ovSummaryF$ovlpNum
      } else {
        ovlapDf[p] = ovSummary$ovlpNum
      }
  }
  
  return(ovlapDf)
  
}