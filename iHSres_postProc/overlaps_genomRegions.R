findOverlaps <- function(posDf){
  # Set column names to the input matrix
  colnames(posDf) = c("startPos", "endPos", "chrom", "pop")
  # Vector with all population IDs
  population = unname(levels(posDf$pop))
  # Vector with all the chromosomes
  chromosome = unique(posDf$chrom)
  
  ## Initiate lists, dataframes etc
  
  ovlapDf = data.frame(matrix(ncol = length(population), nrow = length(population)))
  rownames(ovlapDf) = population
  colnames(ovlapDf) = population
  
  for (p in population){
    iter = 1
    tempLs = list()
    
    for (chr in chromosome){
      # Extract the windows on 'chr' chromosome of 'p' population  
      objDF = posDf %>% filter(pop %in% p, chrom == chr)
      # Extract the windows on the 'chr' chromosome for the rest
      targDf = posDf %>% filter(! pop %in% p, chrom == chr)
      
      # If no windows on the chr chromosome
      if (nrow(objDF) == 0 | nrow(targDf) == 0){
        next
      } else {
        # For each window
        for (ind in 1:nrow(objDF)){
          
          # Store the start and end positions
          startP = objDF$startPos[ind]
          endP = objDF$endPos[ind]
          
          # Check for overlap 
          tempCoord = targDf %>% filter((startPos >= startP & startPos <= endP) |
                                          (endPos >= startP & endPos <= endP) |
                                          (startPos <= startP & endPos >= endP))
          if (nrow(tempCoord) == 0){
            # If no overlap then move to the next iteration
            next
          } else {
            # Store the overlapping co-ordinates into a list
            tempLs[[iter]] = tempCoord 
            iter = iter + 1
          }
        }
      }
        
    }
    # Convert the list into a dataframe
    tempOvDf = do.call(rbind.data.frame, tempLs)
    # Count the overlapping windows for every population
    ovSummary = tempOvDf %>% group_by(pop) %>% summarise(ovlpNum = n())
    # Check for populations with no overlaps on that chromosome
    if (nrow(ovSummary) < length(population)){
      # Find the missing populations
      missPop = population[! population %in% ovSummary$pop]
      # Create a dataframe with 0 counts for the missing populations
      tempPop = data.frame(pop=missPop, ovlpNum=rep(0,length(missPop)))
      # Append it to the main dataframe
      ovSummary_temp = rbind(ovSummary, tempPop)
      # Convert factors to character string
      ovSummary_temp$pop = as.character(ovSummary_temp$pop)
      # Sort the dataframe alphabetically for pop IDs
      ovSummaryF = ovSummary_temp %>% arrange(pop)
      # Store the count columns in the 'p'th column of output dataframe
      ovlapDf[p] = ovSummaryF$ovlpNum
      } else {
        # Store the count columns in the 'p'th column of output dataframe
        ovlapDf[p] = ovSummary$ovlpNum
      }
    # Store the total number of windows for the 'p'th on diagonal
    ovlapDf[p,p] = nrow(posDf %>% filter(pop %in% p)) 
    
  }
  
  return(ovlapDf)
  
}