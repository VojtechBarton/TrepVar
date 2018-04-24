DivRefToPools <- function(ref, pools_coord){
  
  # nalezeni souradnic jednotlivych poolu
  L <- length(ref)
  coords <- c();
  for (i in 1:length(pools_coord)){
    coords[i] = as.numeric(gregexpr(pools_coord[[i]],ref)[[1]])[1]
  }
  coords <- matrix(coords,ncol = 2, byrow = T)
  
  # odstraneni tprC a tprE z reference poolu 1 a 2
  tprC <- ref[coords[5,1]:coords[5,2]]
  tprE <- ref[coords[6,1]:coords[6,2]]
  ref[coords[5,1]:coords[5,2]] <- DNAString('N')
  ref[coords[6,1]:coords[6,2]] <- DNAString('N')
  
  # rozdeleni reference na pooly
  pools <- DNAStringSet()
  for (i in 1:(nrow(coords)-2)){
    if (coords[i,1]>coords[i,2]){
      x <- DNAStringSet(c(ref[coords[i,1]:L], ref[1:coords[i,2]]))
      pools <- c(pools, x)
    } else {
      pools <- c(pools, DNAStringSet(ref[coords[i,1]:coords[i,2]]))
    }
  }
  pools <- c(pools, DNAStringSet(tprC), DNAStringSet(tprE))
  names(pools) <- c('pool1','pool2','pool3','pool4','poolC','poolE')
  
  writeXStringSet(pools[1],'ref_pool1.fa')
  writeXStringSet(pools[2],'ref_pool2.fa')
  writeXStringSet(pools[3],'ref_pool3.fa')
  writeXStringSet(pools[4],'ref_pool4.fa')
  writeXStringSet(pools[5],'ref_poolC.fa')
  writeXStringSet(pools[6],'ref_poolE.fa')
  
}