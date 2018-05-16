alignToRef <- function(data, refSequence, poolsIndices) {
  pools_first_index <-
    c(
      as.numeric(gregexpr(poolsIndices[1], refSequence)[[1]]),
      as.numeric(gregexpr(poolsIndices[3], refSequence)[[1]]),
      as.numeric(gregexpr(poolsIndices[5], refSequence)[[1]]),
      as.numeric(gregexpr(poolsIndices[7], refSequence)[[1]]),
      as.numeric(gregexpr(poolsIndices[9], refSequence)[[1]]),
      as.numeric(gregexpr(poolsIndices[11], refSequence)[[1]])
    )
  genom <- data.frame()
  for (i in 1:length(data)) {
    data[[i]]$Position <- data[[i]]$Position + pools_first_index[i] - 1
    cc <- which(data[[i]]$Position > length(refSequence))
    if (length(cc) > 0) {
      data[[i]]$Position[cc] <-
        data[[i]]$Position[cc] - length(refSequence)
    }
    genom <- rbind(genom, data[[i]])
  }
  genom <- genom[order(genom$Position),]

  return(genom)
}
