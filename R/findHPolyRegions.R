findHPolyRegions <- function(Sequence, n){
  PolyReg <- vector(length = length(Sequence))

  I_start <- c()
  I_length <- c()

  for (i in c("A","C","G","T")){
    I <- gregexpr(paste(c(i,"{", as.character(n), ",}"), collapse = ""), Sequence, ignore.case = T)

    I_start <- c(I_start, I[[1]])
    I_length <- c(I_length, attr(I[[1]],"match.length"))

  }

  for (i in 1:length(I_start)){
    PolyReg[ I_start[i] : (I_start[i]+(I_length[i]-1))] <- T
  }

  return(PolyReg)
}
