exportIndels <- function(data){
  refs <- data$ReferenceAlels
  refs <- separate(data.frame(refs),1,into = c("1","2"), sep = 1)
  cc <- which(refs[,2] != "")

  indels <- data[cc,]
  data <- data[-cc,]

  return(list(data,indels))
}
