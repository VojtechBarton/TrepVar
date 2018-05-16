findVariableSpots <- function(data){

  # coverage depth
  cc <- which(data$Depth > 50)
  data <- data[cc,]

  cc <- which(data$Depth_ref_perc < 0.99)
  data <- data[cc,]

  cc <- which(data$Depth_Alt > 8)
  data <- data[cc,]

  cc <- which(data$Alt_reads_ratio > 0.4)
  data <- data[cc,]

  cc <- which(data$Alt_reads_ratio < 2.3)
  data <- data[cc,]

  return(data)

}
