findFeatures <- function(data){
  Depth <- data$AD
  Depth <- separate(data.frame(Depth),1, into = c("1","2","3","4","5"), sep = ",", convert = T)
  Depth[is.na(Depth)] <- 0
  data$Depth_ref <- Depth[,1]
  data$Depth_Alt <- Depth[,2]
  data$Depth <- rowSums(Depth)
  data$Depth_ref_perc <- data$Depth_ref / data$Depth
  data$Depth_Alt_perc <- data$Depth_Alt / data$Depth

  Depth_Alt_f <- data$ADF
  Depth_Alt_f <- separate(data.frame(Depth_Alt_f),1, into = c("1","2","3","4","5"), sep = ",", convert = T)
  Depth_Alt_f[is.na(Depth_Alt_f)] <- 0
  data$Depth_Alt_f <- Depth_Alt_f[,2]

  Depth_Alt_r <- data$ADR
  Depth_Alt_r <- separate(data.frame(Depth_Alt_r),1, into = c("1","2","3","4","5"), sep = ",", convert = T)
  Depth_Alt_r[is.na(Depth_Alt_r)] <- 0
  data$Depth_Alt_r <- Depth_Alt_r[,2]

  data$Alt_reads_ratio <- data$Depth_Alt_f / data$Depth_Alt_r
  data$Alt_reads_ratio[data$Alt_reads_ratio==Inf] <- 0
  data$Alt_reads_ratio[is.nan(data$Alt_reads_ratio)] <- 0

  return(data)

}
