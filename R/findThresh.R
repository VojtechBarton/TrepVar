findThresh <- function(data){
  # count ratio of non reference alels per position
  data <- 1-data$Depth_ref_perc

  # create cumulative histogram
  y <- c()
  t <- seq(0,0.5,0.0001)
  for (i in t ) {
    a <- length(which(data>i))
    y <- c(y,a)
  }

  # linear model on histogram
  model <- lm(y ~ t)

  r <- as.numeric(model$residuals)

  # hard threshold - first upraising of the alt alleles ratio
  hard_thresh <- t[which.min(r)]

  # soft threshold - add tolerance of mean error from model
  std.error <- mean(abs(r))
  soft_thresh <- t[which(r < (min(r)+std.error))[1]]

  thresh <- c(hard_thresh, soft_thresh)

  return(thresh)

}
