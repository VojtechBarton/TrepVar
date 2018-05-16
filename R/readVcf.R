readVcf <- function(filename) {
  # find first data line (to skip file header)
  n_start <- 0

  con <- file(filename, "r")
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    if (substr(line, 1, 1) != "#") {
      break
    }
    n_start <- n_start + 1
  }
  close(con)

  # import data
  data <- read.table(filename, sep = "\t", skip = n_start)
  data <- data[,-c(1, 3, 6, 7, 9)] #leave only usefull columns
  data[, 3] <-
    gsub(x = as.character(data[, 3]),
         pattern = "<\\*>",
         replacement = "N")

  data <- separate(data, 5, into = c("1", "2", "3", "4"), sep = ":")
  data <- data[, -5]

  colnames(data) <-
    c("Position",
      "ReferenceAlels",
      "AlternativeAlels",
      "INFO",
      "ADF",
      "ADR",
      "AD")

  data <- data[, -4]

  cc <- which(as.character(data$ReferenceAlels) == "N")
  if (length(cc)){
    data <- data[-cc, ]
  }

  return(data)

}
