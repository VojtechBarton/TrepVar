joinRows <- function(data) {
  cond <- c(F, (data[-nrow(data), 1] == data[-1, 1]))
  referenceError <- c()
  for (i in (1:nrow(data))[cond]) {
    if (data$ReferenceAlels[i] == data$ReferenceAlels[i - 1]) {
      if ((data$AlternativeAlels[i] == data$AlternativeAlels[i - 1])) {
        ADF1 <- as.numeric(strsplit(data$ADF[i], split = ",")[[1]])
        ADF2 <-
          as.numeric(strsplit(data$ADF[i - 1], split = ",")[[1]])

        ADR1 <- as.numeric(strsplit(data$ADR[i], split = ",")[[1]])
        ADR2 <-
          as.numeric(strsplit(data$ADR[i - 1], split = ",")[[1]])

        AD1 <- as.numeric(strsplit(data$AD[i], split = ",")[[1]])
        AD2 <-
          as.numeric(strsplit(data$AD[i - 1], split = ",")[[1]])

        data$ADF[i - 1] <-
          paste(as.character(ADF1 + ADF2), collapse = ",")
        data$ADR[i - 1] <-
          paste(as.character(ADR1 + ADR2), collapse = ",")
        data$AD[i - 1] <-
          paste(as.character(AD1 + AD2), collapse = ",")

      } else {
        n_Alts <- c(0, 0, 0, 0)
        f_Alts <- c(0, 0, 0, 0)
        r_Alts <- c(0, 0, 0, 0)



        ADF1 <- as.numeric(strsplit(data$ADF[i], split = ",")[[1]])
        ADF2 <-
          as.numeric(strsplit(data$ADF[i - 1], split = ",")[[1]])

        ADR1 <- as.numeric(strsplit(data$ADR[i], split = ",")[[1]])
        ADR2 <-
          as.numeric(strsplit(data$ADR[i - 1], split = ",")[[1]])

        AD1 <- as.numeric(strsplit(data$AD[i], split = ",")[[1]])
        AD2 <-
          as.numeric(strsplit(data$AD[i - 1], split = ",")[[1]])

        Alt1 <- strsplit(data$AlternativeAlels[i], split = ",")[[1]]
        Alt2 <-
          strsplit(data$AlternativeAlels[i - 1], split = ",")[[1]]

        for (j in 1:length(Alt1)) {
          x <- Alt1[j]
          if (x == "A") {
            n_Alts[1] <- n_Alts[1] + AD1[j + 1]
            r_Alts[1] <- r_Alts[1] + ADR1[j + 1]
            f_Alts[1] <- f_Alts[1] + ADF1[j + 1]
          }
          if (x == "C") {
            n_Alts[2] <- n_Alts[2] + AD1[j + 1]
            r_Alts[2] <- r_Alts[2] + ADR1[j + 1]
            f_Alts[2] <- f_Alts[2] + ADF1[j + 1]
          }
          if (x == "G") {
            n_Alts[3] <- n_Alts[3] + AD1[j + 1]
            r_Alts[3] <- r_Alts[3] + ADR1[j + 1]
            f_Alts[3] <- f_Alts[3] + ADF1[j + 1]
          }
          if (x == "T") {
            n_Alts[4] <- n_Alts[4] + AD1[j + 1]
            r_Alts[4] <- r_Alts[4] + ADR1[j + 1]
            f_Alts[4] <- f_Alts[4] + ADF1[j + 1]
          }

        }

        for (j in 1:length(Alt2)) {
          x <- Alt2[j]
          if (x == "A") {
            n_Alts[1] <- n_Alts[1] + AD2[j + 1]
            r_Alts[1] <- r_Alts[1] + ADR2[j + 1]
            f_Alts[1] <- f_Alts[1] + ADF2[j + 1]
          }
          if (x == "C") {
            n_Alts[2] <- n_Alts[2] + AD2[j + 1]
            r_Alts[2] <- r_Alts[2] + ADR2[j + 1]
            f_Alts[2] <- f_Alts[2] + ADF2[j + 1]
          }
          if (x == "G") {
            n_Alts[3] <- n_Alts[3] + AD2[j + 1]
            r_Alts[3] <- r_Alts[3] + ADR2[j + 1]
            f_Alts[3] <- f_Alts[3] + ADF2[j + 1]
          }
          if (x == "T") {
            n_Alts[4] <- n_Alts[4] + AD2[j + 1]
            r_Alts[4] <- r_Alts[4] + ADR2[j + 1]
            f_Alts[4] <- f_Alts[4] + ADF2[j + 1]
          }

        }

        aleles <- c("A", "C", "G", "T")
        not_zero <- which(n_Alts != 0)
        order <- order(n_Alts, decreasing = TRUE)

        I <- order[which(order %in% not_zero)]

        data$ADF[i - 1] <-
          paste(c(
            as.character(ADF1[1] + ADF2[1]),
            as.character(f_Alts[I]) ,
            "0"
          ), collapse = ",")
        data$ADR[i - 1] <-
          paste(c(
            as.character(ADR1[1] + ADR2[1]),
            as.character(r_Alts[I]) ,
            "0"
          ), collapse = ",")
        data$AD[i - 1] <-
          paste(c(as.character(AD1[1] + AD2[1]), as.character(n_Alts[I]) , "0"), collapse = ",")
        data$AlternativeAlels[i - 1] <-
          paste(c(aleles[I], "N"), collapse = ",")

      }

    } else {
      referenceError <- c(referenceError, i)
    }
  }

  data <- data[-which(cond == T),]
  return(list(data, referenceError))

}
