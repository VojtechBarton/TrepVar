# uprava z csv
readVcf <- function(filename) {
  ## nalezeni zacatku importu
  I_start <- 0
  
  con = file(filename, "r")
  while (TRUE) {
    line = readLines(con, n = 1) #precteni jednoho radku
    if (length(line) == 0) {
      #nalezen konec souboru
      break
    }
    
    if (substr(line, 1, 1) != "#") {
      break
    }
    
    I_start <- I_start + 1
    
  }
  close(con)
  
  ## import dat
  data <- read.table(filename, sep = "\t", skip = I_start)
  data <- data[, -c(3, 6, 7, 8, 9)]
  names(data) <- c("CHROM", "POS", "REF", "ALT", "Depth")
  x <- tstrsplit(data$Depth, split = ":")
  data$Depth <- as.numeric(x[[2]])
  x <- tstrsplit(x[[3]], split = ",")
  data$RC <- as.numeric(x[[1]])
  data$AC1 <- as.numeric(x[[2]])
  data$AC2 <- as.numeric(x[[3]])
  data$AC3 <- as.numeric(x[[4]])
  data$AC4 <- 0
  
  cc <- which(as.character(data$REF) == "N")
  if (length(cc)>0) {
    data <- data[-cc, ]
  }
  
  data$RC[is.na(data$RC)] <- 0
  data$AC1[is.na(data$AC1)] <- 0
  data$AC2[is.na(data$AC2)] <- 0
  data$AC3[is.na(data$AC3)] <- 0
  
  data$ALT <- gsub('<\\*>', "N", as.character(data$ALT))
  
  # Dominant count, frequency; Less count, frequency
  
  # data$DC <- data$RC
  # data$DC[which(data$DC < data$AC1)] <- data$AC1[which(data$DC < data$AC1)]
  # data$DF <- data$DC / data$Depth
  # data$LC <- data$AC1
  # data$LC[which(data$LC > data$RC)] <- data$RC[which(data$LC > data$RC)]
  # data$LF <- data$LC / data$Depth
  
  return(data)
}


## mergeToRef namapuj data z readvcf k referenci
mergeToRef <- function(data, Lref, pos){
  # data  ... list of tables from readvcf()
  # Lref   ... delka reference
  # pos   ... pozice zacatku poolu v referenci
  
  genom <- data.frame()
  
  for(i in 1:length(data)){
      data[[i]]$POS <- data[[i]]$POS + pos[i]-1
      cc <- which(data[[i]]$POS>Lref)
      if (length(cc)>0){
        data[[i]]$POS[cc] <- data[[i]]$POS[cc] - Lref
      }
      genom <- rbind(genom,data[[i]])
  }
  
  genom <- genom[order(genom$POS),]
  
  return(genom)
}


## JoinRows se stejnou pozici v referenci, konsolidace prekryvu
JoinRows <- function(genom){
  cond <- c(F, (genom[-nrow(genom),2]==genom[-1,2]))
  ref_err <- c()
  for (i in (1:nrow(genom))[cond]){
    if (genom$REF[i] == genom$REF[i-1]){
      genom$Depth[i-1] <- genom$Depth[i-1]+genom$Depth[i]
      genom$RC[i-1] <- genom$RC[i-1]+genom$RC[i]
      if (genom$ALT[i-1]=="N"){
        genom[i-1,] <- genom[i,]
      } else if ((genom$ALT[i-1]!="N") & (genom$ALT[i]!="N")){
        counts <- c(0,0,0,0,0) # A,C,G,T,N
        order <- c("A","C","G","T","N")
        for (k in 0:1){
          alts <- c(strsplit(genom$ALT[i-k],split = ",")[[1]], "A","C","G","T","N")
          alts_c <- c(genom[i-k,7],genom[i-k,8],genom[i-k,9],genom[i-k,10],0,0,0,0,0)
          counts[1] <- counts[1]+alts_c[which(alts=="A")[1]]
          counts[2] <- counts[2]+alts_c[which(alts=="C")[1]]
          counts[3] <- counts[3]+alts_c[which(alts=="G")[1]]
          counts[4] <- counts[4]+alts_c[which(alts=="T")[1]]
          counts[5] <- counts[5]+alts_c[which(alts=="N")[1]]
        }
        alts <- order[-which(counts==0)]
        counts <- counts[-which(counts==0)]
        genom$ALT[i-1] <- paste(alts[order(-counts)],collapse = ",")
        counts <- c(sort(counts,decreasing = T),0,0,0,0,0)
        genom[i-1,7] <- counts[1]
        genom[i-1,8] <- counts[2]
        genom[i-1,9] <- counts[3]
        genom[i-1,10] <- counts[4]
        
      }
      
    } else ref_err <- c(ref_err, i)
    
  } 
  genom <- genom[-which(cond==T),]
  return(list(genom,ref_err))
}


## funkce pro nalezeni zasaz het. oblasti

getAreaInfo <- function(pos, AltAl, CDS){
  atrib <- tstrsplit(as.character(names(CDS)), split = " \\[")
  atrib <- data.frame(id = atrib[[1]], gene = atrib[[2]], prot = atrib[[3]], pid = atrib[[4]], loc = atrib[[5]])
  atrib$gene <- tstrsplit(as.character(atrib$gene), split = "=")[[2]]
  atrib$gene <- gsub("\\]", "", atrib$gene)
  
  atrib$prot <- tstrsplit(as.character(atrib$prot), split = "=")[[2]]
  atrib$prot <- gsub("\\]", "", atrib$prot)
  
  atrib$pid <- tstrsplit(as.character(atrib$pid), split = "=")[[2]]
  atrib$pid <- gsub("\\]", "", atrib$pid)
  
  atrib$loc <- tstrsplit(as.character(atrib$loc), split = "=")[[2]]
  atrib$loc <- gsub("\\]", "", atrib$loc)
  
  x <- tstrsplit(atrib$loc, split = "\\.\\.")
  atrib$locA <- x[[1]]
  atrib$locB <- as.numeric(gsub(")","",x[[2]]))
  atrib$compl <- FALSE
  atrib$compl[grep("complement*", atrib$locA)] <- TRUE
  atrib$locA <- as.numeric(gsub("complement\\(", "", atrib$locA))
  
 # data <- data.frame(pos = NA, gene = NA, protein = NA, protein.ID = NA, location = NA, compl = NA,
  #                   codon.ref = NA, codon.alt = NA, AA.ref = NA, AA.alt = NA)
  data <- data.frame()
  sgc11 <- getGeneticCode("11")
  
  
  for (i in 1:length(pos)){
    cc <- which((atrib$locA<=pos[i]) & (atrib$locB>=pos[i]))
    if (length(cc)){
      for (j in cc){
        x <- atrib[j,]
        k <- pos[i]-x$locA+1
        l <- (k-1)%%3
        codon.ref <- CDS[[j]][(k-l):(k-l+2)]
        codon.alt <- codon.ref
        codon.alt[l+1] <- DNAString(AltAl[i])
        
        if (x$compl){
          AA.ref <- translate(reverseComplement(codon.ref) ,genetic.code = sgc11)
          AA.alt <- translate(reverseComplement(codon.alt) ,genetic.code = sgc11)
        } else {
          AA.ref <- translate(codon.ref ,genetic.code = sgc11)
          AA.alt <- translate(codon.alt ,genetic.code = sgc11)
        }
        
        data_add <- data.frame(pos = pos[i], gene = x$gene, protein = x$prot, protein.ID = x$pid, location = x$loc, compl = x$compl, 
                               codon.ref = as.character(codon.ref), codon.alt = as.character(codon.alt), AA.ref = as.character(AA.ref),
                               AA.alt = as.character(AA.alt))
        data <- rbind(data, data_add)
        
      }
      
    } else {
      data_add <- data.frame(pos = pos[i], gene = NA, protein = NA, protein.ID = NA, location = NA, compl = NA,
                         codon.ref = NA, codon.alt = NA, AA.ref = NA, AA.alt = NA)
      data <- rbind(data, data_add)
    }
      
      
    
  }
  
  return(data)
}
