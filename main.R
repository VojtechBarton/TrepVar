setwd("")
### load sequence indicies of the pools
poolsIndicies <- readDNAStringSet("pools_coord.fa")

### files with vcf
pool1 <- "pool1.vcf"
pool2 <- "pool2.vcf"
pool3 <- "pool3.vcf"
pool4 <- "pool4.vcf"
tprC <- "poolC.vcf"
tprE <- "poolE.vcf"

### reference sequence
refSeq <- readDNAStringSet("ref.fasta")[[1]]

### read vcf
pool1 <- readVcf(pool1)
pool2 <- readVcf(pool2)
pool3 <- readVcf(pool3)
pool4 <- readVcf(pool4)
tprC <- readVcf(tprC)
tprE <- readVcf(tprE)

### all pools in one
data <- list(pool1, pool2, pool3, pool4, tprC, tprE)

### create one sorted dataset from pools
genom <- alignToRef(data, refSeq, poolsIndicies)

### export indels
indels <- exportIndels(genom)

genom <- indels[[1]]
indels <- indels[[2]]

### join rows with same position (overlaping pools)
data <- joinRows(genom)

rows_error <- data[[2]]
genom_j <- data[[1]]

### find features
Genom <- findFeatures(genom_j)
Indels <- findFeatures(indels)


### export variable spots
VarSpots <- findVariableSpots(Genom)

thresh <- findThresh(Genom)

save(Genom, Indels, VarSpots, thresh, file = "final.RData")
