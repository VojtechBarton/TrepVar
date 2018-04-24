pool1 <- "pool1_stat.vcf"
pool2 <- "pool2_stat.vcf"
pool3 <- "pool3_stat.vcf"
pool4 <- "pool4_stat.vcf"
TprE <- "poolE_stat.vcf"
TprC <- "poolC_stat.vcf"
ref <- readDNAStringSet("ref.fasta")[[1]]
pos_pool <- c(as.numeric(gregexpr(POOL_ind[1],ref)[[1]]), as.numeric(gregexpr(POOL_ind[3],ref)[[1]]), as.numeric(gregexpr(POOL_ind[5],ref)[[1]]), 
              as.numeric(gregexpr(POOL_ind[7],ref)[[1]]), as.numeric(gregexpr(POOL_ind[9],ref)[[1]]), as.numeric(gregexpr(POOL_ind[11],ref)[[1]]))


library(data.table)

pool1 <- readVcf(pool1)
pool2 <- readVcf(pool2)
pool3 <- readVcf(pool3)
pool4 <- readVcf(pool4)
TprE <- readVcf(TprE)
TprC<- readVcf(TprC)

data <- list(pool1, pool2, pool3, pool4, TprE, TprC)

genom <- mergeToRef(data, length(ref), pos_pool)

genom_j <- JoinRows(genom)
err <- genom_j[[2]] # chybové øádky - èasto tprK a nulové - kontrola
genom_j <- genom_j[[1]]

# Dominant count, frequency; Less count, frequency
data <- genom_j

data$DC <- data$RC
#data$DC[which(data$DC < data$AC1)] <- data$AC1[which(data$DC < data$AC1)]
data$DF <- data$DC / data$Depth
data$LC <- data$AC1
#data$LC[which(data$LC > data$RC)] <- data$RC[which(data$LC > data$RC)]
data$LF <- data$LC / data$Depth

# filtrace
cc <- which(data$Depth < 50)
data <- data[-cc,]

cc <- which(data$DF<0.99)
data <- data[cc,]

cc <- which(data$LC>=6)
data <- data[cc,]

save(genom, genom_j, data, file = "data.RData")

thresh <- data$LF
plot(thresh)

cc <- which(data$LF>= 0.025)
data <- data[cc,]

pos <- data$POS

AltAl <- substr(data$ALT,1,1)

x <- getAreaInfo(pos,AltAl,CDS)

# odstran tprK
cc <- which(x$gene == "tprK")
x <- x[-cc,]
data <- data[-cc,]
AltAl <- AltAl[-cc]

x$AltAl <- AltAl
x$syn <- FALSE
cc <- which(as.character(x$AA.ref) == as.character(x$AA.alt))
x$syn[cc] <- TRUE
 

Soubor <- cbind(data,x)
Soubor$transice <- FALSE
Soubor$DomAl <- Soubor$REF
Soubor$DomAl[which(Soubor$RC<Soubor$AC1)] <- Soubor$AltAl[which(Soubor$RC<Soubor$AC1)]

Soubor$transice[which((Soubor$DomAl=="C") & (Soubor$AltAl=="T"))] <- TRUE
Soubor$transice[which((Soubor$DomAl=="T") & (Soubor$AltAl=="C"))] <- TRUE

Soubor$transice[which((Soubor$DomAl=="A") & (Soubor$AltAl=="G"))] <- TRUE
Soubor$transice[which((Soubor$DomAl=="G") & (Soubor$AltAl=="A"))] <- TRUE



