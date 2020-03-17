# determine taxon importance
# using WA and MAT
# Lukas Jonkers, last edit 05 Dec 2017

library(rioja)
library(foreach)
library(doMC)

# get data
dat <- readRDS('dat.RDS')

# for parallel computing. Perhaps good for randomMAT which is very slow
registerDoMC(cores=detectCores()-1)

# for WA
taxon.imp_func <- function(x, nrep, nTF){
  taxon.imp <- replicate(nrep, randomPTF(x$species, x$SST, fun = WA, nTF = nTF, do.parallel = TRUE), simplify = FALSE)
}

sum_func <- function(j){
  imp.sorted <- as.data.frame(sapply(j, function(x) x$VI[order(row.names(x$VI)), ]))
  row.names(imp.sorted) <- row.names(j[[1]]$VI)[order(row.names(j[[1]]$VI))]
  summarised <- data.frame(
    mean = apply(imp.sorted, 1, 'mean'),
    sd = apply(imp.sorted, 1, 'sd'))
  summarised <- summarised[order(	summarised$mean, decreasing = TRUE), ]
  summarised
}

imp.WA_NAT <- taxon.imp_func(x = dat$NAT, nrep = 10, nTF = 1000)
imp.WA_SAT <- taxon.imp_func(x = dat$SAT, nrep = 10, nTF = 1000)
imp.WA_MDX <- taxon.imp_func(x = dat$MDX, nrep = 10, nTF = 1000)
imp.WA_IND <- taxon.imp_func(x = dat$IND, nrep = 10, nTF = 1000)
imp.WA_PAC <- taxon.imp_func(x = dat$PAC, nrep = 10, nTF = 1000)

imp.WA <- list(NAT = imp.WA_NAT, SAT =  imp.WA_SAT, MDX = imp.WA_MDX, IND = imp.WA_IND, PAC = imp.WA_PAC)
saveRDS(imp.WA, 'imp_WA.RDS')

rank.WA_NAT <- sum_func(imp.WA_NAT)
rank.WA_SAT <- sum_func(imp.WA_SAT)
rank.WA_MDX <- sum_func(imp.WA_MDX)
rank.WA_IND <- sum_func(imp.WA_IND)
rank.WA_PAC <- sum_func(imp.WA_PAC)

rank.WA <- list(NAT = rank.WA_NAT, SAT =rank.WA_SAT, MDX = rank.WA_MDX, IND = rank.WA_IND, PAC = rank.WA_PAC)
saveRDS(rank.WA, 'rank_WA.RDS')

# for MAT, this is really slow.
source('randomMAT.R')
imp.MAT_NAT <- replicate(n = 10, randomMAT(dat$NAT$species, dat$NAT$SST, nTF = 1000, do.parallel = TRUE), simplify  = FALSE)
imp.MAT_SAT <- replicate(n = 10, randomMAT(dat$SAT$species, dat$SAT$SST, nTF = 1000, do.parallel = TRUE), simplify  = FALSE)
imp.MAT_MDX <- replicate(n = 10, randomMAT(dat$MDX$species, dat$MDX$SST, nTF = 1000, do.parallel = TRUE), simplify  = FALSE)
imp.MAT_IND <- replicate(n = 10, randomMAT(dat$IND$species, dat$IND$SST, nTF = 1000, do.parallel = TRUE), simplify  = FALSE)
imp.MAT_PAC <- replicate(n = 10, randomMAT(dat$PAC$species, dat$PAC$SST, nTF = 1000, do.parallel = TRUE), simplify  = FALSE)

imp.MAT <- list(NAT = imp.MAT_NAT, SAT =  imp.MAT_SAT, MDX = imp.MAT_MDX, IND = imp.MAT_IND, PAC = imp.MAT_PAC)
saveRDS(imp.MAT, 'imp_MAT.RDS')

rank.MAT_NAT <- sum_func(imp.MAT_NAT)
rank.MAT_SAT <- sum_func(imp.MAT_SAT)
rank.MAT_MDX <- sum_func(imp.MAT_MDX)
rank.MAT_IND <- sum_func(imp.MAT_IND)
rank.MAT_PAC <- sum_func(imp.MAT_PAC)

rank.MAT <- list(NAT = rank.MAT_NAT, SAT =rank.MAT_SAT, MDX = rank.MAT_MDX, IND = rank.MAT_IND, PAC = rank.MAT_PAC)
saveRDS(rank.MAT, 'rank_MAT.RDS')

# unregister cores
registerDoSEQ()

rm(list = ls())
