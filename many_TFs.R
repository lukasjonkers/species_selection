# perform many TFs with increasing number of species
# species ranked by importance
# Lukas Jonkers, last edit 05 Dec 2017

rank.WA <- readRDS('rank_WA.RDS')
rank.MAT <- readRDS('rank_MAT.RDS')
dat <- readRDS('dat.RDS')

source('many_WA.R')

perf_WA <- lapply(names(rank.WA), manyWA)
names(perf_WA) <- names(rank.WA)
saveRDS(perf_WA, 'performance_WA.RDS')

source('many_MAT.R')

perf_MAT <- lapply(names(rank.MAT), manyMAT)
names(perf_MAT) <- names(rank.MAT)
saveRDS(perf_MAT, 'performance_MAT.RDS')

rm(list = ls())
