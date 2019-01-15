# function to perform many MAT with increasing number of taxa ranked accoring to species importance
# in this function increases the size of the training set when more species are included
# Lukas Jonkers, last edit 14 Dec 2017

manyMAT <- function(x){ # is name of domain
  require(rioja)
  require(fields)
  domain <- dat[[which(names(dat) == x)]]
  SST <- domain$SST
  species <- domain$species
  sample.positions <- cbind(domain$meta$Longitude, domain$meta$Latitude)
  rank <- row.names(rank.MAT[[which(names(rank.MAT) == x)]])
  res <- matrix(NA, ncol = 5, nrow = length(rank)-1)
  resid <- list()
  nsite <- NA
  for(i in 2:length(rank)){
    spec <- species[, which(names(species) %in% rank[1:i])]
    #remove samples where species have zero abundance
    zero.indx <- which(rowSums(spec) < 1e-08)
    if(length(zero.indx)>0){
      spec <- spec[-zero.indx, ]
      env <- SST[-zero.indx]
      xy <- sample.positions[-zero.indx,]
    }else{
      env <- SST
      xy <- sample.positions
    }
    mod <- MAT(spec, env, lean = FALSE)
    h.dist <- rdist.earth(xy, miles = FALSE)
    cv <- crossval(mod, cv.method = 'h-block', h.cutoff = 850, h.dist = h.dist)
    stat <- performance(cv)$crossval[5,]
    res[i-1,] <- stat
    resid[[i-1]] <- cbind.data.frame(SST = env, cv$residuals.cv)
    nsite[i-1] <- nrow(spec)
  }
  res <-as.data.frame(res)
  names(res) <- names(stat)
  list(res = res, residuals.cv = resid, nsite = nsite)
}
