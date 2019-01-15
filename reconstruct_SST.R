# predict SST in downcore datasets

library(rioja)

coredat <- readRDS('fossil_data.RDS')
coredat <- fossil
corename <- 'M35003_4' # change for different core
domain <- 'NAT'

# fossil species data
fos <- as.data.frame(coredat[[which(names(coredat) == corename)]][, -c(1:2)]) # for cores
fos <- as.data.frame(coredat[[which(names(coredat) == corename)]][,-c(1:(which(names(coredat) == 'O_uni'))-1)]) # for MARGO

# training set data
dat <- readRDS('dat.RDS')

# species in training set
spp <- dat[[which(names(dat) == domain)]]$species

# assume species that are in the training set are 0 in core
# species not in down core
sppNotInCore <- names(spp)[!names(spp) %in% names(fos)]
# use Merge to combine spp and fos
species.all <- Merge(fos, spp, split = TRUE)

# SST training set
SST <- dat[[which(names(dat) == domain)]]$SST

# sample positions training set
sample.positions <- cbind(dat[[which(names(dat) == domain)]]$meta$Longitude, dat[[which(names(dat) == domain)]]$meta$Latitude)

# rank training set
method <- 'WA' # change for method
rank <- readRDS(paste0('rank_', method, '.RDS')) # obtained from taxon_raking.R
rank.fos <- row.names(rank[[which(names(rank) == domain)]])

# make TF, crossval and reconstruct SST
manyTFcore <- function(fos, spp, SST, sample.positions, rank.fos, method){ # x is name of domain
  require(rioja)
  require(fields)
  res <- matrix(NA, ncol = 5, nrow = length(rank.fos)-1)
  rec <- matrix(NA, ncol = length(rank.fos)-1, nrow = nrow(fos))
  resid <- list()
  predi <- list()
  nsite <- NA
  # prepare data
  for(i in 2:length(rank.fos)){
    spec <- spp[, which(names(spp) %in% rank.fos[1:i])]
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
    # do TF
    mod <- match.fun(method)(spec, env, lean = FALSE)
    #h.cutoff <- 850 # km, following Trachsel & Telford 2016, CP
    h.dist <- rdist.earth(xy, miles = FALSE)
    cv <- crossval(mod, cv.method = 'h-block', h.cutoff = 850, h.dist = h.dist)
    if(method == 'MAT'){
      stat <- performance(cv)$crossval[5,] # to get stats for the 5 closest analogues
    }else{
      stat <- performance(cv)$crossval[1,]
    }
    pred <- predict(mod, fos, k = 5)
    predi[[i-1]] <- pred
    rec[,i-1] <- pred$fit[,1]
    res[i-1,] <- stat
    resid[[i-1]] <- cbind.data.frame(SST = env, cv$residuals.cv)
    nsite[i-1] <- nrow(spec)
    }
  res <-as.data.frame(res)
  names(res) <- names(stat)
  reconstructed.SST <- as.data.frame(rec)
  list(res = res, residuals.cv = resid, nsite = nsite, reconstructed.SST = reconstructed.SST, out = predi)
}

out <- manyTFcore(species.all$fos, species.all$spp, SST, sample.positions, rank.fos, method)
out$coredat <- coredat
out$rank.fos <- rank.fos
saveRDS(out, file = paste0('corename, '_', method, '.RDS'))
rm(list = ls())
