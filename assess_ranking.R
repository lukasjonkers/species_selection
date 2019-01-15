# explore how well abundance, nichewidth and temperature sensitivity explain species importance

library(plyr)
library(olsrr)
library(ggplot2)
library(reshape2)

dat <- readRDS('dat.RDS')

# maximum abundance (99th percentile)
max.abun <- lapply(dat.regions, function(x)
  apply(x$species, 2, function(y) quantile(y, 0.99))
)

# how well do species reflect temperature?
# use non linear least squares to fit a Gaussian distribution to species data
# and provide goodness of fit statistic (pseudo r2)
gfit <- function(x, y){
  out <- tryCatch({
    if(sum(y>1e-08)>1){
      df <- cbind.data.frame(x = x, y = y)
      # Crudely estimate starting values
      spl <- smooth.spline(x, y, spar = 1)
      m.0 <- spl$x[which.max(spl$y)]
      s.0 <- (max(x[y>1e-8])-min(x[y>1e-8]))/4
      b.0 <- 0  # always zero
      a.0 <- (max(y)-min(y))
      # fit the data
      # define the Gaussian function
      f <- function(theta){ 
        m <- theta[1];
        s <- theta[2];
        a <- theta[3];
        b <- theta[4];
        rhat <- a*exp(-0.5*((x-m)/s)^2) + b
        sum((y - rhat)^2)
      }
      fit <- optim(par = c(m.0, s.0, a.0, 0), fn = f, method="BFGS", control=list(reltol=1e-9))
      fity <- fit$par[3]*exp(-0.5*((df$x-fit$par[1])/fit$par[2])^2) + fit$par[4]
      # total sum of squares
      TSS <- sum((df$y - mean(df$y))^2)
      # residual sum of squares
      residuals <- df$y-fity
      RSS <- sum(residuals^2)
      1- RSS/TSS}else{NA}
      },
      error = function(cond){
        message('error')
        message(cond)
        message()
        return(NA)},
      finally = message('processed')
  )
  return(out)
}

fitr2 <- sapply(dat.regions, function(x){
  apply(x$species, 2, function(y) gfit(x$SST, y))
})

# get estimate of thermal niche width of species
nicheWidth <- lapply(dat.regions, function(x){
  presence <- x$species >= 1e-8
  apply(presence, 2, function(y) quantile(x$SST[y], 0.975) - quantile(x$SST[y], 0.025))
})

# summary df
sumdf <- mapply(cbind.data.frame, max.abun, fitr2, nicheWidth, SIMPLIFY = FALSE)
sumdf <- lapply(sumdf, setNames, nm = c('maxAbundance', 'pseudoR2', 'nicheWidth'))
sumdf <- lapply(sumdf, function(x){
  x$rn = row.names(x)
  x}
  )

rank.WA <- readRDS('rank_WA.RDS') # determined using taxon_ranking.R
rank.WA <- lapply(rank.WA, function(x)
  cbind.data.frame(impWA = x[,1], rn = row.names(x))
)
rank.MAT <- readRDS('rank_MAT.RDS') # determined using taxon_ranking.R
rank.MAT <- lapply(rank.MAT, function(x)
  cbind.data.frame(impMAT = x[,1], rn = row.names(x))
)


joinsum <- mapply(function(x, y, z) join_all(list(x, y, z), by = 'rn'), x = sumdf, y = rank.WA, z = rank.MAT, SIMPLIFY = FALSE)
meltsum <- melt(joinsum, id.vars = c('maxAbundance', 'pseudoR2', 'nicheWidth', 'impWA', 'impMAT'))

saveRDS(meltsum, '~/Dropbox/projects_ongoing/TF_reduction/assess_ranking.RDS')

# correlation between temperature sensitivity and abundance and temperature sensitivity and thermal niche width
round(cor(subset(meltsum, L1 == 'NAT')$maxAbundance, subset(meltsum, L1 == 'NAT')$pseudoR2), 2)
round(cor(subset(meltsum, L1 == 'NAT')$nicheWidth, subset(meltsum, L1 == 'NAT')$pseudoR2), 2)

# (multiple) linear regression to assess species importance (Table 1)
model <- lm(impWA~ maxAbundance + nicheWidth + pseudoR2, data = joinsum$NAT)
summary(model)
ols_step_both_p(model)
