# figures for Sensitivity to species selection indicates the effect of nuisance variables on marine microfossil transfer functions
# Lukas Jonkers and Michal Kucera; https://doi.org/10.5194/cp-2018-107
# only for North Atlantic

library(ggplot2)
library(egg)
library(reshape2)
library(vegan)
library(rgdal)
library(plyr)
library(RColorBrewer)

# Figure 1: species importance and TF performance ####
sppNames <- readRDS('dat.RDS')$species_names
sppNames$mid_name <- sapply(strsplit(as.character(sppNames$original_name), '_'), function(x){
  paste(substr(x[1], 1, 1), x[2],sep =  '. ')
})

rank.WA <- readRDS('rank_WA.RDS')$NAT # obtained from taxon_ranking.R
RMSE.WA <- readRDS('performance_WA.RDS')$NAT$res # obtained from many_TF.R
minSp.WA <- min(which(RMSE.WA$RMSE < min(RMSE.WA$RMSE)*1.1))+1
impWA <- cbind.data.frame(rank.WA, rbind.data.frame(rep(NA, 5), RMSE.WA))
impWA$spp <- sapply(row.names(impWA), function(x) sppNames$mid_name[which(sppNames$short_name == x)])
impWA$spp <- factor(impWA$spp, levels = impWA$spp)
impWA$method <- 'WA'

rank.MAT <- readRDS('rank_MAT.RDS')$NAT
RMSE.MAT <- readRDS('performance_MAT.RDS')$NAT$res
minSp.MAT <- min(which(RMSE.MAT$RMSE < min(RMSE.MAT$RMSE)*1.1))+1
impMAT <- cbind.data.frame(rank.MAT, rbind.data.frame(rep(NA, 5), RMSE.MAT))
impMAT$spp <- sapply(row.names(impMAT), function(x) sppNames$mid_name[which(sppNames$short_name == x)])
impMAT$spp <- factor(impMAT$spp, levels = impMAT$spp)
impMAT$method <- 'MAT'

FigWArank <- ggplot(impWA) +
  geom_point(aes(spp, mean)) +
  geom_vline(xintercept = minSp.WA, colour = 'grey30') +
  geom_errorbar(data = impWA, aes(spp, ymin = mean-sd, ymax=mean+sd), width =.2, position = position_dodge(.9)) +
  geom_point(data = impWA, aes(spp, RMSE), colour = 'firebrick3') +
  ylab(expression('Importance | RMSE ['*degree*C*']')) +
  ggtitle('WA') +
  ylim(c(min(impWA$mean), 6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
        axis.title.x = element_blank(),
        text = element_text(size = 9),
        panel.grid = element_line(size = 0.25))

FigMATrank <- ggplot(impMAT) +
  geom_point(aes(spp, mean/5)) +
  geom_vline(xintercept = minSp.MAT, colour = 'grey30') +
  geom_errorbar(aes(spp, ymin = mean/5-sd, ymax=mean/5+sd), width =.2, position = position_dodge(.9)) +
  geom_point(aes(spp, RMSE), colour = 'firebrick3') +
  ylab(expression('Importance/5 | RMSE ['*degree*C*']')) +
  ggtitle('MAT') +
  ylim(c(min(impMAT$mean/5), 6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
        axis.title.x = element_blank(),
        text = element_text(size = 9),
        panel.grid = element_line(size = 0.25))

Fig1 <- ggarrange(FigMATrank, FigWArank, ncol = 2)

# Figure 2: species ranking and thermal niche ####
ForCenS <- readRDS('dat.RDS')
ForCenS_NAT <- cbind.data.frame(SST = ForCenS$NAT$SST, ForCenS$NAT$species)

# rank species by average abundance
ForCenS_NAT <- cbind.data.frame(SST = ForCenS$NAT$SST,
                                ForCenS$NAT$species[,order(colMeans(ForCenS$NAT$species)*100, decreasing = FALSE)])

# bin species abundance per degree
breaks <- seq(floor(min(ForCenS_NAT$SST)), ceiling(max(ForCenS_NAT$SST)), by = 1)
bin.spec <- ddply(ForCenS_NAT, .(cut(ForCenS_NAT$SST, breaks)), colwise(mean))
midpoints <- gsub('\\(', '', bin.spec[,1])
midpoints <- as.numeric(gsub("(.*),.*", "\\1", midpoints)) +0.5
# scale
binSpp_scaled <- cbind.data.frame(midpoints, sweep(bin.spec[,-1], 2, colSums(bin.spec[,-1]), "/")*100)
names(binSpp_scaled)[-1] <- sapply(names(binSpp_scaled[,-1]), function(x) sppNames$mid_name[which(sppNames$short_name == x)])
binSpp_scaled[binSpp_scaled < 1] <- NA

gplot_binSpp <- melt(binSpp_scaled, id.vars = 'midpoints')

Fig2 <- ggplot(gplot_binSpp, aes(midpoints, variable, fill = value)) +
  geom_tile() +
  xlab(expression('Temperature ['*degree*C*']')) +
  #scale_fill_gradientn(colours = brewer.pal(9, name = 'YlOrRd'), na.value = 'transparent', name = 'relative\nabundance\n[log10(%)]') +
  scale_fill_gradientn(trans = "log",
                       colours = brewer.pal(9, name = 'YlOrRd'),
                       breaks = c(2, 10, 25, 50, 100),
                       na.value = 'transparent', name = 'relative\nabundance\n[%]') +
  theme_bw() +
  theme(axis.text.y = element_text(face = 'italic'),
        axis.title.y = element_blank(),
        text = element_text(size = 9),
        panel.grid = element_line(size = 0.25))

# Figure 3: species abundance, niche width and temperature sensitivity correlation
spp_sum <- readRDS('assess_ranking.RDS') # obtained from assess_ranking.R

ab_sensi <- ggplot(subset(spp_sum, L1 == 'NAT'), aes(maxAbundance, pseudoR2)) +
  geom_smooth(method = 'lm', colour = 'firebrick3') +
  geom_point(alpha = 0.8)+
  ylab(expression('Temperature sensitivity\n [pseudo '*r^{2}*']')) +
  xlab(expression('Max. relative abundance')) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(aspect.ratio=1,
        text = element_text(size = 9),
        panel.grid = element_line(size = 0.25))

width_sensi <- ggplot(subset(meltsum, L1 == 'NAT'), aes(nicheWidth, pseudoR2)) +
  geom_smooth(method = 'lm', colour = 'firebrick3') +
  geom_point(alpha = 0.8)+
  ylab(expression('Temperature sensitivity\n [pseudo '*r^{2}*']')) +
  xlab(expression('Thermal niche width ['*degree*C*']')) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(aspect.ratio=1,
        text = element_text(size = 9),
        panel.grid = element_line(size = 0.25))


Fig3 <- grid.arrange(ab_sensi, width_sensi, layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 2))

# Figure 4: plot different reconstructions ####
# load reconstruction data for cores
MD95_WArec <- readRDS('~/Dropbox/projects_ongoing/TF_reduction/repeat_runs/MD95-2040_WA.RDS')
MD95_MATrec <- readRDS('~/Dropbox/projects_ongoing/TF_reduction/repeat_runs/MD95-2040_MAT.RDS')
M30_WArec <- readRDS('~/Dropbox/projects_ongoing/TF_reduction/repeat_runs/M35003-4_WA.RDS')
M30_MATrec <- readRDS('~/Dropbox/projects_ongoing/TF_reduction/repeat_runs/M35003-4_MAT.RDS')
MARGO_WArec <- readRDS('~/Dropbox/projects_ongoing/TF_reduction/repeat_runs/MARGO_WA.RDS')
MARGO_MATrec <- readRDS('~/Dropbox/projects_ongoing/TF_reduction/repeat_runs/MARGO_MAT.RDS')

# initially for just one core to show that there are different reconstructions
# Iberian Margin core MD95-2040
WAplot <- melt(cbind.data.frame(age = MD95_WArec$coredat$age_ka_BP, MD95_WArec$reconstructed.SST), id.vars = 'age')
WAplot$cat <- c(rep('sub', (minSp.WA-2)*length(MD95_WArec$coredat$age_ka_BP)),
                rep('min', length(MD95_WArec$coredat$age_ka_BP)),
                rep('sup', (length(MD95_WArec$reconstructed.SST)-minSp.WA)*length(MD95_WArec$coredat$age_ka_BP)), 
                rep('all', length(MD95_WArec$coredat$age_ka_BP)))

MATplot <- melt(cbind.data.frame(age = MD95_MATrec$coredat$age_ka_BP, MD95_MATrec$reconstructed.SST), id.vars = 'age')
MATplot$cat <- c(rep('sub', (minSp.MAT-2)*length(MD95_MATrec$coredat$age_ka_BP)),
                rep('min', length(MD95_MATrec$coredat$age_ka_BP)),
                rep('sup', (length(MD95_MATrec$reconstructed.SST)-minSp.MAT)*length(MD95_MATrec$coredat$age_ka_BP)), 
                rep('all', length(MD95_MATrec$coredat$age_ka_BP)))

FigWArec <- ggplot() +
  #geom_line(data = subset(WAplot, cat == 'sub'), aes(age, value, group = variable), colour = 'palevioletred1', size = 0.2) +
  geom_line(data = subset(WAplot, cat == 'sup'), aes(age, value, group = variable), colour = 'grey50') +
  geom_line(data = subset(WAplot, cat == 'min'), aes(age, value, group = variable), colour = 'firebrick3') +
  geom_line(data = subset(WAplot, cat == 'all'), aes(age, value, group = variable), colour = 'black') +
  ylim(0, 18) +
  ylab(expression('SST ['*degree*C*']')) +
  xlab('age [ka BP]') +
  ggtitle('WA') +
  theme_bw() +
  theme(text = element_text(size = 9),
        panel.grid = element_line(size = 0.25))


FigMATrec <- ggplot() +
  #geom_line(data = subset(MATplot, cat == 'sub'), aes(age, value, group = variable), colour = 'palevioletred1', size = 0.2) +
  geom_line(data = subset(MATplot, cat == 'sup'), aes(age, value, group = variable), colour = 'grey50') +
  geom_line(data = subset(MATplot, cat == 'min'), aes(age, value, group = variable), colour = 'firebrick3') +
  geom_line(data = subset(MATplot, cat == 'all'), aes(age, value, group = variable), colour = 'black') +
  ylim(0, 18) +
  ylab(expression('SST ['*degree*C*']')) +
  xlab('age [ka BP]') +
  ggtitle('MAT') +
  #scale_colour_manual(values = c(rep('grey80', ncol(MD95_MATrec$reconstructed.SST) -1), 'black'), guide = 'none') +
  theme_bw() +
  theme(text = element_text(size = 9),
        panel.grid = element_line(size = 0.25))
  
# plot MARGO reconstructions
# only map of difference between reconstruction with min species and with all species
# read the shapefile for the simple worldmap
# wmap <- readOGR(dsn = 'ne_110m_land', layer = 'ne_110m_land')

# function to prepare map data
prep.map <- function(x, nsp, method){
  LGM <- cbind.data.frame(core = x$coredat$Core,
                          lat =  x$coredat$Latitude,
                          lon = x$coredat$Longitude,
                          x$reconstructed.SST)
  # average cores
  mean.LGM <- ddply(LGM, .(core), function(x) apply(x[,-1], 2, function(k) mean(k, na.rm = TRUE)))
  delta <- mean.LGM[, nsp+2] - mean.LGM[, ncol(mean.LGM)]
  cbind(mean.LGM[,c(1:3)], delta, method)
}

plotMARGO <- rbind(prep.map(MARGO_MATrec, minSp.WA-1, 'MAT'), prep.map(MARGO_WArec, minSp.WA-1, 'WA'))

FigMARGO <- ggplot(plotMARGO) +
  #geom_polygon(data = wmap, aes(long,lat, group=group, fill=hole)) + 
  geom_point(aes(lon, lat, colour = delta), size = 2) +
  facet_wrap(~method) +
  scale_colour_gradient2(low = 'dodgerblue', high = 'firebrick3', mid = 'lemonchiffon', name = expression('Temperature offset ['*degree*C*']')) +
  coord_equal() +
  coord_cartesian(xlim = c(-100, 50),ylim = c(0, 85)) +
  scale_fill_manual(values=c("grey80", "white"), guide="none") +
  theme_bw() +
  theme(text = element_text(size = 9),
        panel.grid = element_line(size = 0.25),
        legend.position = 'bottom')

Fig4 <- grid.arrange(FigMATrec, FigWArec, FigMARGO, layout_matrix = matrix(c(rep(c(1, 2), 2), rep(3, 6)), byrow = TRUE, ncol = 2))

# Figure 5: difference from final reconstruction vs species ranked by importance ####
# function to extract mean difference from final reconstruction
getdiff <- function(x, core, method, nsp){
  #spp <- factor(x$rank.fos, levels = x$rank.fos)
  spp <- sapply(x$rank.fos, function(x) sppNames$mid_name[which(sppNames$short_name == x)])
  spp <- factor(spp, levels = spp)
  deltaT <- x$reconstructed.SST - x$reconstructed.SST[, ncol(x$reconstructed.SST)]
  meanAbsDeltaT <- c(NA, apply(deltaT, 2, function(x) mean(abs(x), na.rm = TRUE)))
  medianAbsDeltaT <- c(NA, apply(deltaT, 2, function(x) median(abs(x), na.rm = TRUE)))
  meanDeltaT <- c(NA, apply(deltaT, 2, function(x) mean(x, na.rm = TRUE)))
  medianDeltaT <- c(NA, apply(deltaT, 2, function(x) median(x, na.rm = TRUE)))
  #impor <- importance$mean[row.names(importance) %in% spp]
  err <- c(NA, x$res$RMSE)
  foo <- cbind.data.frame(spp = spp, error = err, meanDeltaT = meanDeltaT, medianDeltaT = medianDeltaT,  medianAbsDeltaT =  medianAbsDeltaT,  meanAbsDeltaT =  meanAbsDeltaT)
  foo$core <- core
  foo$method <- method
  foo$colour <- factor(c(rep(1, nsp-1), 2, rep(3, length(spp)-nsp)))
  foo
}
# apply to core data
diff_MD95_WA <- getdiff(MD95_WArec, 'MD95-2040', 'WA', minSp.WA)
diff_MD95_MAT <- getdiff(MD95_MATrec, 'MD95-2040', 'MAT', minSp.MAT)
diff_M30_WA <- getdiff(M30_WArec, 'M35003-4', 'WA', minSp.WA)
diff_M30_MAT <- getdiff(M30_MATrec, 'M35003-4', 'MAT', minSp.MAT)
diff_MARGO_WA <- getdiff(MARGO_WArec, 'MARGO LGM', 'WA', minSp.WA)
diff_MARGO_MAT <- getdiff(MARGO_MATrec, 'MARGO LGM', 'MAT', minSp.MAT)

gplotdiffMAT <- rbind.data.frame(diff_MD95_MAT, diff_M30_MAT, diff_MARGO_MAT)
gplotdiffMAT$core <- factor(gplotdiffMAT$core, levels = unique(gplotdiffMAT$core))
gplotdiffWA <- rbind.data.frame(diff_MD95_WA, diff_M30_WA, diff_MARGO_WA)
gplotdiffWA$core <- factor(gplotdiffWA$core, levels = unique(gplotdiffWA$core))

Fig5 <- ggarrange(ggplot(gplotdiffMAT) +
                    #geom_point(data = impMAT[-1,], aes(spp, mean, colour = 'firebrick3', alpha = 0.5)) +
                    geom_point(aes(spp, meanDeltaT, colour = colour), alpha = 0.5, shape = 18, size = 2) +
                    geom_point(aes(spp, meanAbsDeltaT, colour = colour)) +
                    facet_grid(factor(core)~.) +
                    theme_bw() +
                    ylim(c(-4.5, 6.5)) +
                    ylab(expression('Temperature offset ['*degree*C*']')) +
                    ggtitle('MAT') +
                    scale_colour_manual(values = c('grey50', 'firebrick3', 'black'), guide = 'none') +
                    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'), 
                          axis.title.x = element_blank(),
                          text = element_text(size = 9),
                          panel.grid = element_line(size = 0.25)),
                ggplot(gplotdiffWA) +
                  geom_point(aes(spp, meanDeltaT, colour = colour), alpha = 0.5, shape = 18, size = 2) +
                  geom_point(aes(spp, meanAbsDeltaT, colour = colour)) +
                  facet_grid(core~.) +
                  theme_bw() +
                  ylim(c(-4.5, 6.5)) +
                  #ylab(expression('Temperature offset ['*degree*C*']')) +
                  ggtitle('WA') +
                  scale_colour_manual(values = c('grey50', 'firebrick3', 'black'), guide = 'none') +
                  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    text = element_text(size = 9),
                    panel.grid = element_line(size = 0.25)),
  ncol = 2
)

# Figure 6: difference from final reconstruction vs transfer function performance ####
# function to plot

gplotdiff <- rbind.data.frame(gplotdiffWA, gplotdiffMAT)

Fig6 <- ggplot(gplotdiff, aes(error, meanAbsDeltaT, colour = core)) +
  geom_path(alpha = 0.5) +
  geom_point(aes(shape = colour)) +
  xlim(c(1, 5)) +
  ylab(expression('Average difference from\nreconstruction with all species ['*degree*C*']')) +
  xlab(expression('Prediction error ['*degree*C*']')) +
  scale_colour_manual(values = c('firebrick3', 'royalblue3', 'darkgoldenrod3')) +
  scale_shape_manual(values = c(1, 8, 16), guide = 'none') +
  scale_alpha_manual(values = c(0.5, 1, 0.5), guide = 'none') +
  theme_bw() +
  coord_fixed() +
  facet_grid(.~method) +
  theme(text = element_text(size = 9),
        panel.grid = element_line(size = 0.25))

# Figure 8: sensitivity to species pruning vs analogue quality and richness ####
minWA <- minSp.WA-1
minMAT <- minSp.MAT-1
sens_MD95_WA <- cbind.data.frame(sd = apply(MD95_WArec$reconstructed.SST[, -(1:minWA)], 1, 'sd'),
                                 richness = rowSums(MD95_WArec$coredat[-c(1, 2)] > 0),
                                 sqcd = rowMeans(MD95_MATrec$out[[length(MD95_MATrec$out)]]$dist.n),
                                 core = 'MD95-2040',
                                 method = 'WA')

sens_MD95_MAT <-cbind.data.frame(sd = apply(MD95_MATrec$reconstructed.SST[, -(1:minMAT)], 1, 'sd'),
                                 richness = rowSums(MD95_MATrec$coredat[-c(1, 2)] > 0),
                                 sqcd = rowMeans(MD95_MATrec$out[[length(MD95_MATrec$out)]]$dist.n),
                                 core = 'MD95-2040',
                                 method = 'MAT')
sens_M30_WA <- cbind.data.frame(sd = apply(M30_WArec$reconstructed.SST[, -(1:minWA)], 1, 'sd'),
                                 richness = rowSums(M30_WArec$coredat[-c(1, 2)] > 0),
                                 sqcd = rowMeans(M30_MATrec$out[[length(M30_MATrec$out)]]$dist.n),
                                 core = 'M35003-4',
                                 method = 'WA')
sens_M30_MAT <-cbind.data.frame(sd = apply(M30_MATrec$reconstructed.SST[, -(1:minMAT)], 1, 'sd'),
                                 richness = rowSums(M30_MATrec$coredat[-c(1, 2)] > 0),
                                 sqcd = rowMeans(M30_MATrec$out[[length(M30_MATrec$out)]]$dist.n),
                                 core = 'M35003-4',
                                 method = 'MAT')
sens_MARGO_WA <- cbind.data.frame(sd = apply(MARGO_WArec$reconstructed.SST[, -(1:minWA)], 1, 'sd'),
                                richness = rowSums(MARGO_WArec$coredat[-c(1:10)] > 0),
                                sqcd = rowMeans(MARGO_MATrec$out[[length(MARGO_MATrec$out)]]$dist.n),
                                core = 'MARGO LGM',
                                method = 'WA')
sens_MARGO_MAT <-cbind.data.frame(sd = apply(MARGO_MATrec$reconstructed.SST[, -(1:minMAT)], 1, 'sd'),
                                richness = rowSums(MARGO_MATrec$coredat[-c(1:10)] > 0),
                                sqcd = rowMeans(MARGO_MATrec$out[[length(MARGO_MATrec$out)]]$dist.n),
                                core = 'MARGO LGM',
                                method = 'MAT')

gplot_sensWA <- rbind.data.frame(sens_MD95_WA, sens_M30_WA, sens_MARGO_WA)
gplot_sensMAT <- rbind.data.frame(sens_MD95_MAT, sens_M30_MAT, sens_MARGO_MAT)

Fig8 <- ggarrange(
  ggplot(gplot_sensMAT, aes(richness, sd)) +
    geom_point(alpha = 0.5, size = 0.5) +
    #geom_smooth(method = 'lm') +
    #scale_color_gradientn(colours = c('navy', 'forestgreen', 'gold', 'firebrick1'), guide = 'none') +
    ggtitle('MAT') +
    ylab('Pruning sensitivity') +
    facet_grid(core~.) +
    theme_bw() +
    theme(text = element_text(size = 9),
          panel.grid = element_line(size = 0.25)),
  ggplot(gplot_sensWA, aes(richness, sd)) +
    geom_point(alpha = 0.5, size = 0.5) + 
    #geom_smooth(method = 'lm') +
    ylim(c(0, 4)) +
    ggtitle('WA') +
    ylab('Pruning sensitivity') +
    facet_grid(core~.) +
    theme_bw() +
    theme(text = element_text(size = 9),
          panel.grid = element_line(size = 0.25)),
  ggplot(gplot_sensMAT, aes(sqcd, sd)) +
    geom_point(alpha = 0.5, size = 0.5) + 
    #geom_smooth(method = 'lm') +
    #ggtitle('MAT') +
    ylab('Pruning sensitivity') +
    xlab('Analogue quality (dissimilarity)') +
    facet_grid(core~.) +
    theme_bw() +
    theme(text = element_text(size = 9),
          panel.grid = element_line(size = 0.25)),
  ggplot(gplot_sensWA, aes(sqcd, sd)) +
    geom_point(alpha = 0.5, size = 0.5) + 
    #geom_smooth(method = 'lm') +
    ylim(c(0, 4)) +
    #ggtitle('WA') +
    ylab('Pruning sensitivity') +
    xlab('Analogue quality (dissimilarity)') +
    facet_grid(core~.) +
    theme_bw() +
    theme(text = element_text(size = 9),
          panel.grid = element_line(size = 0.25)),
  ncol = 2
)

# Figure 8: analogue quality as a function of inferred SST change ####
# indicative of mixing of assemblages across rapid change
# better with 1st axis of PCA?

MD95spp <- sqrt(MD95_WArec$coredat[, -c(1, 2)])
MD95pca <- rda(MD95spp, scale = FALSE)
MD95pcaImp <- summary(MD95pca, scaling = 1)$cont$importance
MD95pcaScore <- MD95pca$CA$u[,1]

M30spp <- sqrt(M30_WArec$coredat[, -c(1, 2)])
M30pca <- rda(M30spp, scale = FALSE)
M30pcaImp <- summary(M30pca, scaling = 1)$cont$importance
M30pcaScore <- M30pca$CA$u[,1]

gplotAnalogueMD95 <- cbind.data.frame(
  sqcd = sens_MD95_WA$sqcd[-length(sens_MD95_WA$sqcd)],
  diffMAT = log10(abs(diff(MD95_MATrec$reconstructed.SST[,33]))/diff(MD95_WArec$coredat$age_ka_BP)),
  diffWA = log10(abs(diff(MD95_WArec$reconstructed.SST[,33]))/diff(MD95_WArec$coredat$age_ka_BP)),
  diffPCA = log(abs(diff(MD95pcaScore))/diff(MD95_WArec$coredat$age_ka_BP))
)
gplotAnalogueMD95$diffMAT[is.infinite(gplotAnalogueMD95$diffMAT)] <- NA

gplotAnalogueM30 <- cbind.data.frame(
  sqcd = sens_M30_WA$sqcd[-length(sens_M30_WA$sqcd)],
  diffMAT = log(abs(diff(M30_MATrec$reconstructed.SST[,33]))/diff(M30_WArec$coredat$age_ka_BP)),
  diffWA = log10(abs(diff(M30_WArec$reconstructed.SST[,33]))/diff(M30_WArec$coredat$age_ka_BP)),
  diffPCA = log10(abs(diff(M30pcaScore))/diff(M30_WArec$coredat$age_ka_BP))
)
gplotAnalogueM30$diffMAT[is.infinite(gplotAnalogueM30$diffMAT)] <- NA

Fig8 <- ggarrange(ggplot(gplotAnalogueMD95, aes(diffWA, sqcd)) +
            geom_smooth(method = 'lm', colour = 'tomato3') + 
            geom_point(alpha = 0.5, size = 0.5) +
            ylab('Dissimilarity') +
            xlab(expression(paste(Delta, 'SST WA [log10('*degree*C*'/kyr)]'))) +
            ggtitle('MD95-2040') +
            ylim(c(0, 0.35)) +
            theme_bw() +
            theme(text = element_text(size = 9),
                    panel.grid = element_line(size = 0.25)),
          ggplot(gplotAnalogueMD95, aes(diffMAT, sqcd)) +
            geom_smooth(method = 'lm',colour = 'tomato3') +   
            geom_point(alpha = 0.5, size = 0.5) +
            ylab('Dissimilarity') +
            xlab(expression(paste(Delta, 'SST MAT [log10('*degree*C*'/kyr)]'))) +
            ylim(c(0, 0.35)) +
            theme_bw() +
            theme(text = element_text(size = 9),
                  panel.grid = element_line(size = 0.25)),
          ggplot(gplotAnalogueMD95, aes(diffPCA, sqcd)) +
            geom_smooth(method = 'lm', colour = 'tomato3') + 
            geom_point(alpha = 0.5, size = 0.5) +
            ylab('Dissimilarity') +
            xlab(expression(paste(Delta, 'PCA1 [log10(1/kyr)]'))) +
            ylim(c(0, 0.35)) +
            theme_bw() +
            theme(text = element_text(size = 9),
                  panel.grid = element_line(size = 0.25)),
          ggplot(gplotAnalogueM30, aes(diffWA, sqcd)) +
            geom_smooth(method = 'lm', colour = 'tomato3') + 
            geom_point(alpha = 0.5, size = 0.5) +
            #ylab('Dissimilarity') +
            xlab(expression(paste(Delta, 'SST WA [log10('*degree*C*'/kyr)]'))) +
            ylim(c(0, 0.35)) +
            ggtitle('M35003-4') +
            theme_bw() +
            theme(text = element_text(size = 9),
                  axis.title.y = element_blank(),
                  panel.grid = element_line(size = 0.25)),
          ggplot(gplotAnalogueM30, aes(diffMAT, sqcd)) +
            geom_smooth(method = 'lm', colour = 'tomato3') + 
            geom_point(alpha = 0.5, size = 0.5) +
            #ylab('Dissimilarity') +
            xlab(expression(paste(Delta, 'SST MAT [log10('*degree*C*'/kyr)]'))) +
            ylim(c(0, 0.35)) +
            theme_bw() +
            theme(text = element_text(size = 9),
                  axis.title.y = element_blank(),
                  panel.grid = element_line(size = 0.25)),
          ggplot(gplotAnalogueM30, aes(diffPCA, sqcd)) +
            geom_smooth(method = 'lm', colour = 'tomato3') + 
            geom_point(alpha = 0.5, size = 0.5) +
            #ylab('Dissimilarity') +
            xlab(expression(paste(Delta, 'PCA1 [log10(1/kyr)]'))) +
            ylim(c(0, 0.35)) +
            theme_bw() +
            theme(text = element_text(size = 9),
                  axis.title.y = element_blank(),
                  panel.grid = element_line(size = 0.25)),
          ncol = 2, byrow = FALSE)
