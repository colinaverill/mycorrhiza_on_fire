#Visualizing N-green interaction based on linear model (no phylo control).
#Summary - trees burn more if N green is higher. 
#They have more BA standing postfire if they were EM, independent of N green effect.
#Ngreen effect holds up to phylo-autocorrelation, em doesn't, but effect size is similar and p=0.098. Probably a power problem.
#Get more species in this phylogeny.
rm(list=ls())
source('paths.r')
library(MCMCglmm)
library(caper)
library(phytools)


#Load data.----
d <-readRDS(spp_level_traits.path)

#fit a model.----
dat <- d[,c('mburnsppba','mnoburnsppba','Ngreen')]
dat <- dat[complete.cases(dat),]
m <- lm(sqrt(mburnsppba) ~ sqrt(mnoburnsppba)*Ngreen, data = d)
#plot(sqrt(dat$mburnsppba) ~ fitted(m))

#get repsonse norms at low, medium and high Ngreen values.----
lo.dat <- data.frame(seq(0, 5, by = 0.1))
colnames(lo.dat) <- 'mnoburnsppba'
lo.dat$mnoburnsppba <- lo.dat$mnoburnsppba^2
lo.dat$Ngreen <- 7.5
md.dat <- lo.dat
md.dat$Ngreen <- 19
hi.dat <- lo.dat
hi.dat$Ngreen <- 30
lo <- predict(m, newdata = lo.dat)
md <- predict(m, newdata = md.dat)
hi <- predict(m, newdata = hi.dat)


#plot - all data.----
par(mfrow = c(1,2))
#Get colors for plotting.
rbPal <- colorRampPalette(c('green','blue'))
d$col <- rbPal(10)[as.numeric(cut(d$Ngreen,breaks = 10))]

#all data plot.
plot(sqrt(mburnsppba) ~ sqrt(mnoburnsppba), data = d, pch = 16, col = d$col,
     bty = 'l', ylab = NA, xlab = NA)
mtext('all data',side = 3, adj = 0.05, line= 0)
mtext('post-fire biomass', side = 2, line = 2)
mtext('pre-fire biomass', side = 1,  line = 2)

#categorical responses norms.
lines(smooth.spline(lo ~ lo.dat$mnoburnsppba), col = 'green', lwd = 2.5)
lines(smooth.spline(md ~ md.dat$mnoburnsppba), col = '#00AA55', lwd = 2.5)
lines(smooth.spline(hi ~ hi.dat$mnoburnsppba), col = 'blue', lwd = 2.5)

#Plot some data.
#plot.
plot(sqrt(mburnsppba) ~ sqrt(mnoburnsppba), data = d, pch = 16, col = d$col,
     xlim = c(0,2), ylim = c(0, 2), 
     bty = 'l', ylab = NA, xlab = NA)
mtext('subsetted data',side = 3, adj = 0.05, line= 0)
mtext('post-fire biomass', side = 2, line = 2)
mtext('pre-fire biomass', side = 1,  line = 2)

#categorical responses norms.
lines(smooth.spline(lo ~ lo.dat$mnoburnsppba), col = 'green', lwd = 2.5)
lines(smooth.spline(md ~ md.dat$mnoburnsppba), col = '#00AA55', lwd = 2.5)
lines(smooth.spline(hi ~ hi.dat$mnoburnsppba), col = 'blue', lwd = 2.5)
