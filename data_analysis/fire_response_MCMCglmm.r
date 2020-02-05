#pgls analysis
#Summary - trees burn more if N green is higher. 
#They have more BA standing postfire if they were EM, independent of N green effect.
#Ngreen effect holds up to phylo-autocorrelation, em doesn't, but effect size is similar and p=0.098. Probably a power problem.
#Get more species in this phylogeny.
rm(list=ls())
source('paths.r')
library(MCMCglmm)
library(caper)
library(phytools)


#Load data and phylogeny.----
d <-readRDS(spp_level_traits.path)
phy <- read.tree(phylogeny.path) #'colin_2018-12--2.tre'

#Some data prep.----
phy$tip.label <- paste0(toupper(substr(phy$tip.label, 1, 1)), substr(phy$tip.label, 2, nchar(phy$tip.label)))
phy$tip.label <- gsub('_',' ',phy$tip.label)
phy$node.label <- NULL
d <- d[d$Species %in% phy$tip.label,]
#Setup phylo random effect.
phy <- force.ultrametric(phy)
rnd <- inverseA(phy)$Ainv


#Run MCMCglmm.----
#Get compelte cases for your analysis.
mcmc.form <- as.formula('sqrt(mburnsppba) ~ sqrt(mnoburnsppba)*Ngreen + em + nfix')
dat <- d[,c('Species','mburnsppba','mnoburnsppba','Ngreen','em','nfix')]
dat <- dat[complete.cases(dat),]
#priors.
priors <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
#fit model.
mod <- MCMCglmm(mcmc.form, random=~Species, ginverse=list(Species=rnd), data=dat, prior = priors)
summary(mod)
lm.mod <- lm(mcmc.form, data = dat)
summary(lm.mod)
