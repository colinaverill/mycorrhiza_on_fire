#assign traits NA tree species.
rm(list = ls())
source('paths.r')
library(data.table)
library(willeerd)
library(picante)
source('functions/congeneric.merge.r')

#set output path.----
output.path <- spp.symbiosis_and_traits.path
#load NA tree species list, phylogeny, myco trait databases and Nfix database.----
d <- read.csv(NA_tree_species.path)
myco.spp <- readRDS(myco_traits.path)
myco.gen <- readRDS(myco_genera_clean.path)
myco.fam <- myco.gen[myco.gen$Level == 'Family',]
myco.gen <- myco.gen[myco.gen$Level == 'Genus' ,]
nodDB <- read.csv(nodDB.path)
colnames(myco.spp)[3] <- 'myco_doi'
trait.data <- readRDS(spp_traits.path)
trait.data$tpl.Species <- as.character(trait.data$tpl.Species)

#load phylogeny.
phy <- read.tree(phylogeny.path)
phy$tip.label <- paste0(toupper(substr(phy$tip.label, 1, 1)), substr(phy$tip.label, 2, nchar(phy$tip.label)))
phy$tip.label <- gsub('_',' ',phy$tip.label)
phy$node.label <- NULL


#separate out genus and species.----
d$GenusSpecies <- as.character(d$GenusSpecies)
d$genus   <- stringr::word(d$GenusSpecies, 1)
d$species <- stringr::word(d$GenusSpecies, 2)

#drop things without a species, save to append at end.----
d$species <- ifelse(d$species == 'NA', NA, d$species)
drop <- d[is.na(d$species),]
d <- d[!(d$GenusSpecies %in% drop$GenusSpecies),]
#Assign mycorrhizal association.----
#Family level match.
d <- merge(d, myco.fam, by.x = 'Family', by.y = 'ID', all.x = T)

#Genus level match.
sub <- d[is.na(d$MYCO_ASSO),]
sub$myco_doi <- NULL
sub$MYCO_ASSO <- NULL
sub$Level <- NULL
sub <- merge(sub, myco.gen, by.x = 'genus', by.y = 'ID', all.x = T)
d <- d[!(d$species %in% sub$species),]
d <- rbind(d, sub)
d$Level <- NULL

#Species level match.
sub <- d[is.na(d$MYCO_ASSO),]
sub$myco_doi <- NULL
sub$MYCO_ASSO <- NULL
sub$Level <- NULL
sub <- merge(sub, myco.spp, by.x = 'GenusSpecies', by.y = 'Species', all.x = T)
#for some reason one species is getting multiple assignments here. This fixes that.
sub <- sub[!duplicated(sub$species),]
#merge assigments back in.
d <- d[!(d$species %in% sub$species),]
d <- rbind(d, sub)

#Assign N-fixation status.----
nodDB <- data.table(nodDB)
yes.nfix  <- c('Rhizobia','likely_Rhizobia','Frankia','Nostocaceae','likely_present','Present')
nodDB[Consensus.estimate %in% yes.nfix , nfix  := 1]
nodDB <- nodDB[nfix == 1,.(genus,nfix)]
nodDB$genus <- as.character(nodDB$genus)
d$nfix <- ifelse(d$genus %in% nodDB$genus, 1, 0)

#assign N and P traits based on Averill et al. 2019 PNAS.----
#9 specices aren't in the phylogeny. Lets get them in there via congeneric merge.
phy.missing <- as.character(d[!(d$GenusSpecies %in% phy$tip.label),]$GenusSpecies)
new.phy <- congeneric.merge(phy.missing, phy, split = " ")

#There are species that dont merge because the genus isn't in the phylogeny. these are added to the drop list. Just 1 species.
no.phy_species <- d$GenusSpecies[!(d$GenusSpecies %in% new.phy$tip.label)]
drop <- plyr::rbind.fill(drop, d[d$GenusSpecies %in% no.phy_species,])
d <- d[!(d$GenusSpecies %in% no.phy_species),]

trait.grab <- c('Ngreen','Nsenes','Nroots','Pgreen','Psenes','Proots')
trait.output <- list()

for(i in 1:length(trait.grab)){
  #This is where our loop will start.
  #grab a trait.
  trait.lab <- trait.grab[i]
  dat <- trait.data[,c('tpl.Species',trait.lab)]
  dat <- dat[complete.cases(dat),]
  #get trait values actually observed.
  observed.trait <- dat[dat$tpl.Species %in% d$GenusSpecies,]
  observed.spp <- observed.trait$tpl.Species
  #infer traits for un-observed species.
  spp.check <- unique(c(dat$tpl.Species, d$GenusSpecies))
  #prune the phylogeny to match.
  to.drop <- new.phy$tip.label[!(new.phy$tip.label %in% spp.check)]
  infer.phy <- ape::drop.tip(new.phy,to.drop)
  infer.phy <- multi2di(infer.phy) #deal with polytomies.
  
  #order traits to match
  test.trait <- dat[,trait.lab]
  names(test.trait) <- dat$tpl.Species
  #infer traits for species without direct observations.
  unobserved.trait <- picante::phyEstimate(infer.phy, test.trait)
  #infer traits for species with direct observations.
  infer.phy <- ape::drop.tip(infer.phy,rownames(unobserved.trait))
  unobserved.trait.2 <- list()
  for(k in 1:length(observed.spp)){
    sub.trait <- test.trait[!(names(test.trait) %in% observed.spp[k])]
    unobserved.trait.2[[k]] <- picante::phyEstimate(infer.phy, sub.trait)
  }
  
  #arrange output.
  unobserved.trait.2 <- do.call(rbind, unobserved.trait.2)
  unobserved.trait <- rbind(unobserved.trait, unobserved.trait.2)
  unobserved.trait$species <- rownames(unobserved.trait)
  unobserved.trait$se <- NULL
  output <- merge(observed.trait, unobserved.trait, by.x = 'tpl.Species', by.y = 'species', all = T)
  labs <- c(paste0(trait.lab,'_observed'), paste0(trait.lab,'_inferred'))
  colnames(output) <- c('species',labs)
  trait.output[[i]] <- output
}

#Merge together observed/inferred trait output.----
trait.grab <- trait.output[[1]]
for(i in 2:length(trait.output)){
  trait.grab <- cbind(trait.grab,trait.output[[i]][,2:3])
}
#merge together final output and save.----
output <- merge(d, trait.grab, by.x = 'GenusSpecies', by.y = 'species')
output <- plyr::rbind.fill(output, drop)
write.csv(output, output.path)
