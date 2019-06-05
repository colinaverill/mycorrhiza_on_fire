#Get species level traits for fire species, as well as environmental covariates.
rm(list=ls())
source('paths.r')
library(Taxonstand)
library(data.table)

#Set output path.----
output.path <- spp_level_traits.path

#Load species names and fire data.-----
d <- read.csv(fire_spp_list.path)
d$X <- NULL
colnames(d)[3] <- 'Species'
#We need to drop rows with no species assignment.
d <- d[!is.na(d$Species),]
d <- d[!(d$genus == 'Unidentified'),]
d <- d[-(grep(' NA',d$Species)),]
d <- d[!is.na(d$Family),]

#load fire data.
fire <- read.csv(response_by_species.path)
fire$X <- NULL
fire$symb <- NULL
colnames(fire)[3] <- 'Species'
#Kill things without species assignment.
fire <- fire[-grep(' NA',fire$Species),]
fire <- fire[-grep('Unidentified',fire$Species),]

#Load data products.----
traits <- readRDS(intra_traits.path)
  myco.spp <- readRDS(myco_traits.path)
  myco.gen <- readRDS(myco_genera_clean.path)
  myco.fam <- myco.gen[myco.gen$Level == 'Family',]
  myco.gen <- myco.gen[myco.gen$Level == 'Genus' ,]
 nodDB <- read.csv(nodDB.path)
colnames(myco.spp)[3] <- 'myco_doi'


#check and standardized species names to the plant list.----
#Looked through the things TPL changed, and it didn't help. Skip this.
#tpl.check <- TPL(d$Species)
#tpl.check$New.Species <- gsub('na','NA',tpl.check$New.Species)
#tpl.check$New.Species <- ifelse(tpl.check$New.Species == 'na','NA',tpl.check$New.Species)
#tpl.check$tpl.Species <- paste0(tpl.check$New.Genus,' ',tpl.check$New.Species)
#d$tpl.Genus <- tpl.check$New.Genus
#d$tpl.Species <- tpl.check$tpl.Species 


#Assign mycorrhizal association.----
#Family and Genus level assignment.
#we know these families and genera are correct, and we know many of these legacy mycorrhizal databases can have errors.
#Therefore, genus level assignment overrides species match from a reference.
#Family level match.
d <- merge(d, myco.fam, by.x = 'Family', by.y = 'ID', all.x = T)
 
#Genus level match.
sub <- d[is.na(d$MYCO_ASSO),]
sub$myco_doi <- NULL
sub$MYCO_ASSO <- NULL
sub$Level <- NULL
sub <- merge(sub, myco.gen, by.x = 'genus', by.y = 'ID', all.x = T)
d <- d[!(d$Species %in% sub$Species),]
d <- rbind(d, sub)
d$Level <- NULL

#Species level match.
sub <- d[is.na(d$MYCO_ASSO),]
sub$myco_doi <- NULL
sub$MYCO_ASSO <- NULL
sub$Level <- NULL
sub <- merge(sub, myco.spp, all.x = T)
#for some reason one species is getting multiple assignments here. This fixes that.
sub <- sub[!duplicated(sub$Species),]
#merge assigments back in.
d <- d[!(d$Species %in% sub$Species),]
d <- rbind(d, sub)

#Assign N-fixation status.----
nodDB <- data.table(nodDB)
yes.nfix  <- c('Rhizobia','likely_Rhizobia','Frankia','Nostocaceae','likely_present','Present')
#yes.nfix2 <- c('Rhizobia','Frankia','Nostocaceae','Present')
nodDB[Consensus.estimate %in% yes.nfix , nfix  := 1]
nodDB <- nodDB[nfix == 1,.(genus,nfix)]
d$nfix <- ifelse(d$genus %in% nodDB$genus, 1, 0)

#Get species level trait means.----
traits.spp <- data.table(traits)
traits.spp <- traits.spp[,.(Species, Ngreen,Pgreen,Nsenes,Psenes,Nroots,Proots,log.LL,root_lifespan)]
traits.spp <- traits.spp[,lapply(.SD,mean,na.rm=TRUE), by=Species]
#convert NaN values to NA.
traits.spp <- traits.spp[,lapply(.SD,function(x){ifelse(is.nan(x),NA,x)})]
#merge into dataframe.
d <- merge(d, traits.spp, all.x = T)

#Merge in fire responses.----
fire$Family <- NULL
fire$genus <- NULL
fire <- fire[!(duplicated(fire$Species)),]
d <- merge(d, fire, all.x = T)
d <- d[d$MYCO_ASSO %in% c('AM','ECM'),]
d$em <- ifelse(d$MYCO_ASSO == 'AM',0,1)
#mod <- lm(sqrt(mburnsppba) ~ sqrt(mnoburnsppba) , data = d)
#mod <- lm(sqrt(mburnsppba) ~ sqrt(mnoburnsppba)*Ngreen + em + nfix, data = d)
#summary(mod)

#Save output.----
saveRDS(d, output.path)