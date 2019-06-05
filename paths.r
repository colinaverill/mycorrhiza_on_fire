#paths.r for myco_fire_traits
data.dir <- '/fs/data3/caverill/myco_fire_trait_data/'

#Raw data products received from Adam.----
fire_spp_list.path <- paste0(data.dir,'treespplist.csv')
response_by_genus.path <- paste0(data.dir,'responsebygenusforcolin.csv')
response_by_species.path <- paste0(data.dir,'responsebysppforcolin.csv')

#Analysis products.----
spp_level_traits.path <- paste0(data.dir,'spp_level_traits.rds')

#Trait databases.----
dir <- '/fs/data3/caverill/myc_traits/'
intra_traits.path <- paste0(dir,'merged_intra_traits_names_hand_checked.rds')
 myco_traits.path <- paste0(dir,'merged_myco_traits.rds')
 myco_genera_clean.path <- paste0(dir,'myco_genera_clean.rds')
       nodDB.path <- paste0(dir,'nodDB_v1.csv') #N fixing traits.

#Phylogeny.----
phylogeny.path <- paste0(dir,'colin_2018-12--2.tre')
       