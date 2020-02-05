congeneric.merge <- function(lookup, tree, split="_"){
  if(!is.data.frame(lookup))
    lookup <- data.frame(clean=lookup, stringsAsFactors=FALSE)
  before <- sum(lookup$clean %in% tree$tip.label)
  for(i in seq(nrow(lookup))){
    #prog.bar(i, nrow(lookup))
    if(!is.null(lookup$clean[i]) & !lookup$clean[i] %in% tree$tip.label){
      genus <- strsplit(lookup$clean[i], split, fixed=TRUE)[[1]][1]
      matches <- unique(grep(genus, tree$tip.label, value=TRUE))
      if(length(matches) > 0){
        tree <- drop.tip(tree, matches[-1])
        tip.length <- find.unique.branch.length(tree, matches[1])
        polytomy <- make.polytomy(unique(c(matches, lookup$clean[i])), (tip.length/2))
        tree <- bind.ultrametric.by.replacement(tree, polytomy, matches[1], tip.length)
      }
    }
  }
  return(tree)
}

#supporting functions.
find.unique.branch.length <- function(tree, tip){	
  which.tip <- which(tree$tip.label == tip)
  which.edge <- which(tree$edge[,2] == which.tip)
  tip.length <- tree$edge.length[which.edge]
  return(tip.length)	
}
make.polytomy <- function(species, tip.length=NA){	
  d.f <- data.frame(spp=factor(species))	
  polytomy <- as.phylo.formula(~spp, data=d.f)	
  if(!is.na(tip.length)) polytomy$edge.length <- rep(tip.length, length(species))	
  return(polytomy)	
}	
bind.ultrametric.by.replacement <- function(backbone, donor, replacing.tip.label, donor.length=NA){	
  bind.point <- which(backbone$tip.label == replacing.tip.label)
  backbone <- bind.tree(backbone, donor, where=bind.point)
  which.tip <- which(backbone$tip.label == donor$tip.label[1])
  which.node <- backbone$edge[which(backbone$edge[,2] == which.tip),1]
  which.edge <- which(backbone$edge[,2] == which.node)
  tip.length <- backbone$edge.length[which.edge]
  if(is.na(donor.length)){
    backbone$edge.length[which.edge] <- tip.length/2
  } else {
    backbone$edge.length[which.edge] <- tip.length - donor.length/2
  }
  return(backbone)	
}