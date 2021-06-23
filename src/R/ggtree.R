library(dplyr)
library(ggtree)
library(here)
library(readr)
library(treeio)

tree <- read.tree(here("data/phylogeny.fr/results/orthologues20210504/Tree-210503-newick.txt"))

nams <- read_tsv(here("data/phylogeny.fr/results/orthologues20210504/Btub-tree-names.txt"),
                 col_types="cc") %>% mutate(NewName=paste(substr(NewName,1,3),substr(NewName,4,9)),
                                            Species=factor(substr(NewName,1,3)))

nams <- nams[match(tree$tip.label,nams$CurrentName),]

tree$tip.label <- nams$NewName

#lst <- sapply(orth$ID,grep,tree$tip.label)
#sel <- ! 1:length(tree$tip.label) %in% sort(unlist(lst))
#orth$treeID <- NA
#orth$treeID[sapply(lst,length)!=0] <- tree$tip.label[unlist(lst)]
#orth$treeID[sapply(lst,length)==0] <- tree$tip.label[sel]

#species <- factor(orth$Species[match(tree$tip.label,orth$treeID)])

cols <- rainbow(n=nlevels(nams$Species))[as.integer(nams$Species)]

d <- data.frame(label=tree$tip.label, species=nams$Species)
tree2 <- full_join(tree, d, by='label')

ggtree(tree2) +
  geom_tiplab() + 
  geom_tippoint(aes(color=species)) +
  geom_treescale(offset=-1)


ggtree(tree2,branch.length='none', layout='daylight',) +
  geom_tiplab(size=2) + 
  geom_tippoint(aes(color=species)) +
 vexpand(0.1,direction=1) + vexpand(0.1,direction=-1) + 
  hexpand(0.1,direction=1) + hexpand(0.1,direction=-1)

ggtree(tree2,branch.length='none', layout='circular',) +
  geom_tiplab(size=2) + 
  geom_tippoint(aes(color=species)) +
  vexpand(0.1,direction=1) + vexpand(0.1,direction=-1) + 
  hexpand(0.1,direction=1) + hexpand(0.1,direction=-1)
