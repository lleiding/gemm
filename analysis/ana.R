#!/usr/bin/env R
##
## code for analyzing data produced by the MIbGxMCmod
## Ludwig Leidinger <ludwig.leidinger@uni-wuerzburg.de>
## Dec., 2017

library(stats)
library(ape)
library(ggtree)
library(dplyr)
library(vegan)

## read main data:
args = commandArgs()
basename = args[length(args)]

allworld = read.table(paste0(basename, ".tsv"), header = T)

## read sequences (fasta format!):
allseqs = read.dna(file=paste0(basename, ".fa"), format="fasta")

## get id of most abundant lineage:
lineage = names(which.max(table(world$lineage)))

## subset
world = allworld[allworld$lineage == lineage, ]
seqs = seqs[grep(lineage, dimnames(seqs)[[1]]),]

## filter sequences according to header:
seqs = seqs[grep("compat", dimnames(seqs)[[1]]),]

## one chromosome copy:
seqs = seqs[c(TRUE,FALSE),]

## compute distances:
dists = dist.dna(seqs, model = "F81") # use JukesCantor distances JC69, alternatively Felsenstein "F81"

## calculate the tree:
tre = hclust(dists, method = "ward.D2") # CAVE: which method? ward.D2 gives nicest results

## cluster tips to create species:
grps = cutree(tre, h = 0.1) # conservative height of 0.1. similarity 0.8 for high tol, 0.95 for low tol
world$species = paste(lineage, grps, sep = ".")
world$specloc = paste(world$id, world$species, sep = ".")
locspecab = table(world$specloc)

## make species table with abundance
species = world[!duplicated(world$specloc),]
species$abundance = table(world$specloc)

l=list()
for(i in unique(none$island)){
    j=ifelse(i==0, 1, 2)
    l[[j]] = names(grps)[none$island==i]
}
tree=groupOTU(as.phylo(tre), rev(l))
ggtree(tree, aes(color=group)) + theme(legend.position="left")

## plot tree:
p = ggtree(as.phylo(ftre))
fdd=data.frame(taxa=as.character(1:nrow(full)),Temperature=as.factor(full$temp))
p <- p %<+% fdd + geom_tippoint(aes(color=Temperature))
p + scale_color_manual(values=terrain.colors(12)[3:6])

## get data from fasta headers:
tseqs119k = seqs119k[grep("qRB3", names(seqs119k))]
tseqs119k = tseqs119k[grep("tempopt", names(tseqs119k))]
temps=unlist(strsplit(names(tseqs119k), "tempopt"))[c(FALSE, TRUE)] %>% strsplit(",") %>% unlist %>% as.numeric


## NMDS:
fnmds=metaMDS(full[,c(10,13,16:22,25)], group = full$island, trymax=100)
plot(fnmds)
plot(envfit(nmds, full[,c(10,13,16:22,25)]))

fnmds0=metaMDS(full0[,c(11,14,15:23,26)], trymax=100)

pdf(file = "/home/ludwig/Desktop/full0.pdf", bg = "white", height = 5, width = 5)
plot(fnmds0, display="sites", type="n")
cols=ifelse(full0$island==1, "cyan", "salmon")
points(fnmds0, display="sites", col=cols, pch = 16)
plot(envfit(fnmds0, full0[,c(11,14,15:23,26)]))
dev.off()
