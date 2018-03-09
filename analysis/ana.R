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
world = read.table("*.tsv", header = T)

## read sequences (fasta format!):
seqs = read.dna(file="*.fa", format="fasta")

## filter sequences according to header:
nseqs = seqs[grep("neutral", dimnames(seqs)[[1]]),]
nseqs = nseqs[grep("sO1E", dimnames(nseqs)[[1]]),]
rawseqs = as.vector(sapply(nseqs, function(x) paste0(as.character(x), collapse=""))) # get a vector of the sequences

## compute distances:
dists = dist.dna(nseqs, model = "JC69") # use JukesCantor distances, alternatively Felsenstein "F81"

## calculate the tree:
tre = hclust(dists, method = "ward.D2") # CAVE: which method? ward.D2 gives nicest results

## cluster tips:
grps = cutree(tre, h = 0.1) # conservative height of 0.1
world$species = as.vector(grps)

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
