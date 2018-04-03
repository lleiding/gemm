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
headers = grep("compat", dimnames(allseqs)[[1]], value = T)[c(TRUE,FALSE)]

rownames(allworld) = headers
allworld$tips = headers
cols = ncol(allworld)
allworld = allworld[,c(cols,1:cols-1)]
names(allworld)[names(allworld) == "id"] = "ID"
allworld$location = paste(allworld$xloc, allworld$yloc, sep = ".")

## get ids of most abundant lineages:
lineages = names(table(allworld$lineage))[table(allworld$lineage) > 10]

for(lineage in lineages){
    ## subset
    world = allworld[allworld$lineage == lineage, ]
    seqs = allseqs[grep(lineage, dimnames(allseqs)[[1]]),]

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
    world$population = paste(world$id, world$species, sep = ".")
    locspecab = table(world$population)

    ## make species table with abundance
    species = world[!duplicated(world$population),]
    species$abundance = as.vector(table(world$population))

    p = ggtree(drop.tip(as.phylo(tre), setdiff(world$tips, species$tips)))
    p = p %<+% species + geom_tippoint(aes(color=location, size=abundance))
    p + theme(legend.position="right")

    ## save phylo plots:
    ggsave(file=paste(basename, lineage, "pdf", sep= "."), height = 10, width = 10)
}

