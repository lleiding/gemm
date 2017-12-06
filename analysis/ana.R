#!/usr/bin/env R
##
## code for analyzing data produced by the MIbGxMCmod
## Ludwig Leidinger <ludwig.leidinger@uni-wuerzburg.de>
## Dec., 2017

library(stats)
library(ape)
library(ggtree)
## library(dplyr)
## library(vegan) # necessary?

## read main data:
world = read.table("*.tsv", header = T)

## read sequences (fasta format!):
seqs = read.dna(file="*.fa", format="fasta")

## compute distances:
dists = dist.dna(seqs, model = "JC69") # use JukesCantor distances

## calculate the tree:
tre = hclust(dists) # CAVE: which method?

## cluster tips:
grps = cutree(tre, h = 0.1) # conservative height of 0.1
world$species = as.vector(grps)

