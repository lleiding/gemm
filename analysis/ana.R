#!/usr/bin/env R
##
## code for analyzing data produced by the MIbGxMCmod
## Ludwig Leidinger <ludwig.leidinger@uni-wuerzburg.de>
## Dec., 2017

library(stats)
library(ape)
library(ggtree)
## library(dplyr)
library(vegan)

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

## plot tree:
p = ggtree(as.phylo(baum))
dd=data.frame(taxa=as.character(1:nrow(world)),Temperature=as.factor(world$temp))
p <- p %<+% dd + geom_tippoint(aes(color=Temperature))
p +scale_color_manual(values=terrain.colors(12)[3:6])

## simulation arena maps - temperature:
row1 = c(288, 293, 293, 298, 303, NA, NA, NA, NA, NA, NA, NA)
row2 = c(288, 293, 293, 298, 303, NA, NA, NA, 303, 298, 303, NA)
row3 = c(288, 293, 293, 298, 303, NA, NA, NA, 298, 293, 298, NA)
row4 = c(288, 293, 293, 298, 303, NA, NA, NA, 303, 298, 303, NA)
row5 = c(288, 293, 293, 298, 303, NA, NA, NA, NA, NA, NA, NA)

maptemp = rbind(row1, row2, row3, row4, row5)
rownames(maptemp) = 0:4
colnames(maptemp) = 0:11

## IDs:
rowa = c(1, 6, 11, 16, 21, NA, NA, NA, NA, NA, NA, NA)
rowb = c(2, 7, 12, 17, 22, NA, NA, NA, 100, 103, 106, NA)
rowc = c(3, 8, 13, 18, 23, NA, NA, NA, 101, 104, 107, NA)
rowd = c(4, 9, 14, 19, 24, NA, NA, NA, 102, 105, 108, NA)
rowe = c(5, 10, 15, 20, 25, NA, NA, NA, NA, NA, NA, NA)

mapids = rbind(rowa, rowb, rowc, rowd, rowe)
rownames(mapids) = 0:4
colnames(mapids) = 0:11

## plot maps:
pdf(file = "mapfile_image_temp.pdf", bg = "white", height = 5, width = 12)
par(mar = c(0, 0, 0, 0))
image(t(maptemp), col = terrain.colors(12)[3:6], xaxt = "n", yaxt = "n")
legend(0.4, 0.55, legend = paste(unique(as.vector(maptemp))[1:4], " K"), col = terrain.colors(12)[3:6], pch = 15, bty ="n", cex=3)
dev.off()

## NMDS:
nmds=metaMDS(world[,c(10,13,16:22,25)], trymax=100)
plot(nmds)
plot(envfit(nmds, world[,c(10,13,16:22,25)]))
