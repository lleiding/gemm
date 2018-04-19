#!/usr/bin/env R

## analyse simulations for testing of co-existing lineages (no/low mutation, high compatibility)

library(stats) # what for?
library(ape)
library(ggtree)
library(dplyr)
library(vegan)

args = commandArgs()
filename = args[length(args)]

elements = unlist(strsplit(filename, split = "_"))

timeidx = grep("000", unlist(strsplit(elements, split = "_")))

basename = paste(elements[-timeidx], collapse = "_")

allfiles = list.files(pattern = basename)

allfiles = grep("000.tsv", allfiles, value = T)

allspecies = c()

for(f in allfiles){
    allworld = read.table(f, header = T)
    elements = unlist(strsplit(f, split = "_"))
    timeidx = grep("000", unlist(strsplit(elements, split = "_")))
    time = elements[timeidx]
    time = sub("t", "", time)
    time = gsub(".tsv", "", time)

    allworld$time = as.numeric(time)
    allworld$location = paste(allworld$xloc, allworld$yloc, sep = ".")

    world = allworld[allworld$size >= allworld$repsize,] # only adults

    lineages = names(table(world$lineage))[table(world$lineage) >= 10] # only lineages with at least ten adults

    world = world[world$lineage %in% lineages,]

    world$population = paste(world$location, world$lineage, sep = ".")

    ## make species table with abundance
    species = world[!duplicated(world$population),]
    species$abundance = as.vector(table(world$population))

    allspecies = rbind(allspecies, species)
}

ggplot(allspecies, aes(time, abundance)) + geom_point(aes(time, abundance, color = lineage, size = temp, alpha = nichea))
ggsave(file=paste0(basename, "_timeabund", ".pdf"), height = 8, width = 12)
