#!/usr/bin/env R
##
## code for analyzing data produced by the MIbGxMCmod
## Ludwig Leidinger <ludwig.leidinger@uni-wuerzburg.de>
## Dec., 2017

library(stats) # what for?
library(ape)
library(ggplot2)
library(ggtree)
library(dplyr)
library(vegan)
library(viridis)
#library(phyloseq)


args = commandArgs()
basename = args[length(args)]
rlyallsp = c()

## read main data:
allworld = read.table(paste0(basename, ".tsv"), header = T)

nameparts = unlist(strsplit(basename, "_"))
nameparts[length(nameparts) - 1] = "t0"
mlname = paste(nameparts, collapse = "_")
mlworld = read.table(paste0(mlname, ".tsv"), header = T)

## read sequences (fasta format!):
allseqs = read.dna(file=paste0(basename, ".fa"), format="fasta")
headers = names(allseqs[sapply(allseqs, length) == 200])[c(TRUE,FALSE)]

mlseqs = read.dna(file=paste0(mlname, ".fa"), format="fasta")
mlheaders = names(mlseqs[sapply(mlseqs, length) == 200])[c(TRUE,FALSE)]

## make sure seqs and data have same length and elements:
allworld = allworld %>% filter(island == 1)

allworld$tips = headers
cols = ncol(allworld)
allworld = allworld[,c(cols,1:(cols-1))]
names(allworld)[names(allworld) == "id"] = "ID"
allworld$location = paste(allworld$xloc, allworld$yloc, sep = ",")
allworld$temp.C = allworld$temp - 273
allworld$habitat = paste0(allworld$temp.C, "C.", allworld$p, "p")
allworld$linkage = "intermediate"
allworld$linkage[allworld$lnkgunits == 2] = "full"
allworld$linkage[allworld$lnkgunits >= 22] = "none"
nameparts = unlist(strsplit(basename, "_"))
time = nameparts[length(nameparts) - 1]
time = strsplit(time, "t")
time = as.numeric(rev(unlist(time))[1])
allworld$time = time

mlworld$tips = mlheaders
cols = ncol(mlworld)
mlworld = mlworld[,c(cols,1:(cols-1))]
names(mlworld)[names(mlworld) == "id"] = "ID"
mlworld$location = paste(mlworld$xloc, mlworld$yloc, sep = ",")
mlworld$temp.C = mlworld$temp - 273
mlworld$linkage = "intermediate"
mlworld$linkage[mlworld$lnkgunits == 2] = "full"
mlworld$linkage[mlworld$lnkgunits >= 22] = "none"
mlworld$time = time

if(mlworld$temp.C > max(allworld$temp.C)){
    mlworld$temp.C = max(allworld$temp.C)
}
if(mlworld$temp.C < max(allworld$temp.C)){
    mlworld$temp.C = min(allworld$temp.C)
}

mlworld$habitat = paste0(mlworld$temp.C, "C.", mlworld$p, "p")

## get ids of lineages with at least four individuals:
lineages = names(table(allworld$lineage))[table(allworld$lineage) >= 5]

allspecies = c()
for(lineage in lineages){
    ## subset
    world = allworld[allworld$lineage == lineage, ]
    if(sum(world$size >= world$repsize) >= 2){ # build trees only for lineages with at least two adult individuals on island
        seqs = allseqs[grep(lineage, names(allseqs))]

        ## filter sequences according to header:
        #seqs = seqs[grep("compat", names(seqs))]
        seqs = seqs[sapply(seqs, length) == 200]
        mlseqs = mlseqs[sapply(mlseqs, length) == 200]
        mlseqs = mlseqs[c(TRUE,FALSE)]
        ## one chromosome copy:
        seqs = seqs[c(TRUE,FALSE)]

        established = world$size > world$seedsize
        world = world[established,]
        seqs = seqs[established]

        if(nrow(world) > 40000){
           inds = sample(1:nrow(world), 40000)
           world = world[inds,]
           seqs = seqs[inds]
        }
        # also get mainland population infos:
        world = rbind(world, mlworld[mlworld$lineage == lineage, ][1,])
        seqs = c(seqs, mlseqs[grep(lineage, names(mlseqs))][1])
        
        ## compute distances:
        dists = dist.dna(seqs, model = "JC69") # use JukesCantor distances JC69: gives most precise results

        ## calculate the tree:
        tre = hclust(dists, method = "average") # CAVE: "average" = UPGMA. NJ?

        ## cluster tips to create species:
        grps = cutree(tre, h = 0.1) # conservative height of 0.1. similarity 0.8 for high tol, 0.95 for low tol

        world$speciesID = NA
        maxgrps = max(grps)
        for(i in unique(grps[tre$order])){
            world$speciesID[grps == i] = maxgrps
            maxgrps = maxgrps - 1
        }
        world$taxon = paste(world$lineage, world$speciesID, sep = "_")
        world$population = paste(world$taxon, world$location, sep = "_")
        locspecab = table(world$population)

        ## make species table with abundance
        species = world[!duplicated(world$population),]
        species$abundance = as.vector(table(world$population))

        phylo = as.phylo(tre)
        phylo = root(phylo, which(phylo$tip.label == mlworld$tips[1]))
        phylo = drop.tip(phylo, setdiff(world$tips, species$tips))

        cladegroups = list()
        for(i in unique(species$speciesID)){
            cladegroups[[i]] = species$tips[species$speciesID == i]
            ##cladegroups[[paste0("c", as.character(i))]] = species$tips[species$speciesID == i]
        }
        nspecs = length(cladegroups)
        phylo = groupOTU(phylo, cladegroups) #, group_name = "speciesID")
  
        p = ggtree(phylo, aes(color=group), size = 1)
        p = p %<+% species + #aes(color=groups) +
            geom_tippoint(colour="transparent", shape = 21, aes(fill=temp.C, size=abundance)) + #, alpha=prec)) +
            scale_color_manual(values = c("black", viridis_pal()(nspecs)),
                               labels=c("", as.character(1:nspecs)),
                               name = "Species ID") +
            scale_fill_viridis(option = "plasma", name = "Temp./°C", breaks = c(19,21,23,25)) +
            scale_alpha(range = c(0.5, 1.0)) +
            guides(color = guide_legend(order = 1),
                   size = guide_legend(title = "Pop. size", order = 2, override.aes = list(fill = "black"))) +
                   ##alpha = guide_legend(title = "Precipitation", override.aes = list(fill = "black"))) +
            geom_treescale(width=0.1) +
            ## scale_x_continuous(breaks=seq(0.0, max(dists)/3, 0.02)) +
            theme(legend.position="right")
        p
        ## old version (worked):
        ##p = ggtree(drop.tip(as.phylo(tre), setdiff(world$tips, species$tips)))                                                 
        ##p = p %<+% species + scale_color_gradient(low="blue", high="red") +
        ##    geom_tippoint(aes(color=prec, size=abundance, alpha=temp.C)) +
        ##    geom_tiplab(aes(subset=!duplicated(speciesID),label=speciesID), geom='text')                                      
        ##p + theme(legend.position="right")
        ## save phylo plots:
        ggsave(file=paste0(basename, "_", lineage, "_tre.pdf"), height = 8, width = 8 * max(dists)/0.5)

        ## store all species:
        allspecies = rbind(allspecies, species)
    }
}

rlyallsp = rbind(rlyallsp, allspecies)

if(length(allspecies) > 1){ # only continue if there were actually phylogenies made
    allspecies = allspecies %>% filter(island == 1)
    allspecies$speciesID = as.factor(allspecies$speciesID)
    m = ggplot(allworld, aes(xloc, yloc))
    m + geom_tile(aes(fill = temp.C, width = 0.98, height = 0.98)) +
        scale_fill_continuous(low="gray60", high="gray90", name = "Temp./°C", breaks = c(19,21,23,25)) +
        scale_color_viridis(name = "Species ID", discrete = TRUE) + #gradientn(colours = rainbow(5)) +
        geom_jitter(data = allspecies, aes(size = abundance, color = speciesID, shape = lineage))+
        guides(color = guide_legend(order = 1),
               size = guide_legend(order = 2, title = "Pop. size"),
               shape = guide_legend(order = 3)) +
        theme_bw()
    ggsave(file=paste0(basename, "_map", ".pdf"), height = 8, width = 10)
    allspecies$rank[order(allspecies$abundance, decreasing =T)] = 1:nrow(allspecies)
    allspecies$localrank = 0
    for(i in unique(allspecies$location)){
        allspecies[allspecies$location == i,]$localrank[order(allspecies[allspecies$location == i,]$abundance, decreasing =T)] = 1:nrow(allspecies[allspecies$location == i,])
        }
    ggplot(allspecies) + geom_col(aes(rank, abundance))
    ggsave(file=paste0(basename, "_rankab", ".pdf"), height = 10, width = 10)
    ggplot(allspecies) + geom_col(aes(localrank, abundance)) + facet_grid(yloc ~ xloc)
    ggsave(file=paste0(basename, "_rankab_cells", ".pdf"), height = 10, width = 10)
    community = with(allspecies, tapply(abundance, list(location, taxon), mean, default = 0))
    sac <- specaccum(community)
    pdf(file=paste0(basename, "_sac", ".pdf"), width = 8, height = 8)
    plot(sac, ci.type="polygon", ci.col="yellow")
    dev.off()
}


## TRAITS:
ggplot(rlyallsp, aes(x = time, y = dispshape, fill = linkage)) + geom_violin(data = filter(rlyallsp, island == 1)) + geom_hline(data = filter(rlyallsp, island == 0), aes(yintercept = dispshape)) + theme_bw() + scale_fill_viridis_d()

cat(c("lineage", "linkage", "n_species"), sep = "\t")
cat("\n")
for(lineage in unique(allspecies$lineage)){
    nspecs = allspecies$species[allspecies$lineage == lineage] %>% unique %>% length
    linkage = allspecies$species[allspecies$lineage == lineage] %>% unique
    cat(c(lineage, linkage, nspecs), sep = "\t")
    cat("\n")
}
corename = unlist(strsplit(basename, "_"))[-length(unlist(strsplit(basename, "_")))]
corename = c(corename, "t1")
corename = paste(corename, collapse = "_")
#mainland = read.table(paste0(corename, ".tsv"), header = T)

traitnames = c("lnkgunits", "ngenes", "temptol", "seedsize", "tempopt", "repsize", "repradius", "dispmean", "prectol", "precopt", "dispshape")

