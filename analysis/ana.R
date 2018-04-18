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
#library(phyloseq)

## read main data:
args = commandArgs()
basename = args[length(args)]

allworld = read.table(paste0(basename, ".tsv"), header = T)

## read sequences (fasta format!):
allseqs = read.dna(file=paste0(basename, ".fa"), format="fasta")
headers = grep("compat", dimnames(allseqs)[[1]], value = T)[c(TRUE,FALSE)]

#rownames(allworld) = headers
allworld$tips = headers
cols = ncol(allworld)
allworld = allworld[,c(cols,1:(cols-1))]
names(allworld)[names(allworld) == "id"] = "ID"
allworld$location = paste(allworld$xloc, allworld$yloc, sep = ".")
allworld$temp.C = allworld$temp - 273
maxtemp = max(allworld$temp.C)
allworld$habitat = paste0(allworld$temp.C, "C.", allworld$nichea, "p")

## get ids of lineages with at least four individuals:
lineages = names(table(allworld$lineage))[table(allworld$lineage) >= 5]

allspecies = c()

for(lineage in lineages){
    ## subset
    world = allworld[allworld$lineage == lineage, ]
    if(sum(world$size >= world$repsize) >= 2){ # build trees only for lineages with at least two adult individuals on island
        seqs = allseqs[grep(lineage, dimnames(allseqs)[[1]]),]

        ## filter sequences according to header:
        seqs = seqs[grep("compat", dimnames(seqs)[[1]]),]

        ## one chromosome copy:
        seqs = seqs[c(TRUE,FALSE),]

        adults = world$size >= world$repsize
        world = world[adults,]
        seqs = seqs[adults,]

        if(nrow(world) > 40000){
           inds = sample(1:nrow(world), 40000)
           world = world[inds,]
           seqs = seqs[inds,]
        }
        
        ## compute distances:
        dists = dist.dna(seqs, model = "F81") # use JukesCantor distances JC69, alternatively Felsenstein "F81"

        ## calculate the tree:
        tre = hclust(dists, method = "ward.D2") # CAVE: which method? ward.D2 gives nicest results

        ## cluster tips to create species:
        grps = cutree(tre, h = 0.2) # conservative height of 0.1. similarity 0.8 for high tol, 0.95 for low tol
        world$grps = grps
        world$speciesID = 0
        maxgrps = max(grps)
        for(i in unique(grps[tre$order])){
            world$speciesID[world$grps == i] = maxgrps
            maxgrps = maxgrps - 1
        }
        #world$speciesID = paste0(grps)
        world$population = paste(world$location, world$speciesID, sep = ".")
        locspecab = table(world$population)

        ## make species table with abundance
        species = world[!duplicated(world$population),]
        species$abundance = as.vector(table(world$population))

        p = ggtree(drop.tip(as.phylo(tre), setdiff(world$tips, species$tips)))
        p = p %<+% species + scale_color_gradient(low="blue", high="red") + geom_tippoint(aes(color=nichea, size=abundance, alpha=temp.C)) + geom_tiplab(aes(subset=!duplicated(speciesID),label=speciesID), geom='text')
        p + theme(legend.position="right")

        ## save phylo plots:
        ggsave(file=paste0(basename, "_", lineage, ".pdf"), height = 10, width = 10)

        ## store all species:
        allspecies = rbind(allspecies, species)
    }
}

if(length(allspecies) > 1){ # only continue if there were actually phylogenies made
    m = ggplot(allworld, aes(xloc, yloc))
    m + geom_tile(aes(fill = temp.C, width = 0.95, height = 0.95)) +
        scale_fill_continuous(low="white", high="black") +
        scale_color_gradientn(colours = rainbow(5)) +
        geom_jitter(data = allspecies, aes(size = abundance, color = speciesID, shape = lineage))
    ggsave(file=paste0(basename, "_map", ".pdf"), height = 8, width = 10)
    allspecies$rank[order(allspecies$abundance, decreasing =T)] = 1:nrow(allspecies)
    ggplot(allspecies) + geom_bar(aes(abundance))
    ggsave(file=paste0(basename, "_rankab", ".pdf"), height = 10, width = 10)
    ggplot(allspecies) + geom_bar(aes(abundance), width=10) + facet_grid(yloc ~ xloc)
    ggsave(file=paste0(basename, "_rankab_cells", ".pdf"), height = 10, width = 10)
    community = with(allspecies, tapply(abundance, list(location, speciesID), mean, default = 0))
    sac <- specaccum(community)
    pdf(file=paste0(basename, "_sac", ".pdf"), width = 8, height = 8)
    plot(sac, ci.type="polygon", ci.col="yellow")
    dev.off()
}

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

save.image(paste0(basename, ".R"))

