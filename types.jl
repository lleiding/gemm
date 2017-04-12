#!/usr/bin/env julia
## contains all custom types necessary for eco-evo-env speciation island model

module GeneInds

## Types:

type Gene
    sequence::String # contains gene base code
    codes::String # function of the gene
    value::Float # numerical effect of function
end

type Chromosome
    genes::Array{Gene,1} # 1D array of genes
    origin::Bool # parental origin of chromosome (paternal/maternal)
end

type Individual
    genome::Array{Chromosome,2} # genome = 2D array of chromosomes (>=1 sets)
    gamete::Array{Chromosome,2} # gamete = 2D array of chromosomes (>=0 sets)
    fitness::Float # reproduction scaling factor representing life history
    stage::String # demographic stage of individual
    isnew::Bool # indicator whether individual is new to a patch
end

type Patch
    community::Array{Individual,2} # holds the population (1D) of present species (2nd D)
    altitude::Float # altitude: corresponds to T
    nichea::Float # additional niches,
    nicheb::Float # e.g. precipitation
end


## methods:

function mutate(gene::Gene, temp::Float, p::Float)
    for i in eachindex(gene.sequence)
        rand < temp*p && gene.sequence[i]=rand(collect("acgt"),1)[1] # for now, but consider indels!
    end
end

function recombinate(chromosome::Chromosome)
end

function 

