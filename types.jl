#!/usr/bin/env julia
## contains all custom types necessary for eco-evo-env speciation island model

module GeneInds

using Distributions

## Types:

type Gene
    sequence::String # contains gene base code
    codes::String # function of the gene
    value::Float64 # numerical effect of function
end

type Chromosome
    genes::Array{Gene,1} # 1D array of genes
    origin::Bool # parental origin of chromosome (paternal/maternal)
end

type Individual
    genome::Array{Chromosome,2} # genome = 2D array of chromosomes (>=1 sets)
    gamete::Array{Chromosome,2} # gamete = 2D array of chromosomes (>=0 sets)
    fitness::Float64 # reproduction scaling factor representing life history
    stage::String # demographic stage of individual
    isnew::Bool # indicator whether individual is new to a patch
end

type Patch
    community::Array{Individual,1} # holds the population (1D) of present species (2nd D)
    altitude::Float64 # altitude: corresponds to T
    nichea::Float64 # additional niches,
    nicheb::Float64 # e.g. precipitation
end

type Island
    patches::Array{Patch,1}
end


## methods:

function mutate(gene::Gene, temp::Float64, p::Float64)
    for i in eachindex(gene.sequence)
        (rand < temp*p) && (gene.sequence[i]=rand(collect("acgt"),1)[1]) # for now, but consider indels!
    end
end

function recombinate(chromosome::Chromosome)
end

function germinate(ind::Individual)
    if ind.stage == "seed"
        ind.stage = "juvenile"
    end
end

function mature(ind::Individual)
    if ind.stage == "juvenile"
        ind.stage = "adult"
    end
end

function reproduce(ind::Individual)
    if ind.stage == "adult"
    end
end

function disperse(ind::Individual)
end

