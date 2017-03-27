#!/usr/bin/env julia
## contains all custom types necessary for eco-evo-env speciation island model

module GeneInds

## Types:

type Gene
    sequence::String # contains gene base code
    codesfor::String # function of the gene
    value::Float # numerical effect of function
end

type Chromosome
    genes::Array{Gene}(1) # 1D array of genes
    origin::Bool # parental origin of chromosome (paternal/maternal)
end

type Individual
    genome::Array{Chromosome} # genome = array of chromosomes
    stage::String # demographic stage of individual
end


## methods:

function mutate(gene::Gene)
    for i in gene
    end
end

function recombinate(chromosome::Chromosome)
end

function 

