#!/usr/bin/env julia
## contains all custom types necessary for eco-evo-env speciation island model

## NOTE: all functions: make methods for ind/pop/island/world level!!!!

module GeneInds

using Distributions

## Global(?) variables/constants:
## sdmut: mutation standard deviation


## Types:

type Trait
    name::String
    value::Float64 # numerical value
##    strength::Float64 # mutation strength
end

type Gene # really need new type?
    sequence::String # contains gene base code
    id::String # gene identifier
    codes::Array{Trait,1}
end

type Chromosome # placeholder, not used for now # maybe implement as just array
    genes::Array{Gene,1} # 1D array of genes
##    origin::Bool # parental origin of chromosome (paternal/maternal)
end
    

type Individual
    genome::Array{Chromosome,1} # genome = 2D array of chromosomes (>=1 sets)
    traits::Array{Trait,1}

    stage::String # demographic stage of individual
    isnew::Bool # indicator whether individual is new to a patch or has already dispersed etc.
    fitness::Float64 # reproduction etc. scaling factor representing life history
    size::Float64 # body size/mass -> may replace stage.
    ## Traits:
    # noff::Gene # mean number of offspring GENE
    # pgerm::Gene # probability of germination GENE
    # pmat::Gene # probability of maturation GENE
    # pdie::Gene # probability of dying GENE

    # pdisp::Gene # probability of dispersal
    # ddisp::Gene # Dispersal distance: radius around origin patch, destination random

    # pmut::Gene # mutation probability
    ## metabolic properties:
    ## bodysizes at different life stages
    ## normalisation coefficients for biological rates
    ## (reproduction, growth,...) <-see above!

    ## maybe necessary stuff/expansions:
    # gamete::Array{Chromosome,1} # gamete = 2D array of chromosomes (>=0 sets)    

    # dead::Bool # is individual dead? if needed...
end

type Patch
    community::Array{Individual,1} # holds the population (1D) of prob not: present species (2nd D)
    altitude::Float64 # altitude: corresponds to T
    nichea::Float64 # additional niches,
    nicheb::Float64 # e.g. precipitation
end

type Island # needs new type? -> rather world
    patches::Array{Patch,2} # 2D grid
end


## methods:
function compete!(patch::Patch)
    while sum(map(x->x.size,patch.community)) > 1 # occupied area larger than available
        victim = rand(1:size(patch.community,1))
        splice!(patch.community, victim)
    end
end

function reproduce!(patch::Patch)
    idx = 1
    while idx <= size(patch.community,1)
        if patch.community[idx].size >= 0.5 * traits[find(x->x.name=="maxsize",traits)][1].value)
            noffs = rand(Poisson(patch.community[idx].traits[find(x->x.name=="noffspring",patch.community[idx].traits)][1].value))
            for i in 1:noffs
                ind = Individual(patch.community[idx].genome,patch.community[idx].traits,"seed",true,
                1.0,patch.community[idx].traits[find(x->x.name=="seedsize",
                                                     patch.community[idx].traits)][1].value)
                mutate!(ind, patch.temperature)
                push!(patch.community,ind)
            end
        end
        idx += 1
    end
end

function mutate!(gene::Gene, temp::Float64, prob::Float64) # or maybe just rare mutation events, where random bp mutates?
    for i in eachindex(gene.sequence)
        if (rand() < temp*prob)
            gene.sequence[i] = rand(collect("acgt"),1)[1] # for now, but consider indels!
            for trait in gene.codes
                trait.value += rand(Normal(0, trait.strength)) # new value for trait
            end
        end
    end
end

# function recombinate!(chromosome::Chromosome)
# end

## following: act @ population/community level?

function evalenv!(patch::Patch,ind::Individual) # necessary or happening all the time? what about changes?
end

function germinate!(ind::Individual)
    ind.dead && return
    ind.isnew = false
    if (ind.stage == "seed") && (rand() <= ind.pgerm) 
        ind.stage = "juvenile"
    end
    ## decide activation/inactivation of alleles here!
    (rand() <= ind.pdie) && (ind.dead = true)
end

function mature!(ind::Individual)
    ## consider alleles for survival etc.!
    ind.dead && return
    if (ind.stage == "juvenile") && (rand() <= ind.pmat) 
        ind.stage = "adult"
    end
end

function reproduce(ind::Individual) # make two functions: one for number of offspring, second to create offspring?!
    ## genetic "fitness"!
    ind.dead && return []
    offspring = []
    if ind.stage == "adult"
        noff = rand(Poisson(ind.noff))
        for i in 1:noff
            child = deepcopy(ind)
            child.isnew = true
            child.stage = "seed"
            push!(offspring, child)
        end
    end
    offspring
end

function disperse!(ind::Individual)
    ## consider genetics + life history! PLUS: investigate effects - which matters more?
    ## direction random, because no information
    ## maybe stop at patch better than origin?
    ## genes for prob. and dist.
    ##
    ## Dispersal modes:
    ## Dispersal kernel for plants (wind dispersed), maybe only in 1 direction?
    ## For insects: ?? behaviour-explicit?
    
    ind.dead && return
    if ind.isnew
        ind.isnew = false
        return
    end
    (rand() > ind.pdisp) && return # abort if no dispersal rolled
    ind.isnew = true
    ydir,xdir = 0,0
    pydist=sqrt(xdir^2+ydir^2)
    while (pydist > ind.ddisp) || (pydist == 0)
        ydir,xdir = rand(-ind.ddisp:ind.ddisp),rand(-ind.ddisp:ind.ddisp)
        pydist=sqrt(xdir^2+ydir^2)
    end
    ydir,xdir # return direction vector: y,x
end

function getDirection(distances::Distribution) ## check what's right!
    ydir,xdir = 0,0
    distance = rand(distances)
    pydist=sqrt(xdir^2+ydir^2)
    while (pydist > distance) || (pydist == 0)
        ydir,xdir = rand(-distance:distance),rand(-distance:distance)
        pydist=sqrt(xdir^2+ydir^2)
    end
    ydir,xdir # return direction vector: y,x
end

## framework for functions:
##(world:World)
function disperse!(world::World)
for p in world.patches # given world contains patches
    i = 0
    ##(p:Patch)
    while i < size(p.community,1) # alternatively length
        i += 1
        ## decide on dispersal:
        p.community[i].isnew && continue
        p.community[i].isnew = true
        (rand > p.community[i].pdisp) && continue
        xdir,ydir = getDirection(p.community[i].ddisp)
        #push ind to new patch
        splice!(p.community,i) # and delete it in current community
        i -= 1 # reset counter
    end
end

function createtraits(traitnames::Array{String,1})
    traits = Trait[]
    for name in traitnames
        push!(traits,Trait(name,rand()))
    end
    traits
end

function creategenes(ngenes::Int64,traits::Array{Trait,1})
    genes = Gene[]
    for gene in 1:ngenes
        sequence = "acgt"^5 # arbitrary start sequence
        id = randstring(8)
        codes = Trait[]
        append!(codes,rand(traits,rand(Poisson(0.5))))
        push!(genes,Gene(sequence,id,codes))
    end
    genes
end

function createchrs(nchrs::Int64,genes::Array{Gene,1})
    ngenes=size(genes,1)
    if nchrs>1
        chrsplits = sort(rand(1:ngenes,nchrs-1))
        chromosomes = Chromosome[]
        for chr in 1:nchrs
            if chr==1 # first chromosome
                push!(chromosomes, Chromosome(genes[1:chrsplits[chr]]))
            elseif chr==nchrs # last chromosome
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):end]))
            else
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):chrsplits[chr]]))
            end
        end
    else # only one chromosome
        chromosomes = [Chromosome(genes)]
    end
    chromosomes
end

function genesis(ninds::Int64=100, maxgenes::Int64=20, maxchrs::Int64=5,
                 traitnames::Array{String,1} = ["pdisp","ddisp","maxsize"]) # arbitrary traitnames for testing
    community = Individual[]
    for ind in 1:ninds
        ngenes = rand(1:maxgenes)
        nchrs = rand(1:maxchrs)
        traits = createtraits(traitnames)
        genes = creategenes(ngenes,traits)
        chromosomes = createchrs(nchrs,genes)
        push!(community, Individual(chromosomes,traits,"adult",false,1.0,
                                    traits[find(x->x.name=="maxsize",traits)][1].value))
    end
    community
end


    
end
