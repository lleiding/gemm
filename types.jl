#!/usr/bin/env julia
## contains all custom types necessary for eco-evo-env speciation island model

## NOTE: all functions: make methods for ind/pop/island/world level!!!!

module GeneInds

using Distributions

## Global(?) variables/constants:
## sdmut: mutation standard deviation


## Types:

mutable struct Trait
    name::String
    value::Float64 # numerical value
    strength::Float64 # mutation strength
    codedby::Array{String,1}
end

mutable struct Gene
    sequence::Array{Char,1} # contains gene base code
    id::String # gene identifier
    codes::Array{Trait,1}
end

struct Chromosome #CAVE! mutable?
    genes::Array{Gene,1} # 1D array of genes
    ##    origin::Bool # parental origin of chromosome (paternal/maternal)
end


mutable struct Individual
    genome::Array{Chromosome,1} # genome = 2D array of chromosomes (>=1 sets)
    traits::Array{Trait,1}
    stage::String # demographic stage of individual
    isnew::Bool # indicator whether individual is new to a patch or has already dispersed etc.
    fitness::Float64 # reproduction etc. scaling factor representing life history
    size::Float64 # body size/mass -> may replace stage.
end

mutable struct Patch
    community::Array{Individual,1} # holds the population (1D) of prob not: present species (2nd D)
    altitude::Float64 # altitude: corresponds to T
    nichea::Float64 # additional niches,
    nicheb::Float64 # e.g. precipitation
end

mutable struct Island # needs new type? -> rather world
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
        if patch.community[idx].size >= 0.5 * patch.community[idx].traits[find(x->x.name=="maxsize",
                                                                               patch.community[idx].traits)][1].value
            noffs = rand(Poisson(patch.community[idx].traits[find(x->x.name=="noffspring",
                                                                  patch.community[idx].traits)][1].value))
            for i in 1:noffs
                ind = deepcopy(patch.community[idx])
                ind.stage = "seed"
                ind.isnew = true
                ind.fitness = 1.0 # TODO: function for this! Consider environment!
                ind.size = ind.traits[find(x->x.name=="seedsize",ind.traits)][1].value
                mutate!(ind, patch.altitude)
                push!(patch.community,ind)
            end
        end
        idx += 1
    end
end

function mutate!(ind::Individual, temp::Float64)
    temp - 293 > e ? (tempdiff = temp - 293) : (tempdiff = e) # difference to standard temperature #CAVE!
    prob = ind.traits[find(x->x.name=="pmut",ind.traits)][1].value
    for chr in ind.genome
        for gene in chr.genes
            for i in eachindex(gene.sequence)
                if rand() < log(tempdiff)*prob
                    newbase = rand(collect("acgt"),1)[1]
                    while newbase == gene.sequence[i]
                        newbase = rand(collect("acgt"),1)[1]
                    end
                    gene.sequence[i] = newbase
                    for trait in gene.codes
                        newvalue = trait.value + rand(Normal(0, trait.strength)) # new value for trait
                        newvalue > 1 && (newvalue=1)
                        newvalue < 0 && (newvalue=0)
                        trait.value = newvalue
                    end
                end
            end
        end
    end
end

# function recombinate!(chromosome::Chromosome)
# end


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
end

function createtraits(traitnames::Array{String,1})
    traits = Trait[]
    for name in traitnames
        push!(traits,Trait(name,rand(),rand(),[]))
    end
    traits
end

function creategenes(ngenes::Int64,traits::Array{Trait,1})
    genes = Gene[]
    for gene in 1:ngenes
        sequence = collect("acgt"^5) # arbitrary start sequence
        id = randstring(8)
        codesfor = Trait[]
        append!(codesfor,rand(traits,rand(Poisson(0.5))))
        for trait in codesfor
            push!(trait.codedby,id) # crossreference!
        end
        push!(genes,Gene(sequence,id,codesfor))
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
                 traitnames::Array{String,1} = ["maxsize","noffspring","pmut","seedsize"]) # minimal required traitnames
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

function checkviability!(patch::Patch)
    idx=1
    while idx <= size(patch.community,1)
        kill = false
        for trait in patch.community[idx].traits
            size(trait.codedby,1) == 0 && (kill = true) # make sure every trait is coded by at least 1 gene
        end
        if kill
            splice!(patch.community,idx) # or else kill it
            idx -= 1
        end
        idx += 1
    end
        
end

    
end

## Test stuff:
##############



function testscenario(timesteps::Int64=1000,npatches::Int64=100)
    world=Patch[]
    for patch = 1:npatches
        push!(world,Patch(genesis(),293,0.5,0.5))
    end
    for i = 1:timesteps
        for patch = 1:npatches
            checkviability!(world[patch])
            #    size(testpatch.community,1)
            compete!(world[patch])
            #    size(testpatch.community,1)
            reproduce!(world[patch]) # TODO: requires certain amount of resource/bodymass dependent on seedsize!
            #    size(testpatch.community,1)
        end
        i == 1 && return world
    end
    world
end

worldstart,worldend = testscnario()
