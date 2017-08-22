#!/usr/bin/env julia
## contains all custom types necessary for eco-evo-env speciation island model

## NOTE: all functions: make methods for ind/pop/island/world level!!!!

using Distributions

const boltz = 1.38064852e-23 #  J/K = m2⋅kg/(s2⋅K)
const act = 1e-19 # activation energy /J, ca. 0.63eV - Brown et al. 2004
const normconst = 1e10 # normalization constant to get biologically realistic orders of magnitude

## Types:

mutable struct Trait
    name::String
    value::Float64 # numerical value
    ##    strength::Float64 # mutation strength
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
    area::Float64
end


## methods:
function compete!(patch::Patch)
    while sum(map(x->x.size,patch.community)) > patch.area # occupied area larger than available
        victim = rand(1:size(patch.community,1))
        splice!(patch.community, victim)
    end
end

function reproduce!(patch::Patch)
    idx = 1
    temp = patch.altitude
    while idx <= size(patch.community,1)
        currentmass = patch.community[idx].size
        seedsize = patch.community[idx].traits[find(x->x.name=="seedsize",patch.community[idx].traits)][1].value
        if currentmass >= 2 * seedsize > 0
            meanoffs = patch.community[idx].traits[find(x->x.name=="reprate",
                                                        patch.community[idx].traits)][1].value
            mass = patch.community[idx].size
            metaboffs =  meanoffs * currentmass^(-1/4) * exp(-act/(boltz*temp)) * normconst
            noffs = rand(Poisson(metaboffs))
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
    prob = ind.traits[find(x->x.name=="mutprob",ind.traits)][1].value
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
                        newvalue = trait.value + rand(Normal(0, mutsd)) # new value for trait
                        (newvalue > 1 && contains(trait.name,"prob")) && (newvalue=1)
                        newvalue < 0 && (newvalue=0)
                        trait.value = newvalue
                    end
                end
            end
        end
    end
end

function grow!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        growthrate = patch.community[idx].traits[find(x->x.name=="growthrate",patch.community[idx].traits)][1].value
        mass = patch.community[idx].size
        newmass = growthrate * mass^(-1/4) * exp(-act/(boltz*temp)) * normconst #CAVE: what to do when negative growth? -> emergent maximum body size!
        if newmass > 0 && mass > 0
            patch.community[idx].size = newmass
        else
            splice!(patch.community, idx)
            idx -= 1
        end
        idx += 1
    end
end

function age!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        ageprob = patch.community[idx].traits[find(x->x.name=="ageprob",patch.community[idx].traits)][1].value
        mass = patch.community[idx].size
        dieprob = ageprob * mass^(-1/4) * exp(-act/(boltz*temp)) * normconst
        if rand() <= dieprob
            splice!(patch.community, idx)
            idx -= 1
        end
        idx += 1
    end
end
  
    
# function establish!(patch::Patch) #TODO!
#     idx = 1
#     while idx <= size(patch.community,1)
#         !patch.community[idx].isnew && continue
#         if patch
#         splice!(patch.community, idx)
#         idx -= 1
#         idx += 1
#     end
# end

#TODO: Dispersal

#TODO: age mortality!



function createtraits(traitnames::Array{String,1})
    traits = Trait[]
    for name in traitnames
        if contains(name,"rate")
            push!(traits,Trait(name,rand()*10,[]))
        else
            push!(traits,Trait(name,rand(),[]))
        end
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
                 traitnames::Array{String,1} = ["ageprob","growthrate","reprate","mutprob","seedsize"]) # minimal required traitnames
    community = Individual[]
    for ind in 1:ninds
        ngenes = rand(1:maxgenes)
        nchrs = rand(1:maxchrs)
        traits = createtraits(traitnames)
        genes = creategenes(ngenes,traits)
        chromosomes = createchrs(nchrs,genes)
        push!(community, Individual(chromosomes,traits,"adult",true,1.0,rand()))#
#                                    traits[find(x->x.name=="maxsize",traits)][1].value))
    end
    community
end

function checkviability!(patch::Patch) # may consider additional rules...
    idx=1
    while idx <= size(patch.community,1)
        kill = false
        for trait in patch.community[idx].traits
            size(trait.codedby,1) == 0 && (kill = true) # make sure every trait is coded by at least 1 gene
        end
        patch.community[idx].size <= 0 && (kill = true)
#        seedsize = patch.community[idx].traits[find(x->x.name=="seedsize",patch.community[idx].traits)][1].value
#        maxsize = patch.community[idx].traits[find(x->x.name=="maxsize",patch.community[idx].traits)][1].value
#        seedsize > maxsize && (kill = true)
        if kill
            splice!(patch.community,idx) # or else kill it
            idx -= 1
        end
        idx += 1
    end
end

## Test stuff:
##############
testpatch=Patch(genesis(),293,0.5,0.5,10)
const timesteps=Int64(round(parse(ARGS[1])))
const mutsd=parse(ARGS[2])
for i = 1:timesteps
    checkviability!(testpatch)
    #    size(testpatch.community,1)
    age!(testpatch)
    grow!(testpatch)
    compete!(testpatch)
#    println(i,"\t",size(testpatch.community,1))
    reproduce!(testpatch) # TODO: requires certain amount of resource/bodymass dependent on seedsize!
end
println(testpatch)

# function testscenario(timesteps::Int64=100,npatches::Int64=10)
#     world=Patch[]
#     for patch = 1:npatches
#         push!(world,Patch(genesis(),293,0.5,0.5,10))
#     end
#     for i = 1:timesteps
#         for patch = 1:npatches
#             checkviability!(world[patch])
#             #    size(testpatch.community,1)
#             grow!(world[patch])
#             compete!(world[patch])
#             #    size(testpatch.community,1)
#             reproduce!(world[patch]) # TODO: requires certain amount of resource/bodymass dependent on seedsize!
#             #    size(testpatch.community,1)
#         end
#         i == 1 && (worldstart=world)
#     end
#     worldstart,world
# end

# worldstart,worldend = testscenario()
