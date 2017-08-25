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
    active::Bool
end

mutable struct Gene
    sequence::Array{Char,1} # contains gene base code
    id::String # gene identifier
    codes::Array{Trait,1}
end

struct Chromosome #CAVE! mutable?
    genes::Array{Gene,1} # 1D array of genes
    maternal::Bool # parental origin of chromosome
end


mutable struct Individual
    genome::Array{Chromosome,1} # genome = 2D array of chromosomes (>=1 sets)
    traits::Dict{String,Float64}
    age::Int64
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
        repprob = patch.community[idx].traits["repprob"]
        if !patch.community[idx].isnew && rand() < repprob
            currentmass = patch.community[idx].size
            seedsize = patch.community[idx].traits["seedsize"]
            if currentmass >= 2 * seedsize > 0
                meanoffs = patch.community[idx].traits["reprate"]
                mass = patch.community[idx].size
                metaboffs =  meanoffs * currentmass^(-1/4) * exp(-act/(boltz*temp)) * normconst
                noffs = rand(Poisson(metaboffs))
                for i in 1:noffs
                    ind = deepcopy(patch.community[idx])
                    ind.age = 0
                    ind.isnew = true
                    ind.fitness = 1.0
                    ind.size = ind.traits["seedsize"]
                    mutate!(ind, patch.altitude)
                    push!(patch.community,ind)
                end
            end
        end
        idx += 1
    end
end

function mutate!(ind::Individual, temp::Float64)
    prob = ind.traits["mutprob"]
    for chr in ind.genome
        for gene in chr.genes
            for i in eachindex(gene.sequence)
                if rand() < prob*exp(-act/(boltz*temp)) * normconst
                    newbase = rand(collect("acgt"),1)[1]
                    while newbase == gene.sequence[i]
                        newbase = rand(collect("acgt"),1)[1]
                    end
                    gene.sequence[i] = newbase
                    for trait in gene.codes
                        trait.value == 0 && (trait.value = rand(Normal(0,0.01)))
                        newvalue = trait.value + rand(Normal(0, trait.value/mutscaling)) # new value for trait
                        (newvalue > 1 && contains(trait.name,"prob")) && (newvalue=1)
                        newvalue < 0 && (newvalue=abs(newvalue))
                        while newvalue == 0 #&& contains(trait.name,"mut")
                            newvalue = trait.value + rand(Normal(0, trait.value/mutscaling))
                        end
                        trait.value = newvalue
                    end
                end
            end
        end
    end
    traitdict = chrms2traits(ind.genome)
    ind.traits = traitdict
end

function grow!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        if !patch.community[idx].isnew
            growthrate = patch.community[idx].traits["growthrate"]
            mass = patch.community[idx].size
            newmass = growthrate * mass^(-1/4) * exp(-act/(boltz*temp)) * normconst #CAVE: what to do when negative growth? -> emergent maximum body size!
            if newmass > 0 && mass > 0
                patch.community[idx].size = newmass
            else
                splice!(patch.community, idx)
                idx -= 1
            end
        end
        idx += 1
    end
end

function age!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        if !patch.community[idx].isnew
            ageprob = patch.community[idx].traits["ageprob"]
            mass = patch.community[idx].size
            dieprob = ageprob * mass^(-1/4) * exp(-act/(boltz*temp)) * normconst
            if rand() > (1-dieprob) * patch.community[idx].fitness
                splice!(patch.community, idx)
                idx -= 1
            else
                patch.community[idx].age += 1
            end
        end
        idx += 1
    end
end


function establish!(patch::Patch) #TODO!
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        if patch.community[idx].isnew
            tempopt = patch.community[idx].traits["tempopt"]
            tempbreadth = patch.community[idx].traits["tempbreadth"]
            if abs(temp-tempopt) > tempbreadth
                splice!(patch.community, idx)
                idx -= 1
            else
                patch.community[idx].isnew = false
                fitness = 1 - (abs(temp-tempopt))/tempbreadth
                fitness > 1 && (fitness = 1)
                fitness < 0 && (fitness = 0)
                patch.community[idx].fitness = fitness
            end
        end
        idx += 1
    end
end

#TODO: Dispersal

function createtraits(traitnames::Array{String,1})
    traits = Trait[]
    for name in traitnames
        if contains(name,"rate")
            push!(traits,Trait(name,rand()*10,[],true))
        elseif contains(name, "temp") && contains(name, "opt")
            push!(traits,Trait(name,rand(Normal(298,5)),[],true)) #CAVE: code values elsewhere?
        elseif contains(name, "breadth")
            push!(traits,Trait(name,abs(rand(Normal(0,5))),[],true)) #CAVE: code values elsewhere?
        elseif contains(name, "mut")
            push!(traits,Trait(name,abs(rand(Normal(0,0.01))),[],true)) #CAVE: code values elsewhere?
        else
            push!(traits,Trait(name,rand(),[],true))
        end
    end
    traits
end

function creategenes(ngenes::Int64,traits::Array{Trait,1})
    genes = Gene[]
    viable = false
    while !viable
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
        viable = true
        for trait in traits
            size(trait.codedby,1) == 0 && (viable = false) # make sure every trait is coded by at least 1 gene
        end
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
                push!(chromosomes, Chromosome(genes[1:chrsplits[chr]],rand([false,true])))
            elseif chr==nchrs # last chromosome
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):end],rand([false,true])))
            else
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):chrsplits[chr]],rand([false,true])))
            end
        end
    else # only one chromosome
        chromosomes = [Chromosome(genes,rand([false,true]))]
    end
    secondset = deepcopy(chromosomes)
    for chrm in secondset
        chrm.maternal = !chrm.maternal
    end
    append!(chromosomes,secondset)
    chromosomes
end

function activategenes!(chrms::Array{Chromosome,1}) #ERROR: ArgumentError: range must be non-empty

    genes = Gene[]
    for chrm in chrms
        append!(genes,chrm.genes)
    end
    traits = Trait[]
    for gene in genes
        append!(traits,gene.codes)
    end
    traits = unique(traits)
    traitnames = map(x->x.name,traits)
    traitnames = unique(traitnames)
    for name in traitnames
        idxs = findin(traits,[name])
        map(x->traits[x].active = false,idxs)
        map(x->traits[x].active = true,rand(idxs))
    end
end

function chrms2traits(chrms::Array{Chromosome,1})
    genes = Gene[]
    for chrm in chrms
        append!(genes,chrm.genes)
    end
    traits = Trait[]
    for gene in genes
        append!(traits,gene.codes)
    end
    traits = unique(traits)
    traitdict = Dict{String,Float64}()
    for trait in traits
        trait.active && (traitdict[trait.name] = trait.value)
    end
    traitdict
end


function genesis(ninds::Int64=1000, maxgenes::Int64=20, maxchrs::Int64=5,
                 traitnames::Array{String,1} = ["ageprob",
                                                "growthrate",
                                                "mutprob",
                                                "repprob",
                                                "reprate",
                                                "seedsize",
                                                "sexprob",
                                                "tempbreadth",
                                                "tempopt"]) # minimal required traitnames
    community = Individual[]
    for ind in 1:ninds
        ngenes = rand(1:maxgenes)
        nchrs = rand(1:maxchrs)
        traits = createtraits(traitnames)
        genes = creategenes(ngenes,traits)
        chromosomes = createchrs(nchrs,genes)
        traitdict = chrms2traits(chromosomes)
        push!(community, Individual(chromosomes,traitdict,0,true,1.0,rand()))
    end
    community
end

function checkviability!(patch::Patch) # may consider additional rules... # maybe obsolete anyhow...
    idx=1
    while idx <= size(patch.community,1)
        kill = false
        patch.community[idx].size <= 0 && (kill = true)
        0 in collect(values(patch.community[idx].traits)) && (kill = true)
        if kill
            splice!(patch.community,idx) # or else kill it
            idx -= 1
        end
        idx += 1
    end
end

## Test stuff:
##############
testpatch=Patch(genesis(),293,0.5,0.5,100)
startpatch=deepcopy(testpatch)
const timesteps=1000 #Int64(round(parse(ARGS[1])))
const mutscaling=50#parse(ARGS[2])
for i = 1:timesteps
    checkviability!(testpatch)
    size(testpatch.community,1)
    establish!(testpatch)
    size(testpatch.community,1)
    age!(testpatch)
    size(testpatch.community,1)
    grow!(testpatch)
    size(testpatch.community,1)
    compete!(testpatch)
    size(testpatch.community,1)
    reproduce!(testpatch) # TODO: requires certain amount of resource/bodymass dependent on seedsize!
    size(testpatch.community,1)
end

mean(map(x->x.traits["mutprob"],testpatch.community))

mean(map(x->x.traits["repprob"],testpatch.community))

mean(map(x->x.traits["ageprob"],testpatch.community))

mean(map(x->x.traits["tempbreadth"],testpatch.community))

mean(map(x->x.traits["reprate"],testpatch.community))

mean(map(x->x.traits["growthrate"],testpatch.community))

mean(map(x->x.traits["tempopt"],testpatch.community))

mean(map(x->x.traits["seedsize"],testpatch.community))

minimum(map(x->x.traits["seedsize"],testpatch.community))

maximum(map(x->x.traits["seedsize"],testpatch.community))

