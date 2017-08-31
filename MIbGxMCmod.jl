#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## contains all custom types and functions necessary for eco-evo-env speciation island model


module MIbGxMCmod


using Distributions


export Patch, # types
    genesis, checkviability!, establish!, age!, grow!, disperse!, compete!, reproduce! # functions


const boltz = 1.38064852e-23 #  J/K = m2⋅kg/(s2⋅K)
const act = 1e-19 # activation energy /J, ca. 0.63eV - Brown et al. 2004
const normconst = 1e10 # normalization constant to get biologically realistic orders of magnitude
const mutscaling=50 #parse(ARGS[2])
const sexualreproduction = true


## Types:
#########

mutable struct Trait
    name::String
    value::Float64
    active::Bool
end

mutable struct Gene
    sequence::Array{Char,1}
    id::String
    codes::Array{Trait,1}
end

mutable struct Chromosome
    genes::Array{Gene,1}
    maternal::Bool # parental origin of chromosome
end

mutable struct Individual
    genome::Array{Chromosome,1}
    traits::Dict{String,Float64}
    age::Int64
    isnew::Bool # indicator whether individual is new to a patch (after dispersal or birth)
    fitness::Float64 # rate scaling factor
    size::Float64 # body mass
end

mutable struct Patch
    id::Int64
    location::Tuple{Float64,Float64}
    altitude::Float64 # corresponds to T
    area::Float64
    isisland::Bool
    nichea::Float64 # additional niches,
    nicheb::Float64 # e.g. precipitation
    community::Array{Individual,1}
end

Patch(id,location,altitude,area,isisland,nichea,nicheb) = Patch(id,location,altitude,area,isisland,nichea,nicheb,Individual[])
Patch(id,location,altitude,area,isisland,nichea) = Patch(id,location,altitude,area,isisland,nichea,0,Individual[])
Patch(id,location,altitude,area,isisland) = Patch(id,location,altitude,area,isisland,0,0,Individual[])
Patch(id,location,altitude,area) = Patch(id,location,altitude,area,false,0,0,Individual[])
Patch(id,location,altitude) = Patch(id,location,altitude,100,false,0,0,Individual[])
Patch(id,location) = Patch(id,location,298,0,0,Individual[])


## Methods/Functions:
#####################

function meiosis(genome::Array{Chromosome,1},maternal::Bool) # TODO: include further dynamics, errors...
    firstset = find(x->x.maternal,genome)
    secondset = find(x->!x.maternal,genome)
    size(firstset,1) != size(secondset,1) && return Chromosome[] # CAVE: more elegant solution...
    gameteidxs = []
    for i in eachindex(firstset)
        push!(gameteidxs,rand([firstset[i],secondset[i]]))
    end
    gamete = deepcopy(genome[gameteidxs]) #TODO somewhere here: crossing over!
    map(x->x.maternal=maternal,gamete)
    gamete
end

function mutate!(ind::Individual, temp::Float64)
    prob = ind.traits["mutprob"]
    for chrm in ind.genome
        for gene in chrm.genes
            for i in eachindex(gene.sequence)
                if rand() <= prob*exp(-act/(boltz*temp)) * normconst
                    newbase = rand(collect("acgt"),1)[1]
                    while newbase == gene.sequence[i]
                        newbase = rand(collect("acgt"),1)[1]
                    end
                    gene.sequence[i] = newbase
                    for trait in gene.codes
                        trait.value == 0 && (trait.value = rand(Normal(0,0.01)))
                        newvalue = trait.value + rand(Normal(0, trait.value/mutscaling)) # new value for trait
                        (newvalue > 1 && contains(trait.name,"prob")) && (newvalue=1)
                        (newvalue > 1 && contains(trait.name,"reptol")) && (newvalue=1)
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

function activategenes!(chrms::Array{Chromosome,1})
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
        idxs = find(x->x.name==name,traits)
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

function checkviability!(world::Array{Patch,1})
    for patch in world
        checkviability!(patch) # pmap(checkviability!,patch) ???
    end
end

function establish!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        hastemptol = haskey(patch.community[idx].traits,"temptol")
        hastempopt = haskey(patch.community[idx].traits,"tempopt")
        if !hastemptol || !hastempopt
            splice!(patch.community, idx) # kill it!
            idx -= 1
        elseif patch.community[idx].isnew
            tempopt = patch.community[idx].traits["tempopt"]
            temptol = patch.community[idx].traits["temptol"]
            if abs(temp-tempopt) > temptol
                splice!(patch.community, idx)
                idx -= 1
            else
                patch.community[idx].isnew = false
                fitness = 1 - (abs(temp-tempopt))/temptol
                fitness > 1 && (fitness = 1)
                fitness < 0 && (fitness = 0)
                patch.community[idx].fitness = fitness
            end
        end
        idx += 1
    end
end


function establish!(world::Array{Patch,1})
    for patch in world
        establish!(patch) # pmap(!,patch) ???
    end
end

function age!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        hasageprob = haskey(patch.community[idx].traits,"ageprob")
        if !hasageprob
            splice!(patch.community, idx)
            idx -= 1
        else
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
        end
        idx += 1
    end
end

function age!(world::Array{Patch,1})
    for patch in world
        age!(patch) # pmap(!,patch) ???
    end
end

function grow!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        hasgrowthrate = haskey(patch.community[idx].traits,"growthrate")
        if !hasgrowthrate
            splice!(patch.community, idx)
            idx -= 1
        else
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
        end
        idx += 1
    end
end

function grow!(world::Array{Patch,1})
    for patch in world
        grow!(patch) # pmap(!,patch) ???
    end
end

function disperse!(world::Array{Patch,1})
    for patch in world
        idx = 1
        while idx <= size(patch.community,1)
            hasdispmean = haskey(patch.community[idx].traits,"dispmean")
            hasdispprob = haskey(patch.community[idx].traits,"dispprob")
            hasdispshape = haskey(patch.community[idx].traits,"dispshape")
            if !hasdispmean || !hasdispprob || !hasdispshape
                splice!(patch.community,idx)
                idx -= 1
            elseif !patch.community[idx].isnew && rand() <= patch.community[idx].traits["dispprob"]
                dispmean = patch.community[idx].traits["dispmean"]
                dispshape = patch.community[idx].traits["dispshape"]
                patch.community[idx].isnew = true
                indleft = splice!(patch.community,idx)
                xdir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # scaling so that geometric mean...
                ydir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # ...follows original distribution
                xdest = patch.location[1]+xdir # CAVE: might be other way around + border conditions
                ydest = patch.location[1]+ydir # CAVE: might be other way around + border conditions
                targets = unique([(floor(xdest),floor(ydest)),(ceil(xdest),floor(ydest)),(ceil(xdest),ceil(ydest)),(floor(xdest),ceil(ydest))])
                possdests = find(x->in(x.location,targets),world)
                if size(possdests,1) > 0 # if no viable target patch, individual dies
                    size(possdests,1) > 1 ? (destination=rand(possdests)) : (destination = possdests[1])
                    push!(world[destination].community,indleft)
                end
                idx -= 1
            end
            idx += 1
        end
    end
end

function compete!(patch::Patch)
    while sum(map(x->x.size,patch.community)) > patch.area # occupied area larger than available
        victim = rand(1:size(patch.community,1))
        splice!(patch.community, victim)
    end
end

function compete!(world::Array{Patch,1})
    for patch in world
        compete!(patch) # pmap(!,patch) ???
    end
end

function reproduce!(patch::Patch) #TODO: refactorize!
    idx = 1
    temp = patch.altitude
    while idx <= size(patch.community,1)
        hasrepprob = haskey(patch.community[idx].traits,"repprob")
        hasreptol = haskey(patch.community[idx].traits,"reptol")
        hasreprate = haskey(patch.community[idx].traits,"reprate")
        hasseedsize = haskey(patch.community[idx].traits,"seedsize")
        hasmutprob = haskey(patch.community[idx].traits,"mutprob")
        if !hasrepprob || !hasreprate || !hasreptol || !hasmutprob || !hasseedsize 
            splice!(patch.community, idx)
            idx -= 1
        elseif !patch.community[idx].isnew && rand() <= patch.community[idx].traits["repprob"]
            currentmass = patch.community[idx].size
            seedsize = patch.community[idx].traits["seedsize"]
            if currentmass >= 2 * seedsize > 0 #CAVE: set rule for this, now arbitrary -> reprod. size?
                meanoffs = patch.community[idx].traits["reprate"]
                reptol = patch.community[idx].traits["reptol"]
                mass = patch.community[idx].size
                metaboffs =  meanoffs * currentmass^(-1/4) * exp(-act/(boltz*temp)) * normconst
                noffs = rand(Poisson(metaboffs))
                if sexualreproduction
                    posspartners = find(x->(1/reptol)*mass>=x.size>=reptol*mass,patch.community)
                    partneridx = rand(posspartners)
                    while size(posspartners,1) > 1 && partneridx == idx #CAVE: unless self-repr. is allowed
                        partneridx = rand(posspartners)
                    end
                    partnergenome = meiosis(patch.community[partneridx].genome, false) #CAVE: maybe move inside offspring loop?
                    mothergenome = meiosis(patch.community[idx].genome, true)
                    for i in 1:noffs
                        (size(partnergenome,1) < 1 || size(mothergenome,1) < 1) && continue
                        genome = deepcopy([partnergenome...,mothergenome...])
                        activategenes!(genome)
                        traits = chrms2traits(genome)
                        age = 0
                        isnew = true
                        fitness = 1.0
                        newsize = seedsize
                        ind = Individual(genome,traits,age,isnew,fitness,newsize)
                        !haskey(ind.traits,"mutprob") && continue # no mutation, no life ;)
                        mutate!(ind, patch.altitude)
                        push!(patch.community,ind)
                    end
                # else #...?    
                end
            end
        end
        idx += 1
    end
end

function reproduce!(world::Array{Patch,1})
    for patch in world
        reproduce!(patch) # pmap(!,patch) ???
    end
end

function createtraits(traitnames::Array{String,1})
    traits = Trait[]
    for name in traitnames
        if contains(name,"rate")
            push!(traits,Trait(name,rand()*10,true))
        elseif contains(name, "temp") && contains(name, "opt")
            push!(traits,Trait(name,rand(Normal(298,5)),true)) #CAVE: code values elsewhere?
        elseif contains(name, "tol") && !contains(name, "rep")
            push!(traits,Trait(name,abs(rand(Normal(0,5))),true)) #CAVE: code values elsewhere?
        elseif contains(name, "mut")
            push!(traits,Trait(name,abs(rand(Normal(0,0.01))),true)) #CAVE: code values elsewhere?
        else
            push!(traits,Trait(name,rand(),true))
        end
    end
    traits
end

function creategenes(ngenes::Int64,traits::Array{Trait,1})
    genes = Gene[]
    for i in 1:ngenes
        sequence = collect("acgt"^5) # arbitrary start sequence
        id = randstring(8)
        codesfor = Trait[]
        push!(genes,Gene(sequence,id,codesfor))
    end
    for trait in traits
        ncodinggenes = rand(Poisson(1))
        while ncodinggenes < 1  # make sure every trait is coded by at least 1 gene
            ncodinggenes = rand(Poisson(1))
        end
        codinggenes = rand(genes,ncodinggenes)
        for gene in codinggenes
            push!(gene.codesfor,trait)
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

function genesis(ninds::Int64=1000, meangenes::Int64=20, meanchrs::Int64=5,
                 traitnames::Array{String,1} = ["ageprob",
                                                "dispmean",
                                                "dispprob",
                                                "dispshape",
                                                "growthrate",
                                                "mutprob",
                                                "repprob",
                                                "reprate",
                                                "reptol",
                                                "seedsize",
                                                "temptol",
                                                "tempopt"]) # minimal required traitnames
    community = Individual[]
    for ind in 1:ninds
        ngenes = rand(Poisson(meangenes))
        nchrs = rand(Poisson(meanchrs))
        traits = createtraits(traitnames)
        genes = creategenes(ngenes,traits)
        chromosomes = createchrs(nchrs,genes)
        traitdict = chrms2traits(chromosomes)
        push!(community, Individual(chromosomes,traitdict,0,true,1.0,rand()))
    end
    community
end

    
end
