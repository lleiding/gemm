#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## contains all custom types and functions necessary for eco-evo-env speciation island model


module MIbGxMCmod


using Distributions#, Plots


export Patch, # types
    setupdatadir, recordcolonizers, readmapfile, writerawdata, writedata, analysis,
    checkviability!, establish!, survive!, grow!, disperse!, compete!, reproduce!, mutate!,
    createworld, updateworld!, endsim, alldead # functions


const boltz = 1.38064852e-23 # J/K = m2⋅kg/(s2⋅K)
const act = 1e-19 # activation energy /J, ca. 0.63eV - Brown et al. 2004
const growthrate = exp(25.2) # global base growth/biomass production from Brown et al. 2004
const mortality = exp(22) # global base mortality from Brown et al. 2004 is 26.3, but competition and dispersal introduce add. mort.
const fertility = exp(30.0) # global base reproduction rate 23.8 from Brown et al. 2004, alternatively 25.0
const phylconstr = 10 #parse(ARGS[2])
# const meangenes = 20 # mean number of genes per individual
const mutationrate = 1e-3 * 0.3e11 # 1 base in 1000, correction factor for metabolic function
const isolationweight = 3 # additional distance to be crossed when dispersing from or to isolated patches
const cellsize = 1e8 # ca. 100 tons carrying cap.
const maxdispmean = 10 # maximum mean dispersal distance
const genelength = 20 # sequence length of genes


## Types:
#########

mutable struct Trait
    name::String
    value::Float64
end

mutable struct Gene
    sequence::String
    codes::Array{Trait, 1}
end

mutable struct Chromosome
    genes::Array{Gene, 1}
    maternal::Bool # parental origin of chromosome
end

mutable struct Individual
    lineage::String
    genome::Array{Chromosome, 1}
    traits::Dict{String, Float64}
    age::Int64
    isnew::Bool # indicator whether individual is new to a patch (after dispersal or birth)
    fitness::Float64 # rate scaling factor
    size::Float64 # body mass
end

mutable struct Patch
    id::Int64
    location::Tuple{Int64, Int64}
    altitude::Float64 # corresponds to T
    area::Float64
    isisland::Bool
    nichea::Float64 # additional niches,
    nicheb::Float64 # e.g. precipitation
    community::Array{Individual, 1}
    isolated::Bool
end

Patch(id,location,altitude,area,isisland,nichea,nicheb,community) = Patch(id,location,altitude,area,isisland,nichea,nicheb,community,false)
Patch(id,location,altitude,area,isisland,nichea,nicheb) = Patch(id,location,altitude,area,isisland,nichea,nicheb,Individual[],false)
Patch(id,location,altitude,area,isisland,nichea) = Patch(id,location,altitude,area,isisland,nichea,0,Individual[],false)
Patch(id,location,altitude,area,isisland) = Patch(id,location,altitude,area,isisland,0,0,Individual[],false)
Patch(id,location,altitude,area) = Patch(id,location,altitude,area,false,0,0,Individual[],false)
Patch(id,location,altitude) = Patch(id,location,altitude,100,false,0,0,Individual[],false)
Patch(id,location) = Patch(id,location,298,100,false,0,0,Individual[],false)


## Methods/Functions:
#####################
"""
    endsim([x], [t])
End simulation for reason `x`. Optionally provide the timestep `t`.
`x` defaults to "alldead".
"""
function endsim(reason::String = "alldead", t::Int64 = 0)
    if reason == "alldead"
        print("All individuals dead. Terminating simulation")
        t > 0 ? println("at t = $t.") : println(".")
        exit()
    end
end

"""
    alldead(w)
Check whether there are any individuals left in the world `w`.
"""
function alldead(world::Array{Patch,1})
    totalcom = map(x -> length(x.community), world)
    sum(totalcom) == 0
end

function meiosis(genome::Array{Chromosome,1},maternal::Bool) # TODO: include further dynamics, errors...
    firstset = find(x->x.maternal,genome)
    secondset = find(x->!x.maternal,genome)
    size(firstset,1) != size(secondset,1) && return Chromosome[] # CAVE: more elegant solution...
    gameteidxs = []
    for i in eachindex(firstset)
        push!(gameteidxs,rand([firstset[i],secondset[i]]))
    end
    gamete = Chromosome[]
    for i in gameteidxs
        push!(gamete, Chromosome(genome[i].genes, maternal))
    end
    deepcopy(gamete)
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
    traitdict = Dict{String,Float64}()
    for traitname in unique(map(x->x.name,traits))
        traitdict[traitname] = mean(map(x->x.value,filter(x->x.name==traitname,traits)))
    end
    traitdict
end

function mutate!(ind::Individual, temp::Float64, settings::Dict{String,Any})
    prob = ind.traits["mutprob"]
    for chrm in ind.genome
        for idx in eachindex(chrm.genes)
            charseq = collect(chrm.genes[idx].sequence)
            for i in eachindex(charseq)
                if rand() <= prob * exp(-act/(boltz*temp))
                    newbase = rand(collect("acgt"),1)[1]
                    while newbase == charseq[i]
                        newbase = rand(collect("acgt"),1)[1]
                    end
                    charseq[i] = newbase
                    for trait in chrm.genes[idx].codes
                        (contains(trait.name, "mutprob") && mutationrate != 0) && continue
                        contains(trait.name, "reptol") && settings["tolerance"] != "evo" && continue # MARK CAVE!
                        oldvalue = trait.value
                        contains(trait.name, "tempopt") && (oldvalue -= 273)
                        while oldvalue <= 0 # make sure sd of Normal dist != 0
                            oldvalue = abs(rand(Normal(0,0.01)))
                        end
                        newvalue = oldvalue + rand(Normal(0, oldvalue/phylconstr)) # CAVE: maybe handle temp + prec separately
                        newvalue < 0 && (newvalue=abs(newvalue))
                        (newvalue > 1 && contains(trait.name,"prob")) && (newvalue=1)
                        while newvalue <= 0 #&& contains(trait.name,"mut")
                            newvalue = trait.value + rand(Normal(0, trait.value/phylconstr))
                        end
                        contains(trait.name, "tempopt") && (newvalue += 273)
                        trait.value = newvalue
                    end
                end
            end
            chrm.genes[idx].sequence = String(charseq)
        end
    end
    ind.traits = chrms2traits(ind.genome)
end

function mutate!(patch::Patch, settings::Dict{String,Any})
    for ind in patch.community
        ind.age == 0 && mutate!(ind, patch.altitude, settings)
    end
end

function mutate!(world::Array{Patch, 1}, settings::Dict{String,Any})
    for patch in world
        patch.isisland && mutate!(patch, settings)
    end
end

function checkviability!(patch::Patch) # may consider additional rules... # maybe obsolete anyhow...
    idx=1
    while idx <= size(patch.community,1)
        dead = false
        patch.community[idx].size <= 0 && (dead = true)
        any(collect(values(patch.community[idx].traits)) .<= 0) && (dead = true)
        patch.community[idx].traits["repsize"] <= patch.community[idx].traits["seedsize"] && (dead = true)
        patch.community[idx].fitness <= 0 && (dead = true)
        if dead
            splice!(patch.community,idx)
            continue
        end
        idx += 1
    end
end

function checkviability!(world::Array{Patch,1})
    for patch in world
        patch.isisland && checkviability!(patch) # pmap(checkviability!,patch) ???
    end
end

function traitsexist(ind::Individual, traitnames::Array{String, 1})
    for trait in traitnames
        !haskey(ind.traits, trait) && return false
    end
    true
end

function gausscurve(b::Float64, c::Float64, x::Float64, a::Float64=1.0)
    if c != 0 && a != 1.0
        a = 1 / (c * sqrt(2 * pi))
        y = a * exp(-(x-b)^2/(2*c^2))
    elseif c != 0
        y = a * exp(-(x-b)^2/(2*c^2))
    else
        y = 0
    end
end

"""
    establish!(p, n)
establishment of individuals in patch `p`. Sets fitness scaling parameter
according to adaptation to number `n` niches of the surrounding environment.
"""
function establish!(patch::Patch, nniches::Int64=1)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        if patch.community[idx].isnew || patch.community[idx].age == 0
            fitness = 1
            if !traitsexist(patch.community[idx], ["temptol", "tempopt"]) # CAVE: should check for more traits!
                splice!(patch.community, idx) # kill it!
                continue
            end
            opt = patch.community[idx].traits["tempopt"]
            tol = patch.community[idx].traits["temptol"]
            fitness *= gausscurve(opt, tol, temp, 0.0)
            fitness > 1 && (fitness = 1) # should be obsolete
            fitness < 0 && (fitness = 0) # should be obsolete
            if nniches >= 2
                if !traitsexist(patch.community[idx], ["prectol", "precopt"]) # CAVE: should check for more traits!
                    splice!(patch.community, idx) # kill it!
                    continue
                end
                opt = patch.community[idx].traits["precopt"]
                tol = patch.community[idx].traits["prectol"]
                fitness *= gausscurve(opt, tol, patch.nichea, 0.0)
                fitness > 1 && (fitness = 1) # should be obsolete
                fitness < 0 && (fitness = 0) # should be obsolete
            end
            if nniches == 3
                opt = patch.community[idx].traits["nicheopt"]
                tol = patch.community[idx].traits["nichetol"]
                fitness *= gausscurve(opt, tol, patch.nicheb, 0.0)
                fitness > 1 && (fitness = 1) # should be obsolete
                fitness < 0 && (fitness = 0) # should be obsolete
            end  
            patch.community[idx].fitness = fitness
            patch.community[idx].isnew = false
        end
        idx += 1
    end
end


function establish!(world::Array{Patch,1}, nniches::Int64=1)
    for patch in world
        patch.isisland && establish!(patch, nniches) # pmap(!,patch) ???
    end
end

"""
    survive!(p)
density independent survival of individuals in patch `p`
"""
function survive!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        if !patch.community[idx].isnew
            mass = patch.community[idx].size
            deathrate = mortality * mass^(-1/4) * exp(-act/(boltz*temp))
            dieprob = 1 - exp(-deathrate)
            if rand() < dieprob
                splice!(patch.community, idx)
                continue
            else
                patch.community[idx].age += 1
            end
        end
        idx += 1
    end
end

function survive!(world::Array{Patch,1})
    for patch in world
        patch.isisland && survive!(patch) # pmap(!,patch) ???
    end
end

"""
    grow!(p)
Growth of individuals in patch `p`
"""
function grow!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        if !traitsexist(patch.community[idx], ["repsize"])
            splice!(patch.community, idx)
            continue
        elseif !patch.community[idx].isnew
            repsize = patch.community[idx].traits["repsize"]
            mass = patch.community[idx].size
            if mass <= repsize # stop growth if reached repsize 
                growth = growthrate * mass^(3/4) * exp(-act/(boltz*temp))
                newmass = mass + growth
                if newmass > 0 && mass > 0
                    patch.community[idx].size = newmass
                else
                    splice!(patch.community, idx)
                    continue
                end
            end
        end
        idx += 1
    end
end

function grow!(world::Array{Patch,1})
    for patch in world
        patch.isisland && grow!(patch) # pmap(!,patch) ???
    end
end

"""
    findisland(w)
within world `w`, find out in which direction from the continent the island(s) lie(s).
"""
function findisland(world::Array{Patch,1})
    xmin = minimum(map(x->x.location[1],world))
    xmax = maximum(map(x->x.location[1],world))
    ymin = minimum(map(x->x.location[2],world))
    ymax = maximum(map(x->x.location[2],world))
    westernborder = filter(x->x.location[1]==xmin,world)
    northernborder = filter(x->x.location[2]==ymax,world)
    easternborder = filter(x->x.location[1]==xmax,world)
    southernborder = filter(x->x.location[2]==ymin,world)
    if all(map(x->x.isisland,westernborder))
        return "west"
    elseif all(map(x->x.isisland,northernborder))
        return "north"
    elseif all(map(x->x.isisland,easternborder))
        return "east"
    elseif all(map(x->x.isisland,southernborder))
        return "south"
    else
        return "none"
    end
end

"""
    checkborderconditions!(w, x, y)
check if coordinates `x` and `y` lie within world `w` and correct if not,
considering defined border conditions.
"""
function checkborderconditions!(world::Array{Patch,1},xdest::Float64,ydest::Float64)
    xmin = minimum(map(x->x.location[1],world))
    xmax = maximum(map(x->x.location[1],world))
    ymin = minimum(map(x->x.location[2],world))
    ymax = maximum(map(x->x.location[2],world))
    xrange = xmax - xmin + 1 # we're counting cells!
    yrange = ymax - ymin + 1 # we're counting cells!
    xshift = xdest - xmin + 1 # 1-based count of cells
    yshift = ydest - ymin + 1 # 1-based count of cells
    xshift > 0 ? outofx = abs(xshift) : outofx = abs(xshift) + 1
    while outofx > xrange
        outofx -= xrange
    end
    outofx -= 1
    yshift > 0 ? outofy = abs(yshift) : outofy = abs(yshift) + 1
    while outofy > yrange
        outofy -= yrange
    end
    outofy -= 1
    islanddirection = findisland(world::Array{Patch,1})
    if islanddirection == "west"
        xdest > xmax && (xdest = xmax - outofx) # east: reflective
        ydest < ymin && (ydest = ymax - outofy) # south: periodic
        ydest > ymax && (ydest = ymin + outofy) # north: periodic
    elseif islanddirection == "north"
        ydest < ymin && (ydest = ymin + outofy) # south: reflective
        xdest < xmin && (xdest = xmax - outofx) # west: periodic
        xdest > xmax && (xdest = xmin + outofx) # east: periodic
    elseif islanddirection == "east"
        xdest < xmin && (xdest = xmin + outofx) # west: reflective
        ydest < ymin && (ydest = ymax - outofy) # south: periodic
        ydest > ymax && (ydest = ymin + outofy) # north: periodic
    elseif islanddirection == "south"
        ydest > ymax && (ydest = ymax + outofy) # north: reflective
        xdest < xmin && (xdest = xmax - outofx) # west: periodic
        xdest > xmax && (xdest = xmin + outofx) # east: periodic
    else
        ydest > ymax && (ydest = ymin + outofy) # north: periodic
        xdest > xmax && (xdest = xmin + outofx) # east: periodic
        ydest < ymin && (ydest = ymax - outofy) # south: periodic
        xdest < xmin && (xdest = xmax - outofx) # west: periodic
    end
    xdest,ydest
end

"""
    disperse!(w)
Dispersal of individuals within world (array of patches) `w`
"""
function disperse!(world::Array{Patch,1}) # TODO: additional border conditions, refoctorize
    colonizers = ""
    for patch in world
        idx = 1
        while idx <= size(patch.community,1)
            if !traitsexist(patch.community[idx], ["dispmean", "dispshape"])
                splice!(patch.community,idx)
                continue
            elseif !patch.community[idx].isnew && patch.community[idx].age == 0
                dispmean = patch.community[idx].traits["dispmean"]
                dispshape = patch.community[idx].traits["dispshape"]
                xdir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # scaling so that geometric mean...
                ydir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # ...follows original distribution
                xdest = patch.location[1]+xdir
                ydest = patch.location[2]+ydir
                # !patch.isisland && checkborderconditions!(world,xdest,ydest)
                targets = unique([(floor(xdest),floor(ydest)),(ceil(xdest),floor(ydest)),(ceil(xdest),ceil(ydest)),(floor(xdest),ceil(ydest))])
                possdests = find(x->in(x.location,targets),world)
                filter!(x -> world[x].isisland, possdests) # disperse only to islands
                if patch.isisland
                    indleft = splice!(patch.community,idx) # only remove individuals from islands!
                end
                if size(possdests,1) > 0 # if no viable target patch, individual dies
                    if !patch.isisland
                        indleft = deepcopy(patch.community[idx])
                    end
                    indleft.isnew = true
                    destination = rand(possdests)
                    originisolated = patch.isolated && rand(Logistic(dispmean,dispshape)) <= isolationweight # additional roll for isolated origin patch
                    targetisolated = world[destination].isolated && rand(Logistic(dispmean,dispshape)) <= isolationweight # additional roll for isolated target patch
                    (!originisolated && !targetisolated) && push!(world[destination].community, indleft) # new independent individual
                    !patch.isisland && world[destination].isisland && (colonizers *= indleft.lineage * ",")
                end
                patch.isisland && (idx -= 1)
            end
            idx += 1
        end
    end
    colonizers
end

function compete!(patch::Patch)
    sort!(patch.community, by = x -> x.fitness)
    while sum(map(x -> x.size, patch.community)) >= patch.area # occupied area larger than available
        victim = rand(Geometric()) + 1 # fitness sorted?
        victim > length(patch.community) && (victim = length(patch.community))
        splice!(patch.community, victim)
    end
end

function compete!(world::Array{Patch,1})
    for patch in world
        patch.isisland && compete!(patch) # pmap(!,patch) ???
    end
end

function iscompatible(mate::Individual, ind::Individual)
    tolerance = ind.traits["reptol"]
    indgene = ""
    for chrm in ind.genome
        for gene in chrm.genes
            any(x -> x.name == "compat", gene.codes) && (indgene *= gene.sequence)
        end
    end
    mategene = ""
    for chrm in mate.genome
        for gene in chrm.genes
            any(x -> x.name == "compat", gene.codes) && (mategene *= gene.sequence)
        end
    end
    basediffs = 0
    for i in eachindex(indgene)
        try
            indgene[i] != mategene[i] && (basediffs += 1) # alternatively use bioinformatics tools
        catch # e.g., in case of differently long genes
            basediffs += 1
        end
    end
    seqidentity = 1 - (basediffs / length(indgene))
    seqidentity < tolerance && return false # CAVE: maybe compare w/ indgenome?
    true
end

function findposspartners(world::Array{Patch,1}, ind::Individual, location::Tuple{Int64, Int64})
    ind.isnew = true
    radius = floor(ind.traits["repradius"] + 0.5) # CAVE: to account for cell width ... or not??
    coordinates = Tuple[]
    for x = -radius:radius, y = -radius:radius
        sqrt(x^2 + y^2) <= radius && push!(coordinates, (x + location[1], y + location[2]))
    end
    coordinates = map(k -> checkborderconditions!(world, k[1], k[2]), coordinates)
    order = map(x -> x.^2, coordinates) |> x -> map(sum, x) |> x -> map(sqrt, x) |> sortperm
    coordinates = coordinates[order]
    posspartners = Individual[]
    idx = 1 # check patches in order of increasing distance
    while length(posspartners) == 0
        idx > length(coordinates) && break
        targetpatch = filter(l -> l.location == coordinates[idx], world)
        length(targetpatch) >= 1 ? community = targetpatch[1].community : community = Individual[]
        append!(posspartners, community)
        filter!(k -> k.age > 0, posspartners)
        filter!(k -> k.lineage == ind.lineage, posspartners)
        filter!(k -> traitsexist(k, ["repsize"]), posspartners)
        filter!(k -> k.size >= k.traits["repsize"], posspartners)
        filter!(k -> !k.isnew, posspartners) # filter out mating individual
        filter!(k -> iscompatible(k, ind), posspartners)
        idx += 1
    end
    ind.isnew = false
    posspartners
end

"""
    reproduce!(w, p)
Reproduction of individuals in a patch `p` whithin a world (array of patches) `w`
"""
function reproduce!(world::Array{Patch,1}, patch::Patch) #TODO: refactorize!
    idx = 1
    temp = patch.altitude
    seedbank = Individual[]
    while idx <= size(patch.community,1)
        if !traitsexist(patch.community[idx], ["repradius", "repsize", "reptol", "seedsize", "mutprob"])
            splice!(patch.community, idx)
            continue
        elseif !patch.community[idx].isnew && patch.community[idx].age > 0
            currentmass = patch.community[idx].size
            seedsize = patch.community[idx].traits["seedsize"]
            if currentmass >= patch.community[idx].traits["repsize"]
                reptol = patch.community[idx].traits["reptol"]
                metaboffs = fertility * currentmass^(-1/4) * exp(-act/(boltz*temp))
                noffs = rand(Poisson(metaboffs * patch.community[idx].fitness)) # add some stochasticity
                posspartners = findposspartners(world, patch.community[idx], patch.location) # this effectively controls frequency of reproduction
                # length(posspartners) == 0 && push!(posspartners, patch.community[idx]) # selfing if no partners # CAVE!
                if length(posspartners) > 0
                    partner = rand(posspartners)
                    parentmass = currentmass - noffs * seedsize # subtract offspring mass from parent
                    if parentmass <= 0
                        idx += 1 #splice!(patch.community, idx)
                        continue
                    else
                        patch.community[idx].size = parentmass
                    end
                    for i in 1:noffs # pmap? this loop could be factorized!
                        partnergenome = meiosis(partner.genome, false) # offspring have different genome!
                        mothergenome = meiosis(patch.community[idx].genome, true)
                        (length(partnergenome) < 1 || length(mothergenome) < 1) && continue
                        genome = vcat(partnergenome,mothergenome)
                        traits = chrms2traits(genome)
                        age = 0
                        isnew = false
                        fitness = 1.0
                        newsize = seedsize
                        ind = Individual(patch.community[idx].lineage, genome,traits,age,isnew,fitness,newsize)
                        push!(seedbank ,ind) # maybe actually deepcopy!?
                    end
                end
            end
        end
        idx += 1
    end
    append!(patch.community, seedbank)
end

function reproduce!(world::Array{Patch,1})
    for patch in world
        patch.isisland && reproduce!(world, patch) # pmap(!,patch) ???
    end
end

function createtraits(traitnames::Array{String,1}, settings::Dict{String,Any}) #TODO: this is all very ugly. (case/switch w/ v. 2.0+?)
    traits = Trait[]
    seedsize = exp(-7 + 20 * rand()) # corresponds to 1mg to 22kg
    repsize = exp(0 + 17 * rand()) # 1g to 24t
    while repsize <= seedsize
        seedsize = exp(-7 + 20 * rand())
        repsize = exp(0 + 17 * rand())
    end
    for name in traitnames
        if contains(name,"rate")
            push!(traits,Trait(name,rand()*100))
        elseif contains(name,"dispshape")
            push!(traits, Trait(name, rand() * maxdispmean))
        elseif contains(name, "tempopt")
            push!(traits,Trait(name, rand() * 40 + 273)) #CAVE: code values elsewhere?
        elseif contains(name, "temptol")
            push!(traits,Trait(name,rand()*5)) #CAVE: code values elsewhere?
        elseif contains(name, "mut")
            mutationrate == 0 ? push!(traits,Trait(name,rand())) : push!(traits,Trait(name,mutationrate)) #CAVE: code values elsewhere?
        elseif contains(name, "repsize")
            push!(traits,Trait(name,repsize)) #CAVE: code values elsewhere?
        elseif contains(name, "seedsize")
            push!(traits,Trait(name,seedsize)) #CAVE: code values elsewhere?
        elseif contains(name, "precopt")
            push!(traits, Trait(name, rand() * 10))
        elseif contains(name, "prectol")
            push!(traits, Trait(name, rand()))
        elseif contains(name, "reptol")
            if settings["tolerance"] == "high"
                push!(traits,Trait(name,0.8))
            elseif settings["tolerance"] == "low"
                push!(traits,Trait(name,0.95)) #CAVE: code values elsewhere?
            else
                push!(traits,Trait(name, 0.5 + rand() * 0.5))
            end
        else
            push!(traits,Trait(name,rand()))
        end
    end
    traits
end

function creategenes(ngenes::Int64,traits::Array{Trait,1})
    genes = Gene[]
    for i in 1:ngenes
        sequence = String(rand(collect("acgt"), genelength)) # arbitrary start sequence
        codesfor = Trait[]
        push!(genes,Gene(sequence, codesfor))
    end
    for trait in traits
        ncodinggenes = rand(Poisson(1))
        while ncodinggenes < 1  # make sure every trait is coded by at least 1 gene
            ncodinggenes = rand(Poisson(1))
        end
        codinggenes = rand(genes,ncodinggenes)
        for gene in codinggenes
            push!(gene.codes,trait)
        end
    end
    if !any(map(x -> length(x.codes) == 0, genes)) # make sure there is a neutral gene!
        push!(genes, Gene(String(rand(collect("acgt"), genelength)), Trait[]))
    end
    push!(genes, Gene(String(rand(collect("acgt"), genelength)), [Trait("compat", 1.0)])) # create extra compatibility gene
    genes
end

function createchrs(nchrs::Int64,genes::Array{Gene,1})
    ngenes=size(genes,1)
    if nchrs>1
        chrsplits = sort(rand(1:ngenes,nchrs-1))
        chromosomes = Chromosome[]
        for chr in 1:nchrs
            if chr==1 # first chromosome
                push!(chromosomes, Chromosome(genes[1:chrsplits[chr]], true))
            elseif chr==nchrs # last chromosome
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):end], true))
            else
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):chrsplits[chr]], true))
            end
        end
    else # only one chromosome
        chromosomes = [Chromosome(genes, true)]
    end
    secondset = deepcopy(chromosomes)
    for chrm in secondset
        chrm.maternal = !chrm.maternal
    end
    append!(chromosomes,secondset)
    chromosomes
end

function genesis(settings::Dict{String,Any},
                 nspecs::Int64=100000000, popsize::Int64 = 0, # about 600-900 species per cell
                 traitnames::Array{String,1} = ["dispmean",
                                                "dispshape",
                                                "mutprob",
                                                "precopt",
                                                "prectol",
                                                "repradius",
                                                "repsize",
                                                "reptol",
                                                "seedsize",
                                                "tempopt",
                                                "temptol"]) # minimal required traitnames
    community = Individual[]
    totalmass = 0.0
    while true 
        lineage = randstring(4)
        meangenes = length(traitnames)
        ngenes = rand(Poisson(meangenes))
        ngenes < 1 && (ngenes = 1)
        if settings["linkage"] == "none"
            nchrms = ngenes
        elseif settings["linkage"] == "full"
            nchrms = 1
        else
            nchrms = rand(1:ngenes)
        end
        traits = createtraits(traitnames, settings)
        genes = creategenes(ngenes,traits)
        chromosomes = createchrs(nchrms,genes)
        traitdict = chrms2traits(chromosomes)
        popsize = round(fertility * traitdict["repsize"]^(-1/4) * exp(-act/(boltz*traitdict["tempopt"]))) # population size determined by adult size and temperature niche optimum
        popmass = popsize * traitdict["seedsize"]
        if totalmass + popmass > nspecs # stop loop if cell is full
            break
        end
        totalmass += popmass
        newind = Individual(lineage, chromosomes, traitdict, 0, false, 1.0, traitdict["seedsize"])
        for i in 1:popsize
            push!(community, newind)
        end
    end
    community
end

function readmapfile(filename::String)
    println("Reading file \"$filename\"...")
    mapstrings = String[]
    open(filename) do file
        mapstrings = readlines(file)
    end
    mapstrings = filter(x->!all(isspace,x),mapstrings) # remove empty lines
    mapstrings = filter(x->x[1]!='#',mapstrings) # remove comment lines
    mapsubstrings = map(split,mapstrings)
    mapentries = map(x->map(String,x),mapsubstrings)
    timesteps = 0
    try
        timesteps = parse(Int,filter(x->size(x,1)==1,mapentries)[1][1])
    catch
        timesteps = 1000
        warn("your mapfile \"$filename\" does not include timestep information. Assuming $timesteps timesteps.")
    end
    mapentries = filter(x->size(x,1)>1,mapentries)
    return timesteps,mapentries
end

function createworld(maptable::Array{Array{String,1},1}, settings::Dict{String,Any})
    println("Creating world...")
    world = Patch[]
    length(maptable) > 1 ? area = cellsize : area = 1
    for entry in maptable
        size(entry,1) < 3 && error("please check your map file for incomplete or faulty entries. \n
Each line must contain patch information with at least \n
\t - a unique integer ID, \n
\t - an integer x coordinate, \n
\t - an integer y coordinate, \n
separated by a whitespace character (<ID> <x> <y>).")
        id = parse(Int64, entry[1])
        xcord = parse(Int64, entry[2])
        ycord = parse(Int64, entry[3])
        size(entry,1) > 3 ? temperature = parse(Float64, entry[4]) : temperature = 298.0
        isisland = false
        if size(entry,1) > 4 && contains(lowercase(entry[5]),"island")
            isisland = true # islands do not receive an initial community
        end
        newpatch = Patch(id,(xcord,ycord),temperature,area,isisland)
        if size(entry,1) > 5 && contains(lowercase(entry[6]),"isolated")
            newpatch.isolated = true
        end
        if size(entry,1) > 6 # || settings["heterogeneity"] == "high"
            try
                newpatch.nichea = parse(Float64, entry[7])
            catch
                newpatch.nichea = rand() * 10
            end
        end
        if size(entry,1) > 7
            newpatch.nicheb = parse(Float64, entry[8])
        end
        !isisland && append!(newpatch.community,genesis(settings))
        push!(world,newpatch)
    end
    world
end

function updateworld!(world::Array{Patch,1},maptable::Array{Array{String,1},1}) #TODO: add functionality to remove patches!
    println("Updating world...")
    area = 1 # CAVE: just for now...
    for entry in maptable
        size(entry,1) < 3 && error("please check your map file for incomplete or faulty entries. \n
                            Each line must contain patch information with at least \n
                            \t - a unique integer ID, \n
                            \t - an integer x coordinate, \n
                            \t - an integer y coordinate, \n
                            separated by a whitespace character (<ID> <x> <y>).")
        id = parse(Int64, entry[1])
        xcord = parse(Int64, entry[2])
        ycord = parse(Int64, entry[3])
        size(entry,1) > 3 ? temperature = parse(Float64, entry[4]) : temperature = 298.0
        isisland = false
        if size(entry,1) > 4 && contains(lowercase(entry[5]),"island")
            isisland = true # islands do not receive an initial community
        end
        isolated = false
        if size(entry,1) > 5 && contains(lowercase(entry[6]),"isolated")
            isolated = true
        end
        try
            p = find(x->x.id==id,world)
            world[p].location = (xcord,ycord)
            world[p].altitude = temperature
            world[p].isisland = isisland
            world[p].isolated = isolated
        catch
            newpatch = Patch(id,(xcord,ycord),temperature,area,isisland)
            newpatch.isolated = isolated
            push!(world,newpatch)
        end
    end
    world
end

"""
    dumpinds(world, io, sep)
Output all data of individuals in `world` as table to `io`. Columns are separated by `sep`.
"""
function dumpinds(world::Array{Patch, 1}, io::IO = STDOUT, sep::String = "\t", onlyisland::Bool = false)
    header = true
    traitkeys = []
    counter = 0
    for patch in world
        (onlyisland && !patch.isisland) && continue # only one individual per species on mainland
        lineage = ""
        for ind in patch.community
            counter += 1
            (!patch.isisland && ind.lineage == lineage) && continue
            if header
                print(io, "id", sep)
                print(io, "xloc", sep)
                print(io, "yloc", sep)
                print(io, "temp", sep)
                print(io, "area", sep)
                print(io, "nichea", sep)
                ## print(io, "nicheb", sep)
                print(io, "island", sep)
                print(io, "isolation", sep)
                print(io, "counter", sep)
                print(io, "lineage", sep)
                print(io, "age", sep)
                print(io, "new", sep)
                print(io, "fitness", sep)
                print(io, "size", sep)
                print(io, "lnkgunits", sep)
                print(io, "ngenes", sep)
                traitkeys = keys(ind.traits)
                for key in traitkeys
                    print(io, key, sep)
                end
                println(io)
                header = false
            end
            print(io, patch.id, sep)
            print(io, patch.location[1], sep)
            print(io, patch.location[2], sep)
            print(io, patch.altitude, sep)
            print(io, patch.area, sep)
            print(io, patch.nichea, sep)
            ## print(io, patch.nicheb, sep)
            patch.isisland ? print(io, 1, sep) : print(io, 0, sep)
            patch.isolated ? print(io, 1, sep) : print(io, 0, sep)
            print(io, counter, sep)
            print(io, ind.lineage, sep)
            print(io, ind.age, sep)
            ind.isnew ? print(io, 1, sep) : print(io, 0, sep)
            print(io, ind.fitness, sep)
            print(io, ind.size, sep)
            print(io, length(ind.genome), sep)
            print(io, sum(map(x -> length(x.genes), ind.genome)), sep)
            for key in traitkeys
                try
                    print(io, ind.traits[key], sep)
                catch
                    print(io, "NA", sep)
                end
            end
            println(io)
            lineage = ind.lineage
        end
    end
end

function makefasta(world::Array{Patch, 1}, io::IO = STDOUT, sep::String = "", onlyisland::Bool = false)
    counter = 0
    for patch in world
        (onlyisland && !patch.isisland) && continue
        lineage = ""
        for ind in patch.community
            counter += 1
            (!patch.isisland && ind.lineage == lineage) && continue # only one individual per species on mainland
            chrmno = 0
            for chrm in ind.genome
                chrmno += 1
                geneno = 0
                for gene in chrm.genes
                    geneno += 1
                    traits = ""
                    if length(gene.codes) == 0
                        traits *= "neutral"
                    else
                        for trait in gene.codes
                            traits *= trait.name * "$(trait.value)" * ","
                        end
                    end
                    header = ">$counter x$(patch.location[1]) y$(patch.location[2]) $(ind.lineage) c$chrmno g$geneno $traits"
                    println(io, header)
                    println(io, gene.sequence)
                end
            end
            lineage = ind.lineage
        end
    end
end

"""
    setupdatadir(dir)
creates the output directory. `dir` is a string containing the directory path.
"""
function setupdatadir(settings::Dict{String, Any})
    if !isdir(settings["dest"])
        try
            mkdir(settings["dest"])
        catch
            warn("could not create directory ", settings["dest"], ". Assuming it already exists...")
        end
    end
end

#TODO: complete docstring!
"""
    writedata(world, path, settings, seed, timestep)
writes simulation output from `world` to separate table and fasta files.
`path`, `seed`, `timestep` and `setting` information is used for file name creation.
"""
function writedata(world::Array{Patch,1}, mappath::String, settings::Dict{String, Any}, seed::Int64, timestep::Int64)
    mapfile = split(mappath, "/")[end]
    length(mapfile) == 0 && return
    basename = "$(settings["dest"])" * "/" * mapfile * "_s" * "$seed" * "_lnk" * settings["linkage"] * "_tol" * settings["tolerance"] * "_t" * "$timestep"
    counter = 0
    extension = ""
    while ispath(basename * extension * ".tsv") || ispath(basename * extension * ".fa")
        extension = "_$counter"
        counter += 1
        if counter > 9
            warn("could not write to ", basename, extension, ": file exists. \n
Continuing anyway - data might be identical.")
            return
        end
    end
    basename *= extension
    filename = basename * ".tsv"
    println("Writing data \"$filename\"")
    open(filename, "w") do file
        dumpinds(world, file, "\t", timestep > 1)
    end
    filename = basename * ".fa"
    println("Writing fasta \"$filename\"")
    open(filename, "w") do file
        makefasta(world, file, "", timestep > 1)
    end
end

#TODO: complete docstring!
"""
    writerawdata(world, path, settings, seed, t)
writes raw julia data of the complete simulation state from `world` to a file in `path`.
`seed`, `setting` and `t` (timestep) information is used for file name creation.
"""
function writerawdata(world::Array{Patch,1}, mappath::String, settings::Dict{String, Any}, seed::Int64, timestep::Int64)
    mapfile = split(mappath, "/")[end]
    filename = "$(settings["dest"])" * "/" * mapfile * "_s" * "$seed" * "_lnk" * settings["linkage"] * "_tol" * settings["tolerance"] * "_t" * "$timestep" * ".jl"
    counter = 0
    extension = ""
    while ispath(filename * extension)
        extension = "_$counter"
        counter += 1
        counter > 9 && error("file \"$filename$extension\" exists. Please clear your directory.")
    end
    filename *= extension
    touch(filename)
    println("Writing raw data to \"$filename\"...")
    if timestep == 1
        open(filename, "w") do file
            println(file, world)
        end
    else
        island = filter(x -> x.isisland, world)
        open(filename, "w") do file
            println(file, island)
        end
    end
end

#TODO: complete docstring!
"""
    recordcolonizers(colos, path, settings, seed, t)
writes raw julia data of the colonizing individuals `colos` at timestep `t` to a file in `path`.
`seed` and `setting` information is used for file name creation.
"""
function recordcolonizers(colonizers::Array{Individual, 1}, mappath::String, settings::Dict{String, Any}, seed::Int64, timestep::Int64)
    mapfile = split(mappath, "/")[end]
    record = (timestep, colonizers)
    filename = "$(settings["dest"])" * "/" * mapfile * "_s" * "$seed" * "_lnk" * settings["linkage"] * "_tol" * settings["tolerance"] * "_colonizers" * ".jl"
    touch(filename)
    println("Colonisation! Writing data to ", filename, "...")
    open(filename, "a") do file
        println(file, record)
    end
end


end
