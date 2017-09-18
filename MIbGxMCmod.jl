#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## contains all custom types and functions necessary for eco-evo-env speciation island model


module MIbGxMCmod


using Distributions, Plots


export Patch, # types
    readmapfile, writerawdata, writedata, analysis,
    checkviability!, establish!, age!, grow!, disperse!, compete!, reproduce!,
    createworld, updateworld! # functions


const boltz = 1.38064852e-23 #  J/K = m2⋅kg/(s2⋅K)
const act = 1e-19 # activation energy /J, ca. 0.63eV - Brown et al. 2004
const normconst = 1e10 # normalization constant to get biologically realistic orders of magnitude
const mutscaling = 50 #parse(ARGS[2])
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
        warn("your mapfile \"$filename\" does not include timestep information. Assumed $timesteps timesteps.")
    end
    mapentries = filter(x->size(x,1)>1,mapentries)
    return timesteps,mapentries
end

function createworld(maptable::Array{Array{String,1},1})
    println("Creating world...")
    world = Patch[]
    area = 100 # CAVE: just for now...
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
        !isisland && append!(newpatch.community,genesis())
        push!(world,newpatch)
    end
    world
end

function updateworld!(world::Array{Patch,1},maptable::Array{Array{String,1},1}) #TODO: add functionality to remove patches!
    println("Updating world...")
    area = 100 # CAVE: just for now...
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

function visualisation(world::Array{Patch,1},firstplot::Bool)
    xcords = map(x->Int(floor(x.location[1])),world)
    ycords = map(x->Int(floor(x.location[2])),world)
    popsizes = map(x->size(x.community,1),world)
    xmin = minimum(xcords)
    xmax = maximum(xcords)
    xshift = 1 - xmin
    ymin = minimum(ycords)
    ymax = maximum(ycords)
    yshift = 1 - ymin
    mat = zeros(length(ymin:ymax),length(xmin:xmax))
    mat[:] = -10
    for i in eachindex(popsizes)
        mat[ycords[i]+yshift,xcords[i]+xshift]=popsizes[i]
    end
    if firstplot
        gr()
        heatmap(mat,aspect_ratio=1,show=true)
    else
        heatmap(mat,aspect_ratio=1,show=true)
    end
end

function analysis(world::Array{Patch,1})
    for patch in world
        print("Patch #",patch.id," @",patch.location," (",patch.isisland,"): \t")
        println(size(patch.community,1))
    end
end

function dumpinds(world::Array{Patch,1},io::IO=STDOUT,sep::String="\t")
    header = true
    traitkeys = []
    for patch in world
        for ind in patch.community
            if header
                print(io, "id", sep)
                print(io, "xloc", sep)
                print(io, "yloc", sep)
                print(io, "temp", sep)
                print(io, "area", sep)
                ## print(io, "nichea", sep)
                ## print(io, "nicheb", sep)
                print(io, "landtype", sep)
                print(io, "isolation", sep)
                print(io, "age", sep)
                print(io, "new", sep)
                print(io, "fitness", sep)
                print(io, "size", sep)
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
            ## print(io, patch.nichea, sep)
            ## print(io, patch.nicheb, sep)
            patch.isisland ? print(io, "island", sep) : print(io, "continent", sep)
            patch.isolated ? print(io, "yes", sep) : print(io, "no", sep)
            print(io, ind.age, sep)
            ind.isnew ? print(io, "yes", sep) : print(io, "no", sep)
            print(io, ind.fitness, sep)
            print(io, ind.size, sep)
            for key in traitkeys
                print(io, ind.traits[key], sep)
            end
            println(io)
        end
    end
end

function writedata(world::Array{Patch,1}, seed::Int64, mapfile::String)
    filename = mapfile * "_seed" * "$seed" * ".out"
    counter = 0
    extension = ""
    while ispath(filename * extension)
        extension = "_$counter"
        counter += 1
        counter > 9 && error("file \"$filename$extension\" exists. Please clear your directory.")
    end
    filename *= extension
    touch(filename)
    println("Writing data to \"$filename\"...")
    open(filename, "w") do file
        dumpinds(world,file)
    end
end

function writerawdata(world::Array{Patch,1}, seed::Int64, mapfile::String, timestep::Int64)
    filename = mapfile * "_seed" * "$seed" * "_t" * "$timestep" * ".jl"
    counter = 0
    extension = ""
    while ispath(filename * extension)
        extension = "_$counter"
        counter += 1
        counter > 9 && error("file \"$filename$extension\" exists. Please clear your directory.")
    end
    filename *= extension
    touch(filename)
    println("Writing data to \"$filename\"...")
    open(filename, "w") do file
        println(file,world)
    end
end

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
    for traitname in unique(map(x->x.name,traits))
        traitdict[traitname] = mean(map(x->x.value,filter(x->x.name==traitname,traits)))
    end
    traitdict
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

function checkviability!(patch::Patch) # may consider additional rules... # maybe obsolete anyhow...
    idx=1
    while idx <= size(patch.community,1)
        dead = false
        patch.community[idx].size <= 0 && (dead = true)
        0 in collect(values(patch.community[idx].traits)) && (dead = true)
        patch.community[idx].traits["repsize"] <= patch.community[idx].traits["seedsize"] && (dead = true)
        if dead
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

function gausscurve(b::Float64, c::Float64, x::Float64, a::Float64=1.0)
    y = a * exp(-(x-b)^2/(2*c^2))
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
        elseif patch.community[idx].isnew || patch.community[idx].age == 0
            tempopt = patch.community[idx].traits["tempopt"]
            temptol = patch.community[idx].traits["temptol"]
            fitness = gausscurve(tempopt, temptol, temp)
            fitness > 1 && (fitness = 1) # should be obsolete
            fitness < 0 && (fitness = 0) # should be obsolete
            patch.community[idx].fitness = fitness
            patch.community[idx].isnew = false
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
                newmass = growthrate * patch.community[idx].fitness * mass^(3/4) * exp(-act/(boltz*temp)) * normconst # -> emergent maximum body size!
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

"""
    findisland(x)
find out in which direction from the continent the island(s) lie(s).
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

"""
function checkborderconditions(world::Array{Patch,1},xdest::Float64,ydest::Float64)
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

function disperse!(world::Array{Patch,1}) # TODO: additional border conditions
    colonisation = false
     for patch in world
        idx = 1
        while idx <= size(patch.community,1)
            hasdispmean = haskey(patch.community[idx].traits,"dispmean")
            hasdispprob = haskey(patch.community[idx].traits,"dispprob")
            hasdispshape = haskey(patch.community[idx].traits,"dispshape")
            if !hasdispmean || !hasdispprob || !hasdispshape
                splice!(patch.community,idx)
                idx -= 1
            elseif patch.community[idx].isnew && patch.community[idx].age == 0 && rand() <= patch.community[idx].traits["dispprob"] * patch.community[idx].fitness
                dispmean = patch.community[idx].traits["dispmean"]
                dispshape = patch.community[idx].traits["dispshape"]
                patch.community[idx].isnew = true
                indleft = splice!(patch.community,idx)
                xdir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # scaling so that geometric mean...
                ydir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # ...follows original distribution
                xdest = patch.location[1]+xdir
                ydest = patch.location[2]+ydir
                !patch.isisland && ((xdest,ydest) = checkborderconditions(world,xdest,ydest))
                targets = unique([(floor(xdest),floor(ydest)),(ceil(xdest),floor(ydest)),(ceil(xdest),ceil(ydest)),(floor(xdest),ceil(ydest))])
                possdests = find(x->in(x.location,targets),world)
                if size(possdests,1) > 0 # if no viable target patch, individual dies
                    destination = rand(possdests)
                    originisolated = patch.isolated && rand() <= indleft.traits["dispprob"]
                    targetisolated = world[destination].isolated && rand() <= indleft.traits["dispprob"]
                    (!originisolated && !targetisolated) && push!(world[destination].community,indleft)
                    !patch.isisland && world[destination].isisland && (colonisation = true)
                end
                idx -= 1
            end
            idx += 1
        end
    end
    colonisation
end

function compete!(patch::Patch)
    sort!(patch.community, by = x -> x.size)
    while sum(map(x->x.size,patch.community)) > patch.area # occupied area larger than available
        victim = rand(Geometric()) + 1
        victim > length(patch.community) && (victim = length(patch.community))
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
        hasreprate = haskey(patch.community[idx].traits,"reprate")
        hasrepsize = haskey(patch.community[idx].traits,"repsize")
        hasreptol = haskey(patch.community[idx].traits,"reptol")
        hasseedsize = haskey(patch.community[idx].traits,"seedsize")
        hasmutprob = haskey(patch.community[idx].traits,"mutprob")
        if !hasrepprob || !hasreprate || !hasrepsize || !hasreptol || !hasmutprob || !hasseedsize
            splice!(patch.community, idx)
            idx -= 1
        elseif !patch.community[idx].isnew && rand() <= patch.community[idx].traits["repprob"] * patch.community[idx].fitness
            currentmass = patch.community[idx].size
            seedsize = patch.community[idx].traits["seedsize"]
            if currentmass >= patch.community[idx].traits["repsize"]
                meanoffs = patch.community[idx].traits["reprate"]
                reptol = patch.community[idx].traits["reptol"]
                mass = patch.community[idx].size
                metaboffs =  meanoffs * currentmass^(-1/4) * exp(-act/(boltz*temp)) * normconst
                noffs = rand(Poisson(metaboffs))
                if sexualreproduction
                    posspartners = find(x->(1/reptol)*mass>=x.size>=reptol*mass,patch.community)
                    partneridx = idx
                    try
                        partneridx = rand(filter(x->x!=idx,posspartners))
                    catch
                        idx += 1
                        continue
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
    seedsize = rand()
    repsize = rand()
    while repsize <= seedsize
        seedsize = rand()
        repsize = rand()
    end
    for name in traitnames
        if contains(name,"rate")
            push!(traits,Trait(name,rand()*100,true))
        elseif contains(name, "temp") && contains(name, "opt")
            push!(traits,Trait(name,rand()*60+263,true)) #CAVE: code values elsewhere?
        elseif contains(name, "tol") && !contains(name, "rep")
            push!(traits,Trait(name,rand()*20,true)) #CAVE: code values elsewhere?
        elseif contains(name, "mut")
            push!(traits,Trait(name,rand(),true)) #CAVE: code values elsewhere?
        elseif contains(name, "repsize")
            push!(traits,Trait(name,repsize,true)) #CAVE: code values elsewhere?
        elseif contains(name, "seedsize")
            push!(traits,Trait(name,seedsize,true)) #CAVE: code values elsewhere?
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
            push!(gene.codes,trait)
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
                                                "repsize",
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
