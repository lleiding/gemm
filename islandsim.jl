#!/usr/bin/env julia

# Island speciation model using GeMM
#
# Ludwig Leidinger 2018
# <l.leidinger@gmx.net>
#
# Daniel Vedder 2018
# <daniel.vedder@stud-mail.uni-wuerzburg.de>
#
# For a list of options, run `julia islandsim.jl --help`
#
# <MAPFILE> is a textfile containing information about the simulation arena
# Every line describes one patch in the following format:
# <ID> <X-COORDINATE> <Y-COORDINATE> [<TYPE>]

thisDir = pwd() * "/src"
any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)

using GeMM

function simulation!(world::Array{Patch,1}, settings::Dict{String,Any}, mapfile::String, seed::Int64, timesteps::Int=1000)
    info("Starting simulation...")
    for t in 1:timesteps
        info("UPDATE $t, population size $(sum(x -> length(x.community), world))")
        (t <= 10 || mod(t, 20) == 0) && writedata(world, mapfile, settings, seed, t)
        establish!(world, settings["nniches"], settings["static"])
        checkviability!(world, settings["static"])
        compete!(world, settings["static"])
        survive!(world, settings["static"])
        disturb!(world, settings["disturbance"], settings["static"])
        grow!(world, settings["static"])
        compete!(world, settings["static"])
        reproduce!(world, settings["static"])
        settings["mutate"] && mutate!(world, settings)
        #TODO invaders = invade!(world, settings["propagule-pressure"])
        #TODO length(colonizers) >= 1 && println("t=$t: colonization by $colonizers")#recordcolonizers(colonizers, mapfile, settings, seed, t)
        colonizers = disperse!(world, settings["static"])
        length(colonizers) >= 1 && println("t=$t: colonization by $colonizers")#recordcolonizers(colonizers, mapfile, settings, seed, t)
    end
end

function runit(settings::Dict{String,Any},seed::Int64=0)
    seed != 0 && srand(seed)
    settings["maps"] = map(x->String(x),split(settings["maps"],","))
    settings["cellsize"] *= 1e6 #convert tonnes to grams
    setupdatadir(settings)
    for i in 1:length(settings["maps"])
        timesteps,maptable = readmapfile(settings["maps"][i])
        i == 1 && (world = createworld(maptable, settings))
        i > 1 && updateworld!(world,maptable,settings["cellsize"])
        simulation!(world, settings, settings["maps"][i], seed, timesteps)
        println("Memory usage:\n\nWorld - $(Base.summarysize(world))")
        println("Individual - $(mean(map(x -> Base.summarysize(x), world[1].community)))")
        writedata(world, settings["maps"][i], settings, seed, -1)
    end
end

## Settings
const allargs = parsecommandline()
if haskey(allargs, "config") && allargs["config"] != nothing
    allargs = parseconfig(allargs["config"], allargs)
end

## Parallel stuff:
const nprocesses = nworkers()
const startseed = allargs["seed"]
const replicates = startseed:startseed+nprocesses-1


@time runit(allargs, startseed)
