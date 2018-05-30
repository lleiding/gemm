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

thisDir = joinpath(pwd(), "src")
any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)

using GeMM

function simulation!(world::Array{Patch,1}, settings::Dict{String,Any}, mapfile::String, timesteps::Int=1000)
    info("Starting simulation")
    for t in 1:timesteps
        info("UPDATE $t, population size $(sum(x -> length(x.community), world))")
        (t == 1 || mod(t, settings["outfreq"]) == 0) && writedata(world, settings, mapfile, t)
        establish!(world, settings["nniches"], settings["static"])
        checkviability!(world, settings["static"])
        compete!(world, settings["static"])
        survive!(world, settings["static"])
        disturb!(world, settings["disturbance"], settings["static"])
        grow!(world, settings["static"])
        compete!(world, settings["static"])
        reproduce!(world, settings)
        settings["mutate"] && mutate!(world, settings)
        invade!(world, settings)
        colonizers = disperse!(world, settings["static"])
        length(colonizers) >= 1 && println("t=$t: colonization by $colonizers")#recordcolonizers(colonizers, settings, t)
    end
end

function runit(settings::Dict{String,Any})
    srand(settings["seed"])
    setupdatadir(settings)
    world = Patch[]
    for i in 1:length(settings["maps"])
        timesteps,maptable = readmapfile(settings["maps"][i])
        i == 1 && (world = createworld(maptable, settings))
        i > 1 && updateworld!(world,maptable,settings["cellsize"])
        simulation!(world, settings, settings["maps"][i], timesteps)
        writedata(world, settings, settings["maps"][i], -1)
        println("WORLD POPULATION: $(sum(x -> length(x.community), world))") #DEBUG
    end
end

## Settings
const settings = getsettings()

@time runit(settings)
