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

function simulation!(world::Array{Patch,1}, settings::Dict{String, Any}, mapfile::String, timesteps::Int=1000)
    simlog("Starting simulation.", settings)
    checkviability!(world, settings)
    for t in 1:timesteps
        simlog("UPDATE $t", settings)
        # ecological processes
        establish!(world, settings["nniches"], settings["static"])
        compete!(world, settings["static"])
        survive!(world, settings["static"])
        disturb!(world, settings["disturbance"], settings["static"])
        grow!(world, settings["static"])
        compete!(world, settings["static"])
        reproduce!(world, settings)
        settings["mutate"] && mutate!(world, settings)
        checkviability!(world, settings)
        settings["burn-in"] < t && invade!(world)
        disperse!(world, settings["static"])
        # model output
        map(p -> simlog("Patch $(p.id): $(length(p.community)) individuals.", settings, 'd'), world)
        if settings["lineages"]
            recordstatistics(world, settings)
            recordlineages(world, settings, t)
        end
        (t == 1 || mod(t, settings["outfreq"]) == 0) && writedata(world, settings, mapfile, t)
    end
end

function runit()
    settings = getsettings()
    srand(settings["seed"])
    setupdatadir(settings)
    world = Patch[]
    for i in 1:length(settings["maps"])
        timesteps,maptable = readmapfile(settings["maps"][i], settings)
        i == 1 && (world = createworld(maptable, settings))
        i > 1 && updateworld!(world,maptable,settings["cellsize"])
        simulation!(world, settings, settings["maps"][i], timesteps)
        writedata(world, settings, settings["maps"][i], -1)
    end
end

@time runit()
@time runit()
