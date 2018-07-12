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

function simulation!(world::Array{Patch,1}, mapfile::String, timesteps::Int=1000)
    simlog("Starting simulation.")
    checkviability!(world, settings["static"])
    for t in 1:timesteps
        simlog("UPDATE $t")
        # ecological processes
        establish!(world, settings["nniches"], settings["static"])
        compete!(world, settings["static"])
        survive!(world, settings["static"])
        disturb!(world, settings["disturbance"], settings["static"])
        grow!(world, settings["static"])
        compete!(world, settings["static"])
        reproduce!(world)
        settings["mutate"] && mutate!(world)
        checkviability!(world, settings["static"])
        settings["burn-in"] < t && invade!(world)
        colonizers = disperse!(world, settings["static"])
        # model output
        map(p -> simlog("Patch $(p.id): $(length(p.community)) individuals.", 'd'), world)
        length(colonizers) >= 1 && simlog("t=$t: colonization by $colonizers", 'd')#recordcolonizers(colonizers, settings, t)
        if settings["lineages"]
            recordstatistics(world)
            recordlineages(world, t)
        end
        (t == 1 || mod(t, settings["outfreq"]) == 0) && writedata(world, mapfile, t)
    end
end

function runit()
    srand(settings["seed"])
    setupdatadir(settings)
    world = Patch[]
    for i in 1:length(settings["maps"])
        timesteps,maptable = readmapfile(settings["maps"][i])
        i == 1 && (world = createworld(maptable))
        i > 1 && updateworld!(world,maptable,settings["cellsize"])
        simulation!(world, settings["maps"][i], timesteps)
        writedata(world, settings["maps"][i], -1)
    end
end

@time runit()
