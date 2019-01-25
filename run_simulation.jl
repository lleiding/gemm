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

# XXX This should be part of GeMM if we're to distribute the model as a package!

using
    GeMM,
    Random

function runsim(config::String = "", seed::Integer = 0, prerun::Bool = false)
    settings = getsettings()
    if prerun
        settings = defaultSettings()
        settings["cellsize"] = 1.0e6
        settings["indsize"] = "adult"
        settings["popsize"] = "metabolic"
        settings["maps"] = [""]
        settings["phylo"] = false
        settings["fasta"] = false
        settings["raw"] = false
        settings["quiet"] = true
        settings["static"] = false
        settings["seed"] = 1
    end
    Random.seed!(settings["seed"])
    !prerun && setupdatadir(settings)
    world = Patch[]
    for i in 1:length(settings["maps"])
        timesteps, maptable = readmapfile(settings["maps"][i], settings)
        i == 1 && (world = createworld(maptable, settings))
        i > 1 && updateworld!(world,maptable, settings)
        !prerun && writedata(world, settings, 0)
        simulate!(world, settings, timesteps)
    end
end


function rungemm(config::String = "", seed::Integer = 0)
    # compilation run:
    runsim("", 1, true)
    # run intended simulation:
    runsim(config, seed)
end

rungemm()
