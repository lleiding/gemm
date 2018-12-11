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

using
    GeMM,
    Random

function runit(prerun::Bool = false)
    settings = getsettings()
    if prerun
        settings = defaultSettings()
        settings["cellsize"] = 1.0e6
        settings["indsize"] = "adult"
        settings["popsize"] = "single"
        settings["maps"] = [getsettings()["maps"][1]]
        settings["phylo"] = false
        settings["quiet"] = true
        settings["static"] = false
    end
    Random.seed!(settings["seed"])
    !prerun && setupdatadir(settings)
    world = Patch[]
    for i in 1:length(settings["maps"])
        timesteps, maptable = readmapfile(settings["maps"][i], settings)
        if prerun
            timesteps = 10
            maptable = [["1", "1", "1", "initpop"], ["2", "2", "1", "isisland"]]
        end
        i == 1 && (world = createworld(maptable, settings))
        i > 1 && updateworld!(world,maptable,settings["cellsize"])
        !prerun && writedata(world, settings, 0)
        simulate!(world, settings, timesteps)
    end
end

runit(true)
@time runit()
