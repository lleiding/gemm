#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## Main simulation
##
## For a list of options, run julia main.jl --help
##
## <MAPFILE> is a textfile containing information about the simulation arena
## Every line describes one patch in the following format:
## <ID> <X-COORDINATE> <Y-COORDINATE> [<TYPE>]

thisDir = pwd()
any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)

using GeMM, ArgParse

function parsecommandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--seed", "-s"
            help = "inital random seed"
            arg_type = Int
            default = 1
        "--maps", "-m"
            help = "list of map files, comma separated"
            arg_type = String
            required = true
        "--linkage", "-l"
            help = "gene linkage (\"none\", \"random\" or \"full\")"
            arg_type = String
            range_tester = x->in(x,["none", "random", "full"])
            required = false
            default = "random"
        "--nniches", "-n"
            help = "number of environmental niche traits (1 -- 3)"
            arg_type = Int
            range_tester = x -> x > 0 && x <= 3
            required = false
            default = 2
        "--tolerance", "-t"
            help = "tolerance of sequence identity when reproducing (\"high\", \"evo\", \"low\" or \"none\")"
            arg_type = String
            range_tester = x->in(x,["high", "evo", "low", "none"])
            required = false
            default = "evo"
        "--dest", "-d"
            help = "output directory. Defaults to current date"
            arg_type = String
            required = false
            default = string(Dates.today())
#        "--flag1"
#            help = "an option without argument, i.e. a flag"
#            action = :store_true
#        "arg1"
#            help = "a positional argument"
#            required = false
        end
    return parse_args(s)
end


function simulation!(world::Array{Patch,1}, settings::Dict{String,Any}, mapfile::String, seed::Int64, timesteps::Int=1000)
    println("Starting simulation...")
    for t in 1:timesteps
        (t == 1 || mod(t, 1000) == 0) && writedata(world, mapfile, settings, seed, t)
        establish!(world, settings["nniches"])
        checkviability!(world)
        compete!(world)
        survive!(world)
        grow!(world)
        compete!(world)
        reproduce!(world)
        mutate!(world, settings)
        colonizers = disperse!(world)
        length(colonizers) >= 1 && println("t=$t: colonization by $colonizers")#recordcolonizers(colonizers, mapfile, settings, seed, t)
    end
end

function runit(firstrun::Bool,settings::Dict{String,Any},seed::Int64=0)
    seed != 0 && srand(seed)
    mapfiles =  map(x->String(x),split(settings["maps"],","))
    if firstrun
        world=createworld([["1","1","1", "298", "continent"]], settings) # initialise as populated continent
        world[1].isisland = true # then set as island to enable processes
        simulation!(world, settings, "", seed, 5)
    else
        for i in 1:length(mapfiles)
            timesteps,maptable = readmapfile(mapfiles[i])
            i == 1 && (world = createworld(maptable, settings))
            i > 1 && updateworld!(world,maptable)
            simulation!(world, settings, mapfiles[i], seed, timesteps)
            writedata(world, mapfiles[i], settings, seed, 0)
        end
    end
end

## Parallel stuff:
const nprocesses = nworkers()
const allargs = parsecommandline()
const startseed = allargs["seed"]
const replicates = startseed:startseed+nprocesses-1

TT = STDOUT # save original STDOUT stream
setupdatadir(allargs)
redirect_stdout()
runit(true, allargs, startseed) # compilation/optimization run
redirect_stdout(TT) # restore STDOUT
@time runit(false, allargs, startseed)
