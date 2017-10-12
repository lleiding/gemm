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

addprocs(0)
thisDir = pwd()
any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)

#@everywhere include("MIbGxMCmod.jl")

# using MIbGxMCmod # needs to be like this somehow...
@everywhere using MIbGxMCmod, ArgParse

@everywhere function parsecommandline()
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
        "--tolerance", "-t"
            help = "tolerance of sequence identity when reproducing (\"high\", \"evo\" or \"low\")"
            arg_type = String
            range_tester = x->in(x,["high", "evo", "low"])
            required = false
            default = "evo"
        "--dest", "-d"
            help = "output directory. Defaults to current date"
            arg_type = String
            required = false
            default = string(Dates.today())
        "--heterogeneity", "-n" # this ...
            help = "island heterogeneity, i.e. environmental niche variability (\"low\" or \"high\")"
            arg_type = String
            range_tester = x->in(x,["low", "high"])
            required = false
        "--complexity", "-c" # ... and that IS CURRENTLY defined in the maps
            help = "topographic island complexity, controlling intra-island isolation (\"low\" or \"high\")"
            arg_type = String
            range_tester = x->in(x,["low", "high"])
            required = false
#        "--flag1"
#            help = "an option without argument, i.e. a flag"
#            action = :store_true
#        "arg1"
#            help = "a positional argument"
#            required = false
        end
    return parse_args(s)
end


@everywhere function simulation(world::Array{Patch,1}, settings::Dict{String,Any}, mapfile::String, timesteps::Int=1000)
    println("Starting simulation...")
    for t in 1:timesteps
        checkviability!(world)
        establish!(world)
        compete!(world)
        grow!(world)
        age!(world)
        reproduce!(world)
        colonizers = disperse!(world)
        length(colonizers) > 1 && recordcolonizers(colonizers, mapfile, settings, t)
        mod(t, 1000) == 0 && writerawdata(world, mapfile, settings, t)
#        (t == 1 || mod(t,100) == 0) && analysis(world)
        #        mod(timesteps,20) == 0 && visualisation(world,t==20)
        #mod(timesteps,10) == 0 && println(t)
    end
end

@everywhere function runit(firstrun::Bool,settings::Dict{String,Any},seed::Int64=0)
    seed != 0 && srand(seed)
    mapfiles =  map(x->String(x),split(allargs["maps"],","))
    if firstrun
        world=createworld([["1","1","1"]],settings)
        simulation(world, settings, "", 10)
    else
        for i in 1:length(mapfiles)
            timesteps,maptable = readmapfile(mapfiles[i])
            i == 1 && (world = createworld(maptable, settings))
            i > 1 && updateworld!(world,maptable)
            simulation(world, settings, mapfiles[i], timesteps)
            writedata(world, seed, mapfiles[i])
        end
    end
end

## Parallel stuff:
const nprocesses = nworkers()
const allargs = parsecommandline()
const startseed = allargs["seed"]
const replicates = startseed:startseed+nprocesses-1

TT = STDOUT # save original STDOUT stream
redirect_stdout()
pmap(x->runit(true,allargs,x),replicates) # compilation/optimization run
setupdatadir(allargs)
redirect_stdout(TT) # restore STDOUT
pmap(x->runit(false,allargs,x),replicates)
