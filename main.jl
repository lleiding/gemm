#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## Main simulation
##
## Synopsis: julia main.jl <MAPFILE> <RANDOMSEED> <TIMESTEPS>
## <MAPFILE> is a textfile containing information about the simulation arena
## Every line describes one patch in the following format:
## <ID> <X-COORDINATE> <Y-COORDINATE> [<TYPE>]

addprocs(0)
thisDir = pwd()
any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)

#@everywhere include("MIbGxMCmod.jl")

# using MIbGxMCmod # needs to be like this somehow...
@everywhere using MIbGxMCmod


@everywhere function simulation(world::Array{Patch,1}, seed::Int64, mapfile::String, timesteps::Int=1000)
    println("Starting simulation...")
    for t in 1:timesteps
        checkviability!(world)
        establish!(world)
        age!(world)
        grow!(world)
        compete!(world)
        reproduce!(world)
        disperse!(world) && writerawdata(world,seed,mapfile,t)
#        (t == 1 || mod(t,100) == 0) && analysis(world)
        #        mod(timesteps,20) == 0 && visualisation(world,t==20)
        #mod(timesteps,10) == 0 && println(t)
    end
end

@everywhere function runit(firstrun::Bool,mapfiles::Array{String,1},seed::Int64=0)
    seed != 0 && srand(seed)
    if firstrun
        world=createworld([["1","1","1"]])
        simulation(world, seed, "", 10)
    else
        for i in 1:length(mapfiles)
            timesteps,maptable = readmapfile(mapfiles[i])
            i == 1 && (world = createworld(maptable))
            i > 1 && updateworld!(world,maptable)
            simulation(world, seed, mapfiles[i], timesteps)
            writedata(world, seed, mapfiles[i])
        end
    end
end

## Parallel stuff:
const nprocesses = nworkers()
length(ARGS) > 0 ? (const startseed = parse(Int,ARGS[1])) : (const startseed = 1)
length(ARGS) > 1 ? (const mapfiles = ARGS[2:end]) : (const mapfiles = ["mapfile"])
const replicates = startseed:startseed+nprocesses-1
pmap(x->runit(true,mapfiles,x),replicates)
pmap(x->runit(false,mapfiles,x),replicates)
