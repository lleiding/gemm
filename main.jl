#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## Main simulation

include("MIbGxMCmod.jl")


using MIbGxMCmod, Plots


function readmapfile(filename::String)
    mapstrings = String[]
    open(filename) do file
        mapstrings = readlines(file)
    end
    mapentries = map(split,mapstrings)
    mapentries = map(x->map(String,x),mapentries)
    filter(x->size(x,1)>0,mapentries)
end

function createworld(maptable::Array{Array{String,1},1})
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
        size(entry,1) > 3 ? temperature = parse(Float64, entry[4]) : temperature = 298
        isisland = false
        if size(entry,1) > 4
            contains(lowercase(entry[5]),"island") && (isisland = true) # islands do not receive an initial community
        end
        newpatch = Patch(id,(xcord,ycord),temperature,area,isisland)
        !isisland && append!(newpatch.community,genesis())
        push!(world,newpatch)
    end
    world
end

function visualisation(world::Array{Patch,1},firstplot::Bool)
    xcords = map(x->Int(floor(x.location[1])),world)
    ycords = map(x->Int(floor(x.location[2])),world)
    values = map(x->size(x.community,1),world)
    xmin = minimum(xcords)
    xmax = maximum(xcords)
    xshift = 1 - xmin
    ymin = minimum(ycords)
    ymax = maximum(ycords)
    yshift = 1 - ymin
    mat = zeros(length(ymin:ymax),length(xmin:xmax))
    mat[:] = -10
    for i in eachindex(values)
        mat[ycords[i]+yshift,xcords[i]+xshift]=values[i]
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

function simulation(world::Array{Patch,1}, timesteps::Int=1000)
    for t in 1:timesteps
        checkviability!(world)
        establish!(world)
        age!(world)
        grow!(world)
        disperse!(world)
        compete!(world)
        reproduce!(world)
        #        mod(timesteps,20) == 0 && analysis(world)
        #        mod(timesteps,20) == 0 && visualisation(world,t==20)
        mod(timesteps,10) == 0 && println(t)
    end
end

function runit(firstrun::Bool)
    length(ARGS) > 0 ? mapfile = ARGS[1] : mapfile = "mapfile"
    length(ARGS) > 1 ? seed = parse(Int,ARGS[2]) : seed = 0
    length(ARGS) > 2 ? timesteps = parse(Int,ARGS[3]) : timesteps = 1000
    seed != 0 && srand(seed)
    if firstrun
        world=createworld([["1","1","1"]])
        simulation(world, 10)
    else
        maptable=readmapfile(mapfile)
        world=createworld(maptable)
        simulation(world, timesteps)
        analysis(world)
    end
end

@time runit(true)
@time runit(false)
@time runit(false)
