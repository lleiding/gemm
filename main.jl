#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## Main simulation

include("MIbGxMCmod.jl")


using MIbGxMCmod


function simulation(world::Array{Patch,1}, timesteps::Int=1000)
    for t in 1:timesteps
        checkviability!(world)
        establish!(world)
        age!(world)
        grow!(world)
        disperse!(world)
        compete!(world)
        reproduce!(world)
    end
end

function createworld(maptable::Array{Array{SubString{String},1},1})
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

function readmapfile(filename::String)
    mapstrings = String[]
    open(filename) do file
        mapstrings = readlines(file)
    end
    mapentries = map(split,mapstrings)
    filter(x->size(x,1)>0,mapentries)
end

function analysis()
end

function runit()
    maptable=readmapfile("mapfile")
    world=createworld(maptable)
    simulation(world)
end

runit()
