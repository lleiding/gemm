#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## Main simulation

include("MIbGxMCmod.jl")


using MIbGxMCmod


const boltz = 1.38064852e-23 #  J/K = m2⋅kg/(s2⋅K)
const act = 1e-19 # activation energy /J, ca. 0.63eV - Brown et al. 2004
const normconst = 1e10 # normalization constant to get biologically realistic orders of magnitude
const mutscaling=50#parse(ARGS[2])
const sexualreproduction = true


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

function createworld(maptable::Array{String,1})
    world = Patch[]
    for entry in maptable
        size(strings,1) < 3 && error("please check your map file for incomplete or faulty entries!")
        id = parse(Int64, strings[1])
        xcord = parse(Int64, strings[2])
        ycord = parse(Int64, strings[3])
        size(strings,1) > 3 ? temperature = parse(Float64, strings[3]) : temperature = 298
        isisland = false
        if size(strings,1) > 3
            contains(lowercase(entry[4]),"island") && isisland = true # islands do not receive an initial community
        end
        newpatch = Patch(id,(xcord,ycord),temperature,area))
    end
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
    world = createworld()
    simulation(world)
    analysis()
end

examplemap = [

exampleworld = 
