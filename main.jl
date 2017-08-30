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

function createworld()
end

function analysis()
end

function runit()
    world = createworld()
    simulation(world)
    analysis()
end

