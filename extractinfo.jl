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

using MIbGxMCmod, ArgParse
