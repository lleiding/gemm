#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## Data extractor for raw julia data
##

thisDir = pwd()
any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)

using MIbGxMCmod, MIbGxMCmod.

function extractinfo()
    inputfile = ARGS[1]
    world = include(inputfile)
    outfile = inputfile * ".tsv"
    if !ispath(outfile)
        open(outfile, "w") do f
            dumpinds(world, f)
        end
    else
        dumpinds(world)
    end
end

extractinfo()
