#!/usr/bin/env julia

## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## Script to extract data in tabular format from raw julia code
##
## For a list of options, run julia getdata.jl --help

thisDir = pwd()
any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)

using MIbGxMCmod, ArgParse

function parsecommandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--dest", "-d"
            help = "output directory. Defaults to working directory"
            arg_type = String
            required = false
            default = "." #string(Dates.today())
        "input"
            help = "name of input file with raw julia data"
            required = true
        end
    return parse_args(s)
end

function extract(settings::Dict{String, Any})
    world = include(settings["input"])
    infile = split(settings["input"], "/")[end]
    open(settings["dest"]*"/"*infile*".out", "w") do file
        MIbGxMCmod.dumpinds(world, file)
    end
    # XXX FASTA not needed for invasion experiments
    # open(settings["dest"]*"/"*infile*".fa", "w") do file
    #     MIbGxMCmod.makefasta(world, file)
    # end
end

allargs = parsecommandline()
extract(allargs)
