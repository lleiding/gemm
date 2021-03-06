#!/usr/bin/env julia
# A wrapper script to launch one or several GeMM simulation runs, parallelising
# them where possible. *Warning:* Because this script includes its own commandline
# parsing, no commandline args are passed on to GeMM!

using # TODO: automatic package installation - ask user?
    Distributed,
    ArgParse,
    Random

@everywhere thisDir = joinpath(pwd(), "src")
@everywhere any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using GeMM

"""
    parseexperiments()

Get configuration file(s), random seed and number of replicates from the commandline.
"""
function parseexperiments()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config", "-c"
            help = "list of configuration files, comma separated"
            arg_type = String
            required = true
        "--seed", "-s"
            help = "initial random seed"
            arg_type = Int
            required = false
            default = 0
        "--replicates", "-n"
            help = "number of replicates to run"
            arg_type = Int
            range_tester = x -> x > 0
            required = false
            default = 1
    end
    args = parse_args(s)
    for a in keys(args)
        (args[a] == nothing) && delete!(args, a)
    end
    args
end

args = parseexperiments()
configs = split(args["config"], ',')

simulations = []

for i in args["seed"]:(args["seed"] + args["replicates"] - 1)
    args["seed"] == 0 ? rep = abs(rand(RandomDevice(), Int32)) : rep = i
    for config in configs
        push!(simulations, [String(config), rep])
    end
end

# Only use `pmap` when we actually have multiple simulations.
# Otherwise, we can save on the overhead (plus, `pmap` obscures
# error messages).
if length(simulations) == 1
    sim = simulations[1]
    rungemm(sim[1], sim[2])
else
    pmap(x -> rungemm(x[1], x[2]), simulations)
end
