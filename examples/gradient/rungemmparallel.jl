#/usr/bin/env julia

using
    Distributed,
    ArgParse

@everywhere thisDir = joinpath(pwd(), "src")
@everywhere any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)

@everywhere using GeMM

"""
    parseexperiments()

Get configuration file(s), random seed and number of replicates from the commandline.
"""
function parseexperiments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--config", "-c"
            help = "list of configuaration files, comma separated"
            arg_type = String
            required = true
        "--seed", "-s"
            help = "inital random seed"
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
    args["seed"] != 0 ? rep = abs(rand(RandomDevice(), Int32)) : rep = i
    for config in configs
        push!(simulations, [config, rep])
    end
end

pmap(x -> rungemm(x[1], x[2]), simulations)


