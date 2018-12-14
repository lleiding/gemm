# Input functions for GeMM

function getsettings()
    defaults = defaultSettings()
    commandline = parsecommandline()
    if isfile(commandline["config"])
        configs = parseconfig(commandline["config"])
    else
        configs = Dict{String, Any}()
    end
    settings = merge(defaults, configs, commandline)
    if settings["seed"] == 0
        settings["seed"] = abs(rand(Int32))
    end
    settings["maps"] = map(x -> String(x), split(settings["maps"], ","))
    if isa(settings["traitnames"], String)
        settings["traitnames"] = map(x -> String(x), split(settings["traitnames"], ","))
    end
    # Flags are automatically set to false by ArgParse if they are not given -
    # this should not override a given config file value
    for s in ["debug", "quiet", "static"]
        if !commandline[s] && s in keys(configs) && configs[s]
            settings[s] = true
        end
    end
    settings
end


function parsecommandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--seed", "-s"
            help = "inital random seed"
            arg_type = Int
        "--maps", "-m"
            help = "list of map files, comma separated"
            arg_type = String
            required = false
        "--config", "-c"
            help = "name of the config file"
            arg_type = String
            required = false
        "--linkage", "-l"
            help = "gene linkage (\"none\", \"random\" or \"full\")"
            arg_type = String
            range_tester = x->in(x,["none", "random", "full"])
            required = false
        "--nniches", "-n"
            help = "number of environmental niche traits (1 -- 3)"
            arg_type = Int
            range_tester = x -> x > 0 && x <= 3
            required = false
        "--tolerance", "-t"
            help = "tolerance of sequence identity when reproducing"
            arg_type = Float64
            range_tester = x -> 0.0 <= x <= 1.0
            required = false
        "--dest", "-d"
            help = "output directory. Defaults to current date"
            arg_type = String
            required = false
        "--debug"
            help = "debug mode. Turns on output of debug statements."
            action = :store_true
        "--quiet"
            help = "quiet mode. Don't print output to screen."
            action = :store_true
        "--static"
            help = "static mainland. Turns off any dynamics on the continent"
            action = :store_true
    end
    args = parse_args(s)
    for a in keys(args)
        (args[a] == nothing) && delete!(args, a)
    end
    args
end

"""
    basicparser(f)
Do elementary parsing on a config or map file.

Reads in the file, strips whole-line and inline comments
and separates lines by whitespace.
Returns a 2d array representing the tokens in each line.
"""
function basicparser(filename::String)
    # Read in the file
    lines = String[]
    open(filename) do file
        lines = readlines(file)
    end
    # Remove comments and tokenize
    lines = map(x -> strip(x), lines)
    filter!(x -> !isempty(x), lines)
    filter!(x -> (x[1] != '#'), lines)
    lines = map(s -> strip(split(s, '#')[1]), lines)
    lines = map(split, lines)
    map(l -> map(s -> convert(String, s), l), lines)
end

function parseconfig(configfilename::String)
    config = basicparser(configfilename)
    settings = Dict{String, Any}()
    defaults = defaultSettings()
    for c in config
        if length(c) != 2
            simlog("Bad config file syntax: $c", settings, 'w')
        elseif c[1] in keys(defaults)
            value = c[2]
            if !(typeof(defaults[c[1]]) <: AbstractString)
                try
                    value = parse(typeof(defaults[c[1]]), c[2]) # or Meta.parse with the old functionality
                catch
                    simlog("$(c[1]) not of type $(typeof(defaults[c[1]])).", settings, 'w')
                end
            end
            settings[c[1]] = value
        else
            simlog(c[1]*" is not a recognized parameter!", settings, 'w') # XXX maybe parse anyway
        end
    end
    settings
end

function readmapfile(mapfilename::String, settings::Dict{String, Any})
    simlog("Reading map file $mapfilename.", settings)
    mapfile = basicparser(mapfilename)
    timesteps = parse(Int, mapfile[1][1])
    if length(mapfile[1]) != 1 || !isa(timesteps, Integer)
        timesteps = 1000
        simlog("Invalid timestep information in the mapfile. Setting timesteps to 1000.", settings, 'w')
    end
    return timesteps,mapfile[2:end]
end

