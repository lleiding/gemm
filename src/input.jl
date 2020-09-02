# Input functions for GeMM

"""
    getsettings()

Combines all configuration options to produce a single settings dict.
Order of precedence: commandline parameters - config file - default values
"""
function getsettings(configfile::String = "", seed::Integer = 0)
    defaults = defaultSettings()
    commandline = parsecommandline() # deprecate?
    if !isempty(configfile) && isfile(configfile)
        configs = parseconfig(configfile)
        commandline["config"] = configfile
    elseif haskey(commandline, "config") && isfile(commandline["config"])
        configs = parseconfig(commandline["config"])
    else
        configs = Dict{String, Any}()
    end
    seed != 0 && (commandline["seed"] = seed)
    settings = merge(defaults, configs, commandline)
    if settings["dest"] == defaults["dest"]
        settings["dest"] = settings["dest"] * "_" * join(split(basename(settings["config"]), ".")) * "_" * string(settings["seed"])
    end
    if settings["seed"] == 0
        settings["seed"] = abs(rand(RandomDevice(), Int32))
    end
    settings["maps"] = map(x -> String(x), split(settings["maps"], ","))
    map!(x -> joinpath(dirname(settings["config"]), x), settings["maps"], settings["maps"])
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

"""
    parsecommandline()

Certain parameters can be set via the commandline.
"""
function parsecommandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
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
    basicparser(filename)

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

"""
    parseconfig(filename)

Parse a configuration file and return a settings dict.

The config syntax is very simple: each line consists of a parameter
name and a value (unquoted), e.g. `nniches 2`. `#` is the comment character.
"""
function parseconfig(configfilename::String)
    config = basicparser(configfilename)
    settings = Dict{String, Any}()
    defaults = defaultSettings()
    for c in config
        if length(c) != 2
            simlog("Bad config file syntax: $c", settings, 'w', "")
        elseif c[1] in keys(defaults)
            value = c[2]
            if !(typeof(defaults[c[1]]) <: AbstractString)
                try
                    value = parse(typeof(defaults[c[1]]), c[2]) # or Meta.parse with the old functionality
                catch
                    simlog("$(c[1]) not of type $(typeof(defaults[c[1]])).",
                           settings, 'w', "")
                end
            end
            settings[c[1]] = value
        else
            # XXX maybe parse anyway
            simlog(c[1]*" is not a recognized parameter!", settings, 'w', "")
        end
    end
    settings
end

"""
    readmapfile(mapfilename, settings)

Parse a map file and return the number of timesteps this map is to be used for
(first line of the file) and the patch definitions. The latter is used by
`createworld` and `updateworld!`.
"""
function readmapfile(mapfilename::String, settings::Dict{String, Any})
    simlog("Reading map file $mapfilename.", settings)
    if isfile(mapfilename)
        maptable = basicparser(mapfilename)
        timesteps = parse(Int, maptable[1][1])
        if length(maptable[1]) != 1 || !isa(timesteps, Integer)
            timesteps = 10
            simlog("Invalid timestep information in the mapfile. Setting timesteps to 10.", settings, 'w')
        end
    else
        simlog("No map definition file provided. Assuming a two-cell world for 10 time steps.", settings, 'w')
        timesteps = 10
        maptable = [["",""], ["1", "1", "1", "initpop"], ["2", "2", "1"]]
    end
    return timesteps,maptable[2:end]
end
