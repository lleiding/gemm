# Input functions for GeMM

"""
    getsettings(configfile, seed)

Combines all configuration options to produce a single settings dict.
Precedence: function arguments - commandline parameters - config file - default values
"""
function getsettings(configfile::String = "", seed::Integer = 0)
    defaults = defaultSettings()
    commandline = parsecommandline()
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
    for s in ["debug", "quiet"]
        if !commandline[s] && s in keys(configs) && configs[s]
            settings[s] = true
        end
    end
    settings
end

"""
    parsecommandline()

Certain software parameters can be set via the commandline.
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

Reads in the file, strips out comments, concatenates lines ending
in '\', and then separates lines by whitespace.
Returns a 2d array representing the tokens in each line.
"""
function basicparser(filename::String)
    # Read in the file
    rawlines = String[]
    try
        open(filename) do file
            rawlines = readlines(file)
        end
    catch
        simlog("Could not open $filename!", 'e')
    end
    # Process each line
    l = 0
    lines = Array{String,1}[]
    while l < length(rawlines)
        l += 1
        ln = strip(split(rawlines[l], '#')[1]) # remove comments
        (isempty(ln)) && continue # remove empty lines
        if ln[end] == '\\' # concatenate lines ending with \
            rawlines[l+1] = ln[1:(end-1)]*strip(rawlines[l+1])
            continue
        end
        tokens = map(s -> convert(String, s), split(ln)) # tokenise
        push!(lines, tokens)
    end
    lines
end

"""
    parseconfig(filename)

Parse a configuration file and return a settings dict.

The config syntax is very simple: each line consists of a parameter name
and a value (unquoted), e.g. `nniches 2`. `#` is the comment character.

NOTE: parameter values must not contain any spaces! (e.g. when passing dicts)
"""
function parseconfig(configfilename::String)
    config = basicparser(configfilename)
    settings = Dict{String, Any}()
    defaults = defaultSettings()
    for c in config
        if length(c) != 2
            simlog("Bad config file syntax: $c", 'w', "")
        elseif c[1] in keys(defaults)
            value = c[2]
            if !(typeof(defaults[c[1]]) <: AbstractString)
                value = eval(Meta.parse(c[2]))
                if !isa(value, typeof(defaults[c[1]]))
                    try
                        value = convert(typeof(defaults[c[1]]), value)
                    catch
                        simlog("$(c[1]) not of type $(typeof(defaults[c[1]])).", 'w', "")
                    end
                end
            end
            settings[c[1]] = value
        else
            simlog(c[1]*" is not a recognized parameter!", 'w', "")
        end
    end
    settings
end

"""
    readmapfile(mapfilename)

Parse a map file and return the number of timesteps this map is to be used for
(first line of the file) and the patch definitions. The latter is used by
`createworld` and `updateworld!`.
"""
function readmapfile(mapfilename::String)
    simlog("Reading map file $mapfilename.")
    if isfile(mapfilename)
        maptable = basicparser(mapfilename)
        timesteps = parse(Int, maptable[1][1])
        if length(maptable[1]) != 1 || !isa(timesteps, Integer)
            timesteps = 10
            simlog("Invalid timestep information in the mapfile. Setting timesteps to 10.", 'w')
        end
    else
        simlog("No map definition file provided. Assuming a two-cell world for 10 time steps.", 'w')
        timesteps = 10
        maptable = [["",""], ["1", "1", "1", "initpop"], ["2", "2", "1"]]
    end
    return timesteps,maptable[2:end]
end
