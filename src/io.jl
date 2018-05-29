# IO functions for GeMM

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
    if isa(settings["cellsize"], Integer)
        settings["cellsize"] *= 1e6 #convert tonnes to grams
    end
    settings
end


function parsecommandline()
    defaults = defaultSettings()
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
            help = "tolerance of sequence identity when reproducing (\"high\", \"evo\", \"low\" or \"none\")"
            arg_type = String
            range_tester = x->in(x,["high", "evo", "low", "none"])
            required = false
        "--dest", "-d"
            help = "output directory. Defaults to current date"
            arg_type = String
            required = false
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
    println("Reading file \"$filename\"...")
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
            warn("Bad config file syntax: $c")
        elseif c[1] in keys(defaults)
            value = parse(c[2])
            if !isa(value, typeof(defaults[c[1]]))
                warn("$(c[1]): expected $(typeof(settings[c[1]])), got $(typeof(value)).")
            end
            settings[c[1]] = value
        else
            warn(c[1]*" is not a recognized parameter!") # XXX maybe parse anyway
        end
    end
    return settings
end

function readmapfile(mapfilename::String)
    mapfile = basicparser(mapfilename)
    timesteps = parse(mapfile[1][1])
    if length(mapfile[1]) != 1 || !isa(timesteps, Integer)
        timesteps = 1000
        warn("Invalid timestep information in the mapfile. Setting timesteps to 1000.")
    end
    return timesteps,mapfile[2:end]
end

"""
    dumpinds(world, io, sep)
Output all data of individuals in `world` as table to `io`. Columns are separated by `sep`.
"""
function dumpinds(world::Array{Patch, 1}, io::IO = STDOUT, sep::String = "\t", onlyisland::Bool = false)
    header = true
    traitkeys = []
    counter = 0
    for patch in world
        (onlyisland && !patch.isisland) && continue # only one individual per species on mainland
        lineage = ""
        for ind in patch.community
            counter += 1
            (!patch.isisland && ind.lineage == lineage) && continue
            if header
                print(io, "id", sep)
                print(io, "xloc", sep)
                print(io, "yloc", sep)
                print(io, "temp", sep)
                print(io, "area", sep)
                print(io, "nichea", sep)
                ## print(io, "nicheb", sep)
                print(io, "island", sep)
                print(io, "isolation", sep)
                print(io, "counter", sep)
                print(io, "lineage", sep)
                print(io, "age", sep)
                print(io, "new", sep)
                print(io, "fitness", sep)
                print(io, "size", sep)
                print(io, "lnkgunits", sep)
                print(io, "ngenes", sep)
                traitkeys = keys(ind.traits)
                for key in traitkeys
                    print(io, key, sep)
                end
                println(io)
                header = false
            end
            print(io, patch.id, sep)
            print(io, patch.location[1], sep)
            print(io, patch.location[2], sep)
            print(io, patch.altitude, sep)
            print(io, patch.area, sep)
            print(io, patch.nichea, sep)
            ## print(io, patch.nicheb, sep)
            patch.isisland ? print(io, 1, sep) : print(io, 0, sep)
            patch.isolated ? print(io, 1, sep) : print(io, 0, sep)
            print(io, counter, sep)
            print(io, ind.lineage, sep)
            print(io, ind.age, sep)
            ind.isnew ? print(io, 1, sep) : print(io, 0, sep)
            print(io, ind.fitness, sep)
            print(io, ind.size, sep)
            print(io, length(ind.genome), sep)
            print(io, sum(map(x -> length(x.genes), ind.genome)), sep)
            for key in traitkeys
                try
                    print(io, ind.traits[key], sep)
                catch
                    print(io, "NA", sep)
                end
            end
            println(io)
            lineage = ind.lineage
        end
    end
end

function makefasta(world::Array{Patch, 1}, io::IO = STDOUT, sep::String = "", onlyisland::Bool = false)
    counter = 0
    for patch in world
        (onlyisland && !patch.isisland) && continue
        lineage = ""
        for ind in patch.community
            counter += 1
            (!patch.isisland && ind.lineage == lineage) && continue # only one individual per species on mainland
            chrmno = 0
            for chrm in ind.genome
                chrmno += 1
                geneno = 0
                for gene in chrm.genes
                    geneno += 1
                    traits = ""
                    if length(gene.codes) == 0
                        traits *= "neutral"
                    else
                        for trait in gene.codes
                            traits *= trait.name * string(trait.value) * ","
                        end
                    end
                    header = ">"*string(counter)*" x"*string(patch.location[1])*" y"*string(patch.location[2])
                    header *= " "*ind.lineage*" c"*string(chrmno)*" g"*string(geneno)*" "*traits
                    println(io, header)
                    println(io, num2seq(gene.sequence))
                end
            end
            lineage = ind.lineage
        end
    end
end

"""
    setupdatadir(dir)
Creates the output directory and copies relevant files into it.
If the output directory already includes files, create a new
directory by appending a counter.
"""
function setupdatadir(settings::Dict{String, Any})
    if isdir(settings["dest"]) && !isempty(readdir(settings["dest"]))
        replicate = split(settings["dest"], "_")[end]
        if all(isnumber, replicate)
            replicate = parse(UInt8, replicate) + 1 # throw an error if replicate > 255
            settings["dest"] = string(split(settings["dest"], "_")[1:(end-1)]...) * "_" * string(replicate)
        else
            replicate = 1
            settings["dest"] = string(settings["dest"], "_", replicate)
        end
        setupdatadir(settings)
    else
        info("Setting up output directory $(settings["dest"])")
        mkpath(settings["dest"])
        writesettings(settings)
        if haskey(settings, "maps")
            for m in settings["maps"]
                cp(m, joinpath(settings["dest"], m), remove_destination = true) # most likely replicates with same parameters
            end
        end
    end
end

function writesettings(settings::Dict{String, Any})
    open(joinpath(settings["dest"], settings["config"]), "w") do f
        println(f, "#\n# --- Island speciation model settings ---")
        println(f, "# This file was generated automatically.")
        println(f, "# Simulation run on $(Dates.format(Dates.now(), "d u Y HH:MM:SS"))\n#\n")
        for k in keys(settings)
            value = settings[k]
            if isa(value, String)
                value = "\"" * value * "\""
            elseif isa(value, Array)
                vstr = "\""  #" this comment is for emacs' syntax highlighting ;)
                for x in value
                    vstr *= string(x) * ","
                end
                value = vstr[1:end-1] * "\"" #" and this one.
            end
            println(f, "$k $value")
        end
    end
end

#TODO: complete docstring!
"""
    writedata(world, settings, timestep)
writes simulation output from `world` to separate table and fasta files.
`timestep` and `setting` information is used for file name creation.
"""
function writedata(world::Array{Patch,1}, settings::Dict{String, Any}, mapfile::String, timestep::Int)
    basename = mapfile * "_t" * string(timestep) * "_s" * string(settings["seed"])
    basename = joinpath(settings["dest"], basename)
    filename = basename * ".tsv"
    println("Writing data \"$filename\"")
    open(filename, "w") do file
        dumpinds(world, file, "\t", settings["static"] || timestep > 1)
    end
    if settings["fasta"]
        filename = basename * ".fa"
        println("Writing fasta \"$filename\"")
        open(filename, "w") do file
            makefasta(world, file, "", settings["static"] || timestep > 1)
        end
    end
end

#TODO: complete docstring!
"""
    writerawdata(world, settings, t)
writes raw julia data of the complete simulation state from `world` to file in the output directory.
`setting` and `t` (timestep) information is used for file name creation.
"""
function writerawdata(world::Array{Patch,1}, settings::Dict{String, Any}, mapfile::String, timestep::Int)
    filename = mapfile * "_t" * string(timestep) * "_s" * string(settings["seed"]) * ".jl"
    filename = joinpath(settings["dest"], filename)
    touch(filename)
    println("Writing raw data to \"$filename\"...")
    if timestep == 1
        open(filename, "w") do file
            println(file, world)
        end
    else
        island = filter(x -> x.isisland, world)
        open(filename, "w") do file
            println(file, island)
        end
    end
end

#TODO: complete docstring!
"""
    recordcolonizers(colos, settings, t)
writes raw julia data of the colonizing individuals `colos` at timestep `t` to file.
`seed` and `setting` information is used for file name creation.
"""
function recordcolonizers(colonizers::Array{Individual, 1}, settings::Dict{String, Any}, mapfile::String, timestep::Int)
    record = (timestep, colonizers)
    filename = mapfile * "_t" * string(timestep) * "_s" * string(settings["seed"] * "_colonizers.jl")
    filename = joinpath(settings["dest"], filename)
    touch(filename)
    println("Colonisation. Writing data to ", filename, "...")
    open(filename, "a") do file
        println(file, record)
    end
end

