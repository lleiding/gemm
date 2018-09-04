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

function printheader(settings::Dict{String, Any}, io::IO = stdout, sep::String = "\t")
        #XXX Transfer to a dynamic system? (As in createworld()?)
        print(io, "patch_no", sep)
        print(io, "xloc", sep)
        print(io, "yloc", sep)
        print(io, "temp", sep)
        print(io, "area", sep)
        print(io, "prec", sep)
        print(io, "island", sep)
        print(io, "isolation", sep)
        print(io, "invasible", sep)
        print(io, "initpop", sep)
        print(io, "id", sep)
        print(io, "lineage", sep)
        print(io, "age", sep)
        print(io, "new", sep)
        print(io, "fitness", sep)
        print(io, "size", sep)
        print(io, "lnkgunits", sep)
        print(io, "ngenes", sep)
        for key in settings["traitnames"]
            print(io, key, sep)
        end
        print(io, "time")
        println(io)
end

"""
    dumpinds(world, io, sep)
Output all data of individuals in `world` as table to `io`. Columns are separated by `sep`.
"""
function dumpinds(world::Array{Patch, 1}, settings::Dict{String, Any}, timestep::Int, io::IO = stdout, sep::String = "\t")
    timestep == 0 && printheader(settings, io, sep)
    onlyisland = settings["static"] && timestep > 1
    for patch in world
        #XXX Sometimes, this only dumps three or four individuals, with a population of >10⁴!
        # (Should be fixed)
        (onlyisland && !patch.isisland) && continue
        lineage = ""
        for ind in patch.community
            # only one individual per species on a static mainland (should be the case anyway)
            (!patch.isisland && settings["static"] && ind.lineage == lineage) && continue
            print(io, patch.id, sep)
            print(io, patch.location[1], sep)
            print(io, patch.location[2], sep)
            print(io, patch.temp, sep)
            print(io, patch.area, sep)
            print(io, patch.prec, sep)
            ## print(io, patch.nicheb, sep)
            patch.isisland ? print(io, 1, sep) : print(io, 0, sep)
            patch.isolated ? print(io, 1, sep) : print(io, 0, sep)
            patch.invasible ? print(io, 1, sep) : print(io, 0, sep)
            patch.initpop ? print(io, 1, sep) : print(io, 0, sep)
            print(io, ind.id, sep)
            print(io, ind.lineage, sep)
            print(io, ind.age, sep)
            ind.marked ? print(io, 1, sep) : print(io, 0, sep)
            print(io, ind.fitness, sep)
            print(io, ind.size, sep)
            print(io, length(ind.genome), sep)
            print(io, sum(map(x -> length(x.genes), ind.genome)), sep)
            for key in settings["traitnames"]
                try
                    print(io, ind.traits[key], sep)
                catch
                    print(io, "NA", sep)
                end
            end
            print(io, timestep)
            println(io)
            lineage = ind.lineage
        end
    end
end

function makefasta(world::Array{Patch, 1}, settings::Dict{String, Any}, io::IO = stdout, onlyisland::Bool = false, sep::String = "_")
    for patch in world
        (onlyisland && !patch.isisland) && continue
        lineage = ""
        for ind in patch.community
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
                            traits *= string(settings["traitnames"][trait.nameindex]) * ":" * string(trait.value) * ","
                        end
                    end
                    header = ">"*string(ind.id)*sep*"x"*string(patch.location[1])*sep*"y"*string(patch.location[2])
                    header *= sep*ind.lineage*sep*"c"*string(chrmno)*sep*"g"*string(geneno)*sep*traits
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
    if isdir(settings["dest"]) # && !isempty(readdir(settings["dest"])) ## prevent mixing up of ouptut data from parallel sims
        replicate = split(settings["dest"], "_")[end]
        if all(isnumeric, replicate)
            replicate = parse(UInt8, replicate) + 1 # throw an error if replicate > 255
            settings["dest"] = string(split(settings["dest"], "_")[1:(end-1)]...) * "_" * string(replicate)
        else
            replicate = 1
            settings["dest"] = string(settings["dest"], "_", replicate)
        end
        setupdatadir(settings)
    else
        mkpath(settings["dest"])
        simlog("Setting up output directory $(settings["dest"])", settings)
        writesettings(settings)
        if haskey(settings, "maps")
            for m in settings["maps"]
                cp(m, joinpath(settings["dest"], m), force = true) # most likely replicates with same parameters
            end
        end
        cp(joinpath("src", "constants.jl"), joinpath(settings["dest"], "constants.jl"), force = true)
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
    basename = mapfile * "_s" * string(settings["seed"])
    basename = joinpath(settings["dest"], basename)
    filename = basename * ".tsv"
    simlog("Writing data \"$filename\"", settings)
    open(filename, "a") do file
        dumpinds(world, settings, timestep, file)
    end
    if settings["fasta"]
        filename = basename * ".fa"
        simlog("Writing fasta \"$filename\"", settings)
        open(filename, "a") do file
            makefasta(world, settings, file, settings["static"] && timestep > 1)
        end
    end
end

"""
    recordstatistics(w)
Write out world properties to the log file for later analysis.
"""
function recordstatistics(world::Array{Patch,1}, settings::Dict{String, Any})
    if !isfile(joinpath(settings["dest"], "diversity.log"))
        simlog("population,freespace,lineages,alpha,beta,gamma", settings,
               'i', "diversity.log", true)
    end
    popsize = sum(x -> length(x.community), world)
    lineages = unique(reduce(vcat, map(p -> collect(keys(p.whoiswho)), world)))
    div = round.(diversity(world), digits = 3)
    space = freespace(world)
    simlog("Population size: $popsize, lineages: $(length(lineages))", settings)
    simlog("$popsize,$space,$(length(lineages)),$(div[1]),$(div[2]),$(div[3])", settings,
           'i', "diversity.log", true)
end

"""
    recordlineages(w)
Save the abundance of each lineage per patch
"""
function recordlineages(world::Array{Patch,1}, settings::Dict{String, Any}, timestep::Int)
    if !isfile(joinpath(settings["dest"], "lineages.log"))
        simlog("t,X,Y,lineage,population", settings, 'i', "lineages.log", true)
    end
    for p in world
        for l in keys(p.whoiswho)
            simlog("$timestep,$(p.location[1]),$(p.location[2]),$l,$(length(p.whoiswho[l]))", settings,
                   'i', "lineages.log", true)
        end
    end    
end

"""
    writephylo(w, settings, t)
Save the phylogeny of individuals currently in world `w`.
"""
function writephylo(world::Array{Patch,1}, settings::Dict{String, Any}, timestep::Int)
    if !isfile(joinpath(settings["dest"], "phylo.tsv"))
        open(joinpath(settings["dest"], "phylo$(settings["seed"]).tsv"), "w") do f
            println(f, "t\tpatch\tid\tparent")
        end
    end
    open(joinpath(settings["dest"], "phylo$(settings["seed"]).tsv"), "a") do f
        for p in world
            for i in 1:size(p.phylo, 1)
                print(f, timestep, "\t")
                print(f, p.id, "\t")
                print(f, p.phylo[i, 1], "\t")
                println(f, p.phylo[i, 2])
            end
        end
    end
end    

"""
    simlog(msg, logfile, category)
Write a log message to STDOUT/STDERR and the specified logfile 
(if logging is turned on in the settings).
Categories: d (debug), i (information, default), w (warn), e (error)
"""
function simlog(msg::String, settings::Dict{String, Any}, category='i', logfile="simulation.log", onlylog=false)
    (isa(category, String) && length(category) == 1) && (category = category[1])
    function logprint(msg::String, settings::Dict{String, Any}, tostderr=false)
        if tostderr || !(settings["quiet"] || onlylog)
            tostderr ? iostr = stderr : iostr = stdout
            println(iostr, msg)
        end
        if settings["logging"]
            open(joinpath(settings["dest"], logfile), "a") do f
                println(f, msg)
            end
        end
    end
    if category == 'i'
        logprint(msg, settings)
    elseif category == 'd'
        settings["debug"] && logprint("DEBUG: "*string(msg), settings)
    elseif category == 'w'
        logprint("WARNING: "*string(msg), settings, true)
    elseif category == 'e'
        logprint("ERROR: "*string(msg), settings, true)
        exit(1)
    else
        simlog("Invalid log category $category.", 'w')
    end
end
