# Output functions for GeMM

function printheader(settings::Dict{String, Any}, io::IO = stdout, sep::String = "\t")
    #XXX Transfer to a dynamic system? (As in createworld()?)
    print(io, "time", sep)
    print(io, "patch_no", sep)
    print(io, "xloc", sep)
    print(io, "yloc", sep)
    print(io, "temp", sep)
    print(io, "area", sep)
    settings["nniches"] > 1 && print(io, "prec", sep)
    settings["nniches"] > 2 && print(io, "nicheb", sep)
    print(io, "island", sep)
    print(io, "isolation", sep)
    print(io, "invasible", sep)
    print(io, "initpop", sep)
    print(io, "id", sep)
    print(io, "lineage", sep)
    print(io, "age", sep)
    print(io, "new", sep)
    print(io, "tempadaption", sep)
    print(io, "precadaption", sep)
    print(io, "size", sep)
    print(io, "lnkgunits", sep)
    print(io, "ngenes", sep)
    for key in settings["traitnames"]
        print(io, key, sep)
    end
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
            print(io, timestep, sep)
            print(io, patch.id, sep)
            print(io, patch.location[1], sep)
            print(io, patch.location[2], sep)
            print(io, patch.temp, sep)
            print(io, patch.area, sep)
            settings["nniches"] > 1 && print(io, patch.prec, sep)
            settings["nniches"] > 2 && print(io, patch.nicheb, sep)
            patch.isisland ? print(io, 1, sep) : print(io, 0, sep)
            patch.isolated ? print(io, 1, sep) : print(io, 0, sep)
            patch.invasible ? print(io, 1, sep) : print(io, 0, sep)
            patch.initpop ? print(io, 1, sep) : print(io, 0, sep)
            print(io, ind.id, sep)
            print(io, ind.lineage, sep)
            print(io, ind.age, sep)
            ind.marked ? print(io, 1, sep) : print(io, 0, sep)
            print(io, ind.tempadaption, sep)
            print(io, ind.precadaption, sep)
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
                vstr = "\""
                        for x in value
                            vstr *= string(x) * ","
                        end
                        value = vstr[1:end-1] * "\""
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
function writedata(world::Array{Patch,1}, settings::Dict{String, Any}, timestep::Int)
    if settings["raw"]
        basename = "s" * string(settings["seed"])
        basename = joinpath(settings["dest"], basename)
        filename = basename * ".tsv"
        simlog("Writing data \"$filename\"", settings)
        open(filename, "a") do file
            dumpinds(world, settings, timestep, file)
        end
    end
    if settings["stats"]
        basename = "stats_s" * string(settings["seed"])
        basename = joinpath(settings["dest"], basename)
        filename = basename * ".tsv"
        simlog("Writing stats to \"$filename\"", settings)
        open(filename, "a") do file
            printpopstats(file, world, settings, timestep)
        end
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
    #CAVE Should not use `whoiswho`
    lineages = unique(reduce(vcat, map(p -> collect(keys(p.whoiswho)), world)))
    div = round.(diversity(world), digits = 3)
    space = freespace(world)
    simlog("Metacommunity size: $popsize, lineages: $(length(lineages))", settings)
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
        #CAVE Should not use `whoiswho`
        for l in keys(p.whoiswho)
            simlog("$timestep,$(p.location[1]),$(p.location[2]),$l,$(length(p.whoiswho[l]))", settings,
                   'i', "lineages.log", true)
        end
    end    
end

"""
    simlog(msg, settings, category)
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
        simlog("Invalid log category $category.", settings, 'w')
    end
end

function printpopheader(io::IO)
    print(io, "time\t", "x\t", "y\t", "temp\t", "prec\t", "area\t", "isisland")
    print(io, "\tlineage", "\tjuveniles", "\tadults", "\tmaxage", "\tmaxsize", "\tmaxtempadaption", "\tmaxprecadaption")
    traitnames =  ["compat", "dispmean", "dispshape", "maxtraitvar",
    "medtraitvar", "mintraitvar", "ngenes", "nlnkgunits", "precopt",
    "prectol", "repsize", "reptol", "seedsize", "tempopt", "temptol"]
    for traitname in traitnames
        print(io, "\t", traitname, "min")
        print(io, "\t", traitname, "max")
        print(io, "\t", traitname, "med")
        print(io, "\t", traitname, "std")
    end
    print(io, "\treplicate", "\tconf")
    println(io)
end
    
function printpopstats(io::IO, world::Array{Patch, 1}, settings::Dict{String, Any}, timestep::Integer)
    timestep == 0 && printpopheader(io)
    traitnames =  ["compat", "dispmean", "dispshape", "maxtraitvar",
    "medtraitvar", "mintraitvar", "ngenes", "nlnkgunits", "precopt",
    "prectol", "repsize", "reptol", "seedsize", "tempopt", "temptol"]
    for patch in world
        lineages = unique(map(i -> i.lineage, patch.community))
        for lineage in lineages
            print(io, timestep, "\t", patch.location[1], "\t", patch.location[2],
                  "\t", patch.temp, "\t", patch.prec, "\t", patch.area, "\t", patch.isisland)
            popidxs = findall(i -> i.lineage == lineage, patch.community)
            population = patch.community[popidxs]
            adultidxs = findall(i -> i.size > i.traits["repsize"], patch.community[popidxs])
            print(io, "\t", population[1].lineage, "\t", length(popidxs) - length(adultidxs), "\t", length(adultidxs), "\t",
                  maximum(map(i -> i.age, population)), "\t", maximum(map(i -> i.size, population)),
                  "\t", maximum(map(i -> i.tempadaption, population)), "\t", maximum(map(i -> i.precadaption, population)))
            poptraitdict = Dict{String, Array{Float64, 1}}()
            for traitname in traitnames
                poptrait = map(i -> i.traits[traitname], population)
                print(io, "\t", minimum(poptrait))
                print(io, "\t", maximum(poptrait))
                print(io, "\t", median(poptrait))
                print(io, "\t", std(skipmissing(poptrait)))
            end
            print(io, "\t", settings["seed"], "\t", settings["config"])
            println(io)
        end
    end
end
