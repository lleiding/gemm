# Output functions for GeMM

"""
    printheader(settings, io, sep)

Print a list of property names to the given IO stream. This is a helper function
for `dumpinds`.
"""
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
    print(io, "new", sep)
    print(io, "tempadaptation", sep)
    print(io, "precadaptation", sep)
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
WARNING: this produces *very* large files!
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
            ind.marked ? print(io, 1, sep) : print(io, 0, sep)
            print(io, ind.tempadaptation, sep)
            print(io, ind.precadaptation, sep)
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

"""
    makefasta(world, settings, io, onlyisland, sep)

Record the genome of every individual currently alive to the given IO stream.
WARNING: this produces *very* large files!
"""
function makefasta(world::Array{Patch, 1}, settings::Dict{String, Any}, io::IO = stdout, onlyisland::Bool = false, sep::String = "_")
    for patch in world
        (onlyisland && !patch.isisland) && continue
        lineage = ""
        for ind in patch.community
            (!patch.isisland && settings["static"] && ind.lineage == lineage) && continue # only one individual per species on mainland
            chrmno = 0
            for chrm in ind.genome
                chrmno += 1
                geneno = 0
                for gene in chrm.genes
                    geneno += 1
                    traits = ""
                    if length(gene.codes) == 0
                        traits *= ","
                    else
                        for trait in gene.codes
                            traits *= string(settings["traitnames"][trait.nameindex]) * ":" * string(trait.value) * ","
                        end
                    end
                    header = ">" * ind.lineage * sep * string(ind.id) * sep * string(chrmno) * sep * string(geneno) * sep * traits
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
                isempty(m) && continue
                cp(m, joinpath(settings["dest"], basename(m)), force = true) # most likely replicates with same parameters
            end
        end
    end
end

"""
    writesettings(settings)

Record the settings actually used for a simulation run (cf. `getsettings`).
Creates a config file that can be used for future replicate runs.
Also records a time stamp and the current git commit.
"""
function writesettings(settings::Dict{String, Any})
    if isempty(basename(settings["config"]))
        settingspath = "settings.conf"
    else
        settingspath = basename(settings["config"])
    end
    open(joinpath(settings["dest"], settingspath), "w") do f
        println(f, "#\n# --- Island speciation model settings ---")
        println(f, "# This file was generated automatically.")
        println(f, "# Simulation run on $(Dates.format(Dates.now(), "d u Y HH:MM:SS"))")
        println(f, "# $(split(read(pipeline(`git log`, `head -1`), String), "\n")[1])\n")
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

"""
    writedata(world, settings, timestep)

Writes simulation output from `world` to separate table and fasta files. (Which
data is recorded depends on the settings.) `timestep` and `setting` information
is used for file name creation.
"""
function writedata(world::Array{Patch,1}, settings::Dict{String, Any}, timestep::Int)
    if settings["raw"]
        filename = "inds_s" * string(settings["seed"])
        filename = joinpath(settings["dest"], filename)
        filename = filename * ".tsv"
        simlog("Writing data \"$filename\"", settings)
        open(filename, "a") do file
            dumpinds(world, settings, timestep, file)
        end
    end
    if settings["stats"]
        filename = "pops_s" * string(settings["seed"])
        filename = joinpath(settings["dest"], filename)
        filename = filename * ".tsv"
        simlog("Writing stats to \"$filename\"", settings)
        open(filename, "a") do file
            printpopstats(file, world, settings, timestep)
        end
    end
    if settings["fasta"]
        filename = "seqs_s" * string(settings["seed"])
        filename = joinpath(settings["dest"], filename)
        filename = filename * ".fa"
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

Save the abundance of each lineage per patch.
"""
function recordlineages(world::Array{Patch,1}, settings::Dict{String, Any}, timestep::Int)
    if !isfile(joinpath(settings["dest"], "lineages.log"))
        simlog("t,X,Y,lineage,abundance,temp,prec", settings, 'i', "lineages.log", true)
    end
    for p in world
        #CAVE Should not use `whoiswho`
        for l in keys(p.whoiswho)
            simlog("$timestep,$(p.location[1]),$(p.location[2]),$l,$(length(p.whoiswho[l])),$(p.temp),$(p.prec)", settings,
                   'i', "lineages.log", true)
        end
    end    
end

"""
    printpopheader(io)

Print a list of property names to the given IO stream. This is a helper function
for `printpopstats`.
"""
function printpopheader(io::IO)
    print(io, "time", "\tx", "\ty", "\ttemp", "\tprec", "\tarea", "\tisisland")
    print(io, "\tlineage", "\tjuveniles", "\tadults", "\tmaxsize", "\ttempadaptationmed", "\tprecadaptationmed")
    traitnames =  ["compat", "compatsd", "dispmean", "dispmeansd", "dispshape", "dispshapesd", 
                   "ngenes", "nlnkgunits", "precopt", "precoptsd", "prectol", "prectolsd",
                   "repsize", "repsizesd", "reptol", "reptolsd", "seedsize", "seedsizesd",
                   "tempopt", "tempoptsd", "temptol", "temptolsd"]
    for traitname in traitnames
        print(io, "\t", traitname, "min")
        print(io, "\t", traitname, "max")
        print(io, "\t", traitname, "med")
        print(io, "\t", traitname, "std")
    end
    print(io, "\treplicate", "\tconf")
    println(io)
end

"""
    printpopstats(io, world, settings, timestep)

Record statistical information (maximum, minimum, median, standard deviation)
for a range of individual properties, as seen over the whole world population.
"""
function printpopstats(io::IO, world::Array{Patch, 1}, settings::Dict{String, Any}, timestep::Integer)
    timestep == 0 && printpopheader(io)
    traitnames =  ["compat", "compatsd", "dispmean", "dispmeansd", "dispshape", "dispshapesd", 
                   "ngenes", "nlnkgunits", "precopt", "precoptsd", "prectol", "prectolsd",
                   "repsize", "repsizesd", "reptol", "reptolsd", "seedsize", "seedsizesd",
                   "tempopt", "tempoptsd", "temptol", "temptolsd"]
    for patch in world
        lineages = unique(map(i -> i.lineage, patch.community))
        for lineage in lineages
            print(io, timestep, "\t", patch.location[1], "\t", patch.location[2],
                  "\t", patch.temp, "\t", patch.prec, "\t", patch.area, "\t", patch.isisland)
            popidxs = findall(i -> i.lineage == lineage, patch.community)
            population = patch.community[popidxs]
            adultidxs = findall(i -> i.size > i.traits["repsize"], patch.community[popidxs])
            print(io, "\t", population[1].lineage, "\t", length(popidxs) - length(adultidxs), "\t", length(adultidxs),
                  "\t", maximum(map(i -> i.size, population)), "\t", median(map(i -> i.tempadaptation, population)),
                  "\t", median(map(i -> i.precadaptation, population)))
            poptraitdict = Dict{String, Array{Float64, 1}}()
            for traitname in traitnames
                poptrait = map(i -> i.traits[traitname], population)
                print(io, "\t", minimum(poptrait))
                print(io, "\t", maximum(poptrait))
                print(io, "\t", median(poptrait))
                print(io, "\t", std(skipmissing(poptrait))) # CAVEAT: this returns NaN if only one individual
            end
            print(io, "\t", settings["seed"], "\t", settings["config"])
            println(io)
        end
    end
end

"""
    simlog(msg, settings, category, logfile, onlylog)

Write a log message to STDOUT/STDERR and the specified logfile 
(if logging is turned on in the settings).

Categories: `d` (debug), `i` (information, default), `w` (warn), `e` (error)

If `logfile` is the empty string (default: "simulation.log"), the message will 
only be printed to the screen. If `onlylog` is true (default: false), the
message is not printed to screen but only to the log.
"""
function simlog(msg::String, settings::Dict{String, Any}, category='i', logfile="simulation.log", onlylog=false)
    (isa(category, String) && length(category) == 1) && (category = category[1])
    function logprint(msg::String, settings::Dict{String, Any}, tostderr=false)
        if tostderr || !(settings["quiet"] || onlylog)
            tostderr ? iostr = stderr : iostr = stdout
            println(iostr, msg)
        end
        if settings["logging"] && length(logfile) > 0
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

"""
    writephylo(phylo, settings, t)

Write the phylogeny of local populations.
Each node represents the entire population of individuals of a certain lineage in the same location.
Output format is a comma-separated file with information on the branching time,
the parent population and the offspring population
"""
function writephylo(phylo, settings, t)
    basename = "phylo"
    basename = joinpath(settings["dest"], basename)
    filename = basename * ".csv"
    simlog("Writing phylo \"$filename\"", settings)
    open(filename, "a") do file
        for i in phylo
            println(file, t, ",", i)
        end
    end
end
