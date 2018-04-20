# IO functions for GeMM

function readmapfile(filename::String)
    println("Reading file \"$filename\"...")
    mapstrings = String[]
    open(filename) do file
        mapstrings = readlines(file)
    end
    mapstrings = filter(x->!all(isspace,x),mapstrings) # remove empty lines
    mapstrings = filter(x->x[1]!='#',mapstrings) # remove comment lines
    mapsubstrings = map(split,mapstrings)
    mapentries = map(x->map(String,x),mapsubstrings)
    timesteps = 0
    try
        timesteps = parse(Int,filter(x->size(x,1)==1,mapentries)[1][1])
    catch
        timesteps = 1000
        warn("your mapfile \"$filename\" does not include timestep information. Assuming $timesteps timesteps.")
    end
    mapentries = filter(x->size(x,1)>1,mapentries)
    return timesteps,mapentries
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
                            traits *= trait.name * "$(trait.value)" * ","
                        end
                    end
                    header = ">$counter x$(patch.location[1]) y$(patch.location[2]) $(ind.lineage) c$chrmno g$geneno $traits"
                    println(io, header)
                    println(io, gene.sequence)
                end
            end
            lineage = ind.lineage
        end
    end
end

"""
    setupdatadir(dir)
creates the output directory. `dir` is a string containing the directory path.
"""
function setupdatadir(settings::Dict{String, Any})
    if !isdir(settings["dest"])
        try
            mkdir(settings["dest"])
        catch
            warn("could not create directory ", settings["dest"], ". Assuming it already exists...")
        end
    end
end

#TODO: complete docstring!
"""
    writedata(world, path, settings, seed, timestep)
writes simulation output from `world` to separate table and fasta files.
`path`, `seed`, `timestep` and `setting` information is used for file name creation.
"""
function writedata(world::Array{Patch,1}, mappath::String, settings::Dict{String, Any}, seed::Int64, timestep::Int64)
    mapfile = split(mappath, "/")[end]
    length(mapfile) == 0 && return
    basename = "$(settings["dest"])" * "/" * mapfile * "_s" * "$seed" * "_lnk" * settings["linkage"] * "_tol" * settings["tolerance"] * "_t" * "$timestep"
    counter = 0
    extension = ""
    while ispath(basename * extension * ".tsv") || ispath(basename * extension * ".fa")
        extension = "_$counter"
        counter += 1
        if counter > 9
            warn("could not write to ", basename, extension, ": file exists. \n
Continuing anyway - data might be identical.")
            return
        end
    end
    basename *= extension
    filename = basename * ".tsv"
    println("Writing data \"$filename\"")
    open(filename, "w") do file
        dumpinds(world, file, "\t", settings["static"] || timestep > 1)
    end
    filename = basename * ".fa"
    println("Writing fasta \"$filename\"")
    open(filename, "w") do file
        makefasta(world, file, "", settings["static"] || timestep > 1)
    end
end

#TODO: complete docstring!
"""
    writerawdata(world, path, settings, seed, t)
writes raw julia data of the complete simulation state from `world` to a file in `path`.
`seed`, `setting` and `t` (timestep) information is used for file name creation.
"""
function writerawdata(world::Array{Patch,1}, mappath::String, settings::Dict{String, Any}, seed::Int64, timestep::Int64)
    mapfile = split(mappath, "/")[end]
    filename = "$(settings["dest"])" * "/" * mapfile * "_s" * "$seed" * "_lnk" * settings["linkage"] * "_tol" * settings["tolerance"] * "_t" * "$timestep" * ".jl"
    counter = 0
    extension = ""
    while ispath(filename * extension)
        extension = "_$counter"
        counter += 1
        counter > 9 && error("file \"$filename$extension\" exists. Please clear your directory.")
    end
    filename *= extension
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
    recordcolonizers(colos, path, settings, seed, t)
writes raw julia data of the colonizing individuals `colos` at timestep `t` to a file in `path`.
`seed` and `setting` information is used for file name creation.
"""
function recordcolonizers(colonizers::Array{Individual, 1}, mappath::String, settings::Dict{String, Any}, seed::Int64, timestep::Int64)
    mapfile = split(mappath, "/")[end]
    record = (timestep, colonizers)
    filename = "$(settings["dest"])" * "/" * mapfile * "_s" * "$seed" * "_lnk" * settings["linkage"] * "_tol" * settings["tolerance"] * "_colonizers" * ".jl"
    touch(filename)
    println("Colonisation! Writing data to ", filename, "...")
    open(filename, "a") do file
        println(file, record)
    end
end

