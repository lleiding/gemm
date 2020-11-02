# initialisation functions for GeMM

"""
    createpop(settings)

Create a new, random individual and replicates it a certain number of times
(depending on metabolic variables) to create a new population of organisms.
Returns an array of individuals.
"""
function createpop(settings::Dict{String, Any})
    traits = createtraits(settings)
    traitdict = gettraitdict(traits, settings["traitnames"])
    popsize = 0
    if occursin("metabolic", settings["popsize"]) || occursin("single", settings["popsize"])
        # population size determined by adult size
        popsize = round(settings["fertility"] * traitdict["repsize"] ^ (-1 / 4) *
                        exp(-act / (boltz * 298.0)))
    elseif occursin("bodysize", settings["popsize"])
        # population size up to 25% of the maximum possible in this cell
        quarterpopsize = Integer(ceil((settings["cellsize"] / traitdict["repsize"]) / 4))
        popsize = rand(1:quarterpopsize) + 1
    elseif occursin("minimal", settings["popsize"]) || popsize == 0
        popsize = 2 #Takes two to tangle ;-)
    else
        simlog("Invalid value for `popsize`: $(settings["popsize"])", settings, 'e')
    end
    lineage = randstring(6)
    ngenes = rand(1:settings["maxloci"]) * length(settings["traitnames"])
    ngenes < 1 && (ngenes = 1)
    genes = creategenes(ngenes, traits, settings)
    if settings["linkage"] == "none"
        nchrms = length(genes)
    elseif settings["linkage"] == "full"
        nchrms = 1
    else
        linkage = rand(1:length(genes))
        nchrms = Integer(round(ngenes/linkage))
    end
    chromosomes = createchrms(nchrms, genes)
    locivar = rand()
    population = Individual[]
    for i in 1:popsize
        id = rand(Int32)
        chromosomes = deepcopy(chromosomes) # XXX this is potentially very slow!
        varyalleles!(chromosomes, settings, locivar)
        traitdict = gettraitdict(chromosomes, settings["traitnames"])
        if settings["indsize"] == "adult"
            indsize = traitdict["repsize"]
        elseif settings["indsize"] == "seed"
            indsize = traitdict["seedsize"]
        else
            # CAVEAT: sizes uniformally distributed?
            indsize = traitdict["seedsize"] + rand() * (traitdict["repsize"] - traitdict["seedsize"])
        end
        push!(population, Individual(lineage, chromosomes, traitdict, true, 1.0, 1.0,
                                     indsize, hermaphrodite, 0, id))
    end
    population
end

"""
    genesis(settings)

Create a new community, composed of random new species populations, for a patch.
Returns an array of individuals.
"""
function genesis(settings::Dict{String, Any})
    community = Individual[]
    totalmass = 0.0
    while true
        population = createpop(settings)
        popsize = length(population)
        # Check the cell capacity
        popmass = sum(map(x -> x.size, population))
        if totalmass + popmass > settings["cellsize"] * settings["overfill"]
            # stop loop if cell is full
            if totalmass >= settings["cellsize"] * 0.9 || occursin("single", settings["popsize"])
                #make sure the cell is full enough
                simlog("Cell is now $(round((totalmass/settings["cellsize"])*100))% full.", settings, 'd') #DEBUG
                break
            else
                continue
            end
        end
        totalmass += popmass
        append!(community, population)
        occursin("single", settings["popsize"]) && break
    end
    simlog("Patch initialized with $(length(community)) individuals.", settings, 'd') #DEBUG
    community
end

"""
    createpatch(patchentry, settings)

Create a new patch and set its parameters by parsing one line of the map file.
"""
function createpatch(patchentry::Array{String,1}, settings::Dict{String, Any})
    size(patchentry,1) < 3 && simlog("please check your map file for incomplete or faulty entries. \n
Each line must contain patch information with at least \n
\t - a unique integer ID, \n
\t - an integer x coordinate, \n
\t - an integer y coordinate, \n
separated by a whitespace character (<ID> <x> <y>).", settings, 'e')
    # create the basic patch
    id = parse(Int, patchentry[1])
    xcord = parse(Int, patchentry[2])
    ycord = parse(Int, patchentry[3])
    area = settings["cellsize"]
    # XXX the 'global' here is a hack so that I can use eval() later on
    # (eval() always works on the global scope)
    global newpatch = Patch(id, (xcord, ycord), area)
    simlog("Creating patch $id at $xcord/$ycord, size $area", settings, 'd') #DEBUG
    # parse other parameter options
    for p in patchentry[4:end]
        varval = split(p, '=')
        var = varval[1]
        if !(var in map(string, fieldnames(Patch)))
            simlog("Unrecognized patch parameter $var.", settings, 'w')
            continue
        elseif length(varval) < 2
            # if no value is specified, assume 'true'
            # (unless the variable is prepended with an exclamation mark, then take 'false')
            if var[1] == '!'
                var = var[2:end]
                val = false
            else
                val = true
            end
        else
            val = Meta.parse(varval[2])
        end
        # check for correct type and modify the new patch
        vartype = typeof(eval(Meta.parse("newpatch."*var)))
        if !isa(val, vartype)
            try
                val = convert(vartype, val)
            catch
                simlog("Invalid patch parameter type $var: $val", settings, 'w')
                continue
            end
        end
        eval(Meta.parse("newpatch."*string(var)*" = $val"))
    end
    newpatch
end

"""
    createworld(maptable, settings)

Use a parsed map file (as returned by `readmapfile`) to create the world. 
Initialises each patch with its parameters and a new community, then returns
an array of patches.
"""
function createworld(maptable::Array{Array{String,1},1}, settings::Dict{String, Any})
    simlog("Creating world...", settings)
    world = Patch[]
    for entry in maptable
        newpatch = createpatch(entry, settings)
        if newpatch.initpop
            if settings["mode"] == "zosterops"
                append!(newpatch.community, zgenesis(newpatch, settings))
            elseif settings["indsize"] != "seed"
                append!(newpatch.community, genesis(settings))
            else
                append!(newpatch.seedbank, genesis(settings))
            end
        end
        push!(world, newpatch)
        global newpatch = nothing #clear memory used in `createpatch()` #XXX does this make sense?
    end
    world
end

"""
    updateworld!(world, maptable, cellsize)

Reinitialise the world from another parsed map file. Works analogously to 
`createworld`. Intended for use in scenarios where the model world changes
during a run (e.g. through global warming or island ontogeny).
"""
#FIXME this function needs to be folded into createworld(), it duplicates way too much code
function updateworld!(world::Array{Patch,1},maptable::Array{Array{String,1},1}, settings::Dict{String, Any})
    cellsize = settings["cellsize"]
    simlog("Updating world...", settings)
    allids = Int[]
    for entry in maptable
        size(entry,1) < 3 && error("please check your map file for incomplete or faulty entries. \n
                            Each line must contain patch information with at least \n
                            \t - a unique integer ID, \n
                            \t - an integer x coordinate, \n
                            \t - an integer y coordinate, \n
                            separated by a whitespace character (<ID> <x> <y>).")
        id = parse(Int, entry[1])
        push!(allids, id)
        xcord = parse(Int, entry[2])
        ycord = parse(Int, entry[3])
        # XXX the 'global' here is a hack so that I can use eval() later on
        # (this always works on the global scope)
        idx = findall(x -> x.id == id, world)
        if length(idx) == 0
            marked = true
            global newpatch = Patch(id, (xcord, ycord), cellsize)
        else
            marked = false
            world[idx[1]].location = (xcord, ycord)
            global newpatch = world[idx[1]]
        end
        # parse other parameter options
        for p in entry[4:end]
            varval = split(p, '=')
            var = varval[1]
            if !(var in map(string, fieldnames(Patch)))
                simlog("Unrecognized patch parameter $var.", settings, 'w')
                continue
            elseif length(varval) < 2
                val = true # if no value is specified, assume 'true'
            else
                val = Meta.parse(varval[2])
            end
            # check for correct type and modify the new patch
            vartype = typeof(eval(Meta.parse("newpatch."*var)))
            if !isa(val, vartype)
                try
                    val = convert(vartype, val)
                catch
                    simlog("Invalid patch parameter type $var: $val", settings, 'w')
                    continue
                end
            end
            eval(Meta.parse("newpatch."*string(var)*" = $val"))
        end
        if marked
            push!(world, newpatch)
            global newpatch = nothing #clear memory
        end
    end
    filter!(x -> x.id in allids, world)
    world
end

