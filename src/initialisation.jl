# initialisation functions for GeMM

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
        quarterpopsize = Integer(floor((settings["cellsize"] / traitdict["repsize"]) / 4))
        popsize = rand(2:quarterpopsize)
    elseif occursin("minimal", settings["popsize"]) || popsize == 0
        popsize = 2 #Takes two to tangle ;-)
    else
        simlog("Invalid value for `popsize`: $(settings["popsize"])", settings, 'e')
    end
    lineage = randstring(4)
    parentid = rand(Int32)
    ngenes = settings["avgnoloci"] * length(settings["traitnames"])
    ngenes < 1 && (ngenes = 1)
    genes = creategenes(ngenes, traits, settings)
    randchrms = rand(1:length(genes))
    if settings["linkage"] == "none"
        nchrms = length(genes)
    elseif settings["linkage"] == "full"
        nchrms = 1
    else
        nchrms = randchrms
    end
    chromosomes = createchrms(nchrms, genes)
    locivar = rand()
    population = Individual[]
    for i in 1:popsize
        id = rand(Int32)
        chromosomes = deepcopy(chromosomes)
        varyalleles!(chromosomes, settings, locivar)
        traitdict = gettraitdict(chromosomes, settings["traitnames"])
        if settings["indsize"] == "adult"
            indsize = traitdict["repsize"]
            age = 1
        elseif settings["indsize"] == "seed"
            indsize = traitdict["seedsize"]
            age = 0
        else
            indsize = traitdict["seedsize"] + rand() * traitdict["repsize"] # XXX: sizes shouldn't be uniformally dist'd
            age = 1
        end
        push!(population, Individual(lineage, chromosomes, traitdict, age, false, 1.0, indsize, id, parentid))
    end
    population
end

function genesis(settings::Dict{String, Any})
    community = Individual[]
    totalmass = 0.0
    while true
        population = createpop(settings)
        popsize = length(population)
        # Check the cell capacity
        popmass = sum(map(x -> x.size, population))
        if totalmass + popmass > settings["cellsize"] # stop loop if cell is full
            if totalmass >= settings["cellsize"] * 0.9 || occursin("single", settings["popsize"]) #make sure the cell is full enough
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

function createworld(maptable::Array{Array{String,1},1}, settings::Dict{String, Any})
    simlog("Creating world...", settings)
    world = Patch[]
    for entry in maptable
        size(entry,1) < 3 && simlog("please check your map file for incomplete or faulty entries. \n
Each line must contain patch information with at least \n
\t - a unique integer ID, \n
\t - an integer x coordinate, \n
\t - an integer y coordinate, \n
separated by a whitespace character (<ID> <x> <y>).", settings, 'e')
        # create the basic patch
        id = parse(Int, entry[1])
        xcord = parse(Int, entry[2])
        ycord = parse(Int, entry[3])
        area = settings["cellsize"]
        simlog("Creating patch $id at $xcord/$ycord, size $area", settings, 'd') #DEBUG
        # XXX the 'global' here is a hack so that I can use eval() later on
        # (this always works on the global scope)
        global newpatch = Patch(id, (xcord, ycord), area)
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
        if newpatch.initpop && settings["indsize"] != "seed"
            append!(newpatch.community, genesis(settings))
        elseif newpatch.initpop && !newpatch.isisland && settings["static"]
            append!(newpatch.seedbank, genesis(settings))
            lineage = ""
            for ind in newpatch.seedbank # store one sample individual for recording purposes
                if ind.lineage != lineage
                    push!(newpatch.community, ind)
                    lineage = ind.lineage
                end
            end
        elseif newpatch.initpop
            append!(newpatch.seedbank, genesis(settings))
        end
        push!(world, newpatch)
        global newpatch = nothing #clear memory
    end
    world
end

function updateworld!(world::Array{Patch,1},maptable::Array{Array{String,1},1},cellsize::Float64)
    #TODO: add functionality to remove patches!
    simlog("Updating world...", settings)
    for entry in maptable
        size(entry,1) < 3 && error("please check your map file for incomplete or faulty entries. \n
                            Each line must contain patch information with at least \n
                            \t - a unique integer ID, \n
                            \t - an integer x coordinate, \n
                            \t - an integer y coordinate, \n
                            separated by a whitespace character (<ID> <x> <y>).")
        id = parse(Int, entry[1])
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
                val = parse(varval[2])
            end
            # check for correct type and modify the new patch
            vartype = typeof(eval(parse("newpatch."*var)))
            if !isa(val, vartype)
                try
                    val = convert(vartype, val)
                catch
                    simlog("Invalid patch parameter type $var: $val", settings, 'w')
                    continue
                end
            end
            eval(parse("newpatch."*string(var)*" = $val"))
        end
        if marked
            push!(world, newpatch)
            global newpatch = nothing #clear memory
        end
    end
    world
end
