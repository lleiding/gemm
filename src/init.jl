# initialisation functions for GeMM

function genesis(settings::Dict{String,Any})
    community = Individual[]
    totalmass = 0.0
    while true 
        lineage = randstring(4)
        meangenes = length(settings["traitnames"])
        ngenes = rand(Poisson(meangenes))
        ngenes < 1 && (ngenes = 1)
        traits = createtraits(settings) ### !!!
        genes = creategenes(ngenes,traits) ### !!!
        if settings["linkage"] == "none"
            nchrms = length(genes)
        elseif settings["linkage"] == "full"
            nchrms = 1
        else
            nchrms = rand(1:length(genes))
        end
        chromosomes = createchrs(nchrms,genes)
        traitdict = chrms2traits(chromosomes, settings["traitnames"]) ### !!!
        popsize = round(fertility * traitdict["repsize"]^(-1/4) * exp(-act/(boltz*traitdict["tempopt"]))) # population size determined by adult size and temperature niche optimum
        settings["initadults"]? indsize = traitdict["repsize"] : indsize = traitdict["seedsize"]
        popmass = popsize * indsize
        if totalmass + popmass > settings["cellsize"] # stop loop if cell is full
            break
        end
        totalmass += popmass
        newind = Individual(lineage, chromosomes, traitdict, 0, false, 1.0, indsize)
        for i in 1:popsize
            !settings["static"] && (newind = deepcopy(newind))
            push!(community, newind)
        end
    end
    community
end

function createworld(maptable::Array{Array{String,1},1}, settings::Dict{String,Any})
    info("Creating world...")
    world = Patch[]
    for entry in maptable
        size(entry,1) < 3 && error("please check your map file for incomplete or faulty entries. \n
Each line must contain patch information with at least \n
\t - a unique integer ID, \n
\t - an integer x coordinate, \n
\t - an integer y coordinate, \n
separated by a whitespace character (<ID> <x> <y>).")
        # create the basic patch
        id = parse(Int, entry[1])
        xcord = parse(Int, entry[2])
        ycord = parse(Int, entry[3])
        area = settings["cellsize"]
        # XXX the 'global' here is a hack so that I can use eval() later on
        # (this always works on the global scope)
        global newpatch = Patch(id, (xcord, ycord))
        # parse other parameter options
        for p in entry[4:end]
            varval = split(p, '=')
            var = varval[1]
            if !(var in map(string, fieldnames(Patch)))
                warn("Unrecognized patch parameter $var.")
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
                    warn("Invalid patch parameter type $var: $val")
                    continue
                end
            end
            eval(parse("newpatch."*string(var)*" = $val"))
        end
        !newpatch.isisland && append!(newpatch.community, genesis(settings))
        push!(world, newpatch)
        global newpatch = nothing #clear memory
    end
    world
end

function updateworld!(world::Array{Patch,1},maptable::Array{Array{String,1},1},cellsize::Float64)
    #FIXME Update to reflect createworld()
    #FIXME quite possibly degraded due to disuse
    #TODO: add functionality to remove patches!
    info("Updating world...")
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
        size(entry,1) > 3 ? temperature = parse(Float64, entry[4]) : temperature = 298.0
        isisland = false
        if size(entry,1) > 4 && contains(lowercase(entry[5]),"island")
            isisland = true # islands do not receive an initial community
        end
        isolated = false
        if size(entry,1) > 5 && contains(lowercase(entry[6]),"isolated")
            isolated = true
        end
        try
            p = find(x->x.id==id,world)
            world[p].location = (xcord,ycord)
            world[p].altitude = temperature
            world[p].isisland = isisland
            world[p].isolated = isolated
        catch
            newpatch = Patch(id,(xcord,ycord),temperature,cellsize,isisland)
            newpatch.isolated = isolated
            push!(world,newpatch)
        end
    end
    world
end
