# initialisation functions for GeMM

function genesis()
    community = Individual[]
    totalmass = 0.0
    while true
        newind = createind()
        # population size determined by adult size and temperature niche optimum
        popsize = round(fertility * newind.traits["repsize"]^(-1/4) * exp(-act/(boltz*newind.traits["tempopt"]))) 
        popmass = popsize * newind.size
        if totalmass + popmass > settings["cellsize"] # stop loop if cell is full
            if popmass >= settings["cellsize"]*0.9 #make sure the cell is full enough
                break
            end
        end
        totalmass += popmass
        for i in 1:popsize
            !settings["static"] && (newind = deepcopy(newind))
            push!(community, newind)
        end
    end
    community
end

function createworld(maptable::Array{Array{String,1},1})
    simlog("Creating world...")
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
        global newpatch = Patch(id, (xcord, ycord), area)
        # parse other parameter options
        for p in entry[4:end]
            varval = split(p, '=')
            var = varval[1]
            if !(var in map(string, fieldnames(Patch)))
                simlog("Unrecognized patch parameter $var.", 'w')
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
                    simlog("Invalid patch parameter type $var: $val", 'w')
                    continue
                end
            end
            eval(parse("newpatch."*string(var)*" = $val"))
        end
        newpatch.initpop && append!(newpatch.community, genesis())
        push!(world, newpatch)
        global newpatch = nothing #clear memory
    end
    world
end

function updateworld!(world::Array{Patch,1},maptable::Array{Array{String,1},1},cellsize::Float64)
    #TODO: add functionality to remove patches!
    simlog("Updating world...")
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
        idx = find(x->x.id == id, world)
        if length(idx) == 0
            isnew = true
            global newpatch = Patch(id, (xcord, ycord), cellsize)
        else
            isnew = false
            global newpatch = world[idx[1]]
        end
        # parse other parameter options
        for p in entry[4:end]
            varval = split(p, '=')
            var = varval[1]
            if !(var in map(string, fieldnames(Patch)))
                simlog("Unrecognized patch parameter $var.", 'w')
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
                    simlog("Invalid patch parameter type $var: $val", 'w')
                    continue
                end
            end
            eval(parse("newpatch."*string(var)*" = $val"))
        end
        if isnew
            push!(world, newpatch)
            global newpatch = nothing #clear memory
        end
    end
    world
end
