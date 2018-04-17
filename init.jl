# initialisation functions for GeMM


function genesis(settings::Dict{String,Any},
                 nspecs::Int64=10000000, popsize::Int64 = 0, # about 600-900 species per cell
                 traitnames::Array{String,1} = ["dispmean",
                                                "dispshape",
                                                "mutprob",
                                                "precopt",
                                                "prectol",
                                                "repradius",
                                                "repsize",
                                                "reptol",
                                                "seedsize",
                                                "tempopt",
                                                "temptol"]) # minimal required traitnames
    community = Individual[]
    totalmass = 0.0
    while true 
        lineage = randstring(4)
        meangenes = length(traitnames)
        ngenes = rand(Poisson(meangenes))
        ngenes < 1 && (ngenes = 1)
        traits = createtraits(traitnames, settings)
        genes = creategenes(ngenes,traits)
        if settings["linkage"] == "none"
            nchrms = length(genes)
        elseif settings["linkage"] == "full"
            nchrms = 1
        else
            nchrms = rand(1:length(genes))
        end
        chromosomes = createchrs(nchrms,genes)
        traitdict = chrms2traits(chromosomes)
        popsize = round(fertility * traitdict["repsize"]^(-1/4) * exp(-act/(boltz*traitdict["tempopt"]))) # population size determined by adult size and temperature niche optimum
        popmass = popsize * traitdict["seedsize"]
        if totalmass + popmass > nspecs # stop loop if cell is full
            break
        end
        totalmass += popmass
        newind = Individual(lineage, chromosomes, traitdict, 0, false, 1.0, traitdict["seedsize"])
        for i in 1:popsize
            push!(community, newind)
        end
    end
    community
end

function createworld(maptable::Array{Array{String,1},1}, settings::Dict{String,Any})
    println("Creating world...")
    world = Patch[]
    length(maptable) > 1 ? area = cellsize : area = 1
    for entry in maptable
        size(entry,1) < 3 && error("please check your map file for incomplete or faulty entries. \n
Each line must contain patch information with at least \n
\t - a unique integer ID, \n
\t - an integer x coordinate, \n
\t - an integer y coordinate, \n
separated by a whitespace character (<ID> <x> <y>).")
        id = parse(Int64, entry[1])
        xcord = parse(Int64, entry[2])
        ycord = parse(Int64, entry[3])
        size(entry,1) > 3 ? temperature = parse(Float64, entry[4]) : temperature = 298.0
        isisland = false
        if size(entry,1) > 4 && contains(lowercase(entry[5]),"island")
            isisland = true # islands do not receive an initial community
        end
        newpatch = Patch(id,(xcord,ycord),temperature,area,isisland)
        if size(entry,1) > 5 && contains(lowercase(entry[6]),"isolated")
            newpatch.isolated = true
        end
        if size(entry,1) > 6 # || settings["heterogeneity"] == "high"
            try
                newpatch.nichea = parse(Float64, entry[7])
            catch
                newpatch.nichea = rand() * 10
            end
        end
        if size(entry,1) > 7
            newpatch.nicheb = parse(Float64, entry[8])
        end
        !isisland && append!(newpatch.community,genesis(settings))
        push!(world,newpatch)
    end
    world
end

function updateworld!(world::Array{Patch,1},maptable::Array{Array{String,1},1}) #TODO: add functionality to remove patches!
    println("Updating world...")
    area = 1 # CAVE: just for now...
    for entry in maptable
        size(entry,1) < 3 && error("please check your map file for incomplete or faulty entries. \n
                            Each line must contain patch information with at least \n
                            \t - a unique integer ID, \n
                            \t - an integer x coordinate, \n
                            \t - an integer y coordinate, \n
                            separated by a whitespace character (<ID> <x> <y>).")
        id = parse(Int64, entry[1])
        xcord = parse(Int64, entry[2])
        ycord = parse(Int64, entry[3])
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
            newpatch = Patch(id,(xcord,ycord),temperature,area,isisland)
            newpatch.isolated = isolated
            push!(world,newpatch)
        end
    end
    world
end
