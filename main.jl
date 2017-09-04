#!/usr/bin/env julia
## Metabolic Individual-based Genetically-explicit Meta-Community Model
## Ludwig Leidinger 2017
## <l.leidinger@gmx.net>
##
## Main simulation
##
## Synopsis: julia main.jl <MAPFILE> <RANDOMSEED> <TIMESTEPS>
## <MAPFILE> is a textfile containing information about the simulation arena
## Every line describes one patch in the following format:
## <ID> <X-COORDINATE> <Y-COORDINATE> [<TYPE>]

include("MIbGxMCmod.jl")


using MIbGxMCmod, Plots


function readmapfile(filename::String)
    mapstrings = String[]
    open(filename) do file
        mapstrings = readlines(file)
    end
    mapstrings = filter(x->!all(isspace,x),mapstrings) # remove empty lines
    mapstrings = filter(x->x[1]!='#',mapstrings) # remove comment lines
    mapentries = map(split,mapstrings)
    mapentries = map(x->map(String,x),mapentries)
    timesteps = 0
    try
        timesteps = parse(Int,filter(x->size(x,1)==1,mapentries)[1][1])
    catch
        timesteps = 1000
        warn("your mapfile \"$filename\" does not include timestep information. Assumed $timesteps timesteps.")
    end
    mapentries = filter(x->size(x,1)>1,mapentries)
    timesteps,mapentries
end

function createworld(maptable::Array{Array{String,1},1})
    world = Patch[]
    area = 100 # CAVE: just for now...
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
        size(entry,1) > 3 ? temperature = parse(Float64, entry[4]) : temperature = 298
        isisland = false
        if size(entry,1) > 4
            contains(lowercase(entry[5]),"island") && (isisland = true) # islands do not receive an initial community
        end
        newpatch = Patch(id,(xcord,ycord),temperature,area,isisland)
        if size(entry,1) > 5
            contains(lowercase(entry[6]),"isolated") && (newpatch.isolated = true)
        end
        !isisland && append!(newpatch.community,genesis())
        push!(world,newpatch)
    end
    world
end

function updateworld!(world::Array{Patch,1},maptable::Array{Array{String,1},1}) #TODO: add functionality to remove patches!
    area = 100 # CAVE: just for now...
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
        size(entry,1) > 3 ? temperature = parse(Float64, entry[4]) : temperature = 298
        isisland = false
        if size(entry,1) > 4
            contains(lowercase(entry[5]),"island") && (isisland = true) # islands do not receive an initial community
        end
        if size(entry,1) > 5
            contains(lowercase(entry[6]),"isolated") && (newpatch.isolated = true)
        end
        try
            p = find(x->x.id==id,world)
            world[p].location = (xcord,ycord)
            world[p].altitude = temperature
            world[p].isisland = isisland
            world[p].isolated = isolated
        catch
            newpatch = Patch(id,(xcord,ycord),temperature,area,isisland)
            push!(world,newpatch)
        end
    end
    world
end

function visualisation(world::Array{Patch,1},firstplot::Bool)
    xcords = map(x->Int(floor(x.location[1])),world)
    ycords = map(x->Int(floor(x.location[2])),world)
    values = map(x->size(x.community,1),world)
    xmin = minimum(xcords)
    xmax = maximum(xcords)
    xshift = 1 - xmin
    ymin = minimum(ycords)
    ymax = maximum(ycords)
    yshift = 1 - ymin
    mat = zeros(length(ymin:ymax),length(xmin:xmax))
    mat[:] = -10
    for i in eachindex(values)
        mat[ycords[i]+yshift,xcords[i]+xshift]=values[i]
    end
    if firstplot
        gr()
        heatmap(mat,aspect_ratio=1,show=true)
    else
        heatmap(mat,aspect_ratio=1,show=true)
    end
end

function analysis(world::Array{Patch,1})
    for patch in world
        print("Patch #",patch.id," @",patch.location," (",patch.isisland,"): \t")
        println(size(patch.community,1))
    end
end

function writedata(world::Array{Patch,1}, seed::Int64, mapfile::String)
    filename = mapfile * "_seed" * "$seed" * ".out"
    counter = 0
    while ispath(filename)
        filename *= "_1"
        counter += 1
        counter >= 3 && error("cannot create file \"$filename\". Please clear your directory.")
    end
    touch(filename)
    open(filename, "w") do file
        print(file,world)
    end
end

function simulation(world::Array{Patch,1}, timesteps::Int=1000)
    for t in 1:timesteps
        checkviability!(world)
        establish!(world)
        age!(world)
        grow!(world)
        disperse!(world)
        compete!(world)
        reproduce!(world)
        #        mod(timesteps,20) == 0 && analysis(world)
        #        mod(timesteps,20) == 0 && visualisation(world,t==20)
        #mod(timesteps,10) == 0 && println(t)
    end
end

function runit(firstrun::Bool)
    length(ARGS) > 0 ? seed = parse(Int,ARGS[1]) : seed = 0
    length(ARGS) > 1 ? mapfiles = ARGS[2:end] : mapfiles = ["mapfile"]   # length(ARGS) > 2 ? timesteps = parse(Int,ARGS[3]) : timesteps = 1000
    seed != 0 && srand(seed)
    if firstrun
        world=createworld([["1","1","1"]])
        simulation(world, 10)
    else
        for i in 1:length(mapfiles)
            timesteps,maptable = readmapfile(mapfiles[i])
            i == 1 && (world = createworld(maptable))
            i > 1 && updateworld!(world,maptable)
            simulation(world, timesteps)
            writedata(world, seed, mapfiles[i])
        end
    end
end

@time runit(true)
@time runit(false)

