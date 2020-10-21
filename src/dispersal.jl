# All functions related to dispersal:
# - dispersal
# - establishment
# - viability checking

"""
    disperse!(world, static)

Dispersal of individuals within the world.
"""
function disperse!(world::Array{Patch,1}, static::Bool = true)
    # TODO: additional border conditions, refactor
    for patch in world
        idx = 1
        while idx <= size(patch.seedbank,1)
            dispmean = patch.seedbank[idx].traits["dispmean"]
            dispshape = patch.seedbank[idx].traits["dispshape"]
            xdir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # scaling so that geometric mean...
            ydir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # ...follows original distribution
            xdest = patch.location[1] + Int(round(xdir))
            ydest = patch.location[2] + Int(round(ydir))
            !patch.isisland ? target = checkborderconditions(world, xdest, ydest) : target = (xdest, ydest)
            possdest = findall(x -> in(x.location, [target]), world)
            static && filter!(x -> world[x].isisland, possdest) # disperse only to islands
            if !static || patch.isisland
                indleft = splice!(patch.seedbank,idx) # only remove individuals from islands!
                idx -= 1
            end
            if length(possdest) > 0 # if no viable target patch, individual dies
                if static && !patch.isisland
                    indleft = patch.seedbank[idx]
                end
                destination = rand(possdest) # currently there is only one possible destination
                push!(world[destination].community, indleft)
            end
            idx += 1
        end
    end
end

"""
    establish!(patch, nniches)

Establishment of individuals in patch `p`: Sets the adaptation parameters (~fitness)
according to an individual's adaptation to the niches of the surrounding environment.

A maximum of two niches (temperature and "precipitation") is currently supported.
"""
function establish!(patch::Patch, nniches::Int=1)
    temp = patch.temp
    idx = 1
    while idx <= size(patch.community,1)
        if patch.community[idx].marked
            opt = patch.community[idx].traits["tempopt"]
            tol = patch.community[idx].traits["temptol"]
            fitness = gausscurve(opt, tol, temp, 0.0)
            fitness > 1 && (fitness = 1) # should be obsolete
            fitness < 0 && (fitness = 0) # should be obsolete
            patch.community[idx].tempadaptation = fitness
            if nniches >= 2
                opt = patch.community[idx].traits["precopt"]
                tol = patch.community[idx].traits["prectol"]
                fitness = gausscurve(opt, tol, patch.prec, 0.0)
                fitness > 1 && (fitness = 1) # should be obsolete
                fitness < 0 && (fitness = 0) # should be obsolete
                patch.community[idx].precadaptation = fitness
            end
            patch.community[idx].marked = false
        end
        idx += 1
    end
end

"""
    establish!(world, nniches, static)

Carry out establishment for each patch in the world.
"""
function establish!(world::Array{Patch,1}, nniches::Int=1, static::Bool = true)
    for patch in world
        (patch.isisland || !static) && establish!(patch, nniches) # pmap(!,patch) ???
    end
end

"""
    checkviability!(community, settings)

Check whether all individuals in the passed community conform to a basic set of
constraints (i.e. all traits are present and certain properties are >= 0).
Individuals that fail the test are removed from the community.
"""
function checkviability!(community::Array{Individual, 1}, settings::Dict{String, Any})
    idx=1
    while idx <= size(community,1)
        reason = ""
        dead = false
        community[idx].size <= 0 && (dead = true) && (reason *= "size ")
        any(collect(values(community[idx].traits)) .< 0) && (dead = true) && (reason *= "traitvalues ")
        community[idx].traits["repsize"] <= community[idx].traits["seedsize"] && (dead = true) && (reason *= "seed/rep ")
        community[idx].tempadaptation < 0 && (dead = true) && (reason *= "fitness ")
        community[idx].precadaptation < 0 && (dead = true) && (reason *= "fitness ")
        community[idx].traits["selfing"] > 1 && (dead = true) && (reason *= "selfing ")
        community[idx].traits["seqsimilarity"] > 1 && (dead = true) && (reason *= "seqsimilarity ")
        !traitsexist(community[idx].traits, settings) && (dead = true) && (reason *= "missingtrait ")
        if dead
            simlog("Individual not viable: $reason. Being killed.", settings, 'w')
            splice!(community,idx)
            continue
        end
        idx += 1
    end
end

"""
    checkviability(world, settings)

Check the viability of all individuals.
"""
function checkviability!(world::Array{Patch,1}, settings::Dict{String, Any})
    for patch in world
        #XXX pmap(checkviability!,patch) ???
        checkviability!(patch.community, settings)
    end
end


"""
    traitsexist(traits, settings)

Check a trait dict to make sure it contains the full set of traitnames required
by the model (as defined in the settings).
"""
function traitsexist(traits::Dict{String, Float64}, settings::Dict{String, Any})
    missingtraits = setdiff(settings["traitnames"], keys(traits))
    if length(missingtraits) > 0
        simlog("Missing trait $missingtraits. Individual might be killed.", settings, 'w')
        return false
    end
    true
end

"""
    traitsexist(individual, settings)

Make sure an individual organism has the full set of traits required by the model
(as defined in the settings).
"""
function traitsexist(ind::Individual, settings::Dict{String, Any})
    traitnames = settings["traitnames"]
    for trait in traitnames
        if !haskey(ind.traits, trait)
            simlog("Individual is missing trait $trait. Might be killed.", settings, 'e')
            return false
        end
    end
    true
end


"""
    gausscurve(b, c, x, a=1.0)

Calculate the value of the Gauss function ("bell curve") at point x; with
a being the maximum height of the curve, b the position of the curve center and
c the standard deviation ("width").
"""
function gausscurve(b, c, x, a = 1.0)
    if c != 0 && a != 1.0
        a = 1 / (c * sqrt(2 * pi))
        y = a * exp(-(x-b)^2/(2*c^2))
    elseif c != 0
        y = a * exp(-(x-b)^2/(2*c^2))
    else
        y = 0.0
    end
end

"""
    findisland(w)

within world `w`, find out in which direction from the continent the island(s) lie(s).
"""
function findisland(world::Array{Patch,1})
    xmin = minimum(map(x->x.location[1],world))
    xmax = maximum(map(x->x.location[1],world))
    ymin = minimum(map(x->x.location[2],world))
    ymax = maximum(map(x->x.location[2],world))
    westernborder = filter(x->x.location[1]==xmin,world)
    northernborder = filter(x->x.location[2]==ymax,world)
    easternborder = filter(x->x.location[1]==xmax,world)
    southernborder = filter(x->x.location[2]==ymin,world)
    if all(map(x->x.isisland,westernborder))
        return "west"
    elseif all(map(x->x.isisland,northernborder))
        return "north"
    elseif all(map(x->x.isisland,easternborder))
        return "east"
    elseif all(map(x->x.isisland,southernborder))
        return "south"
    else
        return "none"
    end
end

"""
    checkborderconditions!(w, x, y)

check if coordinates `x` and `y` lie within world `w` and correct if not,
considering defined border conditions.
"""
function checkborderconditions(world::Array{Patch,1}, xdest::Int, ydest::Int)
    xmin = minimum(map(x->x.location[1],world))
    xmax = maximum(map(x->x.location[1],world))
    ymin = minimum(map(x->x.location[2],world))
    ymax = maximum(map(x->x.location[2],world))
    xrange = xmax - xmin + 1 # we're counting cells!
    yrange = ymax - ymin + 1 # we're counting cells!
    xshift = xdest - xmin + 1 # 1-based count of cells
    yshift = ydest - ymin + 1 # 1-based count of cells
    xshift > 0 ? outofx = abs(xshift) : outofx = abs(xshift) + 1
    while outofx > xrange
        outofx -= xrange
    end
    outofx -= 1
    yshift > 0 ? outofy = abs(yshift) : outofy = abs(yshift) + 1
    while outofy > yrange
        outofy -= yrange
    end
    outofy -= 1
    islanddirection = findisland(world::Array{Patch,1})
    if islanddirection == "west"
        xdest > xmax && (xdest = xmax - outofx) # east: reflective
        ydest < ymin && (ydest = ymax - outofy) # south: periodic
        ydest > ymax && (ydest = ymin + outofy) # north: periodic
    elseif islanddirection == "north"
        ydest < ymin && (ydest = ymin + outofy) # south: reflective
        xdest < xmin && (xdest = xmax - outofx) # west: periodic
        xdest > xmax && (xdest = xmin + outofx) # east: periodic
    elseif islanddirection == "east"
        xdest < xmin && (xdest = xmin + outofx) # west: reflective
        ydest < ymin && (ydest = ymax - outofy) # south: periodic
        ydest > ymax && (ydest = ymin + outofy) # north: periodic
    elseif islanddirection == "south"
        ydest > ymax && (ydest = ymax + outofy) # north: reflective
        xdest < xmin && (xdest = xmax - outofx) # west: periodic
        xdest > xmax && (xdest = xmin + outofx) # east: periodic
    else
        ydest > ymax && (ydest = ymin + outofy) # north: periodic
        xdest > xmax && (xdest = xmin + outofx) # east: periodic
        ydest < ymin && (ydest = ymax - outofy) # south: periodic
        xdest < xmin && (xdest = xmax - outofx) # west: periodic
    end
    xdest, ydest
end
