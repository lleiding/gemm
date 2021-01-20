# All functions needed for the gradient/variable environment experiments.
# Randomly change habitat precipitation and temperature.


"""
    changetemp!(world, sdtemp)

Change the temperature of all patches according to a normal distribution.
"""
function changetemp!(world::Array{Patch,1}, sdtemp::Float64)
    sdtemp == 0 && return
    deltaval = rand(Normal(0.0, sdtemp))
    for patch in world
        patch.temp += deltaval
    end
    markthem!(world)
end

"""
    changeprec!(world, sdprec)

Change the precipitation of all patches according to a normal distribution.
"""
function changeprec!(world::Array{Patch,1}, sdprec::Float64)
    sdprec == 0 && return
    deltaval = rand(Normal(0.0, sdprec))
    for patch in world
        patch.prec += deltaval
    end
    markthem!(world)
end

"""
    changehabitat!(world)

Carry out 'global change' on all patches.
"""
function changehabitat!(world::Array{Patch,1})
    # TODO: record trajectory? input trajectory?
    changetemp!(world, setting("sdtemp"))
    changeprec!(world, setting("sdprec"))
end

"""
    markthem!(community)

Set each individual in the community (= array of individuals) as "marked".
"""
function markthem!(community::Array{Individual, 1})
    for ind in community
        ind.marked = true
    end
end

"""
    markthem!(habitat)

Set each individual in the given patch as "marked".
"""
function markthem!(habitat::Patch)
    markthem!(habitat.community)
end

"""
    markthem!(world)

Set every individual in the world as "marked".
"""
function markthem!(world::Array{Patch, 1})
    for habitat in world
        markthem!(habitat)
    end
end
