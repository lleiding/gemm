# All functions needed for the invasion experiments.
# Create the global species pool, carry out species introductions, and disturbance.

let speciespool = Individual[]
    #XXX should `speciespool` be part of the `world` object holding all grid cells?
    """
        initglobalspeciespool!()

    Initialise the foreign species pool (i.e. an array of random individuals that
    are used as the propagule source for invasion events.
    """
    function initglobalspeciespool!()
        for i in 1:setting("global-species-pool")
            push!(speciespool, createind())
        end
    end

    """
        invade!(patch, pressure)

    Select a given amount of individuals from the global species pool and add
    them to the patch.
    """
    function invade!(patch::Patch, pressure::Int)
        invaders = deepcopy(rand(speciespool, pressure))
        append!(patch.community, invaders)
    end

    """
        invade!(world)

    Introduce non-native species from the species pool to all patches marked
    as invasible.
    """
    global function invade!(world::Array{Patch,1})
        # This function has to be global to escape the let-block of the species pool
        (setting("propagule-pressure") == 0 || setting("global-species-pool") == 0) && return
        (length(speciespool) == 0) && initglobalspeciespool!()
        for patch in world
            patch.invasible && invade!(patch, setting("propagule-pressure"))
        end
    end
end


"""
    disturb!(patch, intensity)

Species-independent mortality due to disturbance on the given patch. The intensity
gives the mortality percentage.
"""
function disturb!(patch::Patch, intensity::Int)
    length(patch.community) <= 0 && return
    deaths = Integer(round(length(patch.community) * (intensity/100)))
    dead = unique(rand(1:length(patch.community), deaths))
    while length(dead) < deaths
        newdead = rand(1:length(patch.community), deaths-length(dead))
        dead = unique(vcat(dead, newdead))
    end
    sort!(dead)
    deleteat!(patch.community, dead)
end

"""
    disturb!(world)

Disturb all patches in the world.
"""
function disturb!(world::Array{Patch,1})
    (setting("disturbance") == 0) && return
    if setting("disturbance") > 100
        simlog("disturbance must be no more than 100%", 'w')
        updatesetting("disturbance", 100)
    end
    for patch in world
        (patch.isisland || !setting("static")) && disturb!(patch, setting("disturbance"))
    end
end
