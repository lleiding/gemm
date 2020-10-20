# All functions needed for the invasion experiments.
# Create the global species pool, carry out species introductions, and disturbance.

let speciespool = Individual[]
    # TODO: `speciespool` should be part of the `world` object holding all grid cells
    """
        initspeciespool!(settings)

    Initialise the foreign species pool (i.e. an array of random individuals that
    are used as the propagule source for invasion events.
    """
    function initspeciespool!(settings::Dict{String, Any})
        #TODO rename initglobalspeciespool!
        for i in 1:settings["global-species-pool"]
            push!(speciespool, createind(settings))
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
        invade!(world, settings)

    Introduce non-native species from the species pool to all patches marked
    as invasible.
    """
    global function invade!(world::Array{Patch,1}, settings::Dict{String, Any})
        # This function has to be global to escape the let-block of the species pool
        (settings["propagule-pressure"] == 0 || settings["global-species-pool"] == 0) && return
        (length(speciespool) == 0) && initspeciespool!(settings)
        for patch in world
            patch.invasible && invade!(patch, settings["propagule-pressure"])
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
    disturb!(world, settings)

Disturb all patches in the world.
"""
function disturb!(world::Array{Patch,1}, settings::Dict{String, Any})
    (settings["disturbance"] == 0) && return
    if settings["disturbance"] > 100
        simlog("disturbance must be no more than 100%", 'w')
        settings["disturbance"] = 100
    end
    for patch in world
        (patch.isisland || !settings["static"]) && disturb!(patch, settings["disturbance"])
    end
end

"""
    createind(settings, marked=false)

Create an individual organism of a new species with a random genome.
"""
function createind(settings::Dict{String, Any}, marked::Bool = false)
    id = rand(Int32)
    parentid = rand(Int32)
    lineage = randstring(4)
    ngenes = settings["maxloci"] * length(settings["traitnames"])
    ngenes < 1 && (ngenes = 1)
    traits = createtraits(settings)
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
    locivar = rand() #XXX Is this legit? Don't quite understand `locivar` yet...
    varyalleles!(chromosomes, settings, locivar)
    traitdict = gettraitdict(chromosomes, settings["traitnames"])
    if settings["indsize"] == "adult"
        indsize = traitdict["repsize"]
    elseif settings["indsize"] == "seed"
        indsize = traitdict["seedsize"]
    else
        indsize = traitdict["seedsize"] + rand() * traitdict["repsize"] # XXX: sizes shouldn't be uniformally dist'd
    end
    Individual(lineage, chromosomes, traitdict, marked, 1.0, 1.0, indsize, hermaphrodite, 0, id)
end
