# All functions specific to the Zosterops experiments.
# Initialise predefined species and redefine dispersal and reproduction

# IMPORTANT NOTES (when using `zosterops` mode)
#
# 1. the second environment niche is now `above-ground carbon` instead of
#    `precipitation` (unfortunately, we can't actually rename it everywhere)
# 2. the `cellsize` setting now determines the patch carrying capacity in
#    individuals, not grams
# 3. the `fertility` setting is the absolute number of offspring per breeding pair
#    (instead of a metabolic coefficient)
# 4. the `tolerance` setting now determines the probability that a mate of
#    another species is accepted
# 5. `degpleiotropy` must be set to 0, otherwise the species initialisation will fail

let zosterops = Individual[]
    """
        initzosteropsspecies(settings)

    Initialise the two relevant Zosterops species archetypes (silvanus/highland, jubaensis/lowland).
    """
    function initzosteropsspecies(settings::Dict{String, Any})
        if settings["degpleiotropy"] != 0 # ensure that "one gene, one trait" is true
            simlog("degpleiotropy must be 0", settings, 'e')
        end
        #XXX magic numbers!
        silvanus = initzosteropsspecies("silvanus", 90, 10, settings) #TODO adjust values
        jubaensis = initzosteropsspecies("jubaensis", 50, 40, settings) #TODO adjust values
        push!(zosterops, silvanus)
        push!(zosterops, jubaensis)
    end

    """
        initzosteropsspecies(name, precopt, prectol, settings)

    Create a new individual then modify its traits to create a Zosterops archetype.
    """
    function initzosteropsspecies(name::String, precopt::Float64, prectol::Float64, settings::Dict{String, Any})
        archetype = createind(settings)
        archetype.lineage = name
        # Find the gene that codes for the relevant trait and change that trait's value
        for chromosome in archetype.genome
            for gene in chromosome
                isempty(gene.codes) && continue
                genetrait = settings["traitnames"][gene.codes[1].nameindex]
                #XXX do we also need to set tempopt/temptol?
                if genetrait == "precopt"
                    gene.codes[1].value = precopt
                elseif genetrait == "prectol"
                    gene.codes[1].value = prectol
                end
            end
        end
        # then recalculate the individual's trait dict and return the finished archetype
        archetype.traits = gettraitdict(archetype.genome, settings["traitnames"])
        return archetype
    end
    
    """
        getzosteropsspecies(name)

    Return a new individual of the named species ("silvanus" or "jubaensis")
    """
    function getzosteropsspecies(name::String, sex::Sex, settings::Dict{String, Any})
        isempty(zosterops) && initzosteropsspecies(settings)
        bird = nothing
        if name == "silvanus"
            bird = deepcopy(zosterops[1])
        elseif name == "jubaensis"
            bird = deepcopy(zosterops[2])
        else
            simlog("Unknown species name: $name", settings, 'e')
        end
        bird.id = rand(Int32)
        varyalleles!(bird.genome, settings, rand())
        bird.traits = gettraitdict(bird.genome, settings["traitnames"])
        bird.sex = sex
        return bird
    end
end

"""
    genesis(settings)

Create a new community of Zosterops breeding pairs (possibly of multiple species).
Returns an array of individuals.
"""
function zgenesis(patch::Patch, settings::Dict{String, Any})
    community = Individual[]
    npairs = rand(0:round(settings["cellsize"]/2))
    (npairs == 0) && return
    for i in 1:npairs
        #TODO
    end
    community
end

"""
    zdisperse!(world, settings)

Dispersal of bird individuals within the world. Males disperse first, looking
for suitable habitats within their dispersal range to establish territories.
Females disperse second, looking for available mates. (Cf. Aben et al. 2016)
"""
function zdisperse!(world::Array{Patch,1}, settings::Dict{String, Any}, sex::Sex=male)
    for patch in world
        newseedbank = Array{Individual,1}()
        while length(patch.seedbank) > 0
            juvenile = pop!(patch.seedbank)
            if juvenile.sex == sex
                zdisperse!(juvenile, world, patch.location,
                           settings["cellsize"], settings["tolerance"])
            else
                push!(newseedbank, juvenile)
            end
        end
        patch.seedbank = newseedbank
    end
    (sex == male) && zdisperse!(world, settings, female)
end

"""
    zdisperse!(bird, world, location, cellsize)

Dispersal of a single bird.
"""
function zdisperse!(bird::Individual, world::Array{Patch,1}, location::Tuple{Int, Int},
                    cellsize::Int, tolerance::Float64)
    # calculate max dispersal distance
    dispmean = patch.seedbank[idx].traits["dispmean"]
    dispshape = patch.seedbank[idx].traits["dispshape"]
    #FIXME we're probably going to need another distribution
    maxdist = rand(Logistic(dispmean,dispshape)) #XXX disperse!() adds `/sqrt(2)`??
    # for each step, calculate the best habitat patch in the surroundings
    x, y = location
    route = [location] # keep track of where we've been
    while maxdist > 0
        target = [(x-1, y-1), (x, y-1), (x+1, y-1),
                  (x-1, y),             (x+1, y),
                  (x-1, y+1), (x, y+1), (x+1, y+1)]
        possdest = findall(p -> in(p.location, target) && !in(p.location, route),  world)
        bestdest = possdest[1]
        for patch in possdest[2:end]
            if abs(patch.prec - bird.traits["precopt"]) < abs(bestdest.prec - bird.traits["precopt"])
                bestdest = patch
            end
        end
        # move there - if it's occupied, repeat, otherwise stay there
        if bird.sex == female #XXX expensive to do here?
            partner = findfirst(b -> ziscompatible(bird, b, tolerance), bestdest.community)
        end
        if (abs(bestdest.prec - bird.traits["precopt"]) <= bird.traits["prectol"] &&
            length(bestdest.community) < cellsize &&
            (bird.sex == male || !isnothing(partner)))
            push!(bestdest.community, bird)
            if bird.sex == female
                bird.partner = partner.id
                partner.partner = bird.id
            end
            return #if we've found a spot, we're done
        end
        x, y = bestdest.location
        push!(route, bestdest.location)
        maxdist -= 1
    end #if the max dispersal distance is reached, the individual simply dies
end

"""
    iscompatible(female, male, tolerance)

Check to see whether two birds are reproductively compatible.
"""
function ziscompatible(f::Individual, m::Individual, tolerance::Float64)
    !(m.sex == male && f.sex == female) && return false
    !(m.size >= m.traits["repsize"] && f.size >= f.traits["repsize"]) && return false
    !(m.partner == 0 && f.partner == 0) && return false
    (m.lineage != f.lineage && rand(Float64) > tolerance) && return false
    #XXX genetic compatibility?
    return true
end
    
"""
    zreproduce!(patch, settings)

Reproduction of Zosterops breeding pairs in a patch.
"""
function zreproduce!(patch::Patch, settings::Dict{String, Any})
    noffs = settings["fertility"]
    for bird in patch.community
        if bird.sex == female && bird.partner != 0
            partner = findfirst(b -> b.id == bird.partner, patch.community)
            (isnothing(partner)) && continue
            append!(patch.seedbank, createoffspring(noffs, bird, partner,
                                                    settings["traitnames"], true))
        end
    end
end

