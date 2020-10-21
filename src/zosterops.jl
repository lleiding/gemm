# All functions specific to the Zosterops experiments.
# Initialise predefined species and redefine dispersal and reproduction

# Note: this experiment uses above-ground carbon (AGC) as the second environmental
# niche instead of precipitation. However, we cannot rename all references to
# `prec*` throughout the source code, so we'll just leave that as it is.


"""
    zdisperse!(world, settings)

Dispersal of bird individuals within the world. Males disperse first, looking
for suitable habitats within their dispersal range to establish territories.
Females disperse second, looking for available mates. (Cf. Aben et al. 2016)
"""
function zdisperse!(world::Array{Patch,1}, settings::Dict{String, Any}, sex::Sex=male)
    for patch in world
        newseedbank = Array{Individual,1}
        while length(patch.seedbank) > 0
            juvenile = pop!(patch.seedbank)
            if juvenile.sex == sex
                zdisperse!(juvenile, world, patch.location)
            else
                push!(newseedbank, juvenile)
            end
        end
        patch.seedbank = newseedbank
    end
    (sex == male) && zdisperse!(world, settings, female)
end

"""
    zdisperse!(bird, world, location)

Dispersal of a single bird.
"""
function zdisperse!(bird::Individual, world::Array{Patch,1}, location::Tuple{Int, Int})
    # calculate max dispersal distance
    dispmean = patch.seedbank[idx].traits["dispmean"]
    dispshape = patch.seedbank[idx].traits["dispshape"]
    #FIXME we're probably going to need another distribution
    maxdist = rand(Logistic(dispmean,dispshape)) #XXX disperse!() add `/sqrt(2)`??
    # for each step, calculate the best habitat patch in the surrounding (excluding the current)
    x, y = location
    route = [location] # keep track of where we've been
    while maxdist > 0
        target = [(x-1,y-1), (x, y-1), (x+1, y-1),
                  (x-1, y), (x+1, y),
                  (x-1, y+1), (x, y+1), (x+1, y+1)]
        possdest = findall(p -> in(p.location, target) && !in(p.location, route),  world)
        bestdest = possdest[1]
        for patch in possdest[2:end]
            if abs(patch.prec - bird.traits["precopt"]) < abs(bestdest.prec - bird.traits["precopt"])
                bestdest = patch
            end
        end
        # move there - if it's occupied, repeat, otherwise stay there
        if true #TODO if the patch still has space, or an available mate
            push!(bestdest.community, bird)
            #TODO find mate
        end
        x, y = bestdest.location
        push!(route, bestdest.location)
        maxdist -= 1
    end
    #if the max dispersal distance is reached, the individual simply dies
end

"""
    zreproduce!(patch, settings)

Reproduction of Zosterops breeding pairs in a patch.
"""
function zreproduce!(patch::Patch}, settings::Dict{String, Any})
    noffs = rand(Poisson(settings["fertility"])) #TODO change to metabolic/species-specific?
    for bird in patch.community
        if bird.sex == female && bird.partner != 0 && bird.size > bird.traits["repsize"]
            partner = missing
            for b in patch.community
                (b.id == bird.partner) && partner = b
            end
            (ismissing(partner)) && continue
            append!(patch.seedbank, createoffspring(noffs, bird, partner,
                                                    settings["traitnames"], true))
        end
    end
end

