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
        for juvenile in patch.seedbank
            (juvenile.sex == sex) && zdisperseindividual!(juvenile, world, patch.location)
        end
    end
    (sex == male) && zdisperse!(world, settings, female)
end

function zdisperseindividual!(bird::Individual, world::Array{Patch,1}, location::Tuple{Int, Int})
    #TODO calculate max dispersal distance
    #TODO for each step, calculate the best habitat patch in the surrounding (excluding the current)
    #TODO move there - if it's occupied, repeat, otherwise stay there
    #TODO if you've reached the max dispersal distance, die
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

