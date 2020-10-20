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
            (juvenile.sex == sex) && zdisperseindividual!(juvenile, patch.location)
        end
    end
    (sex == male) && zdisperse!(world, settings, female)
end

function zdisperseindividual!(bird::Individual, location::Tuple{Int, Int})
    #TODO
end

function zreproduce!(world::Array{Patch,1}, settings::Dict{String, Any})
    #TODO
end

