# All processes related to survival:
# - growth
# - density-independent mortality
# - competition

"""
    survive!(patch, mortality)

Density independent survival of individuals in a patch. The actual mortality
probability is calculated with a metabolic formula, modified by the passed `mortality`
variable and an individual's temperature adaptation.
"""
function survive!(patch::Patch, mortality::Float64)
    temp = patch.temp
    idx = 1
    while idx <= size(patch.community,1)
        if !patch.community[idx].marked
            mass = patch.community[idx].size
            deathrate = mortality * mass^(-1/4) * exp(-act/(boltz*temp))
            dieprob = (1 - exp(-deathrate))
            if rand() * patch.community[idx].tempadaptation < dieprob
                splice!(patch.community, idx)
                continue
            end
        end
        idx += 1
    end
end

"""
    survive!(world, settings)

World-wide mortality. Sounds apocalyptic, but is just a fact of life.
"""
function survive!(world::Array{Patch,1}, settings::Dict{String, Any})
    for patch in world
        (patch.isisland || !settings["static"]) && survive!(patch, settings["mortality"]) # pmap(!,patch) ???
    end
end

"""
    grow!(patch, growthrate, capgrowth)

Growth of individuals in the given patch. The actual growthrate is calculated
with a metabolic formula, modified by the passed `growthrate` variable.
"""
function grow!(patch::Patch, growthrate::Float64, capgrowth::Bool)
    temp = patch.temp
    idx = 1
    while idx <= size(patch.community,1)
        if !patch.community[idx].marked
            repsize = patch.community[idx].traits["repsize"]
            mass = patch.community[idx].size
            if mass <= repsize # stop growth if reached repsize
                growth = growthrate * mass^(3/4) * exp(-act/(boltz*temp))
                newmass = mass + growth
                if capgrowth && newmass > patch.community[idx].traits["repsize"]
                    newmass = patch.community[idx].traits["repsize"]
                end
                if newmass > 0 && mass > 0
                    patch.community[idx].size = newmass
                else
                    splice!(patch.community, idx)
                    continue
                end
            end
        end
        idx += 1
    end
end

"""
    grow!(world, settings)

Carry out growth for all patches.
"""
function grow!(world::Array{Patch,1}, settings::Dict{String, Any})
    for patch in world
        # pmap(!,patch) ???
        if patch.isisland || !settings["static"]
            grow!(patch, settings["growthrate"], settings["capgrowth"])
        end
    end
end

"""
    compete!(patch)

While there are too many organisms in a patch, pick two at random and kill the
one that is less adapted to the local precipitation levels.
"""
function compete!(patch::Patch)
    totalmass = sum(map(x -> x.size, patch.community))
    while totalmass >= patch.capacity # occupied capacity larger than available
        firstind, secondind = rand(eachindex(patch.community), 2)
        firstind == secondind && length(eachindex(patch.community)) > 1 && continue
        if patch.community[firstind].precadaptation < patch.community[secondind].precadaptation
            totalmass -= patch.community[firstind].size
            splice!(patch.community, firstind) # profiling: expensive!
        elseif patch.community[firstind].precadaptation > patch.community[secondind].precadaptation
            totalmass -= patch.community[secondind].size
            splice!(patch.community, secondind) # profiling: expensive!
        else
            victim = rand([firstind, secondind])
            totalmass -= patch.community[victim].size
            splice!(patch.community, victim) # profiling: expensive!
        end
    end
end

"""
    compete!(world, static)

Carry out competition on all patches.
"""
function compete!(world::Array{Patch,1}, static::Bool = true)
    for patch in world
        (patch.isisland || !static) && compete!(patch) # pmap(!,patch) ???
    end
end
