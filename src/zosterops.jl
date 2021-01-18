# All functions specific to the Zosterops experiments.
# Initialise predefined species and redefine dispersal and reproduction

# IMPORTANT NOTES (when using `zosterops` mode)
#
# 1. the second environmental niche is now `above-ground carbon` instead of
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

    Initialise the predefined Zosterops species archetypes (silvanus/highland, jubaensis/lowland).
    """
    function initzosteropsspecies(settings::Dict{String, Any})
        # ensure that "one gene, one trait" is true and that we have species definitions
        (settings["degpleiotropy"] != 0) && simlog("degpleiotropy must be 0", 'e')
        (isnothing(settings["species"])) && simlog("no species defined", 'e')
        # load per-species AGC optimum and tolerance values from settings
        for (s, v) in pairs(settings["species"])
            push!(zosterops, initzosteropsspecies(s, Float64(v[1]), Float64(v[2]), settings))
        end
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
            for gene in chromosome.genes
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
        getzosteropsspecies(name, sex, settings)

    Return a new individual of the named species
    """
    global function getzosteropsspecies(name::String, sex::Sex, settings::Dict{String, Any})
        # This function has to be global to escape the let-block of the species list
        isempty(zosterops) && initzosteropsspecies(settings)
        bird = nothing
        for z in zosterops
            if z.lineage == name
                bird = deepcopy(z)
                break
            end
        end
        (isnothing(bird)) && simlog("Unknown species name: $name", 'e')
        bird.id = rand(Int32)
        varyalleles!(bird.genome, settings, rand()) #XXX do we want mutations?
        bird.traits = gettraitdict(bird.genome, settings["traitnames"])
        bird.size = bird.traits["repsize"] #if we want mutations, we need to fix the size afterwards
        bird.sex = sex
        return bird
    end

    """
        getzosteropsnames(settings)

    Return the names of all defined Zosterops species.
    """
    #XXX utility function - but do I need it?
    global function getzosteropsnames(settings::Dict{String, Any})
        isempty(zosterops) && initzosteropsspecies(settings)
        return map(s->s.lineage, zosterops)
    end
end

"""
    zgenesis(patch, settings)

Create a new community of Zosterops breeding pairs (possibly of multiple species).
Returns an array of individuals.
"""
function zgenesis(patch::Patch, settings::Dict{String, Any})
    community = Array{Individual, 1}()
    (isnothing(settings["species"])) && simlog("no species defined", 'e')
    # check which species can inhabit this patch
    species = Array{String, 1}()
    for s in settings["species"]
        # s[1] = species name, s[2][1] = AGC opt, s[2][2] = AGC tol
        if abs(s[2][1]-patch.prec) <= s[2][2]
            push!(species, s[1])
        end
    end
    (isempty(species)) && return community
    # calculate the number of initial breeding pairs and add a male and a female for each
    npairs = Integer(rand(0:round(settings["cellsize"]/2)))
    for i in 1:npairs
        sp = rand(species)
        m = getzosteropsspecies(sp, male, settings)
        f = getzosteropsspecies(sp, female, settings)
        f.partner = m.id
        m.partner = f.id
        push!(community, m)
        push!(community, f)
        simlog("Adding a pair of Z. $sp", 'd')
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
            #XXX This is probably rather inefficient (creating a new list each time)
            juvenile = pop!(patch.seedbank)
            if juvenile.sex == sex
                zdisperse!(juvenile, world, patch.location,
                           Integer(settings["cellsize"]), settings["tolerance"])
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
    # keep track of where we've been and calculate the max dispersal distance
    x, y = location
    route = [location]
    #XXX disperse!() adds `/sqrt(2)`??
    #FIXME we're probably going to need another distribution
    maxdist = rand(Logistic(bird.traits["dispmean"], bird.traits["dispshape"]))
    while maxdist > 0
        # calculate the best habitat patch in the surroundings (i.e. the closest to AGC optimum)
        target = [(x-1, y-1), (x, y-1), (x+1, y-1),
                  (x-1, y),             (x+1, y),
                  (x-1, y+1), (x, y+1), (x+1, y+1)]
        filter!(c -> !in(c, route), target)
        possdest = map(p -> coordinate(p[1], p[2], world), target)
        filter!(p->!isnothing(p), possdest)
        if iszero(length(possdest))
            simlog("A Z.$(bird.lineage) died after failed dispersal.", 'd')
            return
        end
        bestdest = possdest[1]
        bestfit = abs(bestdest.prec - bird.traits["precopt"])
        for patch in possdest[2:end]
            patchfit = abs(patch.prec - bird.traits["precopt"])
            if patchfit < bestfit
                bestdest = patch
                bestfit = patchfit
            end
        end
        # check if the patch is full, within the bird's AGC range, and (for females) has a mate
        if bird.sex == female
            partner = findfirst(b -> ziscompatible(bird, b, tolerance), bestdest.community)
        end
        if (bestfit <= bird.traits["prectol"] &&
            length(bestdest.community) < cellsize &&
            (bird.sex == male || !isnothing(partner)))
            bird.marked = true
            push!(bestdest.community, bird)
            if bird.sex == female
                bird.partner = partner.id
                partner.partner = bird.id
            end
            simlog("A Z.$(bird.lineage) moved to $(bestdest.location[1])/$(bestdest.location[2]).", 'd')
            return #if we've found a spot, we're done
        end
        x, y = bestdest.location
        push!(route, bestdest.location)
        maxdist -= 1
    end #if the max dispersal distance is reached, the individual simply dies
    simlog("A Z.$(bird.lineage) died after failed dispersal.", 'd')
end

"""
    ziscompatible(female, male, tolerance)

Check to see whether two birds are reproductively compatible.
"""
function ziscompatible(f::Individual, m::Individual, tolerance::Float64)
    !(m.sex == male && f.sex == female) && return false
    !(m.size >= m.traits["repsize"] && f.size >= f.traits["repsize"]) && return false
    !(m.partner == 0 && f.partner == 0) && return false
    (m.lineage != f.lineage && rand(Float64) > tolerance) && return false
    #TODO how about seqsimilarity and genetic compatibility?
    return true
end
    
"""
    zreproduce!(patch, settings)

Reproduction of Zosterops breeding pairs in a patch.
"""
function zreproduce!(patch::Patch, settings::Dict{String, Any})
    noffs = Integer(settings["fertility"])
    for bird in patch.community
        if bird.sex == female && bird.partner != 0
            pt = findfirst(b -> b.id == bird.partner, patch.community)
            if isnothing(pt)
                simlog("A Z.$(bird.lineage) has no partner.", 'd')
                continue
            end
            partner = patch.community[pt]
            simlog("A Z.$(bird.lineage) mated with a Z.$(partner.lineage).")
            #TODO Offspring are assigned the lineage of their mother. Is that what we want?
            append!(patch.seedbank, createoffspring(noffs, bird, partner,
                                                    settings["traitnames"], true))
        end
    end
end


let width = 0, height = 0
    """
        coordinate(x, y, world)

    A utility function to perform a fast look-up for the patch at coordinate x/y.
    Important: this assumes a rectangular world with coordinates in row-major order!
    """
    global function coordinate(x::Int, y::Int, world::Array{Patch,1})
        if iszero(width)
            width = maximum(p -> p.location[1], world)
            height = maximum(p -> p.location[2], world)
        end
        (x <= 0 || y <= 0 || x > width || y > height) && return
        i = ((y-1) * width) + x
        return world[i]
    end
end
