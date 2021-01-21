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

let zosterops = Individual[] #holds the species archetypes
    """
        initzosteropsspecies()

    Initialise the predefined Zosterops species archetypes (silvanus/highland, jubaensis/lowland).
    """
    function initzosteropsspecies()
        # ensure that "one gene, one trait" is true and that we have species definitions
        (setting("degpleiotropy") != 0) && simlog("degpleiotropy must be 0", 'e')
        (isnothing(setting("species"))) && simlog("no species defined", 'e')
        # load per-species AGC optimum and tolerance values from settings
        for (s, v) in pairs(setting("species"))
            push!(zosterops, initzosteropsspecies(s, Float64(v[1]), Float64(v[2])))
        end
    end

    """
        initzosteropsspecies(name, precopt, prectol)

    Create a new individual then modify its traits to create a Zosterops archetype.
    """
    function initzosteropsspecies(name::String, precopt::Float64, prectol::Float64)
        archetype = createind()
        archetype.lineage = name
        # Find the gene that codes for relevant traits and change their values
        for chromosome in archetype.genome
            for gene in chromosome.genes
                isempty(gene.codes) && continue
                genetrait = setting("traitnames")[gene.codes[1].nameindex]
                #XXX do we also need to set tempopt/temptol?
                if genetrait == "precopt"
                    gene.codes[1].value = precopt
                elseif genetrait == "prectol"
                    gene.codes[1].value = prectol
                end
            end
        end
        # then recalculate the individual's trait dict and return the finished archetype
        archetype.traits = gettraitdict(archetype.genome, setting("traitnames"))
        return archetype
    end
    
    """
        getzosteropsspecies(name, sex)

    Return a new individual of the named species
    """
    global function getzosteropsspecies(name::String, sex::Sex)
        # This function has to be global to escape the let-block of the species list
        isempty(zosterops) && initzosteropsspecies()
        bird = nothing
        for z in zosterops
            if z.lineage == name
                bird = deepcopy(z)
                break
            end
        end
        (isnothing(bird)) && simlog("Unknown species name: $name", 'e')
        bird.id = rand(UInt32)
        varyalleles!(bird.genome, rand())
        bird.traits = gettraitdict(bird.genome, setting("traitnames"))
        bird.size = bird.traits["repsize"] # we need to fix the size after again after mutation
        bird.sex = sex
        return bird
    end
end

"""
    zgenesis(patch)

Create a new community of Zosterops breeding pairs (possibly of multiple species).
Returns an array of individuals.
"""
function zgenesis(patch::Patch)
    community = Array{Individual, 1}()
    (isnothing(setting("species"))) && simlog("no species defined", 'e')
    # check which species can inhabit this patch
    species = Array{String, 1}()
    for s in setting("species")
        # s[1] = species name, s[2][1] = AGC opt, s[2][2] = AGC tol
        if abs(s[2][1]-patch.prec) <= s[2][2]
            push!(species, s[1])
        end
    end
    (isempty(species)) && return community
    # calculate the number of initial breeding pairs and add a male and a female for each
    npairs = Integer(rand(0:round(setting("cellsize")/2)))
    for i in 1:npairs
        sp = rand(species)
        m = getzosteropsspecies(sp, male)
        f = getzosteropsspecies(sp, female)
        f.partner = m.id
        m.partner = f.id
        push!(community, m)
        push!(community, f)
        simlog("Adding a pair of Z. $sp", 'd')
    end
    community
end
    
"""
    zdisperse!(world)

Disperse all juvenile individuals within the world. 
"""
function zdisperse!(world::Array{Patch,1})
    for patch in world
        for ind in patch.seedbank
            zdisperse!(ind, world, patch.location)
        end
        patch.seedbank = Array{Individual,1}()
    end
end

"""
    zdisperse!(bird, world, location)

Dispersal of a single bird. Birds look patches with a suitable
habitat and a free territory or available mate. (Cf. Aben et al. 2016)
"""
function zdisperse!(bird::Individual, world::Array{Patch,1}, location::Tuple{Int, Int})
    # keep track of where we've been and calculate the max dispersal distance
    x, y = location
    route = [location]
    #XXX disperse!() adds `/sqrt(2)`??
    #XXX sex-biased maximum dispersal?
    maxdist = rand(Logistic(bird.traits["dispmean"], bird.traits["dispshape"]))
    while maxdist > 0
        # calculate the best habitat patch in the surroundings (i.e. the closest to AGC optimum)
        target = [(x-1, y-1), (x, y-1), (x+1, y-1),
                  (x-1, y),             (x+1, y),
                  (x-1, y+1), (x, y+1), (x+1, y+1)]
        filter!(c -> !in(c, route), target)
        possdest = map(p -> coordinate(p[1], p[2], world), target)
        filter!(p->!isnothing(p), possdest)
        iszero(length(possdest)) && @goto failure
        bestdest = possdest[1]
        bestfit = abs(bestdest.prec - bird.traits["precopt"])
        for patch in possdest[2:end]
            patchfit = abs(patch.prec - bird.traits["precopt"])
            if patchfit < bestfit
                bestdest, bestfit = patch, patchfit
            end
        end
        x, y = bestdest.location
        # check if the patch is within the bird's AGC range and has free space
        if (bestfit <= bird.traits["prectol"] && length(bestdest.community) < setting("cellsize"))
            # look for an available mate, starting with conspecifics
            conspecifics = findall(b -> b.lineage == bird.lineage, bestdest.community)
            partner = findfirst(b -> ziscompatible(bird, b), bestdest.community[conspecifics])
            if isnothing(partner) && setting("speciation") == "off"
                congenerics = setdiff(eachindex(bestdest.community), conspecifics)
                partner = findfirst(b -> ziscompatible(bird, b), bestdest.community[congenerics])
            end
            if !isnothing(partner)
                # if we've found a partner
                partner = bestdest.community[partner]
                bird.partner = partner.id
                partner.partner = bird.id
                @goto success
            elseif length(bestdest.community) <= (setting("cellsize")-2)
                # if there's space for another breeding pair
                @label success
                bird.marked = true
                push!(bestdest.community, bird)
                simlog("$(idstring(bird)) moved to $x/$y.", 'd')
                return #if we've found a spot, we're done
            end
        end
        push!(route, bestdest.location)
        maxdist -= 1
    end #if the max dispersal distance is reached, the individual simply dies
    @label failure
    simlog("A Z.$(bird.lineage) died after failed dispersal.", 'd')
end

"""
    ziscompatible(individual1, individual2)

Check to see whether two birds are reproductively compatible.
"""
function ziscompatible(i1::Individual, i2::Individual)
    (i1.sex == i2.sex) && return false
    !(i1.partner == 0 && i2.partner == 0) && return false
    if setting("speciation") == "off" #check for hybridisation
        (i1.lineage != i2.lineage && rand(Float64) > setting("tolerance")) && return false
    else #check for speciation
        (!iscompatible(i1, i2)) && return false
    end
    simlog("Found a partner: $(idstring(i1)) and $(idstring(i2)).", 'd')
    return true
end
    
"""
    zreproduce!(patch)

Reproduction of Zosterops breeding pairs in a patch.
"""
function zreproduce!(patch::Patch)
    #FIXME reproduction happens too seldom?
    #hypothesis: new partners don't reproduce
    #FIXME partner matching is not constant
    noffs = Integer(setting("fertility"))
    for bird in patch.community
        if bird.sex == female && bird.partner != 0
            pt = findfirst(b -> b.id == bird.partner, patch.community)
            if isnothing(pt)
                simlog("$(idstring(bird)) has no partner.", 'd')
                continue
            end
            partner = patch.community[pt]
            simlog("$(idstring(bird)) mated with $(idstring(partner)), $noffs offspring.", 'd')
            #TODO Offspring are assigned the lineage of their mother. Is that what we want?
            append!(patch.seedbank, createoffspring(noffs, bird, partner,
                                                    setting("traitnames"), true))
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

"""
    idstring(individual)

A small utility function to print a string identifier for a given individual.
"""
function idstring(bird::Individual)
    return "Z."*bird.lineage*" "*string(bird.id)
end
