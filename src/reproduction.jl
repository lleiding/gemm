# All functions related to reproduction.

"""
    reproduce!(patch, settings)

Reproduction of individuals in a patch (default function).
"""
function reproduce!(patch::Patch, settings::Dict{String, Any}) #TODO: refactor!
    #TODO This is one of the most compute-intensive functions in the model - optimise?
    for ind in patch.community
        ind.marked && continue # individual might not have established yet
        ind.size < ind.traits["repsize"] && continue
        metaboffs = settings["fertility"] * ind.size^(-1/4) * exp(-act/(boltz*patch.temp))
        noffs = rand(Poisson(metaboffs))
        noffs < 1 && continue
        partners = findmate(patch.community, ind, settings)
        if length(partners) < 1 && rand() < ind.traits["selfing"]
            partners = [ind]
        elseif length(partners) < 1
            continue
        end
        numpartners = Integer(round(ind.traits["numpollen"]))
        if numpartners == 0
            #FIXME happens regularly!
            # `numpollen` not handled specifically in `createtraits()`?
            simlog("Individual cannot reproduce, `numpollen` too low.", settings, 'd')
            continue
        end
        for ptn in 1:numpartners
            partner = rand(partners)
            # subtract offspring mass from parent #XXX is that sensible, does it make a difference?
            parentmass = ind.size - noffs * ind.traits["seedsize"]
            if parentmass <= 0
                break
            else
                ind.size = parentmass
            end
            append!(patch.seedbank, createoffspring(noffs, ind, partner, settings["traitnames"]))
        end
    end
    simlog("Patch $(patch.id): $(length(patch.seedbank)) offspring", settings, 'd')
end

"""
    greproduce!(patch, settings)

Reproduction of individuals in a patch with global mating.
"""
function greproduce!(patch::Patch, settings::Dict{String, Any})
    for ind in patch.community
        ind.marked && continue # individual might not have established yet
        ind.size < ind.traits["repsize"] && continue
        metaboffs = settings["fertility"] * ind.size^(-1/4) * exp(-act/(boltz*patch.temp))
        noffs = rand(Poisson(metaboffs))
        noffs < 1 && continue
        partners = findmate([(map(x -> x.community, world)...)...], ind, settings)
        if length(partners) < 1 && rand() < ind.traits["selfing"]
            partners = [ind]
        elseif length(partners) < 1
            continue
        end
        numpartners = Integer(round(ind.traits["numpollen"]))
        for ptn in 1:numpartners
            partner = rand(partners, 1)[1]
            parentmass = ind.size - noffs * ind.traits["seedsize"] # subtract offspring mass from parent
            if parentmass <= 0
                continue
            else
                ind.size = parentmass
            end
            append!(patch.seedbank, createoffspring(noffs, ind, partner, settings["traitnames"]))
        end
    end
    simlog("Patch $(patch.id): $(length(patch.seedbank)) offspring", settings, 'd')
end

"""
    reproduce!(world, settings)

Carry out reproduction on all patches.
"""
function reproduce!(world::Array{Patch,1}, settings::Dict{String, Any})
    for patch in world
        if settings["globalmating"]
            greproduce!(patch, settings)
        elseif settings["mode"] == "zosterops"
            zreproduce!(patch, settings)
        else
            (patch.isisland || !settings["static"]) && reproduce!(patch, settings) # pmap(!,patch) ???
        end
    end
end

"""
    findmate(population, individual, settings)

Find a reproduction partner for the given individual in the given population.
"""
function findmate(population::AbstractArray{Individual, 1}, ind::Individual, settings::Dict{String, Any})
    #XXX This function is pretty expensive
    indstate = ind.marked
    ind.marked = true
    mates = Individual[]
    startidx = rand(eachindex(population))
    mateidx = startidx
    while true
        mate = population[mateidx]
        if !mate.marked && iscompatible(mate, ind, settings)
            push!(mates, mate)
            break
        end
        mateidx += 1
        mateidx > length(eachindex(population)) && (mateidx = 1)
        mateidx == startidx && break
    end
    ind.marked = indstate
    mates
end

"""
    createoffspring(noffs, individual, partner, traitnames, dimorphism)

The main reproduction function. Take two organisms and create the given number
of offspring individuals. Returns an array of individuals.
"""
function createoffspring(noffs::Integer, ind::Individual, partner::Individual, traitnames::Array{String, 1}, dimorphism::Bool=false)
    #XXX another expensive function
    offspring = Individual[]
    for i in 1:noffs # pmap? this loop could be factorized!
        partnergenome = meiosis(partner.genome, false) # offspring have different genome!
        mothergenome = meiosis(ind.genome, true)
        (isempty(partnergenome) || isempty(mothergenome)) && continue
        genome = vcat(partnergenome,mothergenome)
        traits = gettraitdict(genome, traitnames)
        marked = true
        fitness = 0.0
        newpartner = 0
        newsize = ind.traits["seedsize"]
        newid = rand(Int32) #XXX Is this collision-safe?
        sex = hermaphrodite
        if dimorphism
            rand(Bool) ? sex = male : sex = female
        end
        ind = Individual(ind.lineage, genome, traits, marked, fitness,
                         fitness, newsize, sex, newpartner, newid)
        push!(offspring, ind)
    end
    offspring
end
