# All functions related to reproduction:
# - mutation
# - reproduction

"""
    mutate!(traits, settings, locivar)

Loop over an array of traits, mutating each value in place along a normal distribution.
`locivar` can be used to scale the variance of the normal distribution used to draw new
trait values (together with `settings[phylconstr]`).
"""
function mutate!(traits::Array{Trait, 1}, settings::Dict{String, Any}, locivar::Float64 = 1.0)
    settings["phylconstr"] * locivar == 0 && return
    for trait in traits
        traitname = settings["traitnames"][trait.nameindex]
        occursin("seqsimilarity", traitname) && settings["fixtol"] && continue
        oldvalue = trait.value
        occursin("tempopt", traitname) && (oldvalue -= 273)
        while oldvalue <= 0 # make sure sd of Normal dist != 0
            oldvalue = abs(rand(Normal(0,0.01)))
        end
        newvalue = rand(Normal(oldvalue, oldvalue * settings["phylconstr"] * locivar))
        (newvalue > 1 && occursin("prob", traitname)) && (newvalue=1.0)
        while newvalue <= 0
            newvalue = rand(Normal(oldvalue, oldvalue * settings["phylconstr"] * locivar))
        end
        occursin("tempopt", traitname) && (newvalue += 273)
        trait.value = newvalue
    end
end

"""
    mutate!(individual, temp, settings)

Mutate an individual's genome (sequence and traits) in place.
"""
function mutate!(ind::Individual, temp::Float64, settings::Dict{String, Any})
    muts = settings["mutationrate"] * exp(-act/(boltz*temp))
    nmuts = rand(Poisson(muts))
    nmuts == 0 && return
    chrmidcs = rand(eachindex(ind.genome), nmuts)
    for c in chrmidcs
        ind.genome[c] = deepcopy(ind.genome[c])
        length(ind.genome[c].genes) == 0 && continue
        g = rand(eachindex(ind.genome[c].genes))
        charseq = collect(num2seq(ind.genome[c].genes[g].sequence))
        i = rand(eachindex(charseq))
        newbase = rand(collect("acgt"),1)[1]
        while newbase == charseq[i]
            newbase = rand(collect("acgt"),1)[1]
        end
        charseq[i] = newbase
        mutate!(ind.genome[c].genes[g].codes, settings)
        ind.genome[c].genes[g].sequence = deepcopy(ind.genome[c].genes[g].sequence)
        if length(charseq) > 21
            ind.genome[c].genes[g].sequence = seq2bignum(String(charseq))
        else
            ind.genome[c].genes[g].sequence = seq2num(String(charseq))
        end
    end
    ind.traits = gettraitdict(ind.genome, settings["traitnames"])
end

"""
    mutate!(patch, setting)

Mutate all seed individuals in a patch.
"""
function mutate!(patch::Patch, settings::Dict{String, Any})
    for ind in patch.seedbank
        mutate!(ind, patch.temp, settings)
    end
end

"""
    mutate!(world, settings)

Mutate the world. (That sounds scary!)
"""
function mutate!(world::Array{Patch, 1}, settings::Dict{String, Any})
    for patch in world
        (patch.isisland || !settings["static"]) && mutate!(patch, settings)
    end
end


"""
    reproduce!(patch, settings)

Reproduction of individuals in a patch.
"""
function reproduce!(patch::Patch, settings::Dict{String, Any}) #TODO: refactor!
    for ind in patch.community
        ind.marked && continue # individual might not have established yet
        ind.size < ind.traits["repsize"] && continue
        metaboffs = settings["fertility"] * ind.size^(-1/4) * exp(-act/(boltz*patch.temp))
        noffs = rand(Poisson(metaboffs))
        noffs < 1 && continue
        partners = findmate(patch.community, ind, settings["traitnames"])
        if length(partners) < 1 && rand() < ind.traits["selfing"]
            partners = [ind]
        else
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
            for ind in patch.community
                ind.marked && continue # individual might not have established yet
                ind.size < ind.traits["repsize"] && continue
                metaboffs = settings["fertility"] * ind.size^(-1/4) * exp(-act/(boltz*patch.temp))
                noffs = rand(Poisson(metaboffs))
                noffs < 1 && continue
                partners = findmate([(map(x -> x.community, world)...)...], ind, settings["traitnames"])
                if length(partners) < 1 && rand() < ind.traits["selfing"]
                    partners = [ind]
                else
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
        else
            (patch.isisland || !settings["static"]) && reproduce!(patch, settings) # pmap(!,patch) ???
        end
    end
end

"""
    findmate(population, individual, traitnames)

Find a reproduction partner for the given individual in the given population.
"""
function findmate(population::AbstractArray{Individual, 1}, ind::Individual, traitnames::Array{String, 1})
    indstate = ind.marked
    ind.marked = true
    mates = Individual[]
    startidx = rand(eachindex(population))
    mateidx = startidx
    while true
        mate = population[mateidx]
        if !mate.marked && iscompatible(mate, ind, traitnames)
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
    createoffspring(noffs, individual, partner, traitnames)

The main reproduction function. Take two organisms and create the given number
of offspring individuals. Returns an array of individuals.
"""
function createoffspring(noffs::Integer, ind::Individual, partner::Individual, traitnames::Array{String, 1})
    offspring = Individual[]
    for i in 1:noffs # pmap? this loop could be factorized!
        partnergenome = meiosis(partner.genome, false) # offspring have different genome!
        mothergenome = meiosis(ind.genome, true)
        (length(partnergenome) < 1 || length(mothergenome) < 1) && continue
        genome = vcat(partnergenome,mothergenome)
        traits = gettraitdict(genome, traitnames)
        marked = true
        fitness = 0.0
        newsize = ind.traits["seedsize"]
        ind = Individual(ind.lineage, genome, traits, marked, fitness,
                         fitness, newsize, rand(Int32))
        push!(offspring, ind)
    end
    offspring
end
