# Main processes for GeMM

function mutate!(ind::Individual, temp::Float64, settings::Dict{String,Any})
    prob = ind.traits["mutprob"]
    for chrm in ind.genome
        for idx in eachindex(chrm.genes)
            charseq = collect(chrm.genes[idx].sequence)
            for i in eachindex(charseq)
                if rand() <= prob * exp(-act/(boltz*temp))
                    newbase = rand(collect("acgt"),1)[1]
                    while newbase == charseq[i]
                        newbase = rand(collect("acgt"),1)[1]
                    end
                    charseq[i] = newbase
                    for trait in chrm.genes[idx].codes
                        (contains(trait.name, "mutprob") && mutationrate != 0) && continue
                        contains(trait.name, "reptol") && settings["tolerance"] != "evo" && continue # MARK CAVE!
                        oldvalue = trait.value
                        contains(trait.name, "tempopt") && (oldvalue -= 273)
                        while oldvalue <= 0 # make sure sd of Normal dist != 0
                            oldvalue = abs(rand(Normal(0,0.01)))
                        end
                        newvalue = oldvalue + rand(Normal(0, oldvalue/phylconstr)) # CAVE: maybe handle temp + prec separately
                        newvalue < 0 && (newvalue=abs(newvalue))
                        (newvalue > 1 && contains(trait.name,"prob")) && (newvalue=1)
                        while newvalue <= 0 #&& contains(trait.name,"mut")
                            newvalue = trait.value + rand(Normal(0, trait.value/phylconstr))
                        end
                        contains(trait.name, "tempopt") && (newvalue += 273)
                        trait.value = newvalue
                    end
                end
            end
            chrm.genes[idx].sequence = String(charseq)
        end
    end
    ind.traits = chrms2traits(ind.genome)
end

function mutate!(patch::Patch, settings::Dict{String,Any})
    for ind in patch.community
        ind.age == 0 && mutate!(ind, patch.altitude, settings)
    end
end

function mutate!(world::Array{Patch, 1}, settings::Dict{String,Any})
    for patch in world
        (patch.isisland || !settings["static"]) && mutate!(patch, settings)
    end
end

function checkviability!(patch::Patch) # may consider additional rules... # maybe obsolete anyhow...
    idx=1
    while idx <= size(patch.community,1)
        dead = false
        patch.community[idx].size <= 0 && (dead = true)
        any(collect(values(patch.community[idx].traits)) .<= 0) && (dead = true)
        patch.community[idx].traits["repsize"] <= patch.community[idx].traits["seedsize"] && (dead = true)
        patch.community[idx].fitness <= 0 && (dead = true)
        if dead
            splice!(patch.community,idx)
            continue
        end
        idx += 1
    end
end

function checkviability!(world::Array{Patch,1}, static::Bool = true)
    for patch in world
        (patch.isisland || !static) && checkviability!(patch) # pmap(checkviability!,patch) ???
    end
end

"""
    establish!(p, n)
establishment of individuals in patch `p`. Sets fitness scaling parameter
according to adaptation to number `n` niches of the surrounding environment.
"""
function establish!(patch::Patch, nniches::Int64=1)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        if patch.community[idx].isnew || patch.community[idx].age == 0
            fitness = 1
            if !traitsexist(patch.community[idx], ["temptol", "tempopt"]) # CAVE: should check for more traits!
                info(STDERR, "Individual killed due to missing trait(s).")
                splice!(patch.community, idx) # kill it!
                continue
            end
            opt = patch.community[idx].traits["tempopt"]
            tol = patch.community[idx].traits["temptol"]
            fitness *= gausscurve(opt, tol, temp, 0.0)
            fitness > 1 && (fitness = 1) # should be obsolete
            fitness < 0 && (fitness = 0) # should be obsolete
            if nniches >= 2
                if !traitsexist(patch.community[idx], ["prectol", "precopt"]) # CAVE: should check for more traits!
                    info(STDERR, "Individual killed due to missing trait(s).")
                    splice!(patch.community, idx) # kill it!
                    continue
                end
                opt = patch.community[idx].traits["precopt"]
                tol = patch.community[idx].traits["prectol"]
                fitness *= gausscurve(opt, tol, patch.nichea, 0.0)
                fitness > 1 && (fitness = 1) # should be obsolete
                fitness < 0 && (fitness = 0) # should be obsolete
            end
            if nniches == 3
                opt = patch.community[idx].traits["nicheopt"]
                tol = patch.community[idx].traits["nichetol"]
                fitness *= gausscurve(opt, tol, patch.nicheb, 0.0)
                fitness > 1 && (fitness = 1) # should be obsolete
                fitness < 0 && (fitness = 0) # should be obsolete
            end  
            patch.community[idx].fitness = fitness
            patch.community[idx].isnew = false
        end
        idx += 1
    end
end


function establish!(world::Array{Patch,1}, nniches::Int64=1, static::Bool = true)
    for patch in world
        (patch.isisland || !static) && establish!(patch, nniches) # pmap(!,patch) ???
    end
end

"""
    survive!(p)
density independent survival of individuals in patch `p`
"""
function survive!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        if !patch.community[idx].isnew
            mass = patch.community[idx].size
            deathrate = mortality * mass^(-1/4) * exp(-act/(boltz*temp))
            dieprob = 1 - exp(-deathrate)
            if rand() < dieprob
                splice!(patch.community, idx)
                continue
            else
                patch.community[idx].age += 1
            end
        end
        idx += 1
    end
end

function survive!(world::Array{Patch,1}, static::Bool = true)
    for patch in world
        (patch.isisland || !static) && survive!(patch) # pmap(!,patch) ???
    end
end

"""
    grow!(p)
Growth of individuals in patch `p`
"""
function grow!(patch::Patch)
    temp = patch.altitude
    idx = 1
    while idx <= size(patch.community,1)
        if !traitsexist(patch.community[idx], ["repsize"])
            info(STDERR, "Individual killed due to missing trait(s).")
            splice!(patch.community, idx)
            continue
        elseif !patch.community[idx].isnew
            repsize = patch.community[idx].traits["repsize"]
            mass = patch.community[idx].size
            if mass <= repsize # stop growth if reached repsize 
                growth = growthrate * mass^(3/4) * exp(-act/(boltz*temp))
                newmass = mass + growth
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

function grow!(world::Array{Patch,1}, static::Bool = true)
    for patch in world
        (patch.isisland || !static) && grow!(patch) # pmap(!,patch) ???
    end
end

"""
    disperse!(w)
Dispersal of individuals within world (array of patches) `w`
"""
function disperse!(world::Array{Patch,1}, static::Bool = true) # TODO: additional border conditions, refoctorize
    colonizers = ""
    for patch in world
        idx = 1
        while idx <= size(patch.community,1)
            if !traitsexist(patch.community[idx], ["dispmean", "dispshape"])
                info(STDERR, "Individual killed due to missing trait(s).")
                splice!(patch.community,idx)
                continue
            elseif !patch.community[idx].isnew && patch.community[idx].age == 0
                dispmean = patch.community[idx].traits["dispmean"]
                dispshape = patch.community[idx].traits["dispshape"]
                xdir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # scaling so that geometric mean...
                ydir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # ...follows original distribution
                xdest = patch.location[1]+xdir
                ydest = patch.location[2]+ydir
                # !patch.isisland && checkborderconditions!(world,xdest,ydest)
                targets = unique([(floor(xdest),floor(ydest)),(ceil(xdest),floor(ydest)),(ceil(xdest),ceil(ydest)),(floor(xdest),ceil(ydest))])
                possdests = find(x->in(x.location,targets),world)
                filter!(x -> world[x].isisland, possdests) # disperse only to islands
                if static && patch.isisland
                    indleft = splice!(patch.community,idx) # only remove individuals from islands!
                end
                if size(possdests,1) > 0 # if no viable target patch, individual dies
                    if static && !patch.isisland
                        indleft = deepcopy(patch.community[idx])
                    end
                    indleft.isnew = true
                    destination = rand(possdests)
                    originisolated = patch.isolated && rand(Logistic(dispmean,dispshape)) <= isolationweight # additional roll for isolated origin patch
                    targetisolated = world[destination].isolated && rand(Logistic(dispmean,dispshape)) <= isolationweight # additional roll for isolated target patch
                    (!originisolated && !targetisolated) && push!(world[destination].community, indleft) # new independent individual
                    !patch.isisland && world[destination].isisland && (colonizers *= indleft.lineage * ",")
                end
                patch.isisland && (idx -= 1)
            end
            idx += 1
        end
    end
    colonizers
end

function compete!(patch::Patch)
    sort!(patch.community, by = x -> x.fitness)
    while sum(map(x -> x.size, patch.community)) >= patch.area # occupied area larger than available
        victim = rand(Geometric()) + 1 # fitness sorted?
        victim > length(patch.community) && (victim = length(patch.community))
        splice!(patch.community, victim)
    end
end

function compete!(world::Array{Patch,1}, static::Bool = true)
    for patch in world
        (patch.isisland || !static) && compete!(patch) # pmap(!,patch) ???
    end
end

"""
    reproduce!(w, p)
Reproduction of individuals in a patch `p` whithin a world (array of patches) `w`
"""
function reproduce!(world::Array{Patch,1}, patch::Patch) #TODO: refactorize!
    idx = 1
    temp = patch.altitude
    seedbank = Individual[]
    while idx <= size(patch.community,1)
        if !traitsexist(patch.community[idx], ["repradius", "repsize", "reptol", "seedsize", "mutprob"])
            info(STDERR, "Individual killed due to missing trait(s).")
            splice!(patch.community, idx)
            continue
        elseif !patch.community[idx].isnew && patch.community[idx].age > 0
            currentmass = patch.community[idx].size
            seedsize = patch.community[idx].traits["seedsize"]
            if currentmass >= patch.community[idx].traits["repsize"]
                reptol = patch.community[idx].traits["reptol"]
                metaboffs = fertility * currentmass^(-1/4) * exp(-act/(boltz*temp))
                noffs = rand(Poisson(metaboffs * patch.community[idx].fitness)) # add some stochasticity
                posspartners = findposspartners(world, patch.community[idx], patch.location) # this effectively controls frequency of reproduction
                # length(posspartners) == 0 && push!(posspartners, patch.community[idx]) # selfing if no partners # CAVE!
                if length(posspartners) > 0
                    partner = rand(posspartners)
                    parentmass = currentmass - noffs * seedsize # subtract offspring mass from parent
                    if parentmass <= 0
                        idx += 1 #splice!(patch.community, idx)
                        continue
                    else
                        patch.community[idx].size = parentmass
                    end
                    for i in 1:noffs # pmap? this loop could be factorized!
                        partnergenome = meiosis(partner.genome, false) # offspring have different genome!
                        mothergenome = meiosis(patch.community[idx].genome, true)
                        (length(partnergenome) < 1 || length(mothergenome) < 1) && continue
                        genome = vcat(partnergenome,mothergenome)
                        traits = chrms2traits(genome)
                        age = 0
                        isnew = false
                        fitness = 0.0
                        newsize = seedsize
                        ind = Individual(patch.community[idx].lineage, newsize, age, fitness, isnew, genome, traits, mtraits, ptraits)
                        push!(seedbank ,ind) # maybe actually deepcopy!?
                    end
                end
            end
        end
        idx += 1
    end
    append!(patch.community, seedbank)
end

function reproduce!(world::Array{Patch,1}, static::Bool = true)
    for patch in world
        (patch.isisland || !static) && reproduce!(world, patch) # pmap(!,patch) ???
    end
end
