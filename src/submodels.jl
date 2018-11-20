# Main processes for GeMM

function mutate!(traits::Array{Trait, 1}, settings::Dict{String, Any})
    for trait in traits
        traitname = settings["traitnames"][trait.nameindex]
        occursin("reptol", traitname) && settings["fixtol"] && continue # MARK CAVE!
        oldvalue = trait.value
        occursin("tempopt", traitname) && (oldvalue -= 273)
        while oldvalue <= 0 # make sure sd of Normal dist != 0
            oldvalue = abs(rand(Normal(0,0.01)))
        end
        newvalue = oldvalue + rand(Normal(0, oldvalue * settings["phylconstr"])) # CAVE: maybe handle temp + prec separately
        newvalue < 0 && (newvalue=abs(newvalue))
        (newvalue > 1 && occursin("prob", traitname)) && (newvalue=1)
        while newvalue <= 0
            newvalue = trait.value + rand(Normal(0, trait.value * settings["phylconstr"]))
        end
        occursin("tempopt", traitname) && (newvalue += 273)
        trait.value = newvalue
    end
end

function mutate!(ind::Individual, temp::Float64, settings::Dict{String, Any})
    muts = settings["mutationrate"] * exp(-act/(boltz*temp)) # settings["mutsperind"]?
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

function mutate!(patch::Patch, settings::Dict{String, Any})
    for ind in patch.seedbank
        mutate!(ind, patch.temp, settings)
    end
end

function mutate!(world::Array{Patch, 1}, settings::Dict{String, Any})
    for patch in world
        (patch.isisland || !settings["static"]) && mutate!(patch, settings)
    end
end

function checkviability!(community::Array{Individual, 1}, settings::Dict{String, Any})
    idx=1
    while idx <= size(community,1)
        reason = ""
        dead = false
        community[idx].size <= 0 && (dead = true) && (reason *= "size ")
        any(collect(values(community[idx].traits)) .< 0) && (dead = true) && (reason *= "traitvalues ")
        community[idx].traits["repsize"] <= community[idx].traits["seedsize"] && (dead = true) && (reason *= "seed/rep ")
        community[idx].fitness < 0 && (dead = true) && (reason *= "fitness ")
        !traitsexist(community[idx].traits, settings) && (dead = true) && (reason *= "missingtrait ")
        if dead
            simlog("Individual not viable: $reason. Being killed.", settings, 'w')
            splice!(community,idx)
            continue
        end
        idx += 1
    end
end

function checkviability!(patch::Patch, settings::Dict{String, Any})
    checkviability!(patch.community, settings)
end

function checkviability!(world::Array{Patch,1}, settings::Dict{String, Any})
    for patch in world
        (patch.isisland || !settings["static"]) && checkviability!(patch, settings) # pmap(checkviability!,patch) ???
    end
end

"""
    establish!(p, n)
establishment of individuals in patch `p`. Sets fitness scaling parameter
according to adaptation to number `n` niches of the surrounding environment.
"""
function establish!(patch::Patch, nniches::Int=1)
    temp = patch.temp
    idx = 1
    while idx <= size(patch.community,1)
        if patch.community[idx].marked
            fitness = 1
            opt = patch.community[idx].traits["tempopt"]
            tol = patch.community[idx].traits["temptol"]
            fitness *= gausscurve(opt, tol, temp, 0.0)
            fitness > 1 && (fitness = 1) # should be obsolete
            fitness < 0 && (fitness = 0) # should be obsolete
            if nniches >= 2
                opt = patch.community[idx].traits["precopt"]
                tol = patch.community[idx].traits["prectol"]
                fitness *= gausscurve(opt, tol, patch.prec, 0.0)
                fitness > 1 && (fitness = 1) # should be obsolete
                fitness < 0 && (fitness = 0) # should be obsolete
            end
            if nniches == 3
                # XXX 'nicheopt' and 'nichetol' don't exist currently!
                opt = patch.community[idx].traits["nicheopt"]
                tol = patch.community[idx].traits["nichetol"]
                fitness *= gausscurve(opt, tol, patch.nicheb, 0.0)
                fitness > 1 && (fitness = 1) # should be obsolete
                fitness < 0 && (fitness = 0) # should be obsolete
            end  
            patch.community[idx].fitness = fitness
            patch.community[idx].marked = false
        end
        idx += 1
    end
end


function establish!(world::Array{Patch,1}, nniches::Int=1, static::Bool = true)
    for patch in world
        (patch.isisland || !static) && establish!(patch, nniches) # pmap(!,patch) ???
    end
end

"""
    survive!(p)
density independent survival of individuals in patch `p`
"""
function survive!(patch::Patch, mortality::Float64)
    temp = patch.temp
    idx = 1
    while idx <= size(patch.community,1)
        if !patch.community[idx].marked
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

function survive!(world::Array{Patch,1}, settings::Dict{String, Any})
    for patch in world
        (patch.isisland || !settings["static"]) && survive!(patch, settings["mortality"]) # pmap(!,patch) ???
    end
end

"""
    disturb!(p,i)
species-independent mortality due to disturbance on patch `p`
"""
function disturb!(patch::Patch, intensity::Int)
    length(patch.community) <= 0 && return
    deaths = Integer(round(length(patch.community) * (intensity/100)))
    dead = unique(rand(1:length(patch.community), deaths))
    while length(dead) < deaths
        newdead = rand(1:length(patch.community), deaths-length(dead))
        dead = unique(vcat(dead, newdead))
    end
    sort!(dead)
    deleteat!(patch.community, dead)
end

function disturb!(world::Array{Patch,1}, settings::Dict{String, Any})
    (settings["disturbance"] == 0) && return
    if settings["disturbance"] > 100
        simlog("disturbance must be no more than 100%", 'w')
        settings["disturbance"] = 100
    end
    for patch in world
        (patch.isisland || !settings["static"]) && disturb!(patch, settings["disturbance"])
    end
end

let speciespool = Individual[] # TODO: `speciespool` should be part of the `world` object holding all grid cells 
    function initspeciespool!(settings::Dict{String, Any})
        for i in 1:settings["global-species-pool"]
            push!(speciespool, createind(settings))
        end
    end
    
    """
        invade!(patch, pressure)
    select `pressure` individuals from the global species pool and add them to the patch.
    """
    function invade!(patch::Patch, pressure::Int)
        invaders = deepcopy(rand(speciespool, pressure))
        append!(patch.community, invaders)
    end

    global function invade!(world::Array{Patch,1}, settings::Dict{String, Any})
        (settings["propagule-pressure"] == 0 || settings["global-species-pool"] == 0) && return
        (length(speciespool) == 0) && initspeciespool!(settings)
        for patch in world
            patch.invasible && invade!(patch, settings["propagule-pressure"])
        end
    end
end

"""
    grow!(p)
Growth of individuals in patch `p`
"""
function grow!(patch::Patch, growthrate::Float64)
    temp = patch.temp
    idx = 1
    while idx <= size(patch.community,1)
        if !patch.community[idx].marked
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

function grow!(world::Array{Patch,1}, settings::Dict{String, Any})
    for patch in world
        (patch.isisland || !settings["static"]) && grow!(patch, settings["growthrate"]) # pmap(!,patch) ???
    end
end

"""
    disperse!(w)
Dispersal of individuals within world (array of patches) `w`
"""
function disperse!(world::Array{Patch,1}, static::Bool = true) # TODO: additional border conditions, refactorize
    for patch in world
        idx = 1
        while idx <= size(patch.seedbank,1)
            dispmean = patch.seedbank[idx].traits["dispmean"]
            dispshape = patch.seedbank[idx].traits["dispshape"]
            xdir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # scaling so that geometric mean...
            ydir = rand([-1,1]) * rand(Logistic(dispmean,dispshape))/sqrt(2) # ...follows original distribution
            xdest = patch.location[1] + Int(round(xdir))
            ydest = patch.location[2] + Int(round(ydir))
            !patch.isisland ? target = checkborderconditions(world, xdest, ydest) : target = (xdest, ydest)
            possdest = findall(x -> in(x.location, [target]), world)
            static && filter!(x -> world[x].isisland, possdest) # disperse only to islands
            if !static || patch.isisland
                indleft = splice!(patch.seedbank,idx) # only remove individuals from islands!
                idx -= 1
            end
            if length(possdest) > 0 # if no viable target patch, individual dies
                if static && !patch.isisland
                    indleft = patch.seedbank[idx]
                end
                destination = rand(possdest) # currently there is only one possible destination
                push!(world[destination].community, indleft)
            end
            idx += 1
        end
    end
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
    reproduce!(p, settings)
Reproduction of individuals in a patch `p`
"""
function reproduce!(patch::Patch, settings::Dict{String, Any}) #TODO: refactorize!
    identifyAdults!(patch)
    patch.phylo = Array{Int}(undef, 0, 2)
    for ind in patch.community
        if !ind.marked && ind.age > 0
            currentmass = ind.size
            seedsize = ind.traits["seedsize"]
            if currentmass >= ind.traits["repsize"]
                reptol = ind.traits["reptol"]
                metaboffs = settings["fertility"] * currentmass^(-1/4) * exp(-act/(boltz*patch.temp))
                noffs = rand(Poisson(metaboffs))# * ind.fitness)) # add some stochasticity
                if noffs < 1
                    #simlog("0 offspring chosen", settings, 'd') #DEBUG - noisy!
                    continue
                end
                partner = findposspartner(patch, ind, settings["traitnames"])
                if partner != nothing
                    parentmass = currentmass - noffs * seedsize # subtract offspring mass from parent
                    if parentmass <= 0
                        continue
                    else
                        ind.size = parentmass
                    end
                    ## save origin of actually reproducing individuals:
                    patch.phylo = vcat(patch.phylo, [ind.id ind.parentid[1]]) # ind.parentid[2]; # for now only use sparse genealogy
                    # partner.id partner.parentid[1] partner.parentid[2]]
                    for i in 1:noffs # pmap? this loop could be factorized!
                        partnergenome = meiosis(partner.genome, false) # offspring have different genome!
                        mothergenome = meiosis(ind.genome, true)
                        (length(partnergenome) < 1 || length(mothergenome) < 1) && continue
                        genome = vcat(partnergenome,mothergenome)
                        traits = gettraitdict(genome, settings["traitnames"])
                        age = 0
                        marked = true
                        fitness = 0.0
                        newsize = seedsize
                        ind = Individual(ind.lineage, genome, traits, age, marked, fitness,
                                         newsize, rand(Int32), ind.id) #, partner.id))
                        push!(patch.seedbank, ind) # maybe actually deepcopy!?
                    end
                end
            end
        end
    end
    checkviability!(patch.seedbank, settings)
    simlog("Patch $(patch.id): $(length(patch.seedbank)) offspring", settings, 'd')
end

function reproduce!(world::Array{Patch,1}, settings::Dict{String, Any})
    for patch in world
        (patch.isisland || !settings["static"]) && reproduce!(patch, settings) # pmap(!,patch) ???
    end
end

function changehabitat!(world::Array{Patch,1}, settings::Dict{String, Any})
    # TODO: record trajectory? input trajectory?
    changetemp!(world, settings["sdtemp"])
    changeprec!(world, settings["sdprec"])
end    
 
function changetemp!(world::Array{Patch,1}, sdtemp::Float64)
    sdtemp == 0 && return
    deltaval = rand(Normal(0.0, sdtemp))
    for patch in world
        patch.temp += deltaval 
    end
    markthem!(world)
end    
 
function changeprec!(world::Array{Patch,1}, sdprec::Float64)
    sdprec == 0 && return
    deltaval = rand(Normal(0.0, sdprec))
    for patch in world
        patch.prec += deltaval 
    end
    markthem!(world)
end    
