# Auxilliary functions for GeMM

"""
    meiosis(genome, maternal)

Carry out meiosis on a genome (marked as maternal or not). Returns a haploid
gamete genome. (genome => array of chromosomes)
"""
function meiosis(genome::Array{Chromosome,1}, maternal::Bool) # TODO: include further dynamics, errors...
    firstset = findall(x -> x.maternal, genome)
    secondset = findall(x -> !x.maternal, genome)
    length(firstset) != length(secondset) && return Chromosome[] # CAVEAT: more elegant solution...
    gameteidxs = []
    for i in eachindex(firstset)
        push!(gameteidxs, rand([firstset[i], secondset[i]]))
    end
    gamete = Chromosome[]
    for i in gameteidxs
        push!(gamete, Chromosome(genome[i].genes, maternal))
    end
    gamete
end

"""
    getmeantraitvalue(traits, traitidx)

Take an array of traits and return the mean value of the indexed trait.
"""
function getmeantraitvalue(traits::Array{Trait, 1}, traitidx::Integer)
    mean(skipmissing(map(x -> x.value, filter(x -> x.nameindex == traitidx, traits))))
end

"""
    getstdtraitvalue(traits, traitidx)

Take an array of traits and return the standard deviation of the indexed trait.
"""
function getstdtraitvalue(traits::Array{Trait, 1}, traitidx::Integer)
    if length(traits) <= 1
        0.0
    else
        std(skipmissing(map(x -> x.value, filter(x -> x.nameindex == traitidx, traits))))
    end
end

"""
    gettraitdict(chromosomes, traitnames)

Convert a genome (an array of chromosomes) into a dict of traits and their values.
"""
function gettraitdict(chrms::Array{Chromosome, 1}, traitnames::Array{String, 1})
    genes = AbstractGene[]
    nchrms = 0
    ngenes = 0
    for chrm in chrms
        append!(genes, chrm.genes)
        nchrms += 1
        ngenes += length(chrm.genes)
    end
    traits = vcat(map(g -> g.codes, genes)...)
    traitdict = Dict{String, Float64}()
    for traitidx in eachindex(traitnames)
        traitdict[traitnames[traitidx]] = getmeantraitvalue(traits, traitidx)
        traitdict[string(traitnames[traitidx]) * "sd"] = getstdtraitvalue(traits, traitidx)
    end
    traitdict["nlnkgunits"] = nchrms
    traitdict["ngenes"] = ngenes
    traitdict
end

"""
    gettraitdict(genes, traitnames)

Calculate the trait dict for an array of genes.
"""
function gettraitdict(genes::Array{AbstractGene, 1}, traitnames::Array{String, 1})
    for gene in genes
        append!(traits, gene.codes)
    end
    traitdict = Dict{String, Float64}()
    for traitidx in unique(map(x -> x.nameindex, traits))
        traitdict[traitnames[traitidx]] = mean(skipmissing(map(x -> x.value, filter(x -> x.nameindex == traitidx, traits))))
    end
    traitdict
end

"""
    gettraitdict(traits, traitnames)

Construct a trait dict from a list of Trait objects.
"""
function gettraitdict(traits::Array{Trait, 1}, traitnames::Array{String, 1})
    traitdict = Dict{String, Float64}()
    for traitidx in unique(map(x -> x.nameindex, traits))
        #XXX Why do we use traitnames here?
        traitdict[traitnames[traitidx]] = mean(skipmissing(map(x -> x.value, filter(x -> x.nameindex == traitidx, traits))))
    end
    traitdict
end

"""
    traitsexist(traits, settings)

Check a trait dict to make sure it contains the full set of traitnames required
by the model (as defined in the settings).
"""
function traitsexist(traits::Dict{String, Float64}, settings::Dict{String, Any})
    missingtraits = setdiff(settings["traitnames"], keys(traits))
    if length(missingtraits) > 0
        simlog("Missing trait $missingtraits. Individual might be killed.", settings, 'w')
        return false
    end
    true
end

"""
    traitsexist(individual, settings)

Make sure an individual organism has the full set of traits required by the model
(as defined in the settings).
"""
function traitsexist(ind::Individual, settings::Dict{String, Any})
    # FIXME If a trait doesn't exist, we need an error
    traitnames = settings["traitnames"]
    for trait in traitnames
        if !haskey(ind.traits, trait)
            simlog("Individual is missing trait $trait. Might be killed.", settings, 'w')
            return false
        end
    end
    true
end

"""
    gausscurve(b, c, x, a=1.0)

Calculate the value of the Gauss function ("bell curve") at point x; with
a being the maximum height of the curve, b the position of the curve center and
c the standard deviation ("width").
"""
function gausscurve(b, c, x, a = 1.0)
    if c != 0 && a != 1.0
        a = 1 / (c * sqrt(2 * pi))
        y = a * exp(-(x-b)^2/(2*c^2))
    elseif c != 0
        y = a * exp(-(x-b)^2/(2*c^2))
    else
        y = 0.0
    end
end

"""
    diversity(w)

Calculate (average) alpha, beta and gamma diversity of the world.
Returns a tuple with the three values (a,b,g).
cf. Veech et al. 2002
"""
function diversity(world::Array{Patch,1})
    # calculate diversity with the Shannon-Wiener index
    function shannon(index::Dict{String,Int})
        isempty(index) && return 0
        total = sum(x -> index[x], keys(index))
        -1 * sum(s -> (index[s]/total)*log(index[s]/total), keys(index))
    end
    alphas = Float64[]
    globalindex = Dict{String, Int}()
    for p in world
        localindex = Dict{String, Int}()
        #TODO Should not use `whoiswho`
        for s in keys(p.whoiswho)
            localindex[s] = length(p.whoiswho[s])
        end
        merge!(+, globalindex, localindex)
        append!(alphas, shannon(localindex))
    end
    alpha = mean(alphas)
    gamma = shannon(globalindex)
    beta = gamma - alpha
    (alpha, beta, gamma)
end

"""
    freespace(world)

Calculate the average amount of free space on each patch.
"""
function freespace(world::Array{Patch,1})
    space = 0
    for p in world
        if length(p.community) != 0
            space += p.area - sum(x -> x.size, p.community)
        else
            space += p.area
        end
    end
    round((space/length(world))/1e6, digits = 3)
end

"""
    findisland(w)

within world `w`, find out in which direction from the continent the island(s) lie(s).
"""
function findisland(world::Array{Patch,1})
    xmin = minimum(map(x->x.location[1],world))
    xmax = maximum(map(x->x.location[1],world))
    ymin = minimum(map(x->x.location[2],world))
    ymax = maximum(map(x->x.location[2],world))
    westernborder = filter(x->x.location[1]==xmin,world)
    northernborder = filter(x->x.location[2]==ymax,world)
    easternborder = filter(x->x.location[1]==xmax,world)
    southernborder = filter(x->x.location[2]==ymin,world)
    if all(map(x->x.isisland,westernborder))
        return "west"
    elseif all(map(x->x.isisland,northernborder))
        return "north"
    elseif all(map(x->x.isisland,easternborder))
        return "east"
    elseif all(map(x->x.isisland,southernborder))
        return "south"
    else
        return "none"
    end
end

"""
    checkborderconditions!(w, x, y)

check if coordinates `x` and `y` lie within world `w` and correct if not,
considering defined border conditions.
"""
function checkborderconditions(world::Array{Patch,1}, xdest::Int, ydest::Int)
    xmin = minimum(map(x->x.location[1],world))
    xmax = maximum(map(x->x.location[1],world))
    ymin = minimum(map(x->x.location[2],world))
    ymax = maximum(map(x->x.location[2],world))
    xrange = xmax - xmin + 1 # we're counting cells!
    yrange = ymax - ymin + 1 # we're counting cells!
    xshift = xdest - xmin + 1 # 1-based count of cells
    yshift = ydest - ymin + 1 # 1-based count of cells
    xshift > 0 ? outofx = abs(xshift) : outofx = abs(xshift) + 1
    while outofx > xrange
        outofx -= xrange
    end
    outofx -= 1
    yshift > 0 ? outofy = abs(yshift) : outofy = abs(yshift) + 1
    while outofy > yrange
        outofy -= yrange
    end
    outofy -= 1
    islanddirection = findisland(world::Array{Patch,1})
    if islanddirection == "west"
        xdest > xmax && (xdest = xmax - outofx) # east: reflective
        ydest < ymin && (ydest = ymax - outofy) # south: periodic
        ydest > ymax && (ydest = ymin + outofy) # north: periodic
    elseif islanddirection == "north"
        ydest < ymin && (ydest = ymin + outofy) # south: reflective
        xdest < xmin && (xdest = xmax - outofx) # west: periodic
        xdest > xmax && (xdest = xmin + outofx) # east: periodic
    elseif islanddirection == "east"
        xdest < xmin && (xdest = xmin + outofx) # west: reflective
        ydest < ymin && (ydest = ymax - outofy) # south: periodic
        ydest > ymax && (ydest = ymin + outofy) # north: periodic
    elseif islanddirection == "south"
        ydest > ymax && (ydest = ymax + outofy) # north: reflective
        xdest < xmin && (xdest = xmax - outofx) # west: periodic
        xdest > xmax && (xdest = xmin + outofx) # east: periodic
    else
        ydest > ymax && (ydest = ymin + outofy) # north: periodic
        xdest > xmax && (xdest = xmin + outofx) # east: periodic
        ydest < ymin && (ydest = ymax - outofy) # south: periodic
        xdest < xmin && (xdest = xmax - outofx) # west: periodic
    end
    xdest, ydest
end

"""
    identifyAdults!(patch)

Build up the `whoiswho` index of individuals and species in a patch.
"""
function identifyAdults!(patch::Patch)
    adultspeciesidx = Dict{String, Array{Int, 1}}()
    for i in eachindex(patch.community)
        patch.community[i].size < patch.community[i].traits["repsize"] && continue
        adultspeciesidx = merge(append!, adultspeciesidx, Dict(patch.community[i].lineage => [i]))
    end
    patch.whoiswho = adultspeciesidx
end

"""
    getseqsimilarity(seqone, seqtwo)

Compare two strings and return similarity.
"""
function getseqsimilarity(indgene::AbstractString, mategene::AbstractString)
    basediffs = 0
    for i in eachindex(indgene) # this is actually faster than `sum(collect(indgene) .== collect(mategene))`
        try
            indgene[i] != mategene[i] && (basediffs += 1) # alternatively use bioinformatics tools
        catch # e.g., in case of differently long genes
            basediffs += 1
        end
    end
    seqidentity = 1 - (basediffs / length(indgene))
end

"""
    getseq(genome, traitidx)

Find and return the sequence of one gene that codes for the given trait `traitidx`.
"""
function getseq(genome::Array{Chromosome, 1}, traitidx::Integer)
    seq = ""
    for chrm in genome
        for gene in chrm.genes
            any(x -> x.nameindex == traitidx, gene.codes) && (seq = num2seq(gene.sequence)) # use one compatibility gene randomly
        end
    end
    seq
end

"""
    iscompatible(mate, individual, traitnames)

Check to see whether two individual organisms are reproductively compatible.
"""
function iscompatible(mate::Individual, ind::Individual, traitnames::Array{String, 1})
    compatidx = findfirst(x -> x == "compat", traitnames)
    indgene = getseq(ind.genome, compatidx)
    mategene = getseq(mate.genome, compatidx)
    seqidentity = getseqsimilarity(indgene, mategene)
    seqidentity >= ind.traits["reptol"]
end

"""
    findmate(patch, individual, traitnames)

Find a reproduction partner for the given individual in the given patch.
"""
function findmate(patch::Patch, ind::Individual, traitnames::Array{String, 1})
    indstate = ind.marked
    ind.marked = true
    mates = Individual[]
    communityidxs = patch.whoiswho[ind.lineage]
    startidx = rand(1:length(communityidxs))
    mateidx = startidx
    while true
        mate = patch.community[communityidxs[mateidx]]
        if !mate.marked && iscompatible(mate, ind, traitnames)
            push!(mates, mate)
            break
        end
        mateidx += 1
        mateidx > length(communityidxs) && (mateidx = 1)
        mateidx == startidx && break
    end
    ind.marked = indstate
    mates
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
        if !mate.marked  && iscompatible(mate, ind, traitnames)
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
                         fitness, newsize, rand(Int32), ind.originid)
        push!(offspring, ind)
    end
    offspring
end

"""
    seq2num(sequence)

Convert a DNA base sequence (a string) into binary and then into an integer.
This saves memory.
"""
function seq2num(sequence::String)
    bases = "acgt"
    binary = ""
    for base in sequence
        binary *= string(findfirst(x -> x == base, bases) + 3, base = 2)
    end
    parse(Int, binary, base = 2) # Int64 allows for max length of 21bp
end

"""
    seq2bignum(sequence)

Convert a DNA base sequence (a string) into binary and then into an BigInt (for
larger genes). This saves memory.
"""
function seq2bignum(sequence::String)
    bases = "acgt"
    binary = ""
    for base in sequence
        binary *= string(findfirst(x -> x == base, bases) + 3, base = 2)
    end
    parse(BigInt, binary, base = 2)
end

"""
    num2seq(n)

Convert an integer into binary and then into a DNA base sequence string.
"""
function num2seq(n::Integer)
    bases = "acgt"
    binary = string(n, base = 2)
    sequence = ""
    for i in 1:3:(length(binary) - 2)
        sequence *= string(bases[parse(Int, binary[i:(i + 2)], base = 2) - 3])
    end
    sequence
end

"""
    createtraits(settings)

Create an array of trait objects generated from the default trait values (with a
random offset).
"""
function createtraits(settings::Dict{String, Any}) #TODO: this is all very ugly. (case/switch w/ v. 2.0+?)
    traitnames = settings["traitnames"]
    traits = Trait[]
    # exponential distributions of body sizes:
    repoffset = settings["maxrepsize"] - settings["minrepsize"]
    seedoffset = settings["maxseedsize"] - settings["minseedsize"]
    tempoffset = settings["maxtemp"] - settings["mintemp"]
    sizes = Vector{Float64}(undef, 2)
    while true
        sizes[1] = exp(settings["minrepsize"] + repoffset * rand())
        sizes[2] = exp(settings["minseedsize"] + seedoffset * rand())
        sizes[1] > sizes[2] && break
    end
    repsize, seedsize = sizes
    for idx in eachindex(traitnames)
        if occursin("dispshape", traitnames[idx])
            push!(traits, Trait(idx, rand() * settings["dispshape"]))
        elseif occursin("dispmean", traitnames[idx])
            push!(traits, Trait(idx, rand() * settings["dispmean"]))
        elseif occursin("precopt", traitnames[idx])
            push!(traits, Trait(idx, rand() * settings["precrange"]))
        elseif occursin("prectol", traitnames[idx])
            push!(traits, Trait(idx, rand() * settings["maxbreadth"]))
        elseif occursin("repsize", traitnames[idx])
            push!(traits, Trait(idx, repsize))
        elseif occursin("reptol", traitnames[idx]) && settings["fixtol"]
            push!(traits, Trait(idx, settings["tolerance"]))
        elseif occursin("reptol", traitnames[idx]) && !settings["fixtol"]
            push!(traits, Trait(idx, rand())) # assortative mating might evolve
        elseif occursin("seedsize", traitnames[idx])
            push!(traits, Trait(idx, seedsize))
        elseif occursin("tempopt", traitnames[idx])
            push!(traits, Trait(idx, settings["mintemp"] + rand() * tempoffset))
        elseif occursin("temptol", traitnames[idx])
            push!(traits, Trait(idx, rand() * settings["maxbreadth"]))
        else
            push!(traits, Trait(idx, rand()))
        end
    end
    traits
end

"""
    creategenes(ngenes, traits, settings)

Randomly create a given number of gene objects, with their base sequence and
associated traits. Returns the result as an array of AbstractGenes.
"""
function creategenes(ngenes::Int, traits::Array{Trait,1}, settings::Dict{String, Any})
    genes = AbstractGene[]
    compatidx = findfirst(x -> x == "compat", settings["traitnames"])
    for i in 1:ngenes
        sequence = String(rand(collect("acgt"), settings["smallgenelength"])) # arbitrary start sequence
        seqint = seq2num(sequence)
        codesfor = Trait[]
        push!(genes,Gene(seqint, codesfor))
    end
    for trait in traits
        if trait.nameindex == compatidx
            continue
        end
        ncodinggenes = rand(Geometric(1 - settings["degpleiotropy"])) + 1
        codinggenes = rand(genes,ncodinggenes)
        for gene in codinggenes
            push!(gene.codes,trait)
        end
    end
    if settings["usebiggenes"]
        push!(genes, BigGene(seq2bignum(String(rand(collect("acgt"), settings["biggenelength"]))), [Trait(compatidx, 0.5)]))
    else
        push!(genes, Gene(seq2num(String(rand(collect("acgt"), settings["smallgenelength"]))), [Trait(compatidx, 0.5)]))
    end
    genes
end

"""
    createchrms(nchrms, genes)

Randomly distribute the passed genes into the given number of chromosomes.
Returns an array of chromosome objects.
"""
function createchrms(nchrms::Int,genes::Array{AbstractGene,1})
    ngenes=size(genes,1)
    if nchrms>1
        chrmsplits = sort(rand(1:ngenes,nchrms-1))
        chromosomes = Chromosome[]
        for chr in 1:nchrms
            if chr==1 # first chromosome
                push!(chromosomes, Chromosome(genes[1:chrmsplits[chr]], true))
            elseif chr==nchrms # last chromosome
                push!(chromosomes, Chromosome(genes[(chrmsplits[chr-1]+1):end], true))
            else
                push!(chromosomes, Chromosome(genes[(chrmsplits[chr-1]+1):chrmsplits[chr]], true))
            end
        end
    else # only one chromosome
        chromosomes = [Chromosome(genes, true)]
    end
    secondset = deepcopy(chromosomes)
    for chrm in secondset
        chrm.maternal = !chrm.maternal
    end
    append!(chromosomes,secondset)
    chromosomes
end

"""
    createind(settings, marked=false)

Create an individual organism with of a new species with a random genome.
"""
function createind(settings::Dict{String, Any}, marked::Bool = false)
    id = rand(Int32)
    parentid = rand(Int32)
    lineage = randstring(4)
    ngenes = settings["maxloci"] * length(settings["traitnames"])
    ngenes < 1 && (ngenes = 1)
    traits = createtraits(settings)
    genes = creategenes(ngenes, traits, settings)
    randchrms = rand(1:length(genes))
    if settings["linkage"] == "none"
        nchrms = length(genes)
    elseif settings["linkage"] == "full"
        nchrms = 1
    else
        nchrms = randchrms
    end
    chromosomes = createchrms(nchrms, genes)
    locivar = rand() #XXX Is this legit? Don't quite understand `locivar` yet...
    varyalleles!(chromosomes, settings, locivar)
    traitdict = gettraitdict(chromosomes, settings["traitnames"])
    if settings["indsize"] == "adult"
        indsize = traitdict["repsize"]
    elseif settings["indsize"] == "seed"
        indsize = traitdict["seedsize"]
    else
        indsize = traitdict["seedsize"] + rand() * traitdict["repsize"] # XXX: sizes shouldn't be uniformally dist'd
    end
    Individual(lineage, chromosomes, traitdict, marked, 1.0, 1.0, indsize, id)
end

"""
    varyalleles!(genes, settings, locivar)

Mutate gene traits in the passed array of genes.
"""
function varyalleles!(genes::Array{AbstractGene, 1}, settings::Dict{String, Any}, locivar::Float64)
    locivar == 0 && return
    for gene in genes
        mutate!(gene.codes, settings, locivar)
    end
end

"""
    varyalleles!(chromosomes, settings, locivar)

Mutate gene traits in the passed array of chromosomes.
"""
function varyalleles!(chrms::Array{Chromosome, 1}, settings::Dict{String, Any}, locivar::Float64)
    locivar == 0 && return
    for chrm in chrms
        varyalleles!(chrm.genes, settings, locivar)
    end
end

"""
    markthem!(community)

Set each individual in the community (= array of individuals) as "marked".
"""
function markthem!(community::Array{Individual, 1})
    for ind in community
        ind.marked = true
    end
end

"""
    markthem!(habitat)

Set each individual in the given patch as "marked".
"""
function markthem!(habitat::Patch)
    markthem!(habitat.community)
end

"""
    markthem!(world)

Set every individual in the world as "marked".
"""
function markthem!(world::Array{Patch, 1})
    for habitat in world
        markthem!(habitat)
    end
end

"""
    updateorigin!(ind)

Update the origin ID for an individual.
"""
function updateorigin!(ind::Individual)
    ind.originid = ind.id
end

"""
    updateorigin!(community)

Update the origin ID for all individuals in a community.
"""
function updateorigin!(community::Array{Individual, 1})
    for ind in community
        updateorigin!(ind)
    end
end

"""
    updateorigin!(patch)

Update the origin ID for all individuals in a patch.
"""
function updateorigin!(habitat::Patch)
    updateorigin!(habitat.community)
end

"""
    updateorigin!(world)

Update the origin ID for all individuals in a community.
"""
function updateorigin!(world::Array{Patch, 1})
    for patch in world
        updateorigin!(patch)
    end
end
    
