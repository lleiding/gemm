# Subsidiary functions for GeMM

function meiosis(genome::Array{Chromosome,1}, maternal::Bool) # TODO: include further dynamics, errors...
    firstset = findall(x -> x.maternal, genome)
    secondset = findall(x -> !x.maternal, genome)
    length(firstset) != length(secondset) && return Chromosome[] # CAVE: more elegant solution...
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

function chrms2traits(chrms::Array{Chromosome, 1}, traitnames::Array{String, 1})
    genes = AbstractGene[]
    for chrm in chrms
        append!(genes, chrm.genes)
    end
    traits = Trait[]
    for gene in genes
        append!(traits, gene.codes)
    end
    traitdict = Dict{String, Float64}()
    for traitidx in unique(map(x -> x.nameindex, traits))
        traitdict[traitnames[traitidx]] = mean(map(x -> x.value, filter(x -> x.nameindex == traitidx, traits)))
    end
    traitdict
end

function traitsexist(traits::Dict{String, Float64}, traitnames::Array{String, 1})
    for trait in traitnames
        if !haskey(traits, trait)
            simlog("Missing trait $trait. Individual might be killed.", 'w')
            return false
        end
    end
    true
end

function traitsexist(ind::Individual, traitnames::Array{String, 1})
    # FIXME If a trait doesn't exist, we need an error
    for trait in traitnames
        if !haskey(ind.traits, trait)
            simlog("Individual is missing trait $trait. Might be killed.", 'w')
            return false
        end
    end
    true
end

function traitsexist(traits::Dict{String, Float64}, traitname::String)
    if !haskey(traits, traitname)
        simlog("Missing trait $traitname. Individual might be killed.", 'w')
        return false
    end
    true
end

function traitsexist(ind::Individual, traitname::String)
    if !haskey(ind.traits, traitname)
        simlog("Individual is missing trait $traitname. Might be killed.", 'w')
        return false
    end
    true
end

function gausscurve(b::Float64, c::Float64, x::Float64, a::Float64=1.0)
    if c != 0 && a != 1.0
        a = 1 / (c * sqrt(2 * pi))
        y = a * exp(-(x-b)^2/(2*c^2))
    elseif c != 0
        y = a * exp(-(x-b)^2/(2*c^2))
    else
        y = 0
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

# Calculate the average amount of free space on each patch
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

function identifyAdults!(patch::Patch)
    adultspeciesidx = Dict{String, Array{Int, 1}}()
    for i in eachindex(patch.community)
        patch.community[i].size < patch.community[i].traits["repsize"] && continue
        adultspeciesidx = merge(append!, adultspeciesidx, Dict(patch.community[i].lineage => [i]))
    end
    patch.whoiswho = adultspeciesidx
end

function iscompatible(mate::Individual, ind::Individual, traitnames::Array{String, 1})
    compatidx = findfirst(x -> x == "compat", traitnames)
    tolerance = ind.traits["reptol"]
    indgene = ""
    for chrm in ind.genome
        for gene in chrm.genes
            any(x -> x.nameindex == compatidx, gene.codes) && (indgene = num2seq(gene.sequence)) # use one compatibility gene randomly
        end
    end
    mategene = ""
    for chrm in mate.genome
        for gene in chrm.genes
            any(x -> x.nameindex == compatidx, gene.codes) && (mategene = num2seq(gene.sequence)) # use one compatibility gene randomly
        end
    end
    basediffs = 0
    for i in eachindex(indgene)
        try
            indgene[i] != mategene[i] && (basediffs += 1) # alternatively use bioinformatics tools
        catch # e.g., in case of differently long genes
            basediffs += 1
        end
    end
    seqidentity = 1 - (basediffs / length(indgene))
    seqidentity < tolerance && return false
    true
end

function findposspartner(patch::Patch, ind::Individual, traitnames::Array{String, 1})
    ind.marked = true
    posspartner = nothing
    communityidxs = patch.whoiswho[ind.lineage]
    startidx = rand(1:length(communityidxs))
    mateidx = startidx
    while true
        mate = patch.community[communityidxs[mateidx]]
        if !mate.marked && iscompatible(mate, ind, traitnames)
            posspartner = mate
            break
        end
        mateidx += 1
        mateidx > length(communityidxs) && (mateidx = 1)
        mateidx == startidx && break
    end
    ind.marked = false
    posspartner 
end

function createtraits(settings::Dict{String, Any}) #TODO: this is all very ugly. (case/switch w/ v. 2.0+?)
    #TODO move traits to defaults.jl -> or the whole function?
    traitnames = settings["traitnames"]
    traits = Trait[]
    # exponential distributions of body sizes:
    seedoffset = settings["maxseedsize"] - settings["minseedsize"]
    repoffset = settings["maxrepsize"] - settings["minrepsize"]
    seedsize = exp(settings["minseedsize"] + seedoffset * rand()) 
    repsize = exp(settings["minrepsize"] + repoffset * rand())
    while repsize <= seedsize
        repsize = exp(settings["minrepsize"] + repoffset * rand())
    end
    for idx in eachindex(traitnames)
        if occursin("rate", traitnames[idx])
            push!(traits, Trait(idx, rand() * 100))
        elseif occursin("dispshape", traitnames[idx])
            push!(traits, Trait(idx, rand() * maxdispmean))
        elseif occursin("tempopt", traitnames[idx])
            push!(traits, Trait(idx, rand() * 25 + 288)) # range 15-40°C
        elseif occursin("temptol", traitnames[idx])
            push!(traits, Trait(idx, (rand() + 0.39) * 5)) #CAVE: code values elsewhere?
        elseif occursin("mut", traitnames[idx])
            mutationrate == 0 ? push!(traits, Trait(idx, rand())) : push!(traits, Trait(idx, mutationrate)) #CAVE: code values elsewhere?
        elseif occursin("repsize", traitnames[idx])
            push!(traits, Trait(idx, repsize)) #CAVE: code values elsewhere?
        elseif occursin("seedsize", traitnames[idx])
            push!(traits, Trait(idx, seedsize)) #CAVE: code values elsewhere?
        elseif occursin("precopt", traitnames[idx])
            push!(traits, Trait(idx, rand() * 10))
        elseif occursin("prectol", traitnames[idx])
            push!(traits, Trait(idx, rand() + 0.39))
        elseif occursin("reptol", traitnames[idx])
            push!(traits, Trait(idx, settings["tolerance"]))
        else
            push!(traits, Trait(idx, rand()))
        end
    end
    traits
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
        ncodinggenes = rand(Poisson(1))
        while ncodinggenes < 1  # make sure every trait is coded by at least 1 gene
            ncodinggenes = rand(Poisson(1))
        end
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

function createchrs(nchrs::Int,genes::Array{AbstractGene,1})
    ngenes=size(genes,1)
    if nchrs>1
        chrsplits = sort(rand(1:ngenes,nchrs-1))
        chromosomes = Chromosome[]
        for chr in 1:nchrs
            if chr==1 # first chromosome
                push!(chromosomes, Chromosome(genes[1:chrsplits[chr]], true))
            elseif chr==nchrs # last chromosome
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):end], true))
            else
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):chrsplits[chr]], true))
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

function createind(settings::Dict{String, Any})
    id = rand(Int32)
    parentid = rand(Int32)
    lineage = randstring(4)
    ngenes = settings["avgnoloci"] * length(settings["traitnames"])
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
    chromosomes = createchrs(nchrms,genes)
    traitdict = chrms2traits(chromosomes, settings["traitnames"])
    if settings["indsize"] == "adult"
        indsize = traitdict["repsize"]
    elseif settings["indsize"] == "seed"
        indsize = traitdict["seedsize"]
    else
        indsize = traitdict["seedsize"] + rand() * traitdict["repsize"] # XXX: sizes shouldn't be uniformally dist'd
    end
    Individual(lineage, chromosomes, traitdict, 0, false, 1.0, indsize, id, parentid)
end

function createpop(settings::Dict{String, Any})
    traits = createtraits(settings)
    if occursin("metabolic", settings["popsize"]) || occursin("single", settings["popsize"])
        # population size determined by adult size and temperature niche optimum
        popsize = round(fertility * traits["repsize"] ^ (-1 / 4) *
                        exp(-act / (boltz * traits["tempopt"])))
    elseif occursin("bodysize", settings["popsize"])
        # population size up to 25% of the maximum possible in this cell
        quarterpopsize = Integer(floor((settings["cellsize"] / traits["repsize"]) / 4))
        popsize = rand(2:quarterpopsize)
    elseif occursin("minimal", settings["popsize"])
        popsize = 2 #Takes two to tangle ;-)
    else
        simlog("Invalid value for `popsize`: $(settings["popsize"])", settings, 'e')
    end
    lineage = randstring(4)
    ngenes = settings["avgnoloci"] * length(settings["traitnames"])
    ngenes < 1 && (ngenes = 1)
    while true
        locivar = rand()
        locivar >= settings["phylconstr"] && break
    end
    population = Individual[]
    for i in 1:popsize
        genes = creategenes(ngenes, traits, settings)
        randchrms = rand(1:length(genes))
        if settings["linkage"] == "none"
            nchrms = length(genes)
        elseif settings["linkage"] == "full"
            nchrms = 1
        else
            nchrms = randchrms
        end
        chromosomes = createchrs(nchrms, genes)
        traitdict = chrms2traits(chromosomes, settings["traitnames"])
        if settings["indsize"] == "adult"
            indsize = traitdict["repsize"]
        elseif settings["indsize"] == "seed"
            indsize = traitdict["seedsize"]
        else
            indsize = traitdict["seedsize"] + rand() * traitdict["repsize"] # XXX: sizes shouldn't be uniformally dist'd
        end
        push!(population, Individual(lineage, chromosomes, traitdict, 0, false, 1.0, indsize, id, parentid))
    end
end

