# Subsidiary functions for GeMM

function meiosis(genome::Array{Chromosome,1},maternal::Bool) # TODO: include further dynamics, errors...
    firstset = find(x->x.maternal,genome)
    secondset = find(x->!x.maternal,genome)
    size(firstset,1) != size(secondset,1) && return Chromosome[] # CAVE: more elegant solution...
    gameteidxs = []
    for i in eachindex(firstset)
        push!(gameteidxs,rand([firstset[i],secondset[i]]))
    end
    gamete = Chromosome[]
    for i in gameteidxs
        push!(gamete, Chromosome(genome[i].genes, maternal))
    end
    deepcopy(gamete)
end

function chrms2traits(chrms::Array{Chromosome, 1}, traitnames::Array{String, 1})
    genes = Gene[]
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
            simlog("Missing trait \"", trait, "\". Individual might be killed.", 'w')
            return false
        end
    end
    true
end

function traitsexist(ind::Individual, traitnames::Array{String, 1})
    # FIXME If a trait doesn't exist, we need an error
    for trait in traitnames
        if !haskey(ind.traits, trait)
            simlog("Individual is missing trait \"", trait, "\". Might be killed.", 'w')
            return false
        end
    end
    true
end

function traitsexist(traits::Dict{String, Float64}, traitname::String)
    if !haskey(traits, traitname)
        simlog("Missing trait \"", traitname, "\". Individual might be killed.", 'w')
        return false
    end
    true
end

function traitsexist(ind::Individual, traitname::String)
    if !haskey(ind.traits, traitname)
        simlog("Individual is missing trait \"", traitname, "\". Might be killed.", 'w')
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
    round((space/length(world))/1e6, 3)
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
function checkborderconditions!(world::Array{Patch,1},xdest::Float64,ydest::Float64)
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
    xdest,ydest
end

function identifyAdults!(patch::Patch)
    adultspeciesidx = Dict{String, Array{Int, 1}}()
    for i in eachindex(patch.community)
        patch.community[i].size < patch.community[i].traits["repsize"] && continue
        adultspeciesidx = merge(append!, adultspeciesidx, Dict(patch.community[i].lineage => [i]))
    end
    patch.whoiswho = adultspeciesidx
end

function iscompatible(mate::Individual, ind::Individual)
    compatidx = findin(settings["traitnames"], ["compat"])[1]
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

function findposspartners(world::Array{Patch,1}, ind::Individual, location::Tuple{Int, Int})
    # TODO This should be rewritten so that it only returns a single individual, not an array
    # (which only ever contains a single organism anyway)
    ind.isnew = true
    radius = floor(ind.traits["repradius"] + 0.5) # CAVE: to account for cell width ... or not??
    coordinates = Tuple[]
    for x = -radius:radius, y = -radius:radius
        sqrt(x^2 + y^2) <= radius && push!(coordinates, (x + location[1], y + location[2]))
    end
    coordinates = map(k -> checkborderconditions!(world, k[1], k[2]), coordinates)
    order = map(x -> x.^2, coordinates) |> x -> map(sum, x) |> x -> map(sqrt, x) |> sortperm
    coordinates = coordinates[order]
    posspartners = Individual[]
    idx = 1 # check patches in order of increasing distance
    while length(posspartners) == 0
        idx > length(coordinates) && break
        targetpatch = filter(l -> l.location == coordinates[idx], world)
        if length(targetpatch) >= 1 && haskey(targetpatch[1].whoiswho, ind.lineage)
            communityidxs = (targetpatch[1].whoiswho[ind.lineage])
        else
            communityidxs = []
        end
        shuffle!(communityidxs)
        for mateidx in communityidxs
            mateidx > length(targetpatch[1].community) && continue #XXX Not a real fix, but should work
            mate = targetpatch[1].community[mateidx] #FIXME Occasionally, this throws a bounds error
            mate.isnew && continue
            !iscompatible(mate, ind) && continue
            push!(posspartners, mate)
            length(posspartners) >= 1 && break
        end
        idx += 1
    end
    ind.isnew = false
    posspartners
end

function createtraits() #TODO: this is all very ugly. (case/switch w/ v. 2.0+?)
    #TODO move traits to constants.jl
    traitnames = settings["traitnames"]
    traits = Trait[]
    #seedsize = exp(-7 + 17 * rand()) # corresponds to 1mg to 22kg
    #repsize = exp(0 + 17 * rand()) # 1g to 24t
    # Larger sizes to reduce community size
    #seedsize = exp(-6 + 16 * rand()) # corresponds to 25mg to 22kg
    #repsize = exp(3 + 14 * rand()) # 20g to 24t
    seedsize = exp(-4 + 14 * rand()) # corresponds to 180mg to 22kg
    repsize = exp(5 + 12 * rand()) # 150g to 24t
    while repsize <= seedsize
        repsize = exp(5 + 12 * rand()) # 20g to 24t
        #seedsize = exp(-7 + 17 * rand())
        #repsize = exp(0 + 17 * rand())
    end
    for idx in eachindex(traitnames)
        if contains(traitnames[idx], "rate")
            push!(traits, Trait(idx, rand() * 100))
        elseif contains(traitnames[idx], "dispshape")
            push!(traits, Trait(idx, rand() * maxdispmean))
        elseif contains(traitnames[idx], "tempopt")
            #push!(traits, Trait(idx, rand() * 40 + 273)) #CAVE: code values elsewhere?
            push!(traits, Trait(idx, rand() * 25 + 288)) # range 15-40°C
        elseif contains(traitnames[idx], "temptol")
            push!(traits, Trait(idx, (rand() + 0.39) * 5)) #CAVE: code values elsewhere?
        elseif contains(traitnames[idx], "mut")
            mutationrate == 0 ? push!(traits, Trait(idx, rand())) : push!(traits, Trait(idx, mutationrate)) #CAVE: code values elsewhere?
        elseif contains(traitnames[idx], "repsize")
            push!(traits, Trait(idx, repsize)) #CAVE: code values elsewhere?
        elseif contains(traitnames[idx], "seedsize")
            push!(traits, Trait(idx, seedsize)) #CAVE: code values elsewhere?
        elseif contains(traitnames[idx], "precopt")
            push!(traits, Trait(idx, rand() * 10))
        elseif contains(traitnames[idx], "prectol")
            push!(traits, Trait(idx, rand() + 0.39))
        elseif contains(traitnames[idx], "reptol")
            if settings["tolerance"] == "high"
                push!(traits, Trait(idx, 0.75))
            elseif settings["tolerance"] == "low"
                push!(traits, Trait(idx, 0.9)) #CAVE: code values elsewhere?
            elseif settings["tolerance"] == "none"
                push!(traits, Trait(idx, 0.01)) #CAVE: code values elsewhere?
            else
                push!(traits, Trait(idx, 0.5 + rand() * 0.5))
            end
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
        binary *= bin(search(bases, base) - 1, 2)
    end
    parse(Int64, binary, 2) # Int64 = max sequence length ~31
end

"""
    num2seq(n)
Convert an integer into binary and then into a DNA base sequence string.
"""
function num2seq(n::Int)
    bases = "acgt"
    binary = bin(n, genelength*2)
    sequence = ""
    for i in 1:2:(length(binary)-1)
        sequence *= string(bases[parse(Int64, binary[i:(i+1)], 2)+1])
    end
    sequence
end

function creategenes(ngenes::Int,traits::Array{Trait,1})
    genes = Gene[]
    for i in 1:ngenes
        sequence = String(rand(collect("acgt"), genelength)) # arbitrary start sequence
        seqint = seq2num(sequence)
        codesfor = Trait[]
        push!(genes,Gene(seqint, codesfor))
    end
    for trait in traits
        ncodinggenes = rand(Poisson(1))
        while ncodinggenes < 1  # make sure every trait is coded by at least 1 gene
            ncodinggenes = rand(Poisson(1))
        end
        codinggenes = rand(genes,ncodinggenes)
        for gene in codinggenes
            push!(gene.codes,trait)
        end
    end
    if !any(map(x -> length(x.codes) == 0, genes)) # make sure there is a neutral gene!
        push!(genes, Gene(seq2num(String(rand(collect("acgt"), genelength))), Trait[]))
    end
    genes
end

function createchrs(nchrs::Int,genes::Array{Gene,1})
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

function createind()
    lineage = randstring(4)
    meangenes = length(settings["traitnames"])
    ngenes = rand(Poisson(meangenes))
    ngenes < 1 && (ngenes = 1)
    traits = createtraits()
    genes = creategenes(ngenes,traits)
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
    settings["initadults"] ? indsize = traitdict["repsize"] : indsize = traitdict["seedsize"]
    Individual(lineage, chromosomes, traitdict, 0, false, 1.0, indsize)
end

