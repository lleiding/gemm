# Subsidiary functions for GeMM

function meiosis(genome::Array{Chromosome,1}, maternal::Bool, mtraits::Dict{String,Float64}, ptraits::Dict{String,Float64}) # TODO: include further dynamics, errors...
    sort!(genome, by = x -> x.number)
    # firstset = find(x->x.maternal,genome) # CAVE: chromosomes could still be mismatched!
    # secondset = find(x->!x.maternal,genome)
    # length(firstset) != length(secondset) && return (genome[firstset], mtraits) # CAVE: more elegant solution...
    # gameteidxs = []
    # for i in eachindex(firstset) # randomly draw indices of maternal and paternal linkage units
    #     push!(gameteidxs,rand([firstset[i],secondset[i]]))
    # end
    gamete = Chromosome[]
    chrmnos = map(x -> x.number, genome)
    for number in unique(chrmnos)
        push!(gamete, rand(filter(x -> x.number == number, genome)))
    end
    traitdict = Dict{String, Float64}() #maybe initialise as copy of mtraits
    for chrm in gamete
        for gene in chrm.genes
            for traitname in gene.codes
                chrm.maternal ? push!(traitdict, traitname => mtraits[traitname]) : push!(traitdict, traitname => ptraits[traitname])
            end
        end
    end
    gamete = deepcopy(gamete)
    for chrm in gamete
        chrm.maternal = maternal
    end
    (gamete, traitdict)
    # mtraitnames = String[]
    # ptraitnames = String[]
    # ## DEBUGGING:
    # alltraits = []
    # allmtraits = []
    # allptraits = []
    # for chrm in genome
    #     for gene in chrm.genes
    #         append!(alltraits, gene.codes)
    #         chrm.maternal && append!(allmtraits, gene.codes)
    #         !chrm.maternal && append!(allptraits, gene.codes)
    #     end
    # end
    # alltraits = unique(alltraits)
    # allmtraits = unique(allmtraits)
    # allptraits = unique(allptraits)
    # ##/DEBUGGING
    # for i in gameteidxs
    #     for gene in genome[i].genes
    #         if genome[i].maternal
    #             append!(mtraitnames, gene.codes) # save all traitnames of genes in linkage unit
    #         else
    #             append!(ptraitnames, gene.codes)
    #         end
    #     end
    #     push!(gamete, Chromosome(genome[i].genes, genome[i].number, maternal))
    # end
    # mtraitnames = unique(mtraitnames)
    # ptraitnames = unique(ptraitnames)
    # length(unique(vcat(mtraitnames, ptraitnames))) != length(mtraits) && println("INCONSISTENCY IN MEIOSIS (alltraits = ", length(alltraits), "(", length(allmtraits), "/", length(allptraits), ")",", haplotraits = ", length(unique(vcat(mtraitnames, ptraitnames))), ")")
    # traitdict = Dict{String, Float64}() #maybe initialise as copy of mtraits
    # for key in keys(mtraits)
    #     in(key, mtraitnames) && (traitdict[key] = mtraits[key])
    # end
    # for key in keys(ptraits)
    #     in(key, ptraitnames) && (traitdict[key] = ptraits[key])
    # end
    # (deepcopy(gamete), traitdict)
end

function chrms2traits(mtraits::Dict{String,Float64}, ptraits::Dict{String,Float64})
    if sort(collect(keys(mtraits))) == sort(collect(keys(ptraits)))
        traitdict = Dict{String,Float64}()
        for traitname in keys(mtraits)
            traitdict[traitname] = mean([mtraits[traitname], ptraits[traitname]])
        end
    else
        length(mtraits) >= length(ptraits) ? traitdict = mtraits : traitdict = ptraits
        warn(STDERR, "parental traits do not match. Returning homozygote.")
    end
    traitdict
end

function traitsexist(ind::Individual, traitnames::Array{String, 1})
    for trait in traitnames
        !haskey(ind.traits, trait) && return false
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
        y = 0.0
    end
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

function iscompatible(mate::Individual, ind::Individual) # TODO: instead of concatenation, choose one compatibility sequence
    length(ind.genome) != length(mate.genome) && return false
    tolerance = ind.traits["reptol"]
    indgene = ""
    for chrm in ind.genome
        for gene in chrm.genes
            any(x -> contains(x, "compat"), gene.codes) && (indgene = gene.sequence)
        end
    end
    indgene == "" && return false
    mategene = ""
    for chrm in mate.genome
        for gene in chrm.genes
            any(x -> contains(x, "compat"), gene.codes) && (mategene = gene.sequence)
        end
    end
    mategene == "" && return false
    basediffs = 0
    for i in eachindex(indgene)
        try
            indgene[i] != mategene[i] && (basediffs += 1) # alternatively use bioinformatics tools
        catch # e.g., in case of differently long genes
            basediffs += 1
        end
    end
    seqidentity = 1 - (basediffs / length(indgene))
    seqidentity < tolerance && return false # CAVE: maybe compare w/ indgenome?
    true
end

function findposspartners(world::Array{Patch,1}, ind::Individual, location::Tuple{Int64, Int64})
    ind.isnew = true
    radius = floor(ind.traits["repradius"] + 0.5) # CAVE: to account for cell width ... or not??
    coordinates = Tuple[]
    for x = -radius:radius, y = -radius:radius
        sqrt(x^2 + y^2) <= radius && push!(coordinates, (x + location[1], y + location[2]))
    end
    coordinates = vcat(map(k -> checkborderconditions!(world, k[1], k[2]), coordinates)...)
    order = map(k -> k.^2, coordinates) |> l -> map(sum, l) |> m -> map(sqrt, m) |> sortperm
    coordinates = coordinates[order]
    posspartners = Individual[]
    idx = 1 # check patches in order of increasing distance
    while length(posspartners) == 0
        idx > length(coordinates) && break
        targetpatch = filter(l -> l.location == coordinates[idx], world)
        length(targetpatch) >= 1 ? community = targetpatch[1].community : community = Individual[]
        shuffle!(community) # randomly select potential mates
        for mate in community
            mate.age == 0 && continue
            mate.lineage != ind.lineage && continue
            !traitsexist(mate, ["repsize"]) && continue
            mate.size < mate.traits["repsize"] && continue
            mate.isnew && continue
            !iscompatible(mate, ind) && continue
            push!(posspartners, mate)
            length(posspartners) > 0 && break # quit searching if one mate is found
        end
        idx += 1
    end
    ind.isnew = false
    posspartners
end

function makeoffspring(noffs::Int64, mate::Individual, partner::Individual)
    seedbank = Individual[]
    for i in 1:noffs # pmap?
        partnergenome, ptraits = meiosis(partner.genome, false, partner.mtraits, partner.ptraits) # offspring have different genome!
        mothergenome, mtraits = meiosis(mate.genome, true, mate.mtraits, mate.ptraits)
        (length(partnergenome) < 1 || length(mothergenome) < 1) && continue
        genome = vcat(partnergenome,mothergenome)
        traits = chrms2traits(mtraits, ptraits)
        age = 0
        isnew = false
        fitness = 0.0
        newsize = traits["seedsize"]
        ind = Individual(mate.lineage, newsize, age, fitness, isnew, genome, traits, mtraits, ptraits)
        push!(seedbank ,ind) # maybe actually deepcopy!?
    end
    seedbank
end

function createtraits(traitnames::Array{String,1}, settings::Dict{String,Any}) #TODO: this is all very ugly. (case/switch w/ v. 2.0+?)
    traits = Dict{String, Float64}()
    seedsize = exp(-7 + 20 * rand()) # corresponds to 1mg to 22kg
    repsize = exp(0 + 17 * rand()) # 1g to 24t
    while repsize <= seedsize
        seedsize = exp(-7 + 20 * rand())
        repsize = exp(0 + 17 * rand())
    end
    for name in traitnames
        if contains(name,"rate")
            traits[name] = rand()*100
        elseif contains(name,"dispshape")
            traits[name] = rand() * maxdispmean
        elseif contains(name, "tempopt")
            traits[name] =  rand() * 40 + 273 #CAVE: code values elsewhere?
        elseif contains(name, "temptol")
            traits[name] = rand() * 5 #CAVE: code values elsewhere?
        elseif contains(name, "mut")
            mutationrate == 0 ? traits[name] = rand() : traits[name] = mutationrate #CAVE: code values elsewhere?
        elseif contains(name, "repsize")
            traits[name] = repsize #CAVE: code values elsewhere?
        elseif contains(name, "seedsize")
            traits[name] = seedsize #CAVE: code values elsewhere?
        elseif contains(name, "precopt")
            traits[name] = rand() * 10
        elseif contains(name, "prectol")
            traits[name] = rand()
        elseif contains(name, "reptol")
            if settings["tolerance"] == "high"
                traits[name] =  0.75
            elseif settings["tolerance"] == "low"
                traits[name] =  0.9 #CAVE: code values elsewhere?
            elseif settings["tolerance"] == "none"
                traits[name] =  0.01 #CAVE: code values elsewhere?
            else
                traits[name] =  0.5 + rand() * 0.5
            end
        else
            traits[name] = rand()
        end
    end
    traits
end

function creategenes(ngenes::Int64,traits::Dict{String, Float64})
    genes = Gene[]
    for i in 1:ngenes
        sequence = String(rand(collect("acgt"), genelength)) # arbitrary start sequence
        codesfor = String[]
        push!(genes,Gene(sequence, codesfor))
    end
    for traitname in collect(keys(traits))
        ncodinggenes = rand(Poisson(1))
        while ncodinggenes < 1  # make sure every trait is coded by at least 1 gene
            ncodinggenes = rand(Poisson(1))
        end
        codinggenes = rand(genes,ncodinggenes)
        for gene in codinggenes
            push!(gene.codes,traitname)
        end
    end
    if !any(map(x -> length(x.codes) == 0, genes)) # make sure there is a neutral gene!
        push!(genes, Gene(String(rand(collect("acgt"), genelength)), String[]))
    end
    genes
end

function createchrs(nchrs::Int64,genes::Array{Gene,1})
    ngenes = length(genes)
    if nchrs > 1
        chrsplits = sort(rand(1:ngenes,nchrs-1))
        chromosomes = Chromosome[]
        for chr in 1:nchrs
            if chr==1 # first chromosome
                push!(chromosomes, Chromosome(genes[1:chrsplits[chr]], chr, true))
            elseif chr==nchrs # last chromosome
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):end], chr, true))
            else
                push!(chromosomes, Chromosome(genes[(chrsplits[chr-1]+1):chrsplits[chr]], chr, true))
            end
        end
    else # only one chromosome
        chromosomes = [Chromosome(genes, 1, true)]
    end
    secondset = deepcopy(chromosomes)
    for chrm in secondset
        chrm.maternal = !chrm.maternal
    end
    append!(chromosomes,secondset)
    chromosomes
end

