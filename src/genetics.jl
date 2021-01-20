# Functions related to an individual's genome:
# - meiosis
# - mutation
# - genome creation
# - auxiliary functions

"""
    meiosis(genome, maternal)

Carry out meiosis on a genome (marked as maternal or not). Returns a haploid
gamete genome. (genome => array of chromosomes)
"""
function meiosis(genome::Array{Chromosome,1}, maternal::Bool)
    gametelength = Int(length(genome)/2)
    gamete = Chromosome[]
    m, p = 1, 1
    while length(gamete) < gametelength
        # find the next maternal/paternal chromosome
        while !genome[m].maternal
            m += 1
        end
        while genome[p].maternal
            p += 1
        end
        # then choose one at random and use it to create a new chromosome for the gamete
        push!(gamete, Chromosome(genome[rand([m,p])].genes, maternal))
        m += 1; p += 1
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
    #XXX can this be made more efficient? It's called really often...
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
    traitdict["ngenes"] = ngenes
    traitdict["nlnkgunits"] = nchrms
    traitdict
end

"""
    gettraitdict(traits, traitnames)

Construct a trait dict from a list of Trait objects.
"""
function gettraitdict(traits::Array{Trait, 1}, traitnames::Array{String, 1})
    traitdict = Dict{String, Float64}()
    for traitidx in unique(map(x -> x.nameindex, traits))
        traitdict[traitnames[traitidx]] = getmeantraitvalue(traits, traitidx)
    end
    traitdict
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
function getseq(genome::Array{Chromosome, 1}, traitidx::Integer, compressgenes::Bool)
    seq = ""
    for chrm in genome
        for gene in chrm.genes
            # use one compatibility gene randomly #XXX what's that comment supposed to mean?
            #TODO have another look at this
            if any(x -> x.nameindex == traitidx, gene.codes)
                compressgenes ? seq = num2seq(gene.sequence) : seq = gene.sequence
            end
        end
    end
    seq
end

"""
    iscompatible(mate, individual)

Check to see whether two individual organisms are reproductively compatible.
"""
function iscompatible(mate::Individual, ind::Individual)
    mate.lineage != ind.lineage && return false ##TODO must be possible under some conditions
    compatidx = findfirst(x -> x == "compat", setting("traitnames"))
    indgene = getseq(ind.genome, compatidx, setting("compressgenes"))
    mategene = getseq(mate.genome, compatidx, setting("compressgenes"))
    seqidentity = getseqsimilarity(indgene, mategene)
    seqidentity >= ind.traits["seqsimilarity"]
end

"""
    seq2num(sequence)

Convert a DNA base sequence (a string) into binary and then into an integer.
This saves memory.
"""
#XXX Actually, this hardly seems to make a difference :-( So using `compressgenes = false`
# seems to be a pretty good choice right now. What we should really do is create gene
# sequences as int right away and save ourselves the intermediate string allocation.
function seq2num(sequence::String)
    # This function effectively uses the same technique as seq2bignum, but it skips the
    # intermediate allocations and is therefore more efficient.
    num::Int64 = 0  # Int64 allows for max length of 21bp
    for b in eachindex(sequence)
        if sequence[end+1-b] == 'a'
            num += 2^(3*(b-1)) * 4 # b'100'
        elseif sequence[end+1-b] == 'c'
            num += 2^(3*(b-1)) * 5 # b'101'
        elseif sequence[end+1-b] == 'g'
            num += 2^(3*(b-1)) * 6 # b'110'
        elseif sequence[end+1-b] == 't'
            num += 2^(3*(b-1)) * 7 # b'111'
        end
    end
    num
end

function intseq(n::Int)
    #XXX doesn't seem to work either? I'm starting to think I'm optimising in the wrong place :-(
    num::Int64 = 0  # Int64 allows for max length of 21bp
    for b in 0:(n-1)
        # basically like seq2bignum(), but skips all intermediate allocations
        num += 2^(3*b) * rand(4:7)
    end
    num
end

"""
    seq2bignum(sequence)

Convert a DNA base sequence (a string) into binary and then into a BigInt (for
larger genes). This saves memory.
"""
function seq2bignum(sequence::String)
    bases = "acgt"
    binary = ""
    for base in sequence
        binary *= string(findfirst(x -> x == base, bases) + 3, base = 2)
    end
    parse(BigInt, binary, base = 2) # BigInt allows arbitrary-length sequences
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
    createtraits()

Create an array of trait objects generated from the default trait values (with a
random offset).
"""
function createtraits()
    #TODO: this is all very ugly. (case/switch w/ v. 2.0+?)
    traitnames = setting("traitnames")
    traits = Trait[]
    # exponential distributions of body sizes:
    repoffset = setting("maxrepsize") - setting("minrepsize")
    seedoffset = setting("maxseedsize") - setting("minseedsize")
    tempoffset = setting("maxtemp") - setting("mintemp")
    sizes = Vector{Float64}(undef, 2)
    while true
        sizes[1] = exp(setting("minrepsize") + repoffset * rand())
        sizes[2] = exp(setting("minseedsize") + seedoffset * rand())
        sizes[1] > sizes[2] && break
    end
    repsize, seedsize = sizes
    for idx in eachindex(traitnames)
        if occursin("dispshape", traitnames[idx])
            push!(traits, Trait(idx, rand() * setting("dispshape")))
        elseif occursin("dispmean", traitnames[idx])
            push!(traits, Trait(idx, rand() * setting("dispmean")))
        elseif occursin("precopt", traitnames[idx])
            push!(traits, Trait(idx, rand() * setting("precrange")))
        elseif occursin("prectol", traitnames[idx])
            push!(traits, Trait(idx, rand() * setting("maxbreadth")))
        elseif occursin("repsize", traitnames[idx])
            push!(traits, Trait(idx, repsize))
        elseif occursin("seqsimilarity", traitnames[idx]) && setting("fixtol")
            push!(traits, Trait(idx, setting("tolerance")))
        elseif occursin("seqsimilarity", traitnames[idx]) && !setting("fixtol")
            push!(traits, Trait(idx, rand())) # assortative mating might evolve
        elseif occursin("seedsize", traitnames[idx])
            push!(traits, Trait(idx, seedsize))
        elseif occursin("tempopt", traitnames[idx])
            push!(traits, Trait(idx, setting("mintemp") + rand() * tempoffset))
        elseif occursin("temptol", traitnames[idx])
            push!(traits, Trait(idx, rand() * setting("maxbreadth")))
        else
            push!(traits, Trait(idx, rand()))
        end
    end
    traits
end

"""
    creategenes(ngenes, traits)

Randomly create a given number of gene objects, with their base sequence and
associated traits. Returns the result as an array of AbstractGenes.
"""
function creategenes(ngenes::Int, traits::Array{Trait,1})
    genes = AbstractGene[]
    bases = ['a','c','g','t']
    compatidx = findfirst(x -> x == "compat", setting("traitnames"))
    # initialise each gene with an arbitrary sequence
    for i in 1:ngenes
        if setting("compressgenes") #default
            push!(genes, Gene(intseq(setting("smallgenelength")), Trait[]))
        else
            sequence = String(rand(bases, setting("smallgenelength")))
            push!(genes, StringGene(sequence, Trait[]))
        end
    end
    # assign traits to the genes (allows for pleiotropy as well as polygenic inheritance)
    for trait in traits
        (trait.nameindex == compatidx) && continue # the compatibility gene is treated separately
        if setting("degpleiotropy") == 0
            # Disable polygenic inheritance and pleiotropy: make sure one gene codes for one trait.
            for gene in genes
                if isempty(gene.codes)
                    push!(gene.codes, trait)
                    break
                end
            end
        else
            # We're actually setting polygenic inheritance here, rather than pleiotropy.
            # Pleiotropy is introduced indirectly, because a higher number of coding genes
            # increases the likelihood of picking a gene that already codes for other traits.
            ncodinggenes = rand(Geometric(1 - setting("degpleiotropy"))) + 1
            codinggenes = rand(genes,ncodinggenes)
            for gene in codinggenes
                push!(gene.codes,trait)
            end
        end
    end
    # append a gene that will be used to determine mating compatibility
    # Note: we only need big genes for the compatibility gene, because we need a longer base
    # sequence than offered by `smallgenelength` if we want to do phylogenetic analyses.
    setting("usebiggenes") ? seql = setting("biggenelength") : seql = setting("smallgenelength")
    cseq = String(rand(bases, seql))
    ctrait = [Trait(compatidx, 0.5)]
    if setting("compressgenes")
        if setting("usebiggenes")
            push!(genes, BigGene(seq2bignum(cseq), ctrait))
        else
            push!(genes, Gene(seq2num(cseq), ctrait))
        end
    else
        push!(genes, StringGene(cseq, ctrait))
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
    varyalleles!(genes, locivar)

Mutate gene traits in the passed array of genes.
"""
function varyalleles!(genes::Array{AbstractGene, 1}, locivar::Float64)
    locivar == 0 && return
    for gene in genes
        mutate!(gene.codes, locivar)
    end
end

"""
    varyalleles!(chromosomes, locivar)

Mutate gene traits in the passed array of chromosomes.
"""
function varyalleles!(chrms::Array{Chromosome, 1}, locivar::Float64)
    locivar == 0 && return
    for chrm in chrms
        varyalleles!(chrm.genes, locivar)
    end
end

"""
    mutate!(traits, locivar)

Loop over an array of traits, mutating each value in place along a normal distribution.
`locivar` can be used to scale the variance of the normal distribution used to draw new
trait values (together with `setting(phylconstr]`).
"""
function mutate!(traits::Array{Trait, 1}, locivar::Float64 = 1.0)
    setting("phylconstr") * locivar == 0 && return
    for trait in traits
        traitname = setting("traitnames")[trait.nameindex]
        occursin("seqsimilarity", traitname) && setting("fixtol") && continue
        oldvalue = trait.value
        occursin("tempopt", traitname) && (oldvalue -= 273)
        while oldvalue <= 0 # make sure sd of Normal dist != 0
            oldvalue = abs(rand(Normal(0,0.01)))
        end
        newvalue = rand(Normal(oldvalue, oldvalue * setting("phylconstr") * locivar))
        (newvalue > 1 && occursin("prob", traitname)) && (newvalue=1.0)
        while newvalue <= 0
            newvalue = rand(Normal(oldvalue, oldvalue * setting("phylconstr") * locivar))
        end
        occursin("tempopt", traitname) && (newvalue += 273)
        trait.value = newvalue
    end
end

"""
    mutate!(individual, temp)

Mutate an individual's genome (sequence and traits) in place.
"""
function mutate!(ind::Individual, temp::Float64)
    muts = setting("mutationrate") * exp(-act/(boltz*temp))
    nmuts = rand(Poisson(muts))
    nmuts == 0 && return
    chrmidcs = rand(eachindex(ind.genome), nmuts)
    for c in chrmidcs
        ind.genome[c] = deepcopy(ind.genome[c])
        length(ind.genome[c].genes) == 0 && continue
        g = rand(eachindex(ind.genome[c].genes))
        gseq = ind.genome[c].genes[g].sequence
        setting("compressgenes") ? charseq = collect(num2seq(gseq)) : charseq = collect(gseq)
        i = rand(eachindex(charseq))
        newbase = rand(collect("acgt"),1)[1]
        while newbase == charseq[i]
            newbase = rand(collect("acgt"),1)[1]
        end
        charseq[i] = newbase
        mutate!(ind.genome[c].genes[g].codes)
        ind.genome[c].genes[g].sequence = deepcopy(ind.genome[c].genes[g].sequence)
        if setting("compressgenes")
            if length(charseq) > 21
                ind.genome[c].genes[g].sequence = seq2bignum(String(charseq))
            else
                ind.genome[c].genes[g].sequence = seq2num(String(charseq))
            end
        else
            ind.genome[c].genes[g].sequence = String(charseq)
        end
    end
    ind.traits = gettraitdict(ind.genome, setting("traitnames"))
end

"""
    mutate!(patch, setting)

Mutate all seed individuals in a patch.
"""
function mutate!(patch::Patch)
    for ind in patch.seedbank
        mutate!(ind, patch.temp)
    end
end

"""
    mutate!(world)

Mutate the world. (That sounds scary!)
"""
function mutate!(world::Array{Patch, 1})
    for patch in world
        (patch.isisland || !setting("static")) && mutate!(patch)
    end
end
