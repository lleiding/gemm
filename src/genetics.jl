# Functions related to an individual's genome:
# - meiosis
# - trait

"""
    meiosis(genome, maternal)

Carry out meiosis on a genome (marked as maternal or not). Returns a haploid
gamete genome. (genome => array of chromosomes)
"""
function meiosis(genome::Array{Chromosome,1}, maternal::Bool) # TODO: include further dynamics, errors...
    firstset = findall(x -> x.maternal, genome)
    secondset = findall(x -> !x.maternal, genome)
    length(firstset) != length(secondset) && return Chromosome[] # CAVEAT: more elegant solution...
    gameteidxs = Int[]
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
    mate.lineage != ind.lineage && return false ##TODO must be possible under some conditions
    compatidx = findfirst(x -> x == "compat", traitnames)
    indgene = getseq(ind.genome, compatidx)
    mategene = getseq(mate.genome, compatidx)
    seqidentity = getseqsimilarity(indgene, mategene)
    seqidentity >= ind.traits["seqsimilarity"]
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
        elseif occursin("seqsimilarity", traitnames[idx]) && settings["fixtol"]
            push!(traits, Trait(idx, settings["tolerance"]))
        elseif occursin("seqsimilarity", traitnames[idx]) && !settings["fixtol"]
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
