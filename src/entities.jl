# Types for GeMM

"""
A simple struct to store organism trait names and values.
"""
mutable struct Trait
    nameindex::UInt8
    value::Float64
end

"""
An abstraction layer to hide the implementation details of genes.
"""
abstract type AbstractGene end

"""
The standard gene struct; including the base sequence (represented in binary and
stored and as an integer), and an array of traits the gene codes for.
"""
mutable struct Gene <: AbstractGene
    sequence::Int
    codes::Array{Trait, 1}
end

"""
A variation of the standard Gene struct that allows for longer base sequence
lengths through its use of BigInt.

Used for the compatibility gene when setting("usebiggenes") == true.
"""
mutable struct BigGene <: AbstractGene
    sequence::BigInt
    codes::Array{Trait, 1}
end

"""
This gene type cuts out the conversion of the sequence to an integer. (This
conversion was introduced to save memory, but with low population sizes
(as in Zosterops mode), saving computation cycles may be more important.)

Used when setting("compressgenes") == false.
"""
mutable struct StringGene <: AbstractGene
    sequence::String
    codes::Array{Trait, 1}
end

"""
An abstraction layer to hide the implementation details of chromosomes.
"""
abstract type Chromosome end

"""
A struct that saves a collection of genes and a tag to denote whether they are
maternal or not.
"""
mutable struct DefaultChromosome <: Chromosome
    genes::Array{AbstractGene, 1}
    maternal::Bool # parental origin of chromosome
end

"""
A struct that saves a collection of genes and tags to denote whether they are
maternal or not, as well as which lineage they originated from (to keep track
of heterozygosity in hybridisation studies).
"""
mutable struct LineageChromosome <: Chromosome
    genes::Array{AbstractGene, 1}
    maternal::Bool # parental origin of chromosome
    lineage::String
end

@enum Sex hermaphrodite male female

"""
One of the core structs of the model, representing an individual organism.
"""
mutable struct Individual
    lineage::String
    genome::Array{Chromosome, 1}
    traits::Dict{String, Float64}
    marked::Bool  # indicator whether individual is new to a patch (after dispersal or birth)
    precadaptation::Float64 # adaption to precipitation
    tempadaptation::Float64 # adaption to temperature
    size::Float64 # body mass
    sex::Sex
    partner::Int  # ID of the partner individual (if applicable, default 0)
    id::Int
end

"""
One of the core structs of the model, representing a one-hectare patch of ground.
"""
mutable struct Patch
    id::Int
    location::Tuple{Int, Int}
    capacity::Float64   # biomass carrying capacity in g ("cellsize")
    temp::Float64   # temperature (physiologically important)
    prec::Float64   # precipitation (no physiological effect, just a generic niche)
    nicheb::Float64 # additional generic niche - currently not used
    community::Array{Individual, 1}
    seedbank::Array{Individual, 1}
    isisland::Bool  # island? (if false -> mainland)  
    invasible::Bool # can exotics land here?
    isolated::Bool  # add a distance penalty when dispersing?
    initpop::Bool   # initialise with a population?
end

# constructors:
# XXX These default values should be defined elsewhere (defaults.jl, to be precise)
Patch(id, location, capacity) =
    Patch(id, location, capacity, 298, 5, 0, Individual[], Individual[],
          false, false, false, false)
