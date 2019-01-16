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
"""
mutable struct BigGene <: AbstractGene
    sequence::BigInt
    codes::Array{Trait, 1}
end

"""
A struct that saves a collection of genes and a tag to denote whether they are
maternal or not.
"""
mutable struct Chromosome
    genes::Array{AbstractGene, 1}
    maternal::Bool # parental origin of chromosome
end

"""
One of the core structs of the model, representing an individual organism.
"""
mutable struct Individual
    lineage::String
    genome::Array{Chromosome, 1}
    traits::Dict{String, Float64}
    age::Int
    marked::Bool # indicator whether individual is new to a patch (after dispersal or birth)
    precadaption::Float64 # adaption to precipitation
    tempadaption::Float64 # adaption to temperature
    size::Float64 # body mass
    id::Int
end

"""
One of the core structs of the model, representing a one-hectare patch of ground.
"""
mutable struct Patch
    id::Int
    location::Tuple{Int, Int}
    area::Float64
    temp::Float64 # temperature (physiologically important)
    isisland::Bool
    prec::Float64 # precipitation (no physiological effect, just a niche)
    nicheb::Float64 # additional niche
    community::Array{Individual, 1}
    seedbank::Array{Individual, 1}
    invasible::Bool
    isolated::Bool
    initpop::Bool # initialise with a population
    whoiswho::Dict{String, Array{Int, 1}}
end

# constructors:
# TODO These default values should be defined elsewhere (defaults.jl, to be precise)
Patch(id, location, area) =
    Patch(id, location, area, 298, false, 5, 0, Individual[], Individual[],
          false, false, false, Dict{String, Array{Int, 1}}())
Patch(id, location) =
    Patch(id, location, 2e7, 298, false, 5, 0, Individual[], Individual[],
          false, false, false, Dict{String, Array{Int, 1}}())


