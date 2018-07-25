# Types for GeMM

mutable struct Trait
    nameindex::UInt8
    value::Float64
end

mutable struct Gene
    sequence::BigInt
    codes::Array{Trait, 1}
end

mutable struct Chromosome
    genes::Array{Gene, 1}
    maternal::Bool # parental origin of chromosome
end

mutable struct Individual
    lineage::String
    genome::Array{Chromosome, 1}
    traits::Dict{String, Float64}
    age::Int
    marked::Bool # indicator whether individual is new to a patch (after dispersal or birth)
    fitness::Float64 # rate scaling factor
    size::Float64 # body mass
    id::Int
    parentid::Tuple{Int, Int}
end

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
    phylo::Array{Int, 2}
end

# constructors:
# XXX These default values should be defined elsewhere (constants.jl, to be precise)
Patch(id,location,area) =
    Patch(id,location,area,298,false,5,0,Individual[], Individual[],
          false, false, false, Dict{String, Array{Int, 1}}(), Array{Int, 2}(0,3))
Patch(id,location) =
    Patch(id,location,2e7,298,false,5,0,Individual[], Individual[],
          false, false, false, Dict{String, Array{Int, 1}}(), Array{Int, 2}(0,3))


