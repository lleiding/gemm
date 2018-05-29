# Types for GeMM

mutable struct Trait
    name::String
    value::Float64
end

mutable struct Gene
    sequence::Int
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
    isnew::Bool # indicator whether individual is new to a patch (after dispersal or birth)
    fitness::Float64 # rate scaling factor
    size::Float64 # body mass
end

mutable struct Patch
    id::Int
    location::Tuple{Int, Int}
    area::Float64
    altitude::Float64 # corresponds to T
    isisland::Bool
    nichea::Float64 # additional niches,
    nicheb::Float64 # e.g. precipitation
    community::Array{Individual, 1}
    invasible::Bool
    isolated::Bool
    whoiswho::Dict{String, Array{Int, 1}}
end

# constructors:
# XXX These default values should be defined elsewhere (constants.jl, to be precise)
Patch(id,location,area) =
    Patch(id,location,area,298,false,0,0,Individual[], false, false, Dict{String, Array{Int, 1}}())
Patch(id,location) =
    Patch(id,location,2e7,298,false,0,0,Individual[], false, false, Dict{String, Array{Int, 1}}())


