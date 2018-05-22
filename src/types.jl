# Types for GeMM

mutable struct Trait
    name::String
    value::Float64
end

mutable struct Gene
    sequence::Int64
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
    age::Int64
    isnew::Bool # indicator whether individual is new to a patch (after dispersal or birth)
    fitness::Float64 # rate scaling factor
    size::Float64 # body mass
end

mutable struct Patch
    id::Int64
    location::Tuple{Int64, Int64}
    altitude::Float64 # corresponds to T
    area::Float64
    isisland::Bool
    nichea::Float64 # additional niches,
    nicheb::Float64 # e.g. precipitation
    community::Array{Individual, 1}
    isolated::Bool
end

# constructors:
Patch(id,location,altitude,area,isisland,nichea,nicheb,community) = Patch(id,location,altitude,area,isisland,nichea,nicheb,community,false)
Patch(id,location,altitude,area,isisland,nichea,nicheb) = Patch(id,location,altitude,area,isisland,nichea,nicheb,Individual[],false)
Patch(id,location,altitude,area,isisland,nichea) = Patch(id,location,altitude,area,isisland,nichea,0,Individual[],false)
Patch(id,location,altitude,area,isisland) = Patch(id,location,altitude,area,isisland,0,0,Individual[],false)
Patch(id,location,altitude,area) = Patch(id,location,altitude,area,false,0,0,Individual[],false)
Patch(id,location,altitude) = Patch(id,location,altitude,100,false,0,0,Individual[],false)
Patch(id,location) = Patch(id,location,298,100,false,0,0,Individual[],false)


