# GeMM files & functions

## defaults.jl

- defaultSettings()

## initialisation.jl

- createpop(settings::Dict{String, Any})
- genesis(settings::Dict{String, Any})
- createworld(maptable::Array{Array{String,1},1}, settings::Dict{String, Any})
- updateworld!(world::Array{Patch,1},maptable::Array{Array{String,1},1}, settings::Dict{String, Any})

## input.jl

- getsettings(configfile::String = "", seed::Integer = 0)
- parsecommandline()
- basicparser(filename::String)
- parseconfig(configfilename::String)
- readmapfile(mapfilename::String, settings::Dict{String, Any})

## output.jl

- printheader(settings::Dict{String, Any}, io::IO = stdout, sep::String = "\t")
- dumpinds(world::Array{Patch, 1}, settings::Dict{String, Any}, timestep::Int, io::IO = stdout, sep::String = "\t")
- makefasta(world::Array{Patch, 1}, settings::Dict{String, Any}, io::IO = stdout, onlyisland::Bool = false, sep::String = "_")
- setupdatadir(settings::Dict{String, Any})
- writesettings(settings::Dict{String, Any})
- writedata(world::Array{Patch,1}, settings::Dict{String, Any}, timestep::Int)
- recordstatistics(world::Array{Patch,1}, settings::Dict{String, Any})
- recordlineages(world::Array{Patch,1}, settings::Dict{String, Any}, timestep::Int)
- printpopheader(io::IO)
- printpopstats(io::IO, world::Array{Patch, 1}, settings::Dict{String, Any}, timestep::Integer)
- simlog(msg::String, settings::Dict{String, Any}, category='i', logfile="simulation.log", onlylog=false)
- logprint(msg::String, settings::Dict{String, Any}, tostderr=false)

## run_simulation.jl

- runsim(config::String = "", seed::Integer = 0)
- rungemm(config::String = "", seed::Integer = 0)

## scheduling.jl

- simulate!(world::Array{Patch,1}, settings::Dict{String, Any}, timesteps::Int=1000, timeoffset::Int = 0)

## submodels.jl

- mutate!(traits::Array{Trait, 1}, settings::Dict{String, Any}, locivar::Float64 = 1.0)
- mutate!(ind::Individual, temp::Float64, settings::Dict{String, Any})
- mutate!(patch::Patch, settings::Dict{String, Any})
- mutate!(world::Array{Patch, 1}, settings::Dict{String, Any})
- checkviability!(community::Array{Individual, 1}, settings::Dict{String, Any})
- checkviability!(patch::Patch, settings::Dict{String, Any})
- checkviability!(world::Array{Patch,1}, settings::Dict{String, Any})
- establish!(patch::Patch, nniches::Int=1)
- establish!(world::Array{Patch,1}, nniches::Int=1, static::Bool = true)
- survive!(patch::Patch, mortality::Float64)
- survive!(world::Array{Patch,1}, settings::Dict{String, Any})
- disturb!(patch::Patch, intensity::Int)
- disturb!(world::Array{Patch,1}, settings::Dict{String, Any})
- initspeciespool!(settings::Dict{String, Any})
- invade!(patch::Patch, pressure::Int)
- global invade!(world::Array{Patch,1}, settings::Dict{String, Any})
- grow!(patch::Patch, growthrate::Float64)
- grow!(world::Array{Patch,1}, settings::Dict{String, Any})
- disperse!(world::Array{Patch,1}, static::Bool = true) 
- compete!(patch::Patch)
- compete!(world::Array{Patch,1}, static::Bool = true)
- reproduce!(patch::Patch, settings::Dict{String, Any}) 
- reproduce!(world::Array{Patch,1}, settings::Dict{String, Any})
- changetemp!(world::Array{Patch,1}, sdtemp::Float64)
- changeprec!(world::Array{Patch,1}, sdprec::Float64)
- changehabitat!(world::Array{Patch,1}, settings::Dict{String, Any})

## auxfuncts.jl

- meiosis(genome::Array{Chromosome,1}, maternal::Bool) 
- getmeantraitvalue(traits::Array{Trait, 1}, traitidx::Integer)
- getstdtraitvalue(traits::Array{Trait, 1}, traitidx::Integer)
- gettraitdict(chrms::Array{Chromosome, 1}, traitnames::Array{String, 1})
- gettraitdict(genes::Array{AbstractGene, 1}, traitnames::Array{String, 1})
- gettraitdict(traits::Array{Trait, 1}, traitnames::Array{String, 1})
- traitsexist(traits::Dict{String, Float64}, settings::Dict{String, Any})
- traitsexist(ind::Individual, settings::Dict{String, Any})
- gausscurve(b, c, x, a = 1.0)
- diversity(world::Array{Patch,1})
- shannon(index::Dict{String,Int})
- freespace(world::Array{Patch,1})
- findisland(world::Array{Patch,1})
- checkborderconditions(world::Array{Patch,1}, xdest::Int, ydest::Int)
- identifyAdults!(patch::Patch)
- getseqsimilarity(indgene::AbstractString, mategene::AbstractString)
- getseq(genome::Array{Chromosome, 1}, traitidx::Integer)
- iscompatible(mate::Individual, ind::Individual, traitnames::Array{String, 1})
- findmate(population::AbstractArray{Individual, 1}, ind::Individual, traitnames::Array{String, 1})
- createoffspring(noffs::Integer, ind::Individual, partner::Individual, traitnames::Array{String, 1})
- seq2num(sequence::String)
- seq2bignum(sequence::String)
- num2seq(n::Integer)
- createtraits(settings::Dict{String, Any}) 
- creategenes(ngenes::Int, traits::Array{Trait,1}, settings::Dict{String, Any})
- createchrms(nchrms::Int,genes::Array{AbstractGene,1})
- createind(settings::Dict{String, Any}, marked::Bool = false)
- varyalleles!(genes::Array{AbstractGene, 1}, settings::Dict{String, Any}, locivar::Float64)
- varyalleles!(chrms::Array{Chromosome, 1}, settings::Dict{String, Any}, locivar::Float64)
- markthem!(community::Array{Individual, 1})
- markthem!(habitat::Patch)
- markthem!(world::Array{Patch, 1})
