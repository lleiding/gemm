# Genome explicit Metacommunity Model - GeMM
#
# Ludwig Leidinger <l.leidinger@gmx.net>, Daniel Vedder <daniel.vedder@stud-mail.uni-wuerzburg.de> 2019
#
# Contains all custom types and functions necessary for an individual based, genome and niche explicit,
# eco-evo-env Metacommunity model. 

module GeMM


using
    ArgParse,
    Dates,
    Distributions,
    Random,
    Statistics


export
    rungemm


include("defaults.jl")

include("entities.jl")

include("input.jl")

include("output.jl")

include("constants.jl")

include("auxfuncts.jl")

include("initialisation.jl")

include("submodels.jl")

include("scheduling.jl")

include("run_simulation.jl")

end
