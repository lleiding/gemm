# Genome explicit Metacommunity Model - GeMM
#
# Ludwig Leidinger <l.leidinger@gmx.net>,
# Daniel Vedder <daniel.vedder@stud-mail.uni-wuerzburg.de>
# (c) 2017-2020, licensed under the terms of the MIT license
#
# Contains all custom types and functions necessary for an individual based,
# genome and niche explicit, eco-evo-env metacommunity model. 

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

include("initialisation.jl")

include("genetics.jl")

include("reproduction.jl")

include("dispersal.jl")

include("survival.jl")

include("habitatchange.jl")

include("invasion.jl")

include("scheduling.jl")

include("run_simulation.jl")

end
