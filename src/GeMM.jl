# Genome explicit Metacommunity Model - GeMM
#
# Ludwig Leidinger 2018
# <l.leidinger@gmx.net>
#
# contains all custom types and functions necessary for eco-evo-env speciation island model

module GeMM


using
    ArgParse,
    Dates,
    Distributions,
    Random,
    Statistics


export
    # types:
    Patch,
    # functions:
    createworld,
    defaultSettings,
    getsettings,
    readmapfile,
    setupdatadir,
    simulate!,
    updateworld!,
    writedata


include("defaults.jl")

include("entities.jl")

include("input.jl")

include("output.jl")

include("constants.jl")

include("auxfuncts.jl")

include("initialisation.jl")

include("submodels.jl")

include("scheduling.jl")


end
