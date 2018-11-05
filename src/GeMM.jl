# Genome explicit Metacommunity Model - GeMM
#
# Ludwig Leidinger 2018
# <l.leidinger@gmx.net>
#
# contains all custom types and functions necessary for eco-evo-env speciation island model

module GeMM


using
    Dates,
    Random,
    Distributions,
    ArgParse


export
    # types:
    Patch,
    # functions:
    defaultSettings,
    parseconfig,
    simlog,
    getsettings,
    setupdatadir,
    recordcolonizers,
    recordstatistics,
    recordlineages,
    readmapfile,
    writephylo,
    writedata,
    checkviability!,
    establish!,
    survive!,
    disturb!,
    invade!,
    grow!,
    disperse!,
    compete!,
    reproduce!,
    mutate!,
    createworld,
    updateworld!

include("defaults.jl")

include("entities.jl")

include("input.jl")

include("output.jl")

include("constants.jl")

include("auxfuncts.jl")

include("initialisation.jl")

include("submodels.jl")

end
