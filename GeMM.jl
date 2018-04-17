# Genome explicit Metacommunity Model - GeMM
#
# Ludwig Leidinger 2017
# <l.leidinger@gmx.net>
#
# contains all custom types and functions necessary for eco-evo-env speciation island model

module GeMM


using Distributions#, Plots


export Patch, # types
    setupdatadir,
    recordcolonizers,
    readmapfile,
    writerawdata,
    writedata,
    analysis,
    checkviability!,
    establish!,
    survive!,
    grow!,
    disperse!,
    compete!,
    reproduce!,
    mutate!,
    createworld,
    updateworld! # functions

include("types.jl")

include("constants.jl")

include("io.jl")

include("init.jl")

include("mainprocs.jl")

include("subprocs.jl")

end
