# Genome explicit Metacommunity Model - GeMM
#
# Ludwig Leidinger 2018
# <l.leidinger@gmx.net>
#
# contains all custom types and functions necessary for eco-evo-env speciation island model

module GeMM


using Distributions


export
    # types:
    Patch,
    # functions:
    parseconfig,
    setupdatadir,
    recordcolonizers,
    readmapfile,
    writerawdata,
    writedata,
    checkviability!,
    establish!,
    survive!,
    disturb!,
    grow!,
    disperse!,
    compete!,
    reproduce!,
    mutate!,
    createworld,
    updateworld!

include("types.jl")

include("constants.jl")

include("io.jl")

include("subprocs.jl")

include("init.jl")

include("mainprocs.jl")

end
