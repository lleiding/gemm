# Genome explicit Metacommunity Model - GeMM
#
# Ludwig Leidinger 2018
# <l.leidinger@gmx.net>
#
# contains all custom types and functions necessary for eco-evo-env speciation island model

module GeMM


using
    Distributions,
    ArgParse


export
    #constants:
    settings,
    # types:
    Patch,
    # functions:
    simlog,
    getsettings,
    setupdatadir,
    recordcolonizers,
    recordstatistics,
    recordlineages,
    readmapfile,
    writerawdata,
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

include("types.jl")

include("io.jl")

include("constants.jl")

include("subprocs.jl")

include("init.jl")

include("mainprocs.jl")

end
