#!/usr/bin/env julia

include("types.jl")
using GeneInds

using Distributions

## main time cycle
for t in 1:maxt
    ## cycle land masses
    for i in islands
        ## cycle patches on land masses
        for p in patches[i]
            ## cycle individuals
            for j in community
            end
            ## /ind cycle
        end
        ## /patch cycle
    end
    ## /island cycle
end
## /time cycle

