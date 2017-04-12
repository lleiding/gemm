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
            for j in community # individuals are sorted according to their phenologies
                germinate(community[j])
                mature(community[j])
                reproduce(community[j])
                disperse(community[j])
            end
            ## /ind cycle
        end
        ## /patch cycle
    end
    ## /island cycle
end
## /time cycle

