#!/usr/bin/env julia

include("types.jl")
using GeneInds

using Distributions

function run_simulation(maxt::Int, islands::Array{Island,1})
    ## main time cycle
    for t in 1:maxt
        ## cycle land masses
        for i in islands
            ## cycle patches on land masses
            for p in i.patches
                ## cycle individuals
                for j in p.community # individuals are sorted according to their phenologies
                    germinate(j)
                    mature(j)
                    reproduce(j)
                    disperse(j)
                end
                ## /ind cycle
            end
            ## /patch cycle
        end
        ## /island cycle
        (t%50 == 0) && println(t)
    end
    ## /time cycle
end

oneisland=[Island([Patch([[],[],1,"seed",true],1,1,1)])]
run_simulation(1000, oneisland)
