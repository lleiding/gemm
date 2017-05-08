#!/usr/bin/env julia

include("types.jl")
using GeneInds: Island, Patch, Individual, evalenv, germinate, mature,
    reproduce, disperse

using Distributions

function run_simulation(maxt::Int, islands::Array{Island,1})
    ## main time cycle
    for t in 1:maxt
        ## cycle land masses
        for i in islands
            ## cycle patches on land masses
            for p in i.patches
                ## cycle individuals
                offspring = []
                counter = 0 # testing
                println(size(p.community, 1)) # testing
                for j in eachindex(p.community) # individuals are sorted according to their phenologies
                    counter += 1 # for testing
                    p.community[j].isnew && evalenv(p, p.community[j])
                    germinate(p.community[j])
                    mature(p.community[j])
                    children = reproduce(p.community[j])
                    (children!=nothing) && (append!(offspring, children))
                    disperse(p.community[j])
                    (counter>20) && (p.community[j].dead = true) # for testing
                end
                deleteat!(p.community, find(x -> x.dead, p.community))
                append!(p.community, offspring)
                ## /ind cycle
            end
            ## /patch cycle
        end
        ## /island cycle
        (t%10 == 0) && println("t=", t, " ##################")
    end
    ## /time cycle
end

oneisland = [Island([Patch([Individual([], [], 1, "seed", true, 5, 0.3, 0.3, 0.1, false)], 1, 1, 1)])]
run_simulation(100, oneisland)
