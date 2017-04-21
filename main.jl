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
                offspring = []
                counter=0
                for j in p.community # individuals are sorted according to their phenologies
                    counter+=1
                    println(counter)
                    println(size(p.community,1))
                    j.isnew && evaluate_environment(p,j)
                    germinate(j)
                    mature(j)
                    children=reproduce(j)
                    (children!=nothing) && (append!(offspring,children))
                    disperse(j)
                    (counter>9) && break
                end
                deleteat!(p.community,find(x -> x.dead,p.community))
                append!(p.community,offspring)
                ## /ind cycle
            end
            ## /patch cycle
        end
        ## /island cycle
        (t%10 == 0) && println(t)
    end
    ## /time cycle
end

oneisland=[Island([Patch([Individual([],[],1,"seed",true,5,0.5,0.5,0.5,false)],1,1,1)])]
run_simulation(100, oneisland)
