#!/usr/bin/env julia
## Module to create and update arenas comprising a static mainland and changing islands
##
## Ludwig Leidinger <l.leidinger@gmx.net>
## 11/08/2017

## TODO: different data type for discriminating mainland and islands

module Archipelago

function createWorld(sizex::Int64,sizey::Int64,mainland::Int64,isolation::Int64)
    world = zeros(Float64,sizey,sizex)
    world[:,1:mainland] = 1 # TODO: sensible rule
end

## create and raise islands
function updateWorld(world::Array{Float64,2},mainland::Int64,isolation::Int64)
    sizex = size(world,2)
    sizey = size(world,1)
    x = rand((mainland+isolation):sizex)
    y = rand(1:sizey)
    world[y,x] += 1 # TODO: different heights/rule
end

end
