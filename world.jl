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
function spawnIsland(world::Array{Float64,2},mainland::Int64,isolation::Int64)
    sizex = size(world,2)
    sizey = size(world,1)
    x = rand((mainland+isolation):sizex)
    y = rand(1:sizey)
    if world[y,x] == 0
        world[y,x] += 1 # TODO: different heights/rule or maybe do it in fact every time, creating more complex islands
    else
        raiseIsland(world,x,y,mainland)
    end
end

function raiseIsland(world::Array{Float64,2},x::Int64,y::Int64,mainland)
    if (x == 0) || (y == 0)
        
end
