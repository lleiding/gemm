#!/usr/bin/env julia
# A very thin wrapper to start a GeMM simulation. If you need something fancier,
# have a look at `rungemmparallel.jl`, or import the GeMM module into your own
# wrapper script.

thisDir = joinpath(pwd(), "src")
any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)
using Pkg
Pkg.activate(".")
using GeMM

rungemm()
