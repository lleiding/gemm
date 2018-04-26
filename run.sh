#!/bin/bash

echo "Starting simulation run."
julia islandsim.jl --maps invasion.map --linkage "none" --nniches 2 --tolerance "none" --dest tests
echo "Done."
