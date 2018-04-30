#!/bin/bash

echo "Starting simulation run."
time julia islandsim.jl --config test.config
echo "Done."
