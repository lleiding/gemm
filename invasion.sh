#!/bin/bash

echo "Starting simulation run."

date
time julia islandsim.jl --config invasion.config
date

echo "Done."
