#!/bin/bash

echo "Starting simulation run."

date
time julia run_simulation.jl --config invasions/invasion.config 
date

echo "Done."
