#!/bin/bash

echo "Copying input files."

cp invasions/invasion.config invasions/invasion.map .

echo "Starting simulation run."

date
time julia run_simulation.jl --config invasion.config 
date

echo "Deleting input files."

rm invasion.config invasion.map

echo "Done."
