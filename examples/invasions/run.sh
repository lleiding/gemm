#!/bin/bash

echo "Copying input files."

cp examples/invasions/invasion.config examples/invasions/invasion.map .

echo "Starting simulation run."

date
time julia rungemmparallel.jl --config invasion.config 
date

echo "Deleting input files."

rm invasion.config invasion.map

echo "Done."
