#!/bin/bash

echo "Starting simulation run."

date
#TODO Test for destination commandline parameter
#rm results/tests/* #TODO Remove this again
time julia islandsim.jl --config invasion.config
#TODO archive results
date

echo "Done."
