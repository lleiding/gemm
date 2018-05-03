#!/bin/bash

echo "Starting simulation run."
#TODO Test for destination commandline parameter
time julia islandsim.jl --config test.config
#TODO archive results

echo "Done."
