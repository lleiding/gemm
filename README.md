# GeMM

*Note: this branch (`invasion-old`) represents the state of the code base in
April 2020, when we finished the experiments for the invasion paper.*

This is the Genome explicit Metacommunity Model, an individual-based model for
eco-evolutionary experiments.

This README covers installation and basic usage of the software. For a detailed 
documentation of the biological side of the model, read `ODD.md`. For a full 
documentation of the software, see `docs/index.html`.

## System requirements

The model has successfully been tested on various Linux machines (64 bit), running Arch Linux and Ubuntu with kernel versions between 4.15.0 and 5.3.7.
You need to have `git` and `julia` 1.0.1 installed (https://julialang.org/downloads/), including packages `Distributions` (v"0.16.4"), and `ArgParse` (v"0.6.1").
Newer versions should work as well.

## Installation guide

Download and extract julia from https://julialang.org/downloads/.
Enter the created directory and run

```
bin/julia
```

to launch the REPL.
Press `]` to enter the Pkg REPL and run

```
add Distributions
add ArgParse
```

Press backspace or ^C to get back to the Julia REPL (https://docs.julialang.org/en/v1/stdlib/Pkg/index.html).

Download GeMM by running

```
git clone https://github.com/lleiding/gemm.git
```

Enter the directory and run

```
julia rungemmparallell.jl -c ""
```

to test if it works.

## Demo

To run an example experiment in a landscape with two environmental gradients to test the effect of temporal environmental variation, run

```
julia -p 40 rungemmparallel.jl -s 1 -n 20 -c examples/gradient/sg,examples/gradient/sgv
```

This executes 40 parallel simulations, 20 for each of the scenarios defined in the configuration files (`examples/gradient/sg` and `examples/gradient/sgv`).

Each of the simulations creates its own output directory of the format date-configfilename-replicate.
These folders configuration files with the exact configuration the respective simulation was run on (`*.conf`),
the definition of the simulation arena (`grad.map`),
and a tab-separated table containing statistical summaries on the characteristics of all populations through time (`stats_s*.tsv`).
Optionally, the model outputs several log files (`diversity.log`, `lineages.log`, `simulation.log`).

Depending on the machine (number of parellel processes, i.e., computing cores), the simulations can take up to several days to complete.

## Instructions for use

To design your own experiments, you need to provide two files.

### Map file(s)

The map file contains the information to define a simulation arena.
A single integer value at the top sets the number of years the arena should run for.
Each subsequent whitespace separated line defines one grid cell in the arena.
The first three integer values represent ID, x-location and y-location respectively.
Additional fields can be used to further characterize grid cell.
Options and their value types are

```
temp::Float64
prec::Float64
isisland::Bool
invasible::Bool
initpop::Bool
```

These set the temperature and precipition values in the grid cell, and define whether the cell is an island, 
can be invaded by alien species and whether the grid cell should receive a
community at simulation initialisation.

An experiment can pass through a series of map files to e.g. simulate environmental change or geomorphological dynamics.

### Configuration file

The configuration contains all custom parameters values that may be tweaked to fit one's experimental setup.
Paramaters are defined in a whitespace separated, key-value fashion: `<key> <value>`.
See `src/defaults.jl` for a list of all available parameters and a brief description.
All parameters not set by the user will revert to the values defined in this list of defaults.
At the least, any configuration file should contain a line to define the map file(s):

```
maps <mapfile[,mapfile2,...]>
```

For passing more than one map file, separate file names with commas and no white space, e.g. `map1.map,map2.map`.

### Run an experiment

To parallely run several replicates of a given experiment definition, run

```
julia -p <NPROCS> rungemmparallel.jl -s <SEED> -n <NREPS> -c <CONFIG>
```

while `<NPROCS>` = number of cores, `<SEED>` = integer value to set a random seed, `<NREPS>` = number of replicates and
`<CONFIG>` the configuration file.

