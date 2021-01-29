**This repository has moved! The new location is now https://github.com/CCTB-Ecomods/gemm**

-----

# GeMM

This is the Genome explicit Metacommunity Model, an individual-based model for
eco-evolutionary experiments.

This README covers installation and basic usage of the software. For a detailed 
documentation of the biological side of the model, read `ODD.md`. For a full 
documentation of the software, see `docs/index.html`.

## System requirements

The model has successfully been tested on various Linux machines (64 bit), running Arch Linux and Ubuntu with kernel versions between 4.15.0 and 5.3.7.
You need to have `git` and `julia` (>=1.5) installed (https://julialang.org/downloads/), including packages `Distributions` (v"0.16.4"), and `ArgParse` (v"0.6.1").
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
julia rungemm.jl
```

to test if it works.

## Instructions for use

To design your own experiments, you need to provide a map and a configuration file. 
(See the folders in `examples/` for our previous work.)

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

**Example map file:**

```
## SAMPLE EXPERIMENT MAP

3 # Timesteps

# Simulation arena:
# <id> <x> <y> <temperature> <precipitation> [parameters]
1	1	1	temp=293	prec=12 initpop
2	2	1	temp=293	prec=59 initpop
3	3	1	temp=293	prec=89 initpop
4	4	1	temp=293	prec=91 initpop
5	5	1	temp=293	prec=58 initpop
6	6	1	temp=293	prec=57 initpop
7	7	1	temp=293	prec=63 initpop
8	8	1	temp=293	prec=44 initpop
9	9	1	temp=293	prec=28 initpop
```

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

**Example config file:**

```
# Sample configuration file for GeMM
# Note: this does not include all settings!
# For a complete parameter list, see `src/defaults.jl`

# input/output settings
seed     0 #random seed
maps     sample.map
dest     results/sample
outfreq  1
logging  true
debug    true
stats    true
lineages true
fasta    off
raw      false

# general model parameters
linkage       none
nniches       2
static        false
mutate        true
usebiggenes   false
compressgenes true
indsize       adult
popsize       metabolic
maxbreadth    5.0
capgrowth     true
```

### Run an experiment

To run a single experiment, execute:

```
./rungemm.jl -c <CONFIG>
```

Other commandline arguments are:

```
usage: rungemm.jl [-s SEED] [-m MAPS] [-c CONFIG] [-d DEST] [--debug]
                  [--quiet] [-h]

optional arguments:
  -s, --seed SEED      inital random seed (type: Int64)
  -m, --maps MAPS      list of map files, comma separated
  -c, --config CONFIG  name of the config file
  -d, --dest DEST      output directory. Defaults to current date
  --debug              debug mode. Turns on output of debug
                       statements.
  --quiet              quiet mode. Don't print output to screen.
  -h, --help           show this help message and exit
```

There's a separate script to run several replicates of a given experiment
in parallel. This is used as follows:

```
julia -p <NPROCS> rungemmparallel.jl -s <SEED> -n <NREPS> -c <CONFIG>
```

where `<NPROCS>` = number of cores, `<SEED>` = integer value to set a random seed, `<NREPS>` = number of replicates and
`<CONFIG>` the configuration file.


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
