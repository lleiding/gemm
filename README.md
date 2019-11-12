# GeMM

This is the Genome explicit Metacommunity Model

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
julia -p 40 rungemmparallel.jl -s 1 -n 20 -c examples/gradient/static.conf,examples/gradient/variable.conf
```

This executes 40 parallel simulations, 20 for each of the scenarios defined in the configuration files (`examples/gradient/static.conf` and `examples/gradient/variable.conf`).

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

## Model description

The  model  description follows  the  ODD  (Overview,  Design  concepts, Details)  protocol  (Grimm  et  al.,  2006; 2010:
https://doi.org/10.1016/j.ecolmodel.2006.04.023 https://doi.org/10.1016/j.ecolmodel.2010.08.019)

### 1. Purpose

This model is designed to simulate a (meta-)community of plant-like individuals.
For this, the model considers factors and processes across genetic, population and ecological levels.
The model is able to produce several patterns across genetic, individual, population and (meta-)community levels,
including adaptation and speciation through divergence of populations.
Thus, the model expands from basic principles to richer representation of real-world scenarios.

### 2. Entities, state variables and scales

Individuals (plant-like) are the basic entity in the model.
Each individual carries a diploid set of one or more linkage units, which in turn are comprised of genes.
Linkage units are always inherited in their entirety during the recombination phase of a reproduction event.
The higher the number of genes per linkage unit, the higher the degree of genetic linkage.
Some of the genes code for one ore more traits (pleiotropy), while a trait can be dependent on more than one gene (polygene). 
The realized trait value is the mean of all the trait alleles (quantitative trait loci).
Traits thus controlled encompass
the initial body mass (size) of offspring, M_s,
the body mass determining onset of maturity and thus reproductive capability, M_r,
mean dispersal distance, D_mean,
the shape of the dispersal kernel, controlling long-distance-dispersal, D_s,
the threshold of sequence similarity between mates determining compatibility, delta_I,
and values representing the optimum and the tolerance (standard deviation) of a physical niche parameter, such as temeperature and precipitation (T_mean and sigma_T or P_mean and sigma_P, resp.).
Alternatively to be controlled by mutable genes, traits can also be fixed.
Additionally, individuals carry attributes which describe their bodymass, M, and their adaptation to the abiotic environmental conditions (fitness), F_T and F_P.
Furthermore, every individual carries a Boolean marker used to store whether a given individual has newly arrived to a grid cell or discriminate individuals from the rest of the community.

The base rates for processes governed by the metabolic theory of ecology (Brown et al. 2004: https://doi.org/10.1890/03-9000) - growth, reproduction, mortality -
are global constants. Mutation rate is also a global constant.

Every individual is placed inside an arena of grid cells, each of which has a unique location (coordinates)
and is characterized by physical properties such as temperature, precipiation and size (carrying capacity).
Over the course of the simulation these properties (location or physical parameters) might change, reflecting
geomorphological dynamics.
All individuals within one grid cell constitute a community.
The characteristics of the grid cells combined with the state of inhabiting individuals constitute the state variables of the model.
Additional patterns or summary statistics may be calculated based on these individual information.

Processes and updates are repeated every timestep, while each timestep can be considered as one year.

### 3. Process overview and scheduling

In each discrete timestep each individual in each grid cell will (in no particular order unless otherwise stated) undergoe
the following processes:
- (1) establishment,
- (2) density independent mortality influenced by adaptation to temperature,
- [disturbance]
- (4) growth,
- (5) competition (individuals are sorted according to their adaptation to precipitation),
- (6) reproduction
- (7) mutation of offspring,
- (8) filtering of unviable individuals,
- [invasion]
- (9) seed dispersal.

Afterwards the physical habitat of a grid cell might change.
If that happens, all individuals within that cell are marked to undergo establishment again.

Updates to individuals and thus the local communities happen instantaneously after a specific process has been executed
(asynchronous updating).

### 4. Design concepts

#### Basic principles.

Metabolic theory of ecology (@submodel level).
Adaptive and non-adaptive radiation/evolution. Submodel level[, but geomorphological change @system level.]
Sexual reproduction.
Niche theory. Both at system and submodel level. Each individual carries unique ecological and/or functional traits,
as well as preferences for their physical environments.
Resource/energy limitation (carrying capacity). System level property, but invoked at submodel level.

#### Emergence.

Species/Populations/ecotypes. Only constrained via genetic properties.
Community trait composition. Interplay of physical properties (environment, geographical properties) and within community
(competition strength via reproduction, growth, etc.).
Species numbers, endemics, speciation rate.

#### Adaption. (+Objectives.)

See entities. Traits follow evolution: a trait changes its value randomly within a given phylogenetic constraint.
The success of the change (fitness) emerges as the result of adaption to the physical environment and the reproductive success of an individual over its competitors.

#### Learning.

#### Prediction.

#### Sensing.

Individuals are directly affected by the properties of their physical environment (e.g., temperature).

#### Interaction.

Individuals directly interact when sexually reproducing. However, they are not affected themself by this interaction.
Instead, the interaction aims solely at determining the genotype of their offspring.
Additionally, competition for resource/energy/space between individuals represents indirect interaction.

#### Stochasticity.

Most of the submodels are carried out by all elegible individuals.
Some submodels (Survival, Competition, Mutation and Dispersal), however, happen with particular probabilities.
In these cases, execution of submodels is decided at random, taking into account individual characteristics, such as body size, fitness, genome size, or dispersal abilities.
All decisions inside all of the submodels are stochastic (e.g., number of offspring) to maintain variability and relax assumptions.

#### Collectives.

#### Observation.

At the start and end of the simulation and at definable regular time intervals, the properties of all individuals
(including the properties of their locations) are recorded and written to files.


### 5. Initialisation

The initialisation step creates lineages with randomly chosen genetic and ecological trait values in each grid cell that is designated to receive an initial community.
This encompasses choosing the number of genes for a lineage, the number of linkage units and the within genome variance of trait values.
Trait with thus distributed trait values are the distributed randomly among the genes.
Population (number of individuals) size of a lineage is determined by the adult body size of individuals from a lineage.
At this point all individuals of a population are identical.
Values for ecological traits are then varied in each gene where a given trait is found, for all individuals of a lineage.
The variation is Normal distributed with the lineage trait value as mean and the product of sigma_l (phylogenetic constraint) and the lineage trait value as standard deviation.
This ensures initial genetic variation within a lineage population.
Thus created populations are added to a grid cell's community until the additional mass of another population would exceed the grid cell's carrying capacity.
Whether a grid cell receives an initial community depends on the map definition.
At the end of initialisation each of the thus populated grid cells holds one or more different populations, each from a separate lineage.

### 6. Input

At the start of a simulation user defined parameters are read, containing also a definition of the simulation arena (map definition).
This definition is provided in a separate plain text file.
Within the text file a line at the top containing a single number defines the number of timesteps the arena definition is valid for.
Every other non-empty line defines one grid cell with a unique identifier (a number), and the location of the grid cell as two coordinates.
Optionally, one can define the type of the grid cell (island or continent), [whether it is isolated, ]the temperature, ~~and the size~~.

Other optional parameters can be set in a separate configuration file and pertain to defining simulation scenarios:

| Name/Function | Default value | Description |
|---------------|---------------|-------------|
|"avgnoloci" | 1 | average number of loci/copies per gene |
|         "biggenelength" | 200,||
|         "burn-in" | 1000, | timesteps before invasion starts|
|         "cellsize" | 20e6, | maximum biomass per hectare in gramm (based on Clark et al. 2001)|
|         "config" | "simulation.conf", | configuration file name|
|         "debug" | false, | write out debug statements|
|         "dest" | string(Dates.today()), | output folder name|
|         "disturbance" | 0, | percentage of individuals killed per update per cell|
|         "fasta" | false, | record fasta data?|
|         "fertility" | exp(28.0), | global base reproduction rate 23.8 from Brown et al. 2004, alternatively 25.0, default 30.0|
|         "fixtol" | true,||
|         "global-species-pool" | 0, | size of the global species pool (invasion source)|
|         "growthrate" | exp(25.2), | global base growth/biomass production from Brown et al. 2004|
|         "indsize" | "seed", | initialize organisms as seed, adult or mixed|
|         "isolationweight" | 3, | additional distance to be crossed when dispersing from or to isolated patches|
|         "lineages" | false, | record lineage and diversity data?|
|         "linkage" | "random", | gene linkage type (random/full/none)|
|         "logging" | false, | write output to logfile|
|         "maps" | "", | comma-separated list of map files|
|         "maxdispmean" | 10, | maximum mean dispersal distance|
|         "maxrepsize" | 14, | maximal repsize in grams calculated as exp(maxrepsize) -> 1.2 t|
|         "maxseedsize" | 10, | maximal seedsize in grams calculated as exp(maxseedsize) -> 22 kg|
|         "maxtemp" | 313, | max optimum temp in K |
|         "minrepsize" | 3, | minimal repsize in grams calculated as exp(minrepsize) -> 20 g|
|         "minseedsize" | -2, | minimal seedsize in grams calculated as exp(minseedsize) -> 0.14 g|
|         "mintemp" | 283, | min optimum temp in K |
|         "mortality" | exp(22), | global base mortality from Brown et al. 2004 is 26.3, but competition and dispersal introduce add. mort.|
|         "mutate" | true, | mutations occur|
|         "mutationrate" | 3.6e10, | one mutation per generation/individual, corrected for metabolic function|
|         "nniches" | 2, | number of environmental niches (max. 3)|
|         "outfreq" | 100, | output frequency|
|         "phylconstr" | 0.1, | phylogenetic constraint during mutation and inter-loci variation. scales trait value as sd.|
|         "phylo" | false, | record phylogeny?|
|         "popsize" | "metabolic", | initialisation algorithm: metabolic/bodysize/minimal/single|
|         "precrange" | 10, | range from 0 for precipitation optimum|
|         "propagule-pressure" | 0, | number of non-native individuals introduced per invasion event|
|         "quiet" | false, | don't write output to screen|
|         "sdtemp" | 0.0, | SD of temperature change per time step|
|         "seed" | 0, | for the RNG, seed = 0 -> random seed|
|         "smallgenelength" | 20,|
|         "static" | true, | mainland sites don't undergo eco-evolutionary| processes|
|         "tolerance" | 0.8, | sequence similarity threshold for reproduction|
|         "traitnames" | ["compat", "dispmean", "dispshape", "precopt", "prectol", "repsize", "reptol", "seedsize", "tempopt", "temptol"], | minimal required traitnames|
|         "usebiggenes" | true||

If a parameter value is not specified by the user, the default value for that parameter set in the simulation code is assumed.


### 7. Submodels

#### Establishment.
Whenever an individual is new to a grid cell (by recent birth, dispersal event or environmental change [TODO]), their physical niche preferences
are compared with the actual niche properties, e.g. the temperature, T, of the present grid cell.
The individual fitness parameter, F, is set according to the deviation from the optimum value considering the niche breadth
as standard deviation of a gaussian curve.
```
F = a * exp(-(T - T_mean)^2 / (2 * sigma_T^2)) )
with a = 1 / (sigma_T * sqrt(2 * pi))
```

#### Competition.
Individuals are sorted according to fitness (low to high).
If the sum of the community's bodymass exceed the available space, individuals will be removed from the local community starting with the least fit individual with high probability
(following a geometric distribution).
Once total bodymass is below carrying capacity, the procedure terminates.

#### Growth.
Given an individual has undergone establishment,
an individual changes its size (M + delta_M)
following the metabolic theory and the global base growth rate, b_0:
```
    delta_M = b_0 * M^(3 / 4) * exp(-E_A / (k_B * T))
```
with E_A as activation energy and k_B the Boltzmann constant.
In case this change results in zero or negative body mass,
the individual is removed from the community.


#### Density independent mortality/Survival.
An individual is removed from the local community with a probability `p_mort` depending on its size `M` and a global base
mortality rate `b_mort`:
```
    p_mort = b_mort * M^(-1 / 4) * exp(-E_A / (k_B * T))
```
If the individual survives, its age increases by one.

#### Reproduction and mutation.
All individuals that have grown to or beyond their individual reproduction sizes may reproduce.
The number of offspring is randomly drawn following a Poisson distribution with mean N determined by the individual's
size `M` and a global base offspring number `N_0`:
```
    N =  N_0 * M^(-1 / 4) * exp(-E_A / (k_B * T))
```
Possible mates are selected within the same grid cell based on whether they belong to the same lineage, have reached maturity (which includes having established on the grid cell) and whether their compatibility sequences are sufficiently similar.
If a suitable partner is found, sets of haploid chromosomes from the diploid sets of both individuals are drawn randomly,
comprising the genome for the offspring.
[In case of a unsuccessful mate search, it is possible to enable self-compatibility, which involves recombination as well.]
At this point, mutations in the offspring's basecode may happen with a given probability.
In the case of mutation all traits associated with the respective gene will randomly change value (normally distributed,
with the standard deviation the product of sigma_l (phylogenetic constraint) and the original.
```
    newvalue = trait.value + rand(Normal(0, trait.value/phylconstr))
```
The new individuals' trait values are then calculated as the means of all alleles and
the individuals added to the community, with their size set to the initial bodymass (seed size).

#### Dispersal.
After reproduction and mutation, each offspring individual may disperse.
For each of these, a new location (i.e. x and y coordinates) is drawn randomly following a logistic distribution with mean and shape parameters (which controls long-distance-dispersal)
taken from the individual's traits. 
If a suitable grid cell is found at the drawn coordinates, the dispersing individual will be placed there and removed from the original community.
The removal happens even when there is no destination grid cell to be found.
Special attention is paid when the destination grid cell is of island type, while the origin is on the mainland and the simulation runs in static mode.
In this case the dispersing individual is copied to the new destination instead of moved.

#### Habitat change
If enabled, both environmental habitat parameters - temperature and precipitation -
change values throughout the simulation arena.
The amount and direction of change is the same for all grid cells across the landscape.
Changes to temperature and precipiation happen independently from one another.
The change is randomly drawn from a Normal distribution with the current value as the mean
and a user defined standard deviation.


### Output/Calculation
The main simulation data output is stored in two separate formats.
The first is a table containg data characterising the individuals.
Each line represents on individual.
The columns describe an individual's current state.
This is characterised by location, environmental conditions, ecological traits
and summary of the genetic architecture.
Additionally or alternatively to the individual level data, the data can be summarized at the population level
(i.e. all individuals of a common lineage within the same grid cell).
The second format is a fasta file containing the entire genome of all individuals.
Association of sequences to individuals, linkage units, genes and coded traits is defined in the
fasta headers.
Output is stored at the beginning and end of a simulation and at user-definable intervals.
The output considers the state of all non-seed individuals at those times.
