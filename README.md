The  model  description follows  the  ODD  (Overview,  Design  concepts, Details)  protocol  (Grimm  et  al.,  2006; 2010)

# 1. Purpose

This model is designed to simulate several patterns across genetic, individual, population and (meta-)community levels.


# 2. Entities, state variables and scales

Individuals (plants or animals) are the basic entity in the model.
Each individual carries a diploid set of chromosomes, which in turn are comprised of genes.
Some of the genes code for one ore more traits (pleiotropy), while a trait can be dependent on more than one gene (polygene).
The realized trait value is the mean of the trait alleles (quantitative trait loci).
Traits thus controlled encompass
the initial body mass (size) of offspring, M_s,
the body mass determining onset of maturity and thus reproductive capability, M_r,
mean dispersal distance, D_mean,
the shape of the dispersal kernel, controlling long-distance-dispersal, D_s,
the threshold of sequence similarity between mates determining compatibility, delta_I,
and values representing the optimum and the tolerance (standard deviation) of a physical niche parameter, such as temeperature and precipitation (T_mean and sigma_T or P_mean and sigma_P, resp.).
Alternatively to be controlled by mutable genes, these traits can also be set to fixed values.
Additionally, individuals carry attributes which describe their bodymass, M, their age, A and their adaptation to the physical environment (fitness), F.
Furthermore, every indivdual carries a Boolean marker used e.g.\ to store wether a given individual has newly arrived to a grid cell.

The base rates for processes governed by the metabolic theory of ecology (Brown et al. 2004) - growth, reproduction, mortality -
are global constants. Mutation rate is also a global constant.

Every individual is placed inside an arena of grid cells, each of which has a unique location (coordinates)
and is characterized by physical properties such as temperature, precipiation and size (carrying capacity).
Over the course of the simulation these properties (location or physical parameters) might change, reflecting
geomorphological dynamics.
All individuals within one grid cell constitute a community.
The characteristics of the grid cells combined with the state of inhabiting individuals constitute the state variables of the model.
Additional patterns or summary statistics may be calculated based on these individual information.

Processes and updates are repeated every timestep, while each timestep represent approximately one year.

# 3. Process overview and scheduling

In each discrete timestep each individual in each grid cell will (in no particular order unless otherwise stated) undergoes
the following processes:
- (1) establishment,
- (2) competition (individuals are sorted according to their fitness),
- (3) density independent mortality,
- [disturbance]
- (4) growth,
- (5) competition (individuals are sorted according to their fitness),
- (6) competition (individuals are sorted according to their fitness),
- (7) reproduction
- (8) mutation of offspring,
- (9) filtering of unviable individuals,
- [invasion]
- (10) seed dispersal.

Updates to individuals and thus the local communities happen instantaneously after a specific process has been executed
(asynchronous updating).

# 4. Design concepts

Basic principles.
-----------------
Metabolic theory of ecology (@submodel level).
Adaptive and non-adaptive radiation/evolution. Submodel level, but geomorphological change @system level.
Sexual reproduction.
Niche theory. Both at system and submodel level. Each individual carries unique ecological and/or functional traits,
as well as preferences for their physical environments.
Resource/energy limitation (carrying capacity). System level property, but invoked at submodel level.

Emergence.
----------
Species/Populations/ecotypes. Only constrained via genetic properties.
Community trait composition. Interplay of physical properties (environment, geographical properties) and within community
(competition strength via reproduction, growth, etc.).
Species numbers, endemics, speciation rate.

Adaption. (+Objectives.)
------------------------
See entities. Traits follow evolution: a trait changes its value randomly within a given phylogenetic constraint.
The success of the change (fitness) emerges as the result of adaption to the physical environment and the reproductive success of an individual over its competitors.

Learning.
---------

Prediction.
-----------

Sensing.
--------
Individuals sense directly the properties of their physical environment (e.g., temperature).

Interaction.
------------
Individuals directly interact when sexually reproducing. However, they are not affected themself by this interaction.
Instead, the interaction aims solely at determining the genotype of their offspring.
Additionally, competition for resource/energy/space between individuals represents indirect interaction.

Stochasticity.
--------------
Most of the submodels are carried out by all elegible individuals.
Some submodels (Survival, Competition, Mutation and Dispersal), however, happen with particular probabilities.
In these cases, execution of submodels is decided at random, taking into account individual characteristics, such as body size, fitness, genome size, or dispersal abilities.
All decisions inside all of the submodels are stochastic (e.g., number of offspring) to maintain variability and relax assumptions.

Collectives.
------------

Observation.
------------
At the start and end of the simulation and at definable regular time intervals, the properties of all individuals
(including the properties of their locations) are recorded and written to files.


# 5. Initialisation

The initialisation step creates lineages with randomly chosen genetic and ecological traits.
Population size of a lineage is determined by the adult body size of individuals from a lineage.
Values for ecological traits are then varied in each gene where a given trait is found, for all individuals of a lineage.
The variation is Normal distributed with the lineage trait value as mean and the product of sigma_l (phylogenetic constraint) and the lineage trait value as standard deviation.
This ensures initial genetic variation within a lineage population.
Thus created populations are added to a grid cell's community until the additional mass of another population would exceed the grid cell's carrying capacity.
Whether a grid cell receives an initial community depends on the map definition.
At the end of initialisation each of the thus populated grid cells holds one or more different populations, each from a separate lineage.

# 6. Input

At the start of a simulation user defined parameters are read, containing also a definition of the simulation arena (map definition).
This definition is provided in a separate plain text file.
Within the text file a line at the top containing a single number defines the number of timesteps the arena definition is valid for.
Every other non-empty line defines one grid cell with a unique identifier (a number), and the location of the grid cell as two coordinates.
Optionally, one can define the type of the grid cell (island or continent), whether it is isolated, the temperature, ~~and the size~~.

Other optional parameters can be set in a separate configuration file and pertain to defining simulation scenarios:
| Name/Function | Default value | Description |
+---------------+---------------+-------------|
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
         the inital random seed,
a comma-separated list of the names of arena (maps) definition files, the degree of genetic linkage (none, random or full),
and the compatibility tolerance of seed size values for reproduction (high, evolving or low).
If a parameter value is not specified by the user, the default value for that parameter set in the simulation code is assumed.


# 7. Submodels

## Establishment.
When an individual is new to a grid cell (by recent birth or dispersal event), their physical niche preferences
are compared with the actual niche properties of the present grid cell.
The individual fitness parameter is set according to the deviation from the optimum value considering the niche breadth
as standard deviation of a gaussian curve.
```
(a = 1 / (tolerance * sqrt(2 * pi)))
fitness = min(1, a * exp(-(value - optimum)^2 / (2 * tolerance^2)) )
```

## Competition.
Individuals are sorted according to fitness (low to high).
Starting with the least fit individuals an individual will be removed from the local community with high probability
(following a geometric distribution), if the sum of the community's bodymasses exceed the available space.
Once total bodymass is below carrying capacity, the procedure is finished.

## Growth.
Given an individual has undergone establishment (new-marker set to "false"), an individual changes its size (mass + delta_mass) following the metabolic theory and the global base growth rate:
```
    delta_mass = growthrate * mass^(3 / 4) * exp(-act / (boltz * temp))
```
In case this change results in negative body mass or the individual's initial body mass was less or equal zero,
the individual is removed from the community.


## Density independent mortality/Survival.
An individual is removed from the local community with a probability `p_mort` depending on its size `mass` and a global base
mortality rate `mort`:
```
    p_mort = mort * mass^(-1 / 4) * exp(-act / (boltz * temp))
```
If it survives instead, its age is increased by one.


## Reproduction and mutation.
All individuals that have grown to or beyond their individual reproduction sizes may reproduce.
The number of offspring is randomly drawn, following a Poisson distribution with a mean `n_offs` determined by the individual's
size `mass` and a global base offspring number `meanoffs`:
```
    n_offs =  meanoffs * currentmass^(-1 / 4) * exp(-act / (boltz * temp))
```
Following the individuals reproductive radius possible partners in the vicinity (grid cells whose distances fall within the
radius) are selected based on whether they belong to the same lineage, have reached maturity (which includes having established on the grid cell) and whether their seed size trait (`seedsize_mate`) falls within the mating individual's tolerance interval (`reptol`), which is determined by the following logical examination
```
(seedsize_mate >= reptol * seedsize) && (reptol * seedsize_mate <= seedsize)
```
If a suitable partner is found, sets of haploid chromosomes from the diploid sets of both individuals are drawn randomly,
comprising the genome for the offspring.
In case of a unsuccessful mate search, it is possible to enable self-compatibility, which involves recombination as well.
At this point every position in the offspring's basecode may mutate with a given probability.
In the case of mutation all traits associated with the respective gene will randomly change value (normally distributed,
with the standard deviation the quotient of the original value over a scaling constant - the phylogenetic constraint `phylconstr`).
```
    newvalue = trait.value + rand(Normal(0, trait.value/phylconstr))
```
The new individuals' trait values are then calculated as the means of all alleles and
the individuals added to the community, marked as new and with their size set to the initial bodymass (seed size).

## Dispersal.
After reproduction and mutation, each offspring individual may disperse.
For each of these, a distance is drawn randomly following a logistic distribution with mean and shape parameters (which controls long-distance-dispersal)
taken from the individual's traits.
Subsequently, a grid cell is chosen randomly that fits the drawn distance.
If such a grid cell is found, the dispersing individual will be placed there and removed from the original community.
The removal happens even when there is no destination grid cell to be found.
In case origin or destination grid cell are marked as isolated, the probability for successful dispersal needs to be
considered again for each isolated grid cell, whose barriers are to be crossed (origin and destination).
The barrier strength represents an additional distance, controlled by a global constant.
Special attention is paid when the destination grid cell is of island type, while the origin is on the mainland.
In this case the dispersing individual is copied to the new destinaction instead of moved.
Additionally, each  of these colonizers is recorded and its properties together with the respective
time step stored for later analysis.

## Output/Calculation
The main simulation data output is stored in two separate formats.
The first is a table containg data characterising the individuals.
Each line represents on individual.
The columns describe an individual's current state.
This is characterised by location, environmental conditions, ecological traits
and summary of the genetic architecture [UPDATE].
The second format is a fasta file containing the entire genome of all individuals.
Association of sequences to individuals, linkage units, genes and coded traits is defined in the
fasta headers [REVIEW].
Output is stored at the beginning and end of a simulation and at user-definable intervals.
The output considers the state of all individuals at those times,
but can be set to ignore individuals of a certain developmental stage or age [TODO].
