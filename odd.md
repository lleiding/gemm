ODD, after Grimm et al. 2010

# 1. Purpose

This model is designed to simulate several patterns across genetic, individual, population and (meta-)community levels.


# 2. Entities, state variable and scales

Individuals (plants or animals) are the basic entity in the model.
Each individual carries a diploid set of chromosomes, which in turn are comprised of genes.
Some of the genes code for on ore more traits (pleiotropy), while a trait can be dependent on more than one gene (polygene).
The realized trait value is the mean of the traits coded by maternal and paternal alleles.
Traits thus controlled encompass
the initial body mass (size) of offspring ("seedsize"),
the body mass determining onset of maturity and thus reproductive capability ("repsize"),
the groth rate ("growthrate"),
the probability for density-independent mortality ("ageprob"),
mean dispersal distance ("dispmean"),
the probability of dispersal ("dispprob"),
the shape of the dispersal kernel, controlling long-distance-dispersal ("dispshape"),
the probability and thus frequence of reproduction ("repprob"),
the radius in which an individual searches for a mate ("repradius"),
the mean number of offspring in a reproduction event (litter size,n "reprate"),
the threshold of sequence identity between mates determining genetic compatibility ("reptol"),
values representing the optimum and the tolerance (standard deviation) of a physical niche parameter ("tempopt", "temptol",
resp.).
Alternatively to be controlled by mutable genes, these traits can also be set to fixed values.
Additionally, individuals carry attributes which describe their bodymass, their age and the combined effects of adaptation
and life history (fitness).
Furthermore, every indivdual carries a marker which denotes whether a given individual has newly arrived to a patch/grid cell.

Every individual is placed inside an arena of grid cells or patches, each of which has a unique location (coordinates)
and is characterized by physical properties such as temperature and size (carrying capacity), but also isolation
(e.g., by barriers).
Over the course of the simulation these properties (location or physical parameters) might change, reflecting
geomorphological dynamics.
All individuals within one patch constitute a community.

Processes and updates are repeated every timestep, while each timestep represent approximately one generation.

# 3. Process overview and scheduling

In each discrete timestep each individual in each patch will (in no particular order unless otherwise stated) undergo
the following processes:
- (a) establishment,
- (b) competition (individuals are sorted according to their body sizes (from small to large)),
- (c) growth,
- (d) Density independent mortality,
- (e) reproduction,
- (f) dispersal.

Updates to individuals and thus the local communities happen instantaneously after a specific process has been employed
(asynchronous updating).

# 4. Design concepts

Basic principles.
-----------------
Metabolic theory of ecology (@submodel level).
Adaptive and non-adaptive radiation/evolution. Submodel level, but geomorpholical change @system level.
Sexual reproduction.
Niche theory. Both at system and submodel level. Each individual carries unique ecological and/or functional traits,
as well as preferences for their physical environments.
Resource/energy limitation (carrying capacity). System level property, but invoked ad submodel level.

Emergence.
----------
Species/Populations/ecotypes. Only constrained via genetic properties.
Community trait composition. Interplay of physical properties (environment, geographical properties) and within community
(competition strength via reproduction, growth, etc.).
Species numbers.

Adaption. (+Objectives.)
------------------------
See entities. Traits follow evolution: a trait changes its value randomly within a given constraint.
The success of the change (fitness) emerges as the reproductive success of an individual over its competitors.

Learning.
---------

Prediction.
-----------

Sensing.
--------
Individuals sense directly the properties of their physical environment (e.g., temperature).

Interaction.
------------
Individuals directly interact when sexually reproducing. However, they themselves are not affected by this interaction but the
interaction aims solely at determining the genotype of their offspring.
Additionally, competition for resource/energy/space between individuals represents an indirect interaction.

Stochasticity.
--------------
Most of the submodels are carried out by all elegible individuals.
Some submodels, however, happen with specific frequencies.
These latter submodels are executed randomly in combination with individual probabilities.
All decisions inside all of the submodels are stochastic (e.g., number of offspring) to maintain variability.

Collectives.
------------

Observation.
------------
At the end of the simulation and whenever the arena definition changes the properties of all individuals
(including the properties of their locations) are recorded and written to files.
Furthermore, all island colonizers are recorded and output in a similar way.


# 5. Initialisation

The initialisation step creates individuals with randomly chosen parameters and traits and deposits clones of one kind
of indivudual ("species") in one patch.
At the end of initialisation each patch holds several populations of clones.

# 6. Input

At the start of a simulation user defined parameters are read, containing also a definition of the simulation arena.
This definition is provided in a separate plain text file.
Within the text file a line at the top containing a single number defines the number of timesteps the arena definition is valid for.
Every other non-empty line defines one patch with a unique identifier (a number),and the location of the patch as two coordinates.
Optionally, one can define the type of the patch (island or continent), whether it is isolated, the temperature, and the size.

Other parameters specified at run time pertain to defining simulation scenarios.
```julia
        "--seed", "-s"
            help = "inital random seed"
        "--maps", "-m"
            help = "list of map files, comma separated"
        "--linkage", "-l"
            help = "gene linkage (\"none\", \"random\" or \"full\")"
            default = "random"
        "--tolerance", "-t"
            help = "tolerance of sequence identity when reproducing (\"high\", \"evo\" or \"low\")"
            default = "evo"
```
[...]

# 7. Submodels

## Establishment.
When an individual is new to a patch (by recent birth or dispersal event), their physical niche preferences
are compared with the actual niche properties of the present patch.
The individual fitness parameter is set according to the deviation from the optimum value considering the niche breadth
as standard deviation of a gaussian curve.
```
    fitness = gausscurve(tempopt, temptol, temp)
```

## Competition.
Individuals are sorted according to their body sizes (from small to large).
Starting with small individuals an individual will be removed from the local community with high probability
(following a geometric distribution), if the sum of the community's bodymasses exceed the available space.
Once total bodymass is below carrying capacity the procedure is finished.

## Growth.
An individual changes its size following the metabolic theory of ecology and an individual growth rate parameter
modulated by the fitness parameter.
```
    newmass = growthrate * patch.community[idx].fitness * mass^(3/4) * exp(-act/(boltz*temp)) * normconst
```

## Density independent mortality.
An individual is removed from the local community with a certain probability that is specified
within that individual's traits.
A lower fitness parameter increases mortality.
```
    dieprob = (ageprob * mass^(-1/4) * exp(-act/(boltz*temp)) * normconst) * (1-fitness)
```

## Reproduction.
Given a individual probability ("repprob", modified by the fitness parameter) an individual decides on whether
to reproduce.
In the case of reproduction and if the individual is larger or equal the individual's reproductive size, first the number
of offspring is randomly drawn, following the individual's trait value and the metabolic theory.
```
    metaboffs =  meanoffs * currentmass^(-1/4) * exp(-act/(boltz*temp)) * normconst
```
Following the individuals reproductive radius possible partners in the vicinity (patches whose distances fall within the
radius) are selected based on whether they share the same chromosome number with the reproducing individual and
whether the sequence identity between both individuals is equal or higher the reproducing individual's tolerance.
If a suitable partner is found sets of haploid chromosomes from the diploid sets of both individuals are drawn randomly,
comprising the genome for the offspring.
At this point every position in the offspring's basecode may mutate with a given probability.
In the case of mutation all traits associated with the respective gene will randomly change value (normally distributed,
with the standard deviation the quotient of the original value over a scaling constant).
```
    newvalue = trait.value + rand(Normal(0, trait.value/mutscaling))
```
The new individuals' trait values are then calculated as the means between maternal and paternal alleles and
the individuals added to the community, marked as new and with their size set to the initial bodymass.

## Dispersal.
An individual disperses with an individual probability.
If dispersal is chosen a distance is drawn randomly following a logistic distribution and mean and shape parameters
taken from the individual's traits.
Subsequently, a patch is chosen randomly that fits the drawn distance.
If such a patch is found the dispersing individual will be placed there and removed from the original community.
The removal happens even when there is no destination patch to be found.
In case origin or destination patch are marked as isolated the probability for successful dispersal needs to be
considered again for each isolated patch, whose barriers are to be crossed (origin and destination).
Special attention is paid when the destination patch is of island type.
In this case the dispersing individual (colonizer) is recorded and its properties together with the respective
time step stored for later analysis.