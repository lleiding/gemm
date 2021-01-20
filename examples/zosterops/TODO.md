# ZOSTEROPS TODO

## General

- [ ] read species & habitat descriptions of *Zosterops* & Taita Hills

- [ ] read up on current literature about evolutionary rescue, introgressive hybridization

- [X] prepare F2 presentation on "Adapting a plant community model to bird populations"

- [X] prepare Ecomods presentation on the design of GeMM

## GeMM adaptation

- [X] introduce `modes` to differentiate between experiment types

- [X] convert `precipitation` niche to `AGC` (above-ground carbon) (isn't that simply a relabelling?)

- [X] setup config file to enable testing

- [X] find a suitable temperature value

- [X] change species initialisation process to read predefined species

- [X] initialise two species
		- highland species *Z. silvanus* (Taita White-eye, a.k.a. subspecies of *Z. poliogaster*)
		- lowland species *Z. jubaensis* (former subspecies of *Z. abyssinicus*)

- [X] introduce sex (birds are not hermaphrodites)

- [X] write an alternate findmate() function
		-> breeding pairs stay faithful for life
		-> occupy & stay in a territory

- [X] make sure growth is capped at `repsize`

- [X] move species properties to settings

- [X] fix seed initialisation bug

- [X] add a StringGene class and `compressgenes` setting

- [X] make sure `map_creator.R` works with Petri's data

- [ ] performance issue (https://docs.julialang.org/en/v1/manual/performance-tips/):
  - [X] profile RAM/CPU trade-off with `compressgenes`
  - [X] circumvent `push!()` calls, preallocate memory instead - doesn't work?
  - [ ] avoid intermediate allocations when compressing gene sequences
  - [ ] avoid string interpolation in output
  - [ ] avoid `deepcopy()`?
  - [X] write a coordinate lookup function to speed up `zdisperse!()`

- [X] turn off mutation except at initialisation

- [X] make `settings` global

- [ ] rewrite dispersal function as discussed on 19/01/21

- [ ] implement lineage tagging for chromosomes

- [ ] adjust body sizes

- [X] juveniles are "born with" adult size

- [ ] hybridisation: only if no conspecific partner available

- [ ] hybridisation: assign offspring to lineage based on phenotypically similarity

- [ ] `createoffspring()` or `zreproduce()` have to deal with lineage labels

- [X] allow for ecological/neutral speciation in `ziscompatible()`

- [ ] rethink species trait definition in config

- [ ] test dispersal and reproduction functions

- [ ] vary hybridization affinity for core experiment

- [ ] prepare map series for SLOSS experiment (constant habitat size, random configuration)
		- single very large habitat (VL)
		- some large habitats (SL)
		- several small habitats (SS)

- [ ] measure & record:
		- degree of heterozygosity
		- genetic diversity
		- population sizes

- [ ] test parameterisation
		- very large *Z. jubaensis* populations?

- [ ] update data output

- [ ] update documentation

## Notes

- how should burn-in work?

- how do we design the SLOSS map series?
  - idea: use USGS forest cover map as "recovery scenario" (all exotic forest transformed to montane)?
  - or just shift species' AGC opt/tol to increase/decrease habitat suitability

- measuring heterozygosity:
  - number of hybrids in the population?
  - tag each chromosome with its original lineage and keep track of each?

- restrict mutations of max & min sizes?
  -> either in `mutate!` or in `checkviability!`

- what determines brood density in lowland species? forest cover?
