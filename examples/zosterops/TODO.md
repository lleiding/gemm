# ZOSTEROPS TODO

## General

- [ ] read species & habitat descriptions of *Zosterops* & Taita Hills

- [ ] read up on current literature about evolutionary rescue, introgressive hybridization

- [X] prepare F2 presentation on "Adapting a plant community model to bird populations"

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

- [ ] make sure `map_creator.R` works with Petri's data

- [ ] vary hybridization affinity for core experiment

- [ ] `createoffspring()` or `zreproduce()` have to deal with lineage labels

- [ ] prepare map series for SLOSS experiment (constant habitat size, random configuration)
		- single very large habitat (VL)
		- some large habitats (SL)
		- several small habitats (SS)

- [ ] measure & record:
		- degree of heterozygosity
		- genetic diversity
		- population sizes

- [ ] test parameterisation

- [ ] update data output

- [ ] update documentation

## Notes

- how should burn-in work?

- how do we design the SLOSS map series?

- restrict mutations of max & min sizes?
  -> either in `mutate!` or in `checkviability!`

- what determines brood density in lowland species? forest cover?
