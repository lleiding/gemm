# ZOSTEROPS TODO

## General

- [ ] read species & habitat descriptions of *Zosterops* & Taita Hills

- [ ] read up on current literature about evolutionary rescue, introgressive hybridization

- [ ] write F2 course report (about what exactly?)

## GeMM adaptation

- [ ] add `elevation`, `forest`, and `habitat` fields to `Patch`

- [ ] convert `precipitation` niche to `elevation` (and keep temperature constant)

- [ ] change species initialisation process to read predefined species

- [ ] initialise two species
		- highland species *Z. silvanus* (Taita White-eye, a.k.a. subspecies of *Z. poliogaster*)
		- lowland species *Z. jubaensis* (former subspecies of *Z. abyssinicus*)

- [ ] vary hybridization affinity for core experiment

- [ ] prepare map series for SLOSS experiment (constant habitat size, random configuration)
		- single very large habitat (VL)
		- some large habitats (SL)
		- several small habitats (SS)

- [ ] measure & record:
		- degree of heterozygosity
		- genetic diversity
		- population sizes

- [ ] small test runs to test parameterisation

## Notes

- do we need the `forest` variable?

- how should burn-in work?

- how do we design the SLOSS map series?

- restrict mutations of max & min sizes?
  -> either in `mutate!` or in `checkviability!`

- can we assume that habitat = "highland species habitat" and
  no habitat = "lowland species habitat"? Or is that only dependant
  on elevation? (But if so, how do we know where exactly the lowland
  species occurs?)

- what determines brood density in lowland species? forest cover?

- Where are the intermediate zones, where lowland and highland species
  can interact?
  -> mate search should not be restricted to the individual's cell,
  but a radius of 3 cells (allows interaction between different species)

- What about sexuality? (Birds are either male or female...)
