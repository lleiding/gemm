*Zosterops* Team Meeting 19/01/2021
===================================

Participants: Juliano, Luc, Jan, Charlotte, Robin, Daniel

### Dispersal

-   in *Zosterops*, females are the dispersive sex: more likely to
    disperse and go further
-   although this sex bias might be weaker in *Z. silvanus* due to
    "island effect" (overall reduced dispersal) → this might be relevant
    in the habitat fragmentation study!
-   post-natal dispersal distance for *Zosterops* is largely unknown
    (max. observed 14.5 km?)
-   for model: don't differentiate between male and female dispersal
    (but perhaps check this assumption in exploratory simulations with
    species-specific dispersal sex bias)

### Mutations

-   turn off mutations for habitat fragmentation study to avoid noisy
    data, but do exploratory simulations to check if there's an effect
-   keep them on for the phylogeny study (*Zosterops* is a "great
    speciator"!)
-   do a single round of mutations when initialising populations from
    species archetype to introduce genetic variability

### Hybridisation

-   two congenerics may hybridise if they're in a common habitat patch
    (above-ground carbon preference overlap) and there is no conspecific
    partner available
-   in habitat fragmentation study: no sequence similarity testing for
    post-zygotic compatibility
-   in phylogenetic study: sequence similarity matching based on neutral
    or functional genes (depending on speciation scenario)
-   assign hybrid individuals to a lineage based on phenotypic
    similarity
-   alternatives: assign each hybrid individual to a random lineage, or
    to the paternal lineage

### Miscellaneous

-   tag chromosomes with archetypical lineage to allow quantification of
    heterozygosity
-   body sizes need to be adjusted downward and made species-specific
-   we don't need juvenile body sizes, growth is too rapid (several
    weeks to adult size)
-   individually configure frequency of each output type?

### Next steps

-   Daniel: adapt model as discussed above
-   Jan: find information on *Zosterops* dispersal distances and body
    sizes
-   Charlotte & Robin: implement configuration parameter to switch
    between ecological/neutral/no speciation
-   all (except Luc): next meetings Tuesday 11 AM

