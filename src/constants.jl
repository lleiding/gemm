# Constants and settings for GeMM

### below variables should become settings elements at some point
const boltz = 1.38064852e-23 # J/K = m2⋅kg/(s2⋅K)
const act = 1e-19 # activation energy /J, ca. 0.63eV - Brown et al. 2004
const growthrate = exp(25.2) # global base growth/biomass production from Brown et al. 2004
const mortality = exp(22) # global base mortality from Brown et al. 2004 is 26.3, but competition and dispersal introduce add. mort.
const fertility = exp(25.0) # global base reproduction rate 23.8 from Brown et al. 2004, alternatively 25.0
const phylconstr = 10 #parse(ARGS[2])
# const meangenes = 20 # mean number of genes per individual
const mutationrate = 1e7 # one mutation per generation with 200bp genes, corrected for metabolic function
const isolationweight = 3 # additional distance to be crossed when dispersing from or to isolated patches
const maxdispmean = 10 # maximum mean dispersal distance
const genelength = 20 # sequence length of genes

# XXX Trait value calculation happens in subprocs.jl -> createtraits()
# XXX Default patch settings can be found in types.jl
