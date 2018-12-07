# Return the default settings. All parameters must be registered here.
function defaultSettings()
    Dict("avgnoloci" => 1, # average number of loci/copies per gene
         "biggenelength" => 200,
         "burn-in" => 1000, # timesteps before invasion starts
         "cellsize" => 20e6, # maximum biomass per hectare in gramm (based on Clark et al. 2001)
         "config" => "simulation.conf", # configuration file name
         "debug" => false, # write out debug statements
         "dest" => string(Dates.today()), # output folder name
         "disturbance" => 0, # percentage of individuals killed per update per cell
         "fasta" => false, # record fasta data?
         "fertility" => exp(28.0), # global base reproduction rate 23.8 from Brown et al. 2004, alternatively 25.0, default 30.0
         "fixtol" => true,
         "global-species-pool" => 0, # size of the global species pool (invasion source)
         "growthrate" => exp(25.2), # global base growth/biomass production from Brown et al. 2004
         "indsize" => "seed", # initialize organisms as seed, adult or mixed
         "isolationweight" => 3, # additional distance to be crossed when dispersing from or to isolated patches
         "lineages" => false, # record lineage and diversity data?
         "linkage" => "random", # gene linkage type (random/full/none)
         "logging" => false, # write output to logfile
         "maps" => "", # comma-separated list of map files
         "dispmean" => 1, # mean dispersal distance in cells
         "dispshape" => 10, # shape parameter fo dispersal kernel - tail fatness
         "maxrepsize" => 14, # maximal repsize in grams calculated as exp(maxrepsize) -> 1.2 t
         "maxseedsize" => 10, # maximal seedsize in grams calculated as exp(maxseedsize) -> 22 kg
         "maxtemp" => 313, # max optimum temp in K 
         "minrepsize" => 3, # minimal repsize in grams calculated as exp(minrepsize) -> 20 g
         "minseedsize" => -2, # minimal seedsize in grams calculated as exp(minseedsize) -> 0.14 g
         "mintemp" => 273, # min optimum temp in K 
         "mortality" => exp(22), # global base mortality from Brown et al. 2004 is 26.3, but competition and dispersal introduce add. mort.
         "mutate" => true, # mutations occur
         "mutationrate" => 3.6e10, # one mutation per generation/individual, corrected for metabolic function
         "nniches" => 2, # number of environmental niches (max. 3)
         "outfreq" => 100, # output frequency
         "phylconstr" => 0.1, # phylogenetic constraint during mutation and inter-loci variation. scales trait value as sd.
         "phylo" => false, # record phylogeny?
         "popsize" => "metabolic", # initialisation algorithm: metabolic/bodysize/minimal/single
         "precrange" => 10, # range from 0 for precipitation optimum
         "propagule-pressure" => 0, # number of non-native individuals introduced per invasion event
         "quiet" => false, # don't write output to screen
         "raw" => true, # record raw data
         "sdtemp" => 0.0, # SD of temperature change per time step
         "sdprec" => 0.0, # SD of precipiation change per time step
         "seed" => 0, # for the RNG, seed = 0 -> random seed
         "smallgenelength" => 20,
         "static" => true, # mainland sites don't undergo eco-evolutionary processes
         "stats" => true, # record population statistics
         "tolerance" => 0.8, # sequence similarity threshold for reproduction
         "traitnames" => ["compat",
                          "dispmean",
                          "dispshape",
                          "precopt",
                          "prectol",
                          "repsize",
                          "reptol",
                          "seedsize",
                          "tempopt",
                          "temptol"], # minimal required traitnames
         "usebiggenes" => true
         ) 
end

