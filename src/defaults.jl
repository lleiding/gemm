# Return the default settings. All parameters must be registered here.
function defaultSettings()
    Dict(# general software settings
         "seed" => 0, # for the RNG, seed = 0 -> random seed
         "maps" => "", # comma-separated list of map files
         "config" => "simulation.conf", # configuration file name
         "dest" => string(Dates.today()), # output folder name
         "outfreq" => 100, # output frequency
         "fasta" => false, # record fasta data?
         "phylo" => false, # record phylogeny?
         "lineages" => false, # record lineage and diversity data?
         "logging" => false, # write output to logfile
         "quiet" => false, # don't write output to screen
         "debug" => false, # write out debug statements
         # main model settings
         "linkage" => "random", # gene linkage type (random/full/none)
         "nniches" => 2, # number of environmental niches (max. 3)
         "tolerance" => 0.8, # sequence similarity threshold for reproduction
         "fixtol" => true,
         "static" => true, # mainland sites don't undergo eco-evolutionary processes
         "mutate" => true, # mutations occur
         "initadults" => false, # initialize organisms as adults
         "initpopsize" => "metabolic", # initialisation algorithm: metabolic/bodysize/minimal/single
         "cellsize" => 20, # maximum biomass per hectare in tonnes (based on Clark et al. 2001) indicate gramm value as Float
         "biggenelength" => 200,
         "smallgenelength" => 20,
         # invasion specific settings
         "propagule-pressure" => 0, # number of non-native individuals introduced per invasion event
         "global-species-pool" => 0, # size of the global species pool (invasion source)
         "disturbance" => 0, # percentage of individuals killed per update per cell
         "burn-in" => 1000, # timesteps before invasion starts
         "minseedsize" => -2, # minimal seedsize in grams calculated as exp(minseedsize) -> 0.14 g
         "maxseedsize" => 10, # maximal seedsize in grams calculated as exp(maxseedsize) -> 22 kg
         "minrepsize" => 3, # minimal repsize in grams calculated as exp(minrepsize) -> 20 g
         "maxrepsize" => 14, # maximal repsize in grams calculated as exp(maxrepsize) -> 1.2 t
         # species trait names
         "traitnames" => ["compat",
                          "dispmean",
                          "dispshape",
                          "mutprob",
                          "precopt",
                          "prectol",
                          "repsize",
                          "reptol",
                          "seedsize",
                          "tempopt",
                          "temptol"] # minimal required traitnames
         ) 
end

