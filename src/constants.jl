# Constants and settings for GeMM

### below variables should become settings elements at some point
const boltz = 1.38064852e-23 # J/K = m2⋅kg/(s2⋅K)
const act = 1e-19 # activation energy /J, ca. 0.63eV - Brown et al. 2004
const growthrate = exp(25.2) # global base growth/biomass production from Brown et al. 2004
const mortality = exp(22) # global base mortality from Brown et al. 2004 is 26.3, but competition and dispersal introduce add. mort.
const fertility = exp(30.0) # global base reproduction rate 23.8 from Brown et al. 2004, alternatively 25.0
const phylconstr = 10 #parse(ARGS[2])
# const meangenes = 20 # mean number of genes per individual
const mutationrate = 1e-3 * 0.3e11 # 1 base in 1000, correction factor for metabolic function
const isolationweight = 3 # additional distance to be crossed when dispersing from or to isolated patches
const maxdispmean = 10 # maximum mean dispersal distance
const genelength = 20 # sequence length of genes

# Return the default settings. All parameters must be registered here.
function defaultSettings()
    Dict(# general software settings
         "seed" => 0, # for the RNG, seed = 0 -> random seed
         "maps" => "", # comma-separated list of map files
         "config" => "simulation.conf", # configuration file name
         "dest" => string(Dates.today()), # output folder name
         "outfreq" => 1000, # output frequency
         "fasta" => true, # record fasta data?
         "lineages" => false, # record lineage and diversity data?
         "logging" => true, # write output to logfile
         "quiet" => false, # don't write output to screen
         "debug" => false, # write out debug statements
         # main model settings
         "linkage" => "random", # gene linkage type
         "nniches" => 2, # number of environmental niches (max. 3)
         "tolerance" => "low", # sequence similarity threshold for reproduction
         "static" => true, # mainland sites don't undergo eco-evolutionary processes
         "mutate" => true, # mutations occur
         "initadults" => false, # initialize organisms as adults
         "initpopsize" => "metabolic", # initialisation algorithm: metabolic/bodysize/minimal
         "cellsize" => 20, # maximum biomass per hectare in tonnes (based on Clark et al. 2001)
         # invasion specific settings
         "propagule-pressure" => 0, # number of non-native individuals introduced per invasion event
         "global-species-pool" => 0, # size of the global species pool (invasion source)
         "disturbance" => 0, # percentage of individuals killed per update per cell
         # species trait names
         "traitnames" => ["compat",
                          "dispmean",
                          "dispshape",
                          "mutprob",
                          "precopt",
                          "prectol",
                          "repradius",
                          "repsize",
                          "reptol",
                          "seedsize",
                          "tempopt",
                          "temptol"] # minimal required traitnames
         ) 
end

const settings = getsettings()

