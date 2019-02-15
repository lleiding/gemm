var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "CurrentModule = GeMM"
},

{
    "location": "#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "This is the documentation for the Island Speciation Model, a.k.a. Genetically explicit Metacommunity Model (GeMM). (Yeah, we\'re still working on the name...)The aim of this model is to create a virtual island ecosystem that can be used to explore ecological and evolutionary hypotheses in silico. It is genetically and spatially explicit, with discrete space and time.This documentation is generated from the source code using Julia\'s inbuilt Documenter module. It is sorted first by functionality, second by source code file."
},

{
    "location": "#GeMM.rungemm",
    "page": "Introduction",
    "title": "GeMM.rungemm",
    "category": "function",
    "text": "rungemm(config, seed)\n\nWrapper for runsim() Runs a simulation using configuration file config, random seed seed and other settings provided via commandline, configuration file or the defaults. Performs a pre-compilation run first.\n\n\n\n\n\n"
},

{
    "location": "#GeMM.runsim",
    "page": "Introduction",
    "title": "GeMM.runsim",
    "category": "function",
    "text": "runsim(config, seed, prerun)\n\nPerforms a simulation run using configuration file config, random seed seed and other settings provided via commandline, configuration file or the defaults. prerun toggles a pre-compilation run\n\n\n\n\n\n"
},

{
    "location": "#Running-the-model-1",
    "page": "Introduction",
    "title": "Running the model",
    "category": "section",
    "text": "The functions in the run_simulation.jl file are used to start a simulation run:Modules = [GeMM]\nPages = [\"run_simulation.jl\"]rungemmparallel.jl bundles these for quick access. From the commandline, call:> julia -p <cores> rungemmparallel.jl -c <configs>where <cores> is the number of processors you want to make available to Julia (one processor per simulation max), and <configs> is a comma-separated list of configuration files that are to be processed.Last updated: 2019-02-15 (commit 27a93d0)  "
},

{
    "location": "model/index.html#",
    "page": "Main Model",
    "title": "Main Model",
    "category": "page",
    "text": ""
},

{
    "location": "model/index.html#Main-Model-1",
    "page": "Main Model",
    "title": "Main Model",
    "category": "section",
    "text": "These are the core functions of the model, implementing the high-level processes of the simulation and the object classes that are involved."
},

{
    "location": "model/index.html#GeMM.AbstractGene",
    "page": "Main Model",
    "title": "GeMM.AbstractGene",
    "category": "type",
    "text": "An abstraction layer to hide the implementation details of genes.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.BigGene",
    "page": "Main Model",
    "title": "GeMM.BigGene",
    "category": "type",
    "text": "A variation of the standard Gene struct that allows for longer base sequence lengths through its use of BigInt.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.Chromosome",
    "page": "Main Model",
    "title": "GeMM.Chromosome",
    "category": "type",
    "text": "A struct that saves a collection of genes and a tag to denote whether they are maternal or not.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.Gene",
    "page": "Main Model",
    "title": "GeMM.Gene",
    "category": "type",
    "text": "The standard gene struct; including the base sequence (represented in binary and stored and as an integer), and an array of traits the gene codes for.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.Individual",
    "page": "Main Model",
    "title": "GeMM.Individual",
    "category": "type",
    "text": "One of the core structs of the model, representing an individual organism.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.Patch",
    "page": "Main Model",
    "title": "GeMM.Patch",
    "category": "type",
    "text": "One of the core structs of the model, representing a one-hectare patch of ground.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.Trait",
    "page": "Main Model",
    "title": "GeMM.Trait",
    "category": "type",
    "text": "A simple struct to store organism trait names and values.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#entities.jl-1",
    "page": "Main Model",
    "title": "entities.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"entities.jl\"]"
},

{
    "location": "model/index.html#GeMM.createpop-Tuple{Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.createpop",
    "category": "method",
    "text": "createpop(settings)\n\nCreate a new, random individual and replicates it a certain number of times (depending on metabolic variables) to create a new population of organisms. Returns an array of individuals.\n\nFIXME Why does this not interface with auxfuncts/createind()?\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.createworld-Tuple{Array{Array{String,1},1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.createworld",
    "category": "method",
    "text": "createworld(maptable, settings)\n\nUse a parsed map file (as returned by readmapfile) to create the world.  Initialises each patch with its parameters and a new community, then returns an array of patches.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.genesis-Tuple{Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.genesis",
    "category": "method",
    "text": "genesis(settings)\n\nCreate a new community, composed of random new species populations, for a patch. Returns an array of individuals.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.updateworld!-Tuple{Array{GeMM.Patch,1},Array{Array{String,1},1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.updateworld!",
    "category": "method",
    "text": "updateworld!(world, maptable, cellsize)\n\nReinitialise the world from another parsed map file. Works analogously to  createworld. Intended for use in scenarios where the model world changes during a run (e.g. through global warming or island growth).\n\n\n\n\n\n"
},

{
    "location": "model/index.html#initialisation.jl-1",
    "page": "Main Model",
    "title": "initialisation.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"initialisation.jl\"]"
},

{
    "location": "model/index.html#GeMM.changehabitat!-Tuple{Array{GeMM.Patch,1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.changehabitat!",
    "category": "method",
    "text": "changehabitat!(world, settings)\n\nCarry out \'global change\' on all patches.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.changeprec!-Tuple{Array{GeMM.Patch,1},Float64}",
    "page": "Main Model",
    "title": "GeMM.changeprec!",
    "category": "method",
    "text": "changeprec!(world, sdprec)\n\nChange the precipitation of all patches according to a normal distribution.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.changetemp!-Tuple{Array{GeMM.Patch,1},Float64}",
    "page": "Main Model",
    "title": "GeMM.changetemp!",
    "category": "method",
    "text": "changetemp!(world, sdtemp)\n\nChange the temperature of all patches according to a normal distribution.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.checkviability!-Tuple{Array{GeMM.Individual,1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.checkviability!",
    "category": "method",
    "text": "checkviability!(community, settings)\n\nCheck whether all individuals in the passed community conform to a basic set of constraints (i.e. all traits are present and certain properties are >= 0). Individuals that fail the test are removed from the community.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.checkviability!-Tuple{Array{GeMM.Patch,1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.checkviability!",
    "category": "method",
    "text": "checkviability(world, settings)\n\nCheck the viability of all individuals.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.checkviability!-Tuple{GeMM.Patch,Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.checkviability!",
    "category": "method",
    "text": "checkviability!(patch, settings)\n\nCheck the viability of the individuals in this patch.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.compete!",
    "page": "Main Model",
    "title": "GeMM.compete!",
    "category": "function",
    "text": "compete!(world, static)\n\nCarry out competition on all patches.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.compete!-Tuple{GeMM.Patch}",
    "page": "Main Model",
    "title": "GeMM.compete!",
    "category": "method",
    "text": "compete!(patch)\n\nWhile there are too many organisms in a patch, pick two at random and kill the one that is less adapted to the local precipitation levels.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.disperse!",
    "page": "Main Model",
    "title": "GeMM.disperse!",
    "category": "function",
    "text": "disperse!(world, static)\n\nDispersal of individuals within the world.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.disturb!-Tuple{Array{GeMM.Patch,1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.disturb!",
    "category": "method",
    "text": "disturb!(world, settings)\n\nDisturb all patches in the world.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.disturb!-Tuple{GeMM.Patch,Int64}",
    "page": "Main Model",
    "title": "GeMM.disturb!",
    "category": "method",
    "text": "disturb!(patch, intensity)\n\nSpecies-independent mortality due to disturbance on the given patch. The intensity gives the mortality percentage.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.establish!",
    "page": "Main Model",
    "title": "GeMM.establish!",
    "category": "function",
    "text": "establish!(world, nniches, static)\n\nCarry out establishment for each patch in the world.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.establish!",
    "page": "Main Model",
    "title": "GeMM.establish!",
    "category": "function",
    "text": "establish!(patch, nniches)\n\nEstablishment of individuals in patch p: Sets the adaption parameters (~fitness) according to an individual\'s adaptation to the niches of the surrounding environment.\n\nA maximum of two niches is currently supported.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.grow!-Tuple{Array{GeMM.Patch,1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.grow!",
    "category": "method",
    "text": "grow!(world, settings)\n\nCarry out growth for all patches.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.grow!-Tuple{GeMM.Patch,Float64}",
    "page": "Main Model",
    "title": "GeMM.grow!",
    "category": "method",
    "text": "grow!(patch, growthrate)\n\nGrowth of individuals in the given patch. The actual growthrate is calculated with a metabolic formula, modified by the passed growthrate variable.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.mutate!",
    "page": "Main Model",
    "title": "GeMM.mutate!",
    "category": "function",
    "text": "mutate!(traits, settings, locivar)\n\nLoop over an array of traits, mutating each value in place along a normal distribution. locivar can be used to scale the variance of the normal distribution used to draw new trait values (together with settings[phylconstr]).\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.mutate!-Tuple{Array{GeMM.Patch,1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.mutate!",
    "category": "method",
    "text": "mutate!(world, settings)\n\nMutate the world. (That sounds scary!)\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.mutate!-Tuple{GeMM.Individual,Float64,Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.mutate!",
    "category": "method",
    "text": "mutate!(individual, temp, settings)\n\nMutate an individual\'s genome (sequence and traits) in place.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.mutate!-Tuple{GeMM.Patch,Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.mutate!",
    "category": "method",
    "text": "mutate!(patch, setting)\n\nMutate all seedling individual\'s in a patch.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.reproduce!-Tuple{Array{GeMM.Patch,1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.reproduce!",
    "category": "method",
    "text": "reproduce!(world, settings)\n\nCarry out reproduction on all patches.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.reproduce!-Tuple{GeMM.Patch,Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.reproduce!",
    "category": "method",
    "text": "reproduce!(patch, settings)\n\nReproduction of individuals in a patch.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.survive!-Tuple{Array{GeMM.Patch,1},Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.survive!",
    "category": "method",
    "text": "survive!(world, settings)\n\nWorld-wide mortality. Sounds apocalyptic, but is just a fact of life.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#GeMM.survive!-Tuple{GeMM.Patch,Float64}",
    "page": "Main Model",
    "title": "GeMM.survive!",
    "category": "method",
    "text": "survive!(patch, mortality)\n\nDensity independent survival of individuals in a patch. The actual mortality  probability is calculated with a metabolic formula, modified by the passed mortality variable and an individual\'s temperature adaptation.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#submodels.jl-1",
    "page": "Main Model",
    "title": "submodels.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"submodels.jl\"]"
},

{
    "location": "model/index.html#GeMM.simulate!",
    "page": "Main Model",
    "title": "GeMM.simulate!",
    "category": "function",
    "text": "simulate!(world, settings, timesteps)\n\nThis is the central function of the model with the main event loop. It defines the scheduling for all submodels and output functions.\n\n\n\n\n\n"
},

{
    "location": "model/index.html#scheduling.jl-1",
    "page": "Main Model",
    "title": "scheduling.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"scheduling.jl\"]"
},

{
    "location": "io/index.html#",
    "page": "Input, Output, and Settings",
    "title": "Input, Output, and Settings",
    "category": "page",
    "text": ""
},

{
    "location": "io/index.html#Input,-Output,-and-Settings-1",
    "page": "Input, Output, and Settings",
    "title": "Input, Output, and Settings",
    "category": "section",
    "text": "These functions are responsible for reading in all model configurations (passed by config file or commandline), administrating them during a run, and printing any output."
},

{
    "location": "io/index.html#constants.jl-1",
    "page": "Input, Output, and Settings",
    "title": "constants.jl",
    "category": "section",
    "text": "This file defines constants that are needed for calculations involving the metabolic theory of ecology.Currently includes boltz (the Boltzmann constant) and act (the activation energy)."
},

{
    "location": "io/index.html#GeMM.defaultSettings-Tuple{}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.defaultSettings",
    "category": "method",
    "text": "defaultSettings()\n\nDefines the list of configuration variables and returns their default values in a dict.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#defaults.jl-1",
    "page": "Input, Output, and Settings",
    "title": "defaults.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"defaults.jl\"]"
},

{
    "location": "io/index.html#GeMM.basicparser-Tuple{String}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.basicparser",
    "category": "method",
    "text": "basicparser(filename)\n\nDo elementary parsing on a config or map file.\n\nReads in the file, strips whole-line and inline comments and separates lines by whitespace. Returns a 2d array representing the tokens in each line.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.getsettings",
    "page": "Input, Output, and Settings",
    "title": "GeMM.getsettings",
    "category": "function",
    "text": "getsettings()\n\nCombines all configuration options to produce a single settings dict. Order of precedence: commandline parameters - config file - default values\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.parsecommandline-Tuple{}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.parsecommandline",
    "category": "method",
    "text": "parsecommandline()\n\nCertain parameters can be set via the commandline.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.parseconfig-Tuple{String}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.parseconfig",
    "category": "method",
    "text": "parseconfig(filename)\n\nParse a configuration file and return a settings dict.\n\nThe config syntax is very simple: each line consists of a parameter name and a value (unquoted), e.g. nniches 2. # is the comment character.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.readmapfile-Tuple{String,Dict{String,Any}}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.readmapfile",
    "category": "method",
    "text": "readmapfile(mapfilename, settings)\n\nParse a map file and return the number of timesteps this map is to be used for (first line of the file) and the patch definitions. The latter is used by createworld and updateworld!.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#input.jl-1",
    "page": "Input, Output, and Settings",
    "title": "input.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"input.jl\"]"
},

{
    "location": "io/index.html#GeMM.dumpinds",
    "page": "Input, Output, and Settings",
    "title": "GeMM.dumpinds",
    "category": "function",
    "text": "dumpinds(world, io, sep)\n\nOutput all data of individuals in world as table to io. Columns are separated by sep. WARNING: this produces very large files!\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.makefasta",
    "page": "Input, Output, and Settings",
    "title": "GeMM.makefasta",
    "category": "function",
    "text": "makefasta(world, settings, io, onlyisland, sep)\n\nRecord the genome of every individual currently alive to the given IO stream. WARNING: this produces very large files!\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.printheader",
    "page": "Input, Output, and Settings",
    "title": "GeMM.printheader",
    "category": "function",
    "text": "printheader(settings, io, sep)\n\nPrint a list of property names to the given IO stream. This is a helper function for dumpinds.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.printpopheader-Tuple{IO}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.printpopheader",
    "category": "method",
    "text": "printpopheader(io)\n\nPrint a list of property names to the given IO stream. This is a helper function for printpopstats.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.printpopstats-Tuple{IO,Array{GeMM.Patch,1},Dict{String,Any},Integer}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.printpopstats",
    "category": "method",
    "text": "printpopstats(io, world, settings, timestep)\n\nRecord statistical information (maximum, minimum, median, standard deviation) for a range of individual properties, as seen over the whole world population.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.recordlineages-Tuple{Array{GeMM.Patch,1},Dict{String,Any},Int64}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.recordlineages",
    "category": "method",
    "text": "recordlineages(w)\n\nSave the abundance of each lineage per patch.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.recordstatistics-Tuple{Array{GeMM.Patch,1},Dict{String,Any}}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.recordstatistics",
    "category": "method",
    "text": "recordstatistics(w)\n\nWrite out world properties to the log file for later analysis.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.setupdatadir-Tuple{Dict{String,Any}}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.setupdatadir",
    "category": "method",
    "text": "setupdatadir(dir)\n\nCreates the output directory and copies relevant files into it. If the output directory already includes files, create a new directory by appending a counter.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.simlog",
    "page": "Input, Output, and Settings",
    "title": "GeMM.simlog",
    "category": "function",
    "text": "simlog(msg, settings, category, logfile, onlylog)\n\nWrite a log message to STDOUT/STDERR and the specified logfile  (if logging is turned on in the settings).\n\nCategories: d (debug), i (information, default), w (warn), e (error)\n\nIf logfile is the empty string (default: \"simulation.log\"), the message will  only be printed to the screen. If onlylog is true (default: false), the message is not printed to screen but only to the log.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.writedata-Tuple{Array{GeMM.Patch,1},Dict{String,Any},Int64}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.writedata",
    "category": "method",
    "text": "writedata(world, settings, timestep)\n\nWrites simulation output from world to separate table and fasta files. (Which data is recorded depends on the settings.) timestep and setting information is used for file name creation.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#GeMM.writesettings-Tuple{Dict{String,Any}}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.writesettings",
    "category": "method",
    "text": "writesettings(settings)\n\nRecord the settings actually used for a simulation run (cf. getsettings). Creates a config file that can be used for future replicate runs. Also records a time stamp and the current git commit.\n\n\n\n\n\n"
},

{
    "location": "io/index.html#output.jl-1",
    "page": "Input, Output, and Settings",
    "title": "output.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"output.jl\"]"
},

{
    "location": "aux/index.html#",
    "page": "Internal functions",
    "title": "Internal functions",
    "category": "page",
    "text": ""
},

{
    "location": "aux/index.html#Internal-functions-1",
    "page": "Internal functions",
    "title": "Internal functions",
    "category": "section",
    "text": "These are functions that the model uses in the background. They include a range of things: sub-submodels such findposspartner, utility functions like  gausscurve, or performance hacks like seq2num."
},

{
    "location": "aux/index.html#GeMM.checkborderconditions-Tuple{Array{GeMM.Patch,1},Int64,Int64}",
    "page": "Internal functions",
    "title": "GeMM.checkborderconditions",
    "category": "method",
    "text": "checkborderconditions!(w, x, y)\n\ncheck if coordinates x and y lie within world w and correct if not, considering defined border conditions.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.createchrms-Tuple{Int64,Array{GeMM.AbstractGene,1}}",
    "page": "Internal functions",
    "title": "GeMM.createchrms",
    "category": "method",
    "text": "createchrms(nchrms, genes)\n\nRandomly distribute the passed genes into the given number of chromosomes. Returns an array of chromosome objects.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.creategenes-Tuple{Int64,Array{GeMM.Trait,1},Dict{String,Any}}",
    "page": "Internal functions",
    "title": "GeMM.creategenes",
    "category": "method",
    "text": "creategenes(ngenes, traits, settings)\n\nRandomly create a given number of gene objects, with their base sequence and associated traits. Returns the result as an array of AbstractGenes.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.createind",
    "page": "Internal functions",
    "title": "GeMM.createind",
    "category": "function",
    "text": "createind(settings, marked=false)\n\nCreate an individual organism with of a new species with a random genome.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.createoffspring-Tuple{Integer,GeMM.Individual,GeMM.Individual,Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.createoffspring",
    "category": "method",
    "text": "createoffspring(noffs, individual, partner, traitnames)\n\nThe main reproduction function. Take two organisms and create the given number of offspring individuals. Returns an array of individuals.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.createtraits-Tuple{Dict{String,Any}}",
    "page": "Internal functions",
    "title": "GeMM.createtraits",
    "category": "method",
    "text": "createtraits(settings)\n\nCreate an array of trait objects generated from the default trait values (with a random offset).\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.diversity-Tuple{Array{GeMM.Patch,1}}",
    "page": "Internal functions",
    "title": "GeMM.diversity",
    "category": "method",
    "text": "diversity(w)\n\nCalculate (average) alpha, beta and gamma diversity of the world. Returns a tuple with the three values (a,b,g). cf. Veech et al. 2002\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.findisland-Tuple{Array{GeMM.Patch,1}}",
    "page": "Internal functions",
    "title": "GeMM.findisland",
    "category": "method",
    "text": "findisland(w)\n\nwithin world w, find out in which direction from the continent the island(s) lie(s).\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.findposspartner-Tuple{GeMM.Patch,GeMM.Individual,Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.findposspartner",
    "category": "method",
    "text": "findposspartner(patch, individual, traitnames)\n\nFind a reproduction partner for the given individual in the given patch.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.freespace-Tuple{Array{GeMM.Patch,1}}",
    "page": "Internal functions",
    "title": "GeMM.freespace",
    "category": "method",
    "text": "freespace(world)\n\nCalculate the average amount of free space on each patch.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.gausscurve",
    "page": "Internal functions",
    "title": "GeMM.gausscurve",
    "category": "function",
    "text": "gausscurve(b, c, x, a=1.0)\n\nCalculate the value of the Gauss function (\"bell curve\") at point x; with a being the maximum height of the curve, b the position of the curve center and c the standard deviation (\"width\").\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.getmeantraitvalue-Tuple{Array{GeMM.Trait,1},Integer}",
    "page": "Internal functions",
    "title": "GeMM.getmeantraitvalue",
    "category": "method",
    "text": "getmeantraitvalue(traits, traitidx)\n\nTake an array of traits and return the mean value of the indexed trait.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.getstdtraitvalue-Tuple{Array{GeMM.Trait,1},Integer}",
    "page": "Internal functions",
    "title": "GeMM.getstdtraitvalue",
    "category": "method",
    "text": "getstdtraitvalue(traits, traitidx)\n\nTake an array of traits and return the standard deviation of the indexed trait.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.gettraitdict-Tuple{Array{GeMM.AbstractGene,1},Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.gettraitdict",
    "category": "method",
    "text": "gettraitdict(genes, traitnames)\n\nCalculate the trait dict for an array of genes.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.gettraitdict-Tuple{Array{GeMM.Chromosome,1},Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.gettraitdict",
    "category": "method",
    "text": "gettraitdict(chromosomes, traitnames)\n\nConvert a genome (an array of chromosomes) into a dict of traits and their values.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.gettraitdict-Tuple{Array{GeMM.Trait,1},Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.gettraitdict",
    "category": "method",
    "text": "gettraitdict(traits, traitnames)\n\nConstruct a trait dict from a list of Trait objects.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.identifyAdults!-Tuple{GeMM.Patch}",
    "page": "Internal functions",
    "title": "GeMM.identifyAdults!",
    "category": "method",
    "text": "identifyAdults!(patch)\n\nBuild up the whoiswho index of individuals and species in a patch.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.iscompatible-Tuple{GeMM.Individual,GeMM.Individual,Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.iscompatible",
    "category": "method",
    "text": "iscompatible(mate, individual, traitnames)\n\nCheck to see whether two individual organisms are reproductively compatible.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.markthem!-Tuple{Array{GeMM.Individual,1}}",
    "page": "Internal functions",
    "title": "GeMM.markthem!",
    "category": "method",
    "text": "markthem!(community)\n\nSet each individual in the community (= array of individuals) as \"marked\".\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.markthem!-Tuple{Array{GeMM.Patch,1}}",
    "page": "Internal functions",
    "title": "GeMM.markthem!",
    "category": "method",
    "text": "markthem!(world)\n\nSet every individual in the world as \"marked\".\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.markthem!-Tuple{GeMM.Patch}",
    "page": "Internal functions",
    "title": "GeMM.markthem!",
    "category": "method",
    "text": "markthem!(habitat)\n\nSet each individual in the given patch as \"marked\".\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.meiosis-Tuple{Array{GeMM.Chromosome,1},Bool}",
    "page": "Internal functions",
    "title": "GeMM.meiosis",
    "category": "method",
    "text": "meiosis(genome, maternal)\n\nCarry out meiosis on a genome (marked as maternal or not). Returns a haploid gamete genome. (genome => array of chromosomes)\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.num2seq-Tuple{Integer}",
    "page": "Internal functions",
    "title": "GeMM.num2seq",
    "category": "method",
    "text": "num2seq(n)\n\nConvert an integer into binary and then into a DNA base sequence string.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.seq2bignum-Tuple{String}",
    "page": "Internal functions",
    "title": "GeMM.seq2bignum",
    "category": "method",
    "text": "seq2bignum(sequence)\n\nConvert a DNA base sequence (a string) into binary and then into an BigInt (for larger genes). This saves memory.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.seq2num-Tuple{String}",
    "page": "Internal functions",
    "title": "GeMM.seq2num",
    "category": "method",
    "text": "seq2num(sequence)\n\nConvert a DNA base sequence (a string) into binary and then into an integer. This saves memory.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.traitsexist-Tuple{Dict{String,Float64},Dict{String,Any}}",
    "page": "Internal functions",
    "title": "GeMM.traitsexist",
    "category": "method",
    "text": "traitsexist(traits, settings)\n\nCheck a trait dict to make sure it contains the full set of traitnames required by the model (as defined in the settings).\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.traitsexist-Tuple{GeMM.Individual,Dict{String,Any}}",
    "page": "Internal functions",
    "title": "GeMM.traitsexist",
    "category": "method",
    "text": "traitsexist(individual, settings)\n\nMake sure an individual organism has the full set of traits required by the model (as defined in the settings).\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.varyalleles!-Tuple{Array{GeMM.AbstractGene,1},Dict{String,Any},Float64}",
    "page": "Internal functions",
    "title": "GeMM.varyalleles!",
    "category": "method",
    "text": "varyalleles!(genes, settings, locivar)\n\nMutate gene traits in the passed array of genes.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#GeMM.varyalleles!-Tuple{Array{GeMM.Chromosome,1},Dict{String,Any},Float64}",
    "page": "Internal functions",
    "title": "GeMM.varyalleles!",
    "category": "method",
    "text": "varyalleles!(chromosomes, settings, locivar)\n\nMutate gene traits in the passed array of chromosomes.\n\n\n\n\n\n"
},

{
    "location": "aux/index.html#auxfuncts.jl-1",
    "page": "Internal functions",
    "title": "auxfuncts.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"auxfuncts.jl\"]"
},

]}
