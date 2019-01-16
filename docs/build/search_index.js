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
    "text": "This is the documentation for the Island Speciation Model, a.k.a. Genetically explicit Metacommunity Model (GeMM). (Yeah, we\'re still working on the name...)The aim of this model is to create a virtual island ecosystem that can be used to explore ecological and evolutionary hypotheses in silico. It is genetically and spatially explicit, with discrete space and time.This documentation is generated from the source code using Julia\'s inbuilt Documenter module. It is sorted first by functionality, second by source code file.(c) 2016-2018 Ludwig Leidinger, Daniel Vedder"
},

{
    "location": "model/#",
    "page": "Main Model",
    "title": "Main Model",
    "category": "page",
    "text": ""
},

{
    "location": "model/#Main-Model-1",
    "page": "Main Model",
    "title": "Main Model",
    "category": "section",
    "text": "These are the core functions of the model, implementing the high-level processes of the simulation and the object classes that are involved."
},

{
    "location": "model/#GeMM.Patch",
    "page": "Main Model",
    "title": "GeMM.Patch",
    "category": "type",
    "text": "One of the core structs of the model, representing a one-hectare patch of ground.\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.AbstractGene",
    "page": "Main Model",
    "title": "GeMM.AbstractGene",
    "category": "type",
    "text": "An abstraction layer to hide the implementation details of genes.\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.BigGene",
    "page": "Main Model",
    "title": "GeMM.BigGene",
    "category": "type",
    "text": "A variation of the standard Gene struct that allows for longer base sequence lengths through its use of BigInt.\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.Chromosome",
    "page": "Main Model",
    "title": "GeMM.Chromosome",
    "category": "type",
    "text": "A struct that saves a collection of genes and a tag to denote whether they are maternal or not.\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.Gene",
    "page": "Main Model",
    "title": "GeMM.Gene",
    "category": "type",
    "text": "The standard gene struct; including the base sequence (represented in binary and stored and as an integer), and an array of traits the gene codes for.\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.Individual",
    "page": "Main Model",
    "title": "GeMM.Individual",
    "category": "type",
    "text": "One of the core structs of the model, representing an individual organism.\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.Trait",
    "page": "Main Model",
    "title": "GeMM.Trait",
    "category": "type",
    "text": "A simple struct to store organism trait names and values.\n\n\n\n\n\n"
},

{
    "location": "model/#entities.jl-1",
    "page": "Main Model",
    "title": "entities.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"entities.jl\"]"
},

{
    "location": "model/#initialisation.jl-1",
    "page": "Main Model",
    "title": "initialisation.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"initalisation.jl\"]"
},

{
    "location": "model/#GeMM.disperse!",
    "page": "Main Model",
    "title": "GeMM.disperse!",
    "category": "function",
    "text": "disperse!(w)\n\nDispersal of individuals within world (array of patches) w\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.disturb!-Tuple{Patch,Int64}",
    "page": "Main Model",
    "title": "GeMM.disturb!",
    "category": "method",
    "text": "disturb!(p,i)\n\nspecies-independent mortality due to disturbance on patch p\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.establish!",
    "page": "Main Model",
    "title": "GeMM.establish!",
    "category": "function",
    "text": "establish!(p, n)\n\nestablishment of individuals in patch p. Sets adaption parameter according to adaptation to number n niches of the surrounding environment.\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.grow!-Tuple{Patch,Float64}",
    "page": "Main Model",
    "title": "GeMM.grow!",
    "category": "method",
    "text": "grow!(p)\n\nGrowth of individuals in patch p\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.reproduce!-Tuple{Patch,Dict{String,Any}}",
    "page": "Main Model",
    "title": "GeMM.reproduce!",
    "category": "method",
    "text": "reproduce!(p, settings)\n\nReproduction of individuals in a patch p\n\n\n\n\n\n"
},

{
    "location": "model/#GeMM.survive!-Tuple{Patch,Float64}",
    "page": "Main Model",
    "title": "GeMM.survive!",
    "category": "method",
    "text": "survive!(p)\n\ndensity independent survival of individuals in patch p\n\n\n\n\n\n"
},

{
    "location": "model/#submodels.jl-1",
    "page": "Main Model",
    "title": "submodels.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"submodels.jl\"]"
},

{
    "location": "model/#scheduling.jl-1",
    "page": "Main Model",
    "title": "scheduling.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"scheduling.jl\"]"
},

{
    "location": "io/#",
    "page": "Input, Output, and Settings",
    "title": "Input, Output, and Settings",
    "category": "page",
    "text": ""
},

{
    "location": "io/#Input,-Output,-and-Settings-1",
    "page": "Input, Output, and Settings",
    "title": "Input, Output, and Settings",
    "category": "section",
    "text": "These functions are responsible for reading in all model configurations (passed by config file or commandline), administrating them during a run, and printing any output."
},

{
    "location": "io/#constants.jl-1",
    "page": "Input, Output, and Settings",
    "title": "constants.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"constants.jl\"]"
},

{
    "location": "io/#defaults.jl-1",
    "page": "Input, Output, and Settings",
    "title": "defaults.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"defaults.jl\"]"
},

{
    "location": "io/#GeMM.basicparser-Tuple{String}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.basicparser",
    "category": "method",
    "text": "basicparser(f)\n\nDo elementary parsing on a config or map file.\n\nReads in the file, strips whole-line and inline comments and separates lines by whitespace. Returns a 2d array representing the tokens in each line.\n\n\n\n\n\n"
},

{
    "location": "io/#input.jl-1",
    "page": "Input, Output, and Settings",
    "title": "input.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"input.jl\"]"
},

{
    "location": "io/#GeMM.setupdatadir-Tuple{Dict{String,Any}}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.setupdatadir",
    "category": "method",
    "text": "setupdatadir(dir)\n\nCreates the output directory and copies relevant files into it. If the output directory already includes files, create a new directory by appending a counter.\n\n\n\n\n\n"
},

{
    "location": "io/#GeMM.writedata-Tuple{Array{Patch,1},Dict{String,Any},Int64}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.writedata",
    "category": "method",
    "text": "writedata(world, settings, timestep)\n\nwrites simulation output from world to separate table and fasta files. timestep and setting information is used for file name creation.\n\n\n\n\n\n"
},

{
    "location": "io/#GeMM.dumpinds",
    "page": "Input, Output, and Settings",
    "title": "GeMM.dumpinds",
    "category": "function",
    "text": "dumpinds(world, io, sep)\n\nOutput all data of individuals in world as table to io. Columns are separated by sep.\n\n\n\n\n\n"
},

{
    "location": "io/#GeMM.recordlineages-Tuple{Array{Patch,1},Dict{String,Any},Int64}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.recordlineages",
    "category": "method",
    "text": "recordlineages(w)\n\nSave the abundance of each lineage per patch\n\n\n\n\n\n"
},

{
    "location": "io/#GeMM.recordstatistics-Tuple{Array{Patch,1},Dict{String,Any}}",
    "page": "Input, Output, and Settings",
    "title": "GeMM.recordstatistics",
    "category": "method",
    "text": "recordstatistics(w)\n\nWrite out world properties to the log file for later analysis.\n\n\n\n\n\n"
},

{
    "location": "io/#GeMM.simlog",
    "page": "Input, Output, and Settings",
    "title": "GeMM.simlog",
    "category": "function",
    "text": "simlog(msg, settings, category)\n\nWrite a log message to STDOUT/STDERR and the specified logfile  (if logging is turned on in the settings). Categories: d (debug), i (information, default), w (warn), e (error)\n\n\n\n\n\n"
},

{
    "location": "io/#output.jl-1",
    "page": "Input, Output, and Settings",
    "title": "output.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"output.jl\"]"
},

{
    "location": "aux/#",
    "page": "Internal functions",
    "title": "Internal functions",
    "category": "page",
    "text": ""
},

{
    "location": "aux/#Internal-functions-1",
    "page": "Internal functions",
    "title": "Internal functions",
    "category": "section",
    "text": "These are functions that the model uses in the background. They include a range of things: sub-submodels such findposspartner, utility functions like  gausscurve, or performance hacks like seq2num."
},

{
    "location": "aux/#GeMM.checkborderconditions-Tuple{Array{Patch,1},Int64,Int64}",
    "page": "Internal functions",
    "title": "GeMM.checkborderconditions",
    "category": "method",
    "text": "checkborderconditions!(w, x, y)\n\ncheck if coordinates x and y lie within world w and correct if not, considering defined border conditions.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.createchrms-Tuple{Int64,Array{GeMM.AbstractGene,1}}",
    "page": "Internal functions",
    "title": "GeMM.createchrms",
    "category": "method",
    "text": "createchrms(nchrms, genes)\n\nRandomly distribute the passed genes into the given number of chromosomes. Returns an array of chromosome objects.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.creategenes-Tuple{Int64,Array{GeMM.Trait,1},Dict{String,Any}}",
    "page": "Internal functions",
    "title": "GeMM.creategenes",
    "category": "method",
    "text": "creategenes(ngenes, traits, settings)\n\nRandomly create a given number of gene objects, with their base sequence and associated traits. Returns the result as an array of AbstractGenes.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.createind",
    "page": "Internal functions",
    "title": "GeMM.createind",
    "category": "function",
    "text": "createind(settings, marked=false)\n\nCreate an individual organism with of a new species with a random genome.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.createoffspring-Tuple{Integer,GeMM.Individual,GeMM.Individual,Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.createoffspring",
    "category": "method",
    "text": "createoffspring(noffs, individual, partner, traitnames)\n\nThe main reproduction function. Take two organisms and create the given number of offspring individuals. Returns an array of individuals.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.createtraits-Tuple{Dict{String,Any}}",
    "page": "Internal functions",
    "title": "GeMM.createtraits",
    "category": "method",
    "text": "createtraits(settings)\n\nCreate an array of trait objects generated from the default trait values (with a random offset).\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.diversity-Tuple{Array{Patch,1}}",
    "page": "Internal functions",
    "title": "GeMM.diversity",
    "category": "method",
    "text": "diversity(w)\n\nCalculate (average) alpha, beta and gamma diversity of the world. Returns a tuple with the three values (a,b,g). cf. Veech et al. 2002\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.findisland-Tuple{Array{Patch,1}}",
    "page": "Internal functions",
    "title": "GeMM.findisland",
    "category": "method",
    "text": "findisland(w)\n\nwithin world w, find out in which direction from the continent the island(s) lie(s).\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.findposspartner-Tuple{Patch,GeMM.Individual,Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.findposspartner",
    "category": "method",
    "text": "findposspartner(patch, individual, traitnames)\n\nFind a reproduction partner for the given individual in the given patch.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.freespace-Tuple{Array{Patch,1}}",
    "page": "Internal functions",
    "title": "GeMM.freespace",
    "category": "method",
    "text": "freespace(world)\n\nCalculate the average amount of free space on each patch.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.gausscurve",
    "page": "Internal functions",
    "title": "GeMM.gausscurve",
    "category": "function",
    "text": "gausscurve(b, c, x, a=1.0)\n\nCalculate the value of the Gauss function (\"bell curve\") at point x; with a being the maximum height of the curve, b the position of the curve center and c the standard deviation (\"width\").\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.getmeantraitvalue-Tuple{Array{GeMM.Trait,1},Integer}",
    "page": "Internal functions",
    "title": "GeMM.getmeantraitvalue",
    "category": "method",
    "text": "getmeantraitvalue(traits, traitidx)\n\nTake an array of traits and return the mean value of the indexed trait.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.getstdtraitvalue-Tuple{Array{GeMM.Trait,1},Integer}",
    "page": "Internal functions",
    "title": "GeMM.getstdtraitvalue",
    "category": "method",
    "text": "getstdtraitvalue(traits, traitidx)\n\nTake an array of traits and return the standard deviation of the indexed trait.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.gettraitdict-Tuple{Array{GeMM.AbstractGene,1},Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.gettraitdict",
    "category": "method",
    "text": "gettraitdict(genes, traitnames)\n\nCalculate the trait dict for an array of genes.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.gettraitdict-Tuple{Array{GeMM.Chromosome,1},Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.gettraitdict",
    "category": "method",
    "text": "gettraitdict(chromosomes, traitnames)\n\nConvert a genome (an array of chromosomes) into a dict of traits and their values.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.gettraitdict-Tuple{Array{GeMM.Trait,1},Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.gettraitdict",
    "category": "method",
    "text": "gettraitdict(traits, traitnames)\n\nConstruct a trait dict from a list of Trait objects.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.identifyAdults!-Tuple{Patch}",
    "page": "Internal functions",
    "title": "GeMM.identifyAdults!",
    "category": "method",
    "text": "identifyAdults!(patch)\n\nBuild up the whoiswho index of individuals and species in a patch.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.iscompatible-Tuple{GeMM.Individual,GeMM.Individual,Array{String,1}}",
    "page": "Internal functions",
    "title": "GeMM.iscompatible",
    "category": "method",
    "text": "iscompatible(mate, individual, traitnames)\n\nCheck to see whether two individual organisms are reproductively compatible.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.markthem!-Tuple{Array{GeMM.Individual,1}}",
    "page": "Internal functions",
    "title": "GeMM.markthem!",
    "category": "method",
    "text": "markthem!(community)\n\nSet each individual in the community (= array of individuals) as \"marked\".\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.markthem!-Tuple{Array{Patch,1}}",
    "page": "Internal functions",
    "title": "GeMM.markthem!",
    "category": "method",
    "text": "markthem!(world)\n\nSet every individual in the world as \"marked\".\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.markthem!-Tuple{Patch}",
    "page": "Internal functions",
    "title": "GeMM.markthem!",
    "category": "method",
    "text": "markthem!(habitat)\n\nSet each individual in the given patch as \"marked\".\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.meiosis-Tuple{Array{GeMM.Chromosome,1},Bool}",
    "page": "Internal functions",
    "title": "GeMM.meiosis",
    "category": "method",
    "text": "meiosis(genome, maternal)\n\nCarry out meiosis on a genome (marked as maternal or not). Returns a haploid gamete genome. (genome => array of chromosomes)\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.num2seq-Tuple{Integer}",
    "page": "Internal functions",
    "title": "GeMM.num2seq",
    "category": "method",
    "text": "num2seq(n)\n\nConvert an integer into binary and then into a DNA base sequence string.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.seq2bignum-Tuple{String}",
    "page": "Internal functions",
    "title": "GeMM.seq2bignum",
    "category": "method",
    "text": "seq2bignum(sequence)\n\nConvert a DNA base sequence (a string) into binary and then into an BigInt (for larger genes). This saves memory.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.seq2num-Tuple{String}",
    "page": "Internal functions",
    "title": "GeMM.seq2num",
    "category": "method",
    "text": "seq2num(sequence)\n\nConvert a DNA base sequence (a string) into binary and then into an integer. This saves memory.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.traitsexist-Tuple{Dict{String,Float64},Dict{String,Any}}",
    "page": "Internal functions",
    "title": "GeMM.traitsexist",
    "category": "method",
    "text": "traitsexist(traits, settings)\n\nCheck a trait dict to make sure it contains the full set of traitnames required by the model (as defined in the settings).\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.traitsexist-Tuple{GeMM.Individual,Dict{String,Any}}",
    "page": "Internal functions",
    "title": "GeMM.traitsexist",
    "category": "method",
    "text": "traitsexist(individual, settings)\n\nMake sure an individual organism has the full set of traits required by the model (as defined in the settings).\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.varyalleles!-Tuple{Array{GeMM.AbstractGene,1},Dict{String,Any},Float64}",
    "page": "Internal functions",
    "title": "GeMM.varyalleles!",
    "category": "method",
    "text": "varyalleles!(genes, settings, locivar)\n\nMutate gene traits in the passed array of genes.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.varyalleles!-Tuple{Array{GeMM.Chromosome,1},Dict{String,Any},Float64}",
    "page": "Internal functions",
    "title": "GeMM.varyalleles!",
    "category": "method",
    "text": "varyalleles!(chromosomes, settings, locivar)\n\nMutate gene traits in the passed array of chromosomes.\n\n\n\n\n\n"
},

{
    "location": "aux/#auxfuncts.jl-1",
    "page": "Internal functions",
    "title": "auxfuncts.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"auxfuncts.jl\"]"
},

]}
