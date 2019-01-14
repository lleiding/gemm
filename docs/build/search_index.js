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
    "text": "This is the documentation for the Genetically explicit Metacommunity Model.(c) 2016-2018 Ludwig Leidinger, Daniel Vedder"
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
    "text": ""
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
    "text": ""
},

{
    "location": "aux/#GeMM.checkborderconditions-Tuple{Array{Patch,1},Int64,Int64}",
    "page": "Internal functions",
    "title": "GeMM.checkborderconditions",
    "category": "method",
    "text": "checkborderconditions!(w, x, y)\n\ncheck if coordinates x and y lie within world w and correct if not, considering defined border conditions.\n\n\n\n\n\n"
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
    "location": "aux/#GeMM.num2seq-Tuple{Integer}",
    "page": "Internal functions",
    "title": "GeMM.num2seq",
    "category": "method",
    "text": "num2seq(n)\n\nConvert an integer into binary and then into a DNA base sequence string.\n\n\n\n\n\n"
},

{
    "location": "aux/#GeMM.seq2num-Tuple{String}",
    "page": "Internal functions",
    "title": "GeMM.seq2num",
    "category": "method",
    "text": "seq2num(sequence)\n\nConvert a DNA base sequence (a string) into binary and then into an integer. This saves memory.\n\n\n\n\n\n"
},

{
    "location": "aux/#auxfuncts.jl-1",
    "page": "Internal functions",
    "title": "auxfuncts.jl",
    "category": "section",
    "text": "Modules = [GeMM]\nPages = [\"auxfuncts.jl\"]"
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
    "text": ""
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

]}
