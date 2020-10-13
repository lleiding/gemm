# Scheduling of processes in GeMM

"""
    simulate!(world, settings, timesteps)

This is the central function of the model with the main event loop. It defines
the scheduling for all submodels and output functions.
"""
function simulate!(world::Array{Patch,1}, settings::Dict{String, Any}, timesteps::Int=1000, timeoffset::Int = 0)
    simlog("Starting simulation.", settings)
    checkviability!(world, settings)
    for t in (timeoffset + 1):(timeoffset + timesteps)
        simlog("UPDATE $t", settings)
        # ecological processes are outsourced to specialised methods below
        if settings["mode"] == "default"]
            defaultexperiment(world, settings)
        elseif settings["mode"] == "invasion"
            invasionexperiment(world, settings)
        elseif settings["mode"] == "zosterops"
            zosteropsexperiment(world, settings)
        else
            simlog("Mode setting not recognised: $(settings["mode"])", settings, 'e')
        end
        # map(p -> simlog("Patch $(p.id): $(length(p.community)) individuals.", settings, 'd'), world) # disabled - performance concerns
        if settings["lineages"]
            recordstatistics(world, settings)
            recordlineages(world, settings, t)
        end
        if mod(t, settings["outfreq"]) == 0 && any([settings["fasta"] != "off", settings["raw"], settings["stats"]])
            writedata(world, settings, t)
        end
    end
end

"""
    defaultexperiment(world, settings)

The standard annual update procedure, designed primarily for plant communities.
"""
function defaultexperiment(world::Array{Patch,1}, settings::Dict{String, Any})
    establish!(world, settings["nniches"], settings["static"])
    survive!(world, settings)
    grow!(world, settings)
    compete!(world, settings["static"])
    reproduce!(world, settings)
    if settings["mutate"]
        mutate!(world, settings)
    end
    disperse!(world, settings["static"])
    checkviability!(world, settings)
    changehabitat!(world, settings) # model output
end

"""
    invasionexperiment(world, settings)

The annual update procedure for the invasion experiments.
"""
function invasionexperiment(world::Array{Patch,1}, settings::Dict{String, Any})
    establish!(world, settings["nniches"], settings["static"])
    survive!(world, settings)
    grow!(world, settings)
    compete!(world, settings["static"])
    reproduce!(world, settings)
    if 0 < settings["burn-in"] < t
        disturb!(world, settings)
        invade!(world, settings)
    end
    disperse!(world, settings["static"])
    checkviability!(world, settings)
end

"""
    zosteropsexperiment(world, settings)

The annual update procedure for the Zosterops experiments, this time for bird populations.
"""
function zosteropsexperiment(world::Array{Patch,1}, settings::Dict{String, Any})
    #TODO
    establish!(world, settings["nniches"], settings["static"])
    survive!(world, settings)
    grow!(world, settings)
    compete!(world, settings["static"])
    reproduce!(world, settings)
    if settings["mutate"]
        mutate!(world, settings)
    end
    disperse!(world, settings["static"])
    checkviability!(world, settings)
    changehabitat!(world, settings) # model output
end

# FIXME add more modes for gradient/habitat change experiments?
