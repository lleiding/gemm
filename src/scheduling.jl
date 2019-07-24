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
        # ecological processes
        establish!(world, settings["nniches"], settings["static"])
        survive!(world, settings)
        grow!(world, settings)
        compete!(world, settings["static"])
        reproduce!(world, settings)
        if settings["mutate"]
            mutate!(world, settings)
        end
        if 0 < settings["burn-in"] < t
            disturb!(world, settings)
            invade!(world, settings)
        end
        disperse!(world, settings["static"])
        checkviability!(world, settings)
        changehabitat!(world, settings) # model output
        # map(p -> simlog("Patch $(p.id): $(length(p.community)) individuals.", settings, 'd'), world) # disabled - performance concerns
        if settings["lineages"]
            recordstatistics(world, settings)
            recordlineages(world, settings, t)
        end
        if mod(t, settings["outfreq"]) == 0 && any([settings["fasta"], settings["raw"], settings["stats"]])
            writedata(world, settings, t)
            updateorigin!(world)
        end
    end
end

