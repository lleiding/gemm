# Scheduling of processes in GeMM

"""
    simulate!(world, settings, timesteps)

This is the central function of the model with the main event loop. It defines
the scheduling for all submodels and output functions.
"""
function simulate!(world::Array{Patch,1}, settings::Dict{String, Any}, timesteps::Int=1000)
    simlog("Starting simulation.", settings)
    checkviability!(world, settings)
    for t in 1:timesteps
        simlog("UPDATE $t", settings)
        # ecological processes
        establish!(world, settings["nniches"], settings["static"])
        survive!(world, settings)
        grow!(world, settings)
        compete!(world, settings["static"])
        phylo = reproduce!(world, settings)
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
        if settings["writephylo"]
            writephylo(phylo, settings, t)
        end
        if mod(t, settings["outfreq"]) == 0
            writedata(world, settings, t)
        end
    end
end

