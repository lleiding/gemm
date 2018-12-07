# Scheduling of processes in GeMM

function simulate!(world::Array{Patch,1}, settings::Dict{String, Any}, timesteps::Int=1000)
    simlog("Starting simulation.", settings)
    checkviability!(world, settings)
    for t in 1:timesteps
        simlog("UPDATE $t", settings)
        # ecological processes
        establish!(world, settings["nniches"], settings["static"])
        survive!(world, settings)
        disturb!(world, settings)
        grow!(world, settings)
        compete!(world, settings["static"])
        reproduce!(world, settings)
        if settings["mutate"]
            mutate!(world, settings)
        end
        checkviability!(world, settings)
        if 0 < settings["burn-in"] < t
            invade!(world, settings)
        end
        disperse!(world, settings["static"])
        changehabitat!(world, settings) # model output
        # map(p -> simlog("Patch $(p.id): $(length(p.community)) individuals.", settings, 'd'), world) # disabled - performance concerns
        if settings["lineages"]
            recordstatistics(world, settings)
            recordlineages(world, settings, t)
        end
        if mod(t, settings["outfreq"]) == 0
            writedata(world, settings, t)
        end
        if settings["phylo"]
            writephylo(world, settings, t)
        end
    end
end

