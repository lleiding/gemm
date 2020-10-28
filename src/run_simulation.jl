"""
    runsim(config, seed, prerun)

Performs a simulation run using configuration file `config`, random seed `seed`
and other settings provided via commandline, configuration file or the defaults.
"""
function runsim(config::String = "", seed::Integer = 0)
    settings = getsettings(config, seed)
    settings["seed"] = seed
    Random.seed!(settings["seed"])
    setupdatadir(settings)
    world = Patch[]
    timesteps = 0
    timeoffset = 0
    correctmode!(settings)
    for i in 1:length(settings["maps"])
        timeoffset += timesteps
        timesteps, maptable = readmapfile(settings["maps"][i], settings)
        i == 1 && (world = createworld(maptable, settings))
        i > 1 && updateworld!(world, maptable, settings)
        i == 1 && writedata(world, settings, timeoffset)
        simulate!(world, settings, timesteps, timeoffset)
    end
    world
end


"""
    rungemm(config, seed)

Wrapper for `runsim()`
Runs a simulation using configuration file `config`, random seed `seed`
and other settings provided via commandline, configuration file or the defaults.
"""
function rungemm(config::String = "", seed::Integer = 0)
    @time world = runsim(config, seed)
end

"""
    correctmode!(settings)

A heuristic function to keep backward compatibility. (The `mode` setting was
only introduced after several different experiments had already been carried
out with GeMM, so old config files don't include it.)
"""
function correctmode!(settings::Dict{String, Any})
    if settings["mode"] == "default"
        # detect invasion mode
        if settings["global-species-pool"] > 0 || settings["propagule-pressure"] > 0
            simlog("Detected implicit invasion mode, updated settings.", settings, 'w')
            settings["mode"] = "invasion"
        end
        # XXX further modes?
    end
end
