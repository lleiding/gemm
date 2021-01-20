"""
    runsim(config, seed, prerun)

Performs a simulation run using configuration file `config`, random seed `seed`
and other settings provided via commandline, configuration file or the defaults.
"""
function runsim(config::String = "", seed::Integer = 0)
    initsettings(defaultSettings()) #needed for log calls during `getsettings()`
    initsettings(getsettings(config, seed))
    Random.seed!(setting("seed"))
    setupdatadir()
    world = Patch[]
    timesteps = 0
    timeoffset = 0
    correctmode!()
    for i in 1:length(setting("maps"))
        timeoffset += timesteps
        timesteps, maptable = readmapfile(setting("maps")[i])
        i == 1 && (world = createworld(maptable))
        i > 1 && updateworld!(world, maptable)
        i == 1 && writedata(world, timeoffset)
        simulate!(world, timesteps, timeoffset)
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
    correctmode!()

A heuristic function to keep backward compatibility. (The `mode` setting was
only introduced after several different experiments had already been carried
out with GeMM, so old config files don't include it.)
"""
function correctmode!()
    if setting("mode") == "default"
        # detect invasion mode
        if setting("global-species-pool") > 0 || setting("propagule-pressure") > 0
            simlog("Detected implicit invasion mode, updated settings.", 'w')
            updatesetting("mode", "invasion")
        end
        # XXX further modes?
    end
end

let settings::Dict{String, Any}
    """
        initsettings(newsettings)

    Define a new settings dict.
    """
    global function initsettings(newsettings::Dict{String, Any})
        settings = newsettings
    end

    """
        setting(param)

    Return a configuration parameter from the settings.
    """
    global function setting(param::String)
        settings[param]
    end

    """
        updatesetting(param, value)

    Change the value of an individual config parameter.
    Use with caution!
    """
    global function updatesetting(param::String, value::Any)
        settings[param] = value
    end

    """
        settingkeys()

    Return all keys in the settings dict.
    """
    global function settingkeys()
        keys(settings)
    end
end
