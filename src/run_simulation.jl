"""
    runsim(config, seed, prerun)

Performs a simulation run using configuration file `config`, random seed `seed`
and other settings provided via commandline, configuration file or the defaults.
`prerun` toggles a pre-compilation run
"""
function runsim(config::String = "", seed::Integer = 0, prerun::Bool = false)
    settings = getsettings(config, seed)
    settings["seed"] = seed
    if prerun
        settings = defaultSettings()
        settings["cellsize"] = 1.0e6
        settings["indsize"] = "adult"
        settings["popsize"] = "metabolic"
        settings["maps"] = [""]
        settings["phylo"] = false
        settings["fasta"] = false
        settings["raw"] = false
        settings["quiet"] = true
        settings["static"] = false
        settings["seed"] = 1
    end
    Random.seed!(settings["seed"])
    !prerun && setupdatadir(settings)
    world = Patch[]
    for i in 1:length(settings["maps"])
        i == 1 ? (timeoffset = 0) : (timeoffset = timesteps)
        timesteps, maptable = readmapfile(settings["maps"][i], settings)
        i == 1 && (world = createworld(maptable, settings))
        i > 1 && updateworld!(world, maptable, settings)
        !prerun && writedata(world, settings, timeoffset)
        simulate!(world, settings, timesteps, timeoffset)
    end
    world
end


"""
    rungemm(config, seed)

Wrapper for `runsim()`
Runs a simulation using configuration file `config`, random seed `seed`
and other settings provided via commandline, configuration file or the defaults.
Performs a pre-compilation run first.
"""
function rungemm(config::String = "", seed::Integer = 0)
    # compilation run:
    runsim("", 1, true)
    # run intended simulation:
    @time world = runsim(config, seed)
end
