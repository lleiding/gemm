# Scheduling of processes in GeMM

"""
    simulate!(world, timesteps)

This is the central function of the model with the main event loop. It defines
the scheduling for all submodels and output functions.
"""
function simulate!(world::Array{Patch,1}, timesteps::Int=1000, timeoffset::Int = 0)
    simlog("Starting simulation.")
    checkviability!(world)
    for t in (timeoffset + 1):(timeoffset + timesteps)
        simlog("UPDATE $t")
        # ecological processes are outsourced to specialised methods below
        if setting("mode") == "default"
            defaultexperiment(world)
        elseif setting("mode") == "invasion"
            invasionexperiment(world, t)
        elseif setting("mode") == "zosterops"
            zosteropsexperiment(world)
        else
            simlog("Mode setting not recognised: $(setting("mode"))", 'e')
        end
        if setting("lineages")
            recordstatistics(world)
            recordlineages(world, t)
        end
        if mod(t, setting("outfreq")) == 0 && any([setting("fasta") != "off", setting("raw"), setting("stats")])
            writedata(world, t)
        end
    end
end

"""
    defaultexperiment(world)

The standard annual update procedure, designed primarily for plant communities.
"""
function defaultexperiment(world::Array{Patch,1})
    establish!(world, setting("nniches"), setting("static"))
    survive!(world)
    grow!(world)
    compete!(world, setting("static"))
    reproduce!(world)
    if setting("mutate")
        mutate!(world)
    end
    disperse!(world, setting("static"))
    checkviability!(world)
    changehabitat!(world) # model output
end

"""
    invasionexperiment(world)

The annual update procedure for the invasion experiments.
"""
function invasionexperiment(world::Array{Patch,1}, t::Int)
    establish!(world, setting("nniches"), setting("static"))
    survive!(world)
    grow!(world)
    compete!(world, setting("static"))
    reproduce!(world)
    if 0 < setting("burn-in") < t
        disturb!(world)
        invade!(world)
    end
    disperse!(world, setting("static"))
    checkviability!(world)
end

"""
    zosteropsexperiment(world)

The annual update procedure for the Zosterops experiments, this time for bird populations.
"""
function zosteropsexperiment(world::Array{Patch,1})
    #TODO
    establish!(world, setting("nniches"), setting("static"))
    survive!(world)
    reproduce!(world)
    if setting("mutate")
        mutate!(world)
    end
    zdisperse!(world)
    checkviability!(world)
    #changehabitat!(world) # model output
end

# TODO add more modes for gradient/habitat change experiments?
