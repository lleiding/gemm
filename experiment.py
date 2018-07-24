#!/usr/bin/python3
#
# Kick off my island experiments
# Daniel Vedder, 19/6/2018
#

import os, sys, time, random

global simname, replicates

simname = "experiment"
replicates = 1

# First commandline arg gives the simulation name, the second the number of replicates.
# If the simname contains the string "default", or "default" is appended as a fourth
# argument, the default simulation is run. Otherwise, an invasion experiment is set up.
if len(sys.argv) >= 2:
    simname = sys.argv[1]
if len(sys.argv) >= 3:
    replicates = int(sys.argv[2])

# These settings stay constant throughout all simulation runs
constant_settings = {"cellsize":1,
                     "outfreq":50,
                     "logging":"true",
                     "debug":"true", #default: false
                     "fasta":"false",
                     "lineages":"true",
                     "linkage":'"none"',
                     "static":"false",
                     "nniches":2,
                     "tolerance":'"none"',
                     "static":"false",
                     "mutate":"false",
                     "initadults":"true",
                     "initpopsize":'"metabolic"',
                     "burn-in": 1000,
                     "global-species-pool":100}

# These settings are varied (the first value is the default,
# every combination of the rest is tested)
varying_settings = {"maps":["invasion.map",
                            "invasion_hot.map",
                            "invasion_cold.map"],
                    "propagule-pressure":[0,1,10],
                    "disturbance":[0,1,10]}

def slurm(config):
    "Send a job to slurm"
    if "gaia" in os.uname().nodename:
        cmd = "sbatch -c 2 --mem 50GB ./islandsim.jl --config "+str(config)
        os.system(cmd)
        os.system("rm "+config) #cleanup
        time.sleep(3) #prevent output folder merges
    else:
        print("Not on gaia, slurm is probably not available. Retaining job "+config)

def write_config(config, maps, prop_pressure, disturbance, seed):
    "Write out a config file with the given values"
    cf = open(config, 'w')
    cf.write("# Island speciation model for invasion experiments\n")
    cf.write("# This config file was generated automatically.\n")
    cf.write("# "+time.asctime()+"\n")
    cf.write('\ndest "results/'+config.replace(".conf", "")+'"\n')
    cf.write("\n# Constant settings:\n")
    for k in constant_settings.keys():
        cf.write(k + " " + str(constant_settings[k]) + "\n")
    cf.write("\n# Variable settings:\n")
    cf.write("seed "+str(seed)+"\n")
    cf.write("maps "+str(maps)+"\n")
    cf.write("propagule-pressure "+str(prop_pressure)+"\n")
    cf.write("disturbance "+str(disturbance)+"\n")
    cf.close()
    
def run_defaults():
    "Create a series of runs with the default values (no invasion events)"
    global simname, replicates
    print("Running default simulation with "+str(replicates)+" replicates.")
    write_config(simname+".conf",
                 varying_settings["maps"][0],
                 varying_settings["propagule-pressure"][0],
                 varying_settings["disturbance"][0], 0)
    i = 0
    while i < replicates:
        slurm(simname+".conf")
        i = i+1

def run_experiment():
    "Create a full experiment with all parameter combinations"
    global simname, replicates
    i = 0
    while i < replicates:
        for tm in varying_settings["maps"]:
            for pp in varying_settings["propagule-pressure"][1:]:
                for db in varying_settings["disturbance"][1:]:
                    if '_' in tm: temp = tm.split("_")[1].split(".")[0]
                    else: temp = "default"
                    spec = temp+"_"+str(pp)+"PP_"+str(db)+"DB"
                    print("Running simulation with specification "+spec+" for "
                          +str(replicates)+" replicates.")
                    seed = random.randint(0,10000)
                    runname = simname+"_r"+str(i+1)+"_"+spec
                    write_config(runname+".conf", tm, pp, db, seed)
                    write_config(runname+"_control.conf", tm, 0, db, seed)
                    slurm(runname+".conf")
                    slurm(runname+"_control.conf")
        i = i + 1
    print("Done.")

if __name__ == '__main__':
    if "default" in simname or (len(sys.argv) >= 4 and sys.argv[3] == "default"):
        run_defaults()
    else:
        run_experiment()
