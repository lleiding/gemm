#!/usr/bin/python3
#
# Kick off my island experiments
# Daniel Vedder, 19/6/2018
#

import os, sys, time

global simname, replicates

simname = "experiment"
replicates = 1

# First commandline arg gives the simulation name, the second the number of replicates
if len(sys.argv) >= 2:
    simname = sys.argv[1]
if len(sys.argv) >= 3:
    replicates = int(sys.argv[2])

# These settings stay constant throughout all simulation runs
constant_settings = {"seed":0,
                     "maps":'"invasion.map"',
                     "outfreq":10,
                     "logging":"true",
                     "debug":"false",
                     "fasta":"false",
                     "lineages":"true",
                     "linkage":'"none"',
                     "static":"false",
                     "nniches":2,
                     "tolerance":'"none"',
                     "static":"false",
                     "mutate":"false",
                     "initadults":"true",
                     "metabolicpopsize":"false",
                     "global-species-pool":100}

# These settings are varied (the first value is the default,
# every combination of the rest is tested)
varying_settings = {"cellsize":[20, 10, 100],
                    "propagule-pressure":[0,1,10],
                    "disturbance":[0,1,10]}

def slurm(config):
    "Send a job to slurm"
    cmd = "sbatch -c 2 --mem 50GB ./islandsim.jl --config "+str(config)
    os.system(cmd)
    time.sleep(3) #prevent output folder merges

def write_config(config, cellsize, prop_pressure, disturbance):
    "Write out a config file with the given values"
    global simname
    cf = open(config, 'w')
    cf.write("# Island speciation model for invasion experiments\n")
    cf.write("# This config file was generated automatically.\n")
    cf.write("# "+time.asctime()+"\n")
    cf.write('\ndest "results/'+simname+'"\n')
    cf.write("\n# Constant settings:\n")
    for k in constant_settings.keys():
        cf.write(k + " " + str(constant_settings[k]) + "\n")
    cf.write("\n# Variable settings:\n")
    cf.write("cellsize "+str(cellsize)+"\n")
    cf.write("propagule-pressure "+str(prop_pressure)+"\n")
    cf.write("disturbance "+str(disturbance)+"\n")
    cf.close()
    
def run_defaults():
    "Create a series of runs with the default values (no invasion events)"
    global simname, replicates
    print("Running default simulation with "+str(replicates)+" replicates.")
    write_config(simname+".conf",
                 varying_settings["cellsize"][0],
                 varying_settings["propagule-pressure"][0],
                 varying_settings["disturbance"][0])
    i = 0
    while i < replicates:
        slurm(simname+".conf")
        i = i+1

def run_experiment():
    "Create a full experiment with all parameter combinations"
    global simname, replicates
    for cs in varying_settings["cellsize"][1:]:
        for pp in varying_settings["propagule-pressure"][1:]:
            for db in varying_settings["disturbance"][1:]:
                spec = str(cs)+"CS_"+str(pp)+"PP_"+str(db)+"DB"
                print("Running simulation with specification "+spec+" for "
                      +str(replicates)+" replicates.")
                simname = simname+"_"+spec
                write_config(simname+".conf", cs, pp, db)
                i = 0
                while i < replicates:
                    slurm(simname+".conf")
                    i = i + 1
                os.system("rm "+simname+".conf") #cleanup
    print("Done.")

if __name__ == '__main__':
    if "default" in simname:
        run_defaults()
    else:
        run_experiment()
