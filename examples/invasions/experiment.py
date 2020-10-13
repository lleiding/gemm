#!/usr/bin/python3
#
# Kick off my island experiments
# Daniel Vedder, 19/6/2018
#        updated 22/1/2019
#

# NOTE: make sure to copy/symlink this to the model root folder before running

import os, sys, shutil, time, random

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
constant_settings = {"cellsize":2e6, # 2 tonnes/hectar -> default 20?
                     "outfreq":50,
                     "raw":"false",
                     "logging":"true",
                     "debug":"false", #default: false
                     "fasta":"false",
                     "lineages":"true",
                     "linkage":"none",
                     "static":"false",
                     "nniches":2,
                     "tolerance":1.0,
                     "static":"false",
                     "mutate":"false",
                     "usebiggenes":"false",
                     "indsize":"mixed",
                     "popsize":"metabolic",
                     "minseedsize":0, # 1g
                     "minrepsize":5,  # 148g
                     "mode":"invasion",
                     "burn-in": 500,
                     "global-species-pool":100}

# These settings are varied (the first value is the default,
# every combination of the rest is tested)
varying_settings = {"maps":["invasion.map",
                            "invasion_hot.map",
                            "invasion_cold.map"],
                    "propagule-pressure":[0,1,10],
                    "disturbance":[0,1,10]}

def archive_code():
    "Save the current codebase in a tar archive."
    tarname = time.strftime("codebase_%d%b%y.tar.gz")
    print("Archiving codebase in "+tarname)
    os.system("git log -1 > commit.txt")
    cmd = "tar czfh "+tarname+" README.md commit.txt rungemmparallel.jl experiment.py analyse.R src/*"
    os.system(cmd)
    os.remove("commit.txt")

def slurm(config):
    "Send a job to slurm"
    if "gaia" in os.uname().nodename:
        cmd = "sbatch -c 2 --mem 50GB ./run_simulation.jl --config "+str(config)
        os.system(cmd)
        os.system("rm "+config) #cleanup
        time.sleep(3) #prevent output folder merges
    else:
        print("Not on gaia, slurm is probably not available. Retaining job "+config)

def write_config(config, maps, mintemp, prop_pressure, disturbance, seed):
    "Write out a config file with the given values"
    cf = open(config, 'w')
    cf.write("# Island speciation model for invasion experiments\n")
    cf.write("# This config file was generated automatically.\n")
    cf.write("# "+time.asctime()+"\n")
    cf.write("\ndest results/"+config.replace(".conf", "")+"\n")
    cf.write("\n# Constant settings:\n")
    for k in constant_settings.keys():
        cf.write(k + " " + str(constant_settings[k]) + "\n")
    cf.write("\n# Variable settings:\n")
    cf.write("seed "+str(seed)+"\n")
    cf.write("maps "+str(maps)+"\n")
    if mintemp > 0:
        cf.write("mintemp "+str(mintemp)+"\n")
        cf.write("maxtemp "+str(mintemp+15)+"\n")
    cf.write("propagule-pressure "+str(prop_pressure)+"\n")
    cf.write("disturbance "+str(disturbance)+"\n")
    cf.close()
    
def run_defaults():
    "Create a series of runs with the default values (no invasion events)"
    global simname, replicates
    print("Running default simulation with "+str(replicates)+" replicates.")
    if varying_settings["maps"][0] not in os.listdir():
        shutil.copy("examples/invasions/"+varying_settings["maps"][0], ".")
    write_config(simname+".conf",
                 varying_settings["maps"][0],
                 varying_settings["propagule-pressure"][0],
                 varying_settings["disturbance"][0], 0)
    i = 0
    while i < replicates:
        slurm(simname+".conf")
        i = i+1

def run_experiment(control=False):
    "Create a full experiment with all parameter combinations"
    global simname, replicates
    i = 0
    while i < replicates:
        seed = random.randint(0,100000)
        for tm in varying_settings["maps"][1:]:
            if tm not in os.listdir():
                shutil.copy("examples/invasions/"+tm, ".")
            # figure out the range of optimum temperature
            if "hot" in tm: mt = 298
            elif "cold" in tm: mt = 278
            elif "default" in tm: mt = 288
            else: mt = -1
            for pp in varying_settings["propagule-pressure"][1:]:
                for db in varying_settings["disturbance"][1:]:
                    if '_' in tm: temp = tm.split("_")[1].split(".")[0]
                    else: temp = "default"
                    spec = temp+"_"+str(pp)+"PP_"+str(db)+"DB"
                    print("Running simulation with specification "+spec+" for "
                          +str(replicates)+" replicates.")
                    runname = simname+"_r"+str(i+1)+"_"+spec
                    write_config(runname+".conf", tm, mt, pp, db, seed)
                    slurm(runname+".conf")
                    if control:
                        write_config(runname+"_control.conf", tm, mt, 0, db, seed)
                        slurm(runname+"_control.conf")
        i = i + 1
    print("Done.")

if __name__ == '__main__':
    archive_code()
    if simname == "archive":
        pass #only archive the code
    elif "default" in simname or (len(sys.argv) >= 4 and sys.argv[3] == "default"):
        run_defaults()
    else:
        run_experiment()
