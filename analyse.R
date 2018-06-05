#!/usr/bin/Rscript
# This is an analysis script for the island invasion model that creates
# graphs of the total and per species populations over time.

# The working directory may be specified via the commandline, otherwise it
# defaults to results/tests
outdir = paste("results/", commandArgs()[length(commandArgs())], sep="")
if (!(file.exists(outdir) && file.info(outdir)$isdir)) {
    outdir = "results/tests"
}

# Read the population sizes from a logfile
readPopSizes = function(logfile) {
    log = readLines(logfile)
    popsizes = c()
    for (l in log) {
        if (grepl("UPDATE", l)[1]) {
            w = strsplit(l, " ")[[1]]
            popsizes = c(popsizes, as.numeric(w[length(w)]))
        }
    }
    return(popsizes)
}

# Plot the population sizes over time
plotPopSizes = function(logfile, toFile=TRUE) {
    popsizes = readPopSizes(logfile)
    if (toFile) {
        jpeg(paste(outdir, "population.jpg", sep="/"), height=max(popsizes)*0.25,
             width=length(popsizes)*20)
    }
    plot(popsizes, xlab="Time", ylab="Population size", ylim=c(0, max(popsizes)),
         col="red", type='l')
    if (toFile) dev.off()
}

plotPopSizes(paste(outdir, "simulation.log", sep="/"))
