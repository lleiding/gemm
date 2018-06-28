#!/usr/bin/Rscript
# This is an analysis script for the island invasion model that creates
# graphs of the total and per species populations over time.

# The working directory may be specified via the commandline, otherwise it
# defaults to results/tests
simname = commandArgs()[length(commandArgs())]
outdir = paste("results/", simname, sep="")
if (!(file.exists(outdir) && file.info(outdir)$isdir)) {
    simname = "tests"
    outdir = paste("results/", simname, sep="")
}


# Plot the distribution of traits in the population at the given timestep (-1 => END)
plotTraits = function(timestep=-1, toFile=TRUE) {
    print("Plotting traits...")
    traitfile = grep(paste("t", timestep, sep=""), list.files(outdir), value=T)
    if (length(traitfile) == 0) {
        # If the desired timestep doesn't exist, take the newest timestep we have
        timestep = (length(grep(".tsv", list.files(outdir), value=T))-1) * 10
        traitfile = grep(paste("t", timestep, sep=""), list.files(outdir), value=T)
    }
    traitfilepath = paste(outdir, traitfile, sep="/")
    ts = read.table(traitfilepath, header=T)
    jpeg(paste(outdir, "/", simname, "_traits.jpg", sep=""), height=720, width=1800)
    boxplot(ts$fitness*10, log(ts$size), log(ts$seedsize), log(ts$repsize), ts$lnkgunits/10,
            ts$ngenes/10, ts$temptol, ts$tempopt/100, ts$prectol, ts$precopt, ts$compat*10,
            ts$repradius*10, ts$dispmean, ts$dispshape,
            names=c("Fitness (x10)", "log(Size)", "log(Seed size)", "log(Reprod. size)",
                    "Chromosome (x0.1)", "Genes (x0.1)", "Temp. tolerance",
                    "Temp. opt. (x0.01)", "Precip. tolerance",
                    "Precip. optimum", "Compatibility (x10)", "Rep. radius (x10)",
                    "Dispersal mean", "Dispersal shape"))
    legend("top", c(paste("N = ", length(ts$counter))))
    dev.off()
}

# Plot population size and diversity indices over time
plotDiversity = function(logfile="diversity.log") {
    logfile = paste(outdir, logfile, sep="/")
    data = read.csv(logfile)
    data$lineages = data$lineages/10 #otherwise the Y axis is too big
    # Plot population sizes
    print("Plotting population development...")
    jpeg(paste(outdir, "/", simname, "_population.jpg", sep=""), height=720,
         width=length(data$population)*(2000/length(data$population)))
    plot(data$population, xlab="Time", ylab="Population size", ylim=c(0, max(data$population)),
         col="red", type='l')
    dev.off()
    # Plot diversity development
    print("Plotting diversity...")
    ymax = max(data$lineages, data$alpha, data$beta, data$gamma, data$freespace)
    jpeg(paste(outdir, "/", simname, "_diversity.jpg", sep=""), height=720,
         width=length(data$population)*(2000/length(data$population)))
    plot(data$lineages, col="orange", type='l', lty=2, ylim=c(0,ymax),
         xlab="Time", ylab="Diversity")
    lines(data$freespace, col="cyan", type='l', lty=2)
    lines(data$alpha, col="blue", type='l')
    lines(data$beta, col="green", type='l')
    lines(data$gamma, col="red", type='l')
    legend("topright", c("Lineages (x 0.1)", "Free space per tile", "Alpha diversity",
                         "Beta diversity", "Gamma diversity"),
           col=c("orange", "cyan", "blue", "green", "red"), lwd=2)
    dev.off()
}

visualize = function(toFile=TRUE) {
    plotDiversity()
    plotTraits()
    print("Done.")
}

visualize()
