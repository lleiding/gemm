#!/usr/bin/Rscript
# This is an analysis script for the island invasion model that creates
# graphs of the total and per species populations over time.

library(ggplot2)

# The working directory may be specified via the commandline, otherwise it
# defaults to results/tests
simname = commandArgs()[length(commandArgs())]
outdir = paste0("results/", simname)


# Plot the distribution of traits in the population at the given timestep (-1 => END)
plotTraits = function(timestep=-1, toFile=TRUE) {
    print("Plotting traits...")
    traitfile = grep(paste0("t", timestep), list.files(outdir), value=T)
    if (length(traitfile) == 0) {
        # If the desired timestep doesn't exist, take the newest timestep we have
        timestep = (length(grep(".tsv", list.files(outdir), value=T))-1) * 10
        traitfile = grep(paste("t", timestep, sep=""), list.files(outdir), value=T)
    }
    traitfilepath = paste(outdir, traitfile, sep="/")
    if (!file.exists(traitfilepath)) {
        print(paste("WARNING: traitfile not found", traitfilepath))
        return()
    }
    ts = read.table(traitfilepath, header=T)
    jpeg(paste0(outdir, "/", simname, "_traits.jpg"), height=720, width=1800)
    boxplot(ts$fitness*100, log(ts$size), log(ts$seedsize), log(ts$repsize), ts$lnkgunits/10,
            ts$ngenes/10, ts$temptol, ts$tempopt/100, ts$prectol, ts$precopt, ts$compat*10,
            ts$repradius*10, ts$dispmean, ts$dispshape,
            names=c("Fitness (x100)", "log(Size)", "log(Seed size)", "log(Reprod. size)",
                    "Chromosome (x0.1)", "Genes (x0.1)", "Temp. tolerance",
                    "Temp. opt. (x0.01)", "Precip. tolerance",
                    "Precip. optimum", "Compatibility (x10)", "Rep. radius (x10)",
                    "Dispersal mean", "Dispersal shape"))
    legend("top", c(paste("Individuals:", length(ts$counter)),
                    paste("Lineages:", length(unique(ts$lineage)))))
    dev.off()
}

# Plot population size and diversity indices over time
plotDiversity = function(logfile="diversity.log") {
    logfile = paste(outdir, logfile, sep="/")
    if (!file.exists(logfile)) {
        print(paste("WARNING: logfile not found", logfile))
        return()
    }
    data = read.csv(logfile)
    data$lineages = data$lineages/10 #otherwise the Y axis is too big
    # Plot population sizes
    print("Plotting population development...")
    jpeg(paste0(outdir, "/", simname, "_population.jpg"), height=720,
         width=length(data$population)*(2000/length(data$population)))
    plot(data$population, xlab="Time", ylab="Population size", ylim=c(0, max(data$population)),
         col="red", type='l')
    dev.off()
    # Plot diversity development
    print("Plotting diversity...")
    ymax = max(data$lineages, data$alpha, data$beta, data$gamma, data$freespace)
    jpeg(paste0(outdir, "/", simname, "_diversity.jpg"), height=720,
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

plotMap = function(timestep=-1, compensate=TRUE) {
    print(paste0("Plotting map at timestep ", timestep, "..."))
    tsvfile = grep(paste0("t", timestep), list.files(outdir), value=T)
    if (length(tsvfile) == 0 && compensate) {
        # If the desired timestep doesn't exist, take the newest timestep we have
        timestep = (length(grep(".tsv", list.files(outdir), value=T))-1) * 10
        tsvfile = grep(paste("t", timestep, sep=""), list.files(outdir), value=T)
    }
    tsvfilepath = paste(outdir, tsvfile, sep="/")
    if (!file.exists(tsvfilepath)) {
        print(paste("WARNING: tsvfile not found", tsvfilepath))
        return()
    }
    ts = read.table(tsvfilepath, header=T)
    jpeg(paste0(outdir, "/", simname, "_map_", timestep, ".jpg"), height=720, width=720)
    m = ggplot(ts, aes(xloc, yloc))
    m + geom <- tile(aes(fill = temp.C, width = 0.95, height = 0.95)) +
        scale <- fill <- continuous(low="white", high="black") +
            scale <- color <- gradientn(colours = rainbow(5)) +
                geom <- jitter(data = allspecies, aes(size = abundance, color = lineage))
    #FIXME What is `allspecies?`
    dev.off()
}

plotTimeSeries = function(step=1) {
    files = grep("_t", list.files(outdir), value=T)
    for (f in files) {
        timestep = as.numeric(substring(strsplit(f, "_")[[1]][2], 2))
        if (timestep %% step == 0)
            plotMap(timestep, FALSE)
    }
}

visualize = function(toFile=TRUE) {
    plotDiversity()
    plotTraits()
    plotTimeSeries(500)
}

# If the simname is given as 'all', process every folder in 'results'
if (simname == "all") {
    for (f in list.files("results")) {
        print(paste("Processing", f))
        simname = f
        outdir = paste0("results/", simname)
        if (file.info(outdir)$isdir) visualize()
    }
    print("Done.")
} else {
    # Otherwise, just look at the specified directory (or the default, if
    # the specified doesn't exist)
    if (!(file.exists(outdir) && file.info(outdir)$isdir)) {
        simname = "tests"
        outdir = paste0("results/", simname)
    }
    visualize()
    print("Done.")
}
