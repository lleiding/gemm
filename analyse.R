#!/usr/bin/Rscript
# This is an analysis script for the island invasion model that creates
# graphs of the total and per species populations over time.

library(ggplot2)

# The working directory may be specified via the commandline, otherwise it
# defaults to results/tests
resultdir = "results"
simname = commandArgs()[length(commandArgs())]
outdir = paste0(resultdir, "/", simname)


### ANALYSE THE WHOLE EXPERIMENT

collateSpeciesTable = function(rundir, timestep, compensate=TRUE, showinvaders=TRUE) {
    ## load raw data
    tsvfile = grep(".tsv", grep(paste0("t", timestep, "_"), list.files(rundir), value=T), value=T)
    if (length(tsvfile) == 0 && compensate) {
        # If the desired timestep doesn't exist, take the newest timestep we have
        timestep = (length(grep(".tsv", list.files(rundir), value=T))-1) * 10
        tsvfile = grep(paste("t", timestep, sep=""), list.files(rundir), value=T)
    }
    tsvfilepath = paste(rundir, tsvfile, sep="/")
    if (!file.exists(tsvfilepath) || file.info(tsvfilepath)$isdir || file.info(tsvfilepath)$size == 0) {
        print(paste("WARNING: tsvfile not found", tsvfilepath))
        return()
    }
    print(paste("Collating data from", tsvfilepath)) #DEBUG
    ts = read.table(tsvfilepath, header=T)
    ts$temp.C = ts$temp - 273
    ## load comparison data (to show invasive species)
    if (showinvaders && timestep < 1000 && timestep != -1) showinvaders = FALSE
    if (showinvaders) {
        tsvfile2 = grep(".tsv", grep("t1000_", list.files(rundir), value=T), value=T)
        tsvfilepath2 = paste(rundir, tsvfile2, sep="/")
        ts2 = read.table(tsvfilepath2, header=T)
        nativespecs = unique(ts2$lineage)
    }
    ## create an index of species abundance
    allspecies = c()
    for (p in unique(ts$id)) {
        patch = subset(ts, id == p)
        x = patch$xloc[1]
        y = patch$yloc[1]
        for (s in unique(patch$lineage)) {
            n = length(subset(patch, lineage == s)$lineage)
            if (showinvaders && !(s %in% nativespecs)) i = TRUE
            else i = FALSE
            allspecies = rbind(allspecies, c(x,y,s,n,i))
        }
    }
    colnames(allspecies) = c("xloc", "yloc", "lineage", "abundance", "alien")
    allspecies = as.data.frame(allspecies)
    allspecies$abundance = as.numeric(as.character(allspecies$abundance))
    return(allspecies)
}

analyseEstablishment = function() {
    print("Analysing invasion success")
    ## create the results table
    results = array(dim=c(2,2,2,6,4), dimnames=list(temperature=c("T35", "T25"),
                                                    disturbance=c("1DB", "10DB"),
                                                    propagules=c("1PP", "10PP"),
                                                    replicates=c("r1", "r2", "r3", "r4",
                                                                 "r5", "avg"),
                                                    diversity=c("natives", "aliens",
                                                                "invasives", "ratio")))
    for (d in list.files(resultdir)) {
        dir = paste0(resultdir, "/", d)
        if (file.info(dir)$isdir && !grepl("control",d)) {
            ## figure out the scenario
            specs = collateSpeciesTable(dir, -1, FALSE, TRUE)
            if (is.null(specs)) next
            repl = strsplit(d, "_")[[1]][2]
            dist = strsplit(d, "_")[[1]][5]
            prop = strsplit(d, "_")[[1]][4]
            if (grepl("default", d)) temp = "T25"
            else temp = "T35"
            ## calculate diversity
            natives = length(unique(subset(specs, alien==FALSE)$lineage))
            aliens = length(unique(subset(specs, alien==TRUE)$lineage))
            invasives = 0
            for (a in unique(subset(specs, alien==TRUE)$lineage)) {
                if (length(subset(specs, lineage==a)$lineage) > 6) invasives = invasives+1
            }
            nnative = sum(subset(specs, alien==FALSE)$abundance)
            nalien = sum(subset(specs, alien==TRUE)$abundance)
            ratio = nalien / nnative # ratio of total abundances
            results[temp, dist, prop, repl,] = c(natives, aliens, invasives, ratio)
        }
    }
    ##take the averages
    for (s in c("natives", "aliens", "invasives", "ratio")) {
        results["T35","1DB","1PP","avg",s] = mean(results["T35","1DB","1PP",1:5,s], na.rm=TRUE)
        results["T35","1DB","10PP","avg",s] = mean(results["T35","1DB","10PP",1:5,s], na.rm=TRUE)
        results["T35","10DB","1PP","avg",s] = mean(results["T35","10DB","1PP",1:5,s], na.rm=TRUE)
        results["T35","10DB","10PP","avg",s] = mean(results["T35","10DB","10PP",1:5,s], na.rm=TRUE)
        results["T25","1DB","1PP","avg",s] = mean(results["T25","1DB","1PP",1:5,s], na.rm=TRUE)
        results["T25","1DB","10PP","avg",s] = mean(results["T25","1DB","10PP",1:5,s], na.rm=TRUE)
        results["T25","10DB","1PP","avg",s] = mean(results["T25","10DB","1PP",1:5,s], na.rm=TRUE)
        results["T25","10DB","10PP","avg",s] = mean(results["T25","10DB","10PP",1:5,s], na.rm=TRUE)
    }
    return(results)
}

plotEstablishment = function() {
    results = analyseEstablishment()
    ##TODO save to file
    ##TODO create graphics - one per category in 'diversity'
}
        

### ANALYSE AN INDIVIDUAL RUN

## Plot the distribution of traits in the population at the given timestep (-1 => END)
plotTraits = function(timestep=-1, toFile=TRUE) {
    print("Plotting traits...")
    traitfile = grep(paste0("t", timestep, "_"), list.files(outdir), value=T)
    if (length(traitfile) == 0) {
        # If the desired timestep doesn't exist, take the newest timestep we have
        timestep = (length(grep(".tsv", list.files(outdir), value=T))-1) * 10
        traitfile = grep(paste("t", timestep, sep=""), list.files(outdir), value=T)
    }
    traitfilepath = paste(outdir, traitfile, sep="/")
    if (!file.exists(traitfilepath) || file.info(traitfilepath)$size == 0) {
        print(paste("WARNING: traitfile not found", traitfilepath))
        return()
    }
    ts = read.table(traitfilepath, header=T)
    jpeg(paste0(outdir, "/", simname, "_traits_t", timestep, ".jpg"), height=720, width=1800)
    boxplot(ts$fitness*100, log(ts$size), log(ts$seedsize), log(ts$repsize), ts$lnkgunits/10,
            ts$ngenes/10, ts$temptol, ts$tempopt-293, ts$prectol, ts$precopt, ts$compat*10,
            ts$repradius*10, ts$dispmean, ts$dispshape,
            names=c("Fitness (x100)", "log(Size)", "log(Seed size)", "log(Reprod. size)",
                    "Chromosome (x0.1)", "Genes (x0.1)", "Temp. tolerance",
                    "Temp. opt. (-293)", "Precip. tolerance",
                    "Precip. optimum", "Compatibility (x10)", "Rep. radius (x10)",
                    "Dispersal mean", "Dispersal shape"))
    legend("top", c(paste("Individuals:", length(ts$counter)),
                    paste("Lineages:", length(unique(ts$lineage)))))
    dev.off()
}

# Plot population size and diversity indices over time
plotDiversity = function(logfile="diversity.log") {
    logfile = paste(outdir, logfile, sep="/")
    if (!file.exists(logfile) || file.info(logfile)$size == 0) {
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

plotMap = function(timestep=-1, compensate=TRUE, showinvaders=TRUE) {
    #FIXME Make sure `alien` status is displayed the same each time, adjust legend, create custom mapping
    print(paste0("Plotting map at timestep ", timestep, "..."))
    allspecies = collateSpeciesTable(outdir, timestep, compensate, showinvaders)
    m = ggplot(ts, aes(xloc, yloc))
    m + geom_tile(aes(fill = temp.C)) + labs(x="Longitude", y="Latitude") +
        scale_fill_continuous(low="lightgrey", high="darkgrey") +
        scale_color_gradientn(colours = rainbow(5)) +
        geom_jitter(data = allspecies, aes(size = abundance, color = lineage, shape = alien))
    ggsave(file=paste0(outdir, "/", simname, "_map_t", timestep, ".jpg"), height=8, width=10)
}

plotTimeSeries = function(step=1) {
    files = grep(".tsv", grep("_t", list.files(outdir), value=T), value=T)
    for (f in files) {
        timestep = as.numeric(substring(strsplit(f, "_")[[1]][2], 2))
        if (timestep %% step == 0 || abs(timestep) == 1)
            plotMap(timestep, FALSE)
    }
}

visualizeRun = function() {
    plotDiversity()
    plotTraits()
    plotTraits(1000)
    plotTimeSeries(500)
}


### CALL THE APPROPRIATE FUNCTIONS
    
# If the simname is given as 'all', process every folder in 'results'
if (simname == "all") {
    plotEstablishment()
    for (f in list.files(resultdir)) {
        print(paste("Processing", f))
        simname = f
        outdir = paste0(resultdir, "/", simname)
        if (file.info(outdir)$isdir) visualizeRun()
    }
    print("Done.")
} else {
    # Otherwise, just look at the specified directory (or the default, if
    # the specified doesn't exist)
    if (!(file.exists(outdir) && file.info(outdir)$isdir)) {
        simname = "tests"
        outdir = paste0(resultdir, "/", simname)
    }
    visualizeRun()
    print("Done.")
}
