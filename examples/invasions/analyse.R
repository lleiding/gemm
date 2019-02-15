#!/usr/bin/Rscript
### This is an analysis script for the island invasion model that creates
### graphs of population and diversity development and invasion success

library(ggplot2)

## The working directory may be specified via the commandline, otherwise it
## defaults to results/tests
resultdir = "results"
outdir = paste0(resultdir, "/", commandArgs()[length(commandArgs())])

## D-Day and apocalypse
invasionstart = 500
worldend = 1500

## Read in the output file and figure out which lineages are native
collateSpeciesTable = function(rundir)
{
    st = read.csv(paste0(rundir, "/lineages.log"))
    natives = unique(subset(st, t==invasionstart)$lineage)
    alien = rep(c(FALSE), length(subset(st, t <= invasionstart)[,1]))
    alien = c(alien, !subset(st, t > invasionstart)$lineage %in% natives)
    st = cbind(st, alien)
    return(st)
}

### ANALYSE THE WHOLE EXPERIMENT

analyseEstablishment = function(timestep=worldend) {
    print("Analysing invasion success")
    ## figure out how many replicates there are
    nrep = length(grep("control", grep("invasion", list.files(resultdir), value=T), value=T, invert=T)) / 8
    reps = c()
    for (i in 1:nrep) reps = c(reps, paste0("r", i))
    ## create the results table
    results = array(dim=c(2,2,2,nrep+1,3), dimnames=list(temperature=c("T15", "T35"),
                                                         disturbance=c("1DB", "10DB"),
                                                         propagules=c("1PP", "10PP"),
                                                         replicates=c(reps, "avg"),
                                                         diversity=c("natives", "aliens", "invasives")))
    for (d in list.files(resultdir)) {
        dir = paste0(resultdir, "/", d)
        if (file.info(dir)$isdir && !grepl("control",d)) {
            print(d)
            ## figure out the scenario
            specs = subset(collateSpeciesTable(dir), t==timestep)
            if (is.null(specs)) next
            repl = strsplit(d, "_")[[1]][2]
            dist = strsplit(d, "_")[[1]][5]
            prop = strsplit(d, "_")[[1]][4]
            if (grepl("cold", d)) temp = "T15"
            else if (grepl("default", d)) temp = "T25" #currently not used, will give an error
            else if (grepl("hot", d)) temp = "T35"
            else temp = "T_error" # The terror of an unknown value...
            ## calculate diversity
            natives = length(unique(subset(specs, alien==FALSE)$lineage))
            aliens = length(unique(subset(specs, alien==TRUE)$lineage))
            invasives = 0
            for (a in unique(subset(specs, alien==TRUE)$lineage)) {
                if (length(subset(specs, lineage==a)$lineage) > 6) invasives = invasives+1
            }
            nnative = sum(subset(specs, alien==FALSE)$abundance)
            nalien = sum(subset(specs, alien==TRUE)$abundance)
            results[temp, dist, prop, repl,] = c(natives, aliens, invasives)
        }
    }
    ##take the averages
    for (s in dimnames(results)$diversity) {
        results["T35","1DB","1PP","avg",s] = mean(results["T35","1DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T35","1DB","10PP","avg",s] = mean(results["T35","1DB","10PP",1:nrep,s], na.rm=TRUE)
        results["T35","10DB","1PP","avg",s] = mean(results["T35","10DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T35","10DB","10PP","avg",s] = mean(results["T35","10DB","10PP",1:nrep,s], na.rm=TRUE)
        results["T15","1DB","1PP","avg",s] = mean(results["T15","1DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T15","1DB","10PP","avg",s] = mean(results["T15","1DB","10PP",1:nrep,s], na.rm=TRUE)
        results["T15","10DB","1PP","avg",s] = mean(results["T15","10DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T15","10DB","10PP","avg",s] = mean(results["T15","10DB","10PP",1:nrep,s], na.rm=TRUE)
    }
    return(results)
}

plotEstablishment = function(results) {
    print("Plotting establishment matrices")
    for (d in dimnames(results)$diversity) {
        jpeg(paste0(d,".jpg"), height=480, width=960, quality=100)
        mppl = results[,,"1PP","avg",d]
        mpph = results[,,"10PP","avg",d]
        par(mfrow=c(1,2), cex=1.2)
        maxv = max(mppl,mpph,na.rm=TRUE)
        ##XXX Does the following still give problems when 0 < maxv < 1?
        if (maxv <= 0) cols = grey(c(0))
        else if (maxv >= 1) cols = grey(c(1))
        else cols = grey(rev((0:maxv)/maxv))
        image(c(15,35),c(1,10),mppl, col=cols,xaxp=c(15,35,1),yaxp=c(1,10,1),
              xlab="Temperature (°C)",ylab="Disturbance (%)", main="Propagule pressure: 1")
        text(15,1,round(mppl[1,1],2),col="blue",cex=3)
        text(15,10,round(mppl[1,2],2),col="blue",cex=3)
        text(35,1,round(mppl[2,1],2),col="blue",cex=3)
        text(35,10,round(mppl[2,2],2),col="blue",cex=3)
        image(c(15,35),c(1,10),mpph, col=cols,xaxp=c(15,35,1),yaxp=c(1,10,1),
              xlab="Temperature (°C)",ylab="Disturbance (%)",
              main="Propagule pressure: 10")
        text(15,1,round(mpph[1,1],2),col="blue",cex=3)
        text(15,10,round(mpph[1,2],2),col="blue",cex=3)
        text(35,1,round(mpph[2,1],2),col="blue",cex=3)
        text(35,10,round(mpph[2,2],2),col="blue",cex=3)
        dev.off()
    }
}

plotFactors = function(results, var="invasives") {
    ## Correlate each experimental factor with invasion success
    ## (`var` specifies whether to consider aliens (sensu lato) or invasives (sensu stricto))
    print("Plotting factor boxplots")
    jpeg("factors.jpg",height=400, width=1200, quality=100)
    par(mfrow=c(1,3),cex=1.3)
    nrep = length(dimnames(results)$replicates) - 1
    templ = as.vector(results["T15",,,1:nrep,var])
    temph = as.vector(results["T35",,,1:nrep,var])
    boxplot(templ, temph, names=c("25°C","35°C"), col="lightblue",
            main="Temperature", ylab=paste("Number of", substr(var,1,nchar(var)-1), "species"))
    distl = as.vector(results[,"1DB",,1:nrep,"aliens"])
    disth = as.vector(results[,"10DB",,1:nrep,"aliens"])
    boxplot(distl, disth, names=c("1% mortality", "10% mortality"), col="lightblue",
            main="Disturbance")
    propl = as.vector(results[,,"1PP",1:nrep,"aliens"])
    proph = as.vector(results[,,"10PP",1:nrep,"aliens"])
    boxplot(propl, proph, names=c("1 propagule","10 propagules"), col="lightblue",
            main="Propagule pressure")
    dev.off()
}

### ANALYSE AN INDIVIDUAL RUN

# Plot population size and diversity indices over time
plotDiversity = function(outdir, maxt=3000, logfile="diversity.log") {
    simname = strsplit(outdir, "/")[[1]][2]
    logfile = paste(outdir, logfile, sep="/")
    if (!file.exists(logfile) || file.info(logfile)$size == 0) {
        print(paste("WARNING: logfile not found", logfile))
        return()
    }
    data = read.csv(logfile)
    if (maxt < 0 || maxt > length(data$population)) maxt = length(data$population)
    data$lineages = data$lineages/10 #otherwise the Y axis is too big
    # Plot population sizes
    print("Plotting population development...")
    jpeg(paste0(outdir, "/", simname, "_population.jpg"), quality=100, height=720,
         width=maxt*(1300/maxt))
    par(cex=1.6)
    plot(data$population[1:maxt], xlab="Time", ylab="Population size",
         ylim=c(0, max(data$population)), col="red", type='l')
    abline(v=invasionstart,lty=2,col="darkgreen")
    dev.off()
    # Plot diversity development
    print("Plotting diversity...")
    ymax = 3 #max(data$alpha, data$beta, data$gamma, data$freespace)
    jpeg(paste0(outdir, "/", simname, "_diversity.jpg"), quality=100, height=720,
         width=maxt*(1300/maxt))
    par(cex=1.6)
    plot(data$lineages[1:maxt], col="orange", type='l', lty=2, ylim=c(0,ymax),
         xlab="Time", ylab="Diversity")
    lines(data$freespace[1:maxt], col="cyan", type='l', lty=2)
    lines(data$alpha[1:maxt], col="blue", type='l')
    lines(data$beta[1:maxt], col="green", type='l')
    lines(data$gamma[1:maxt], col="red", type='l')
    abline(v=invasionstart,lty=2,col="darkgreen")
    legend("topright", c("Lineages (x 0.1)", "Free space per tile", "Alpha diversity",
                         "Beta diversity", "Gamma diversity"),
           col=c("orange", "cyan", "blue", "green", "red"), lwd=2)
    dev.off()
}

plotMap = function(data, timestep, simname) {
    ##FIXME Make sure `alien` status is displayed the same each time, adjust legend, create custom mapping
    print(paste0("Plotting map ", simname, " at timestep ", timestep, "..."))
    if (is.null(data)) return()
    data$X = data$X + 1
    data$Y = data$Y + 1
    data$temp = data$temp - 273
    m = ggplot(data, aes(X, Y))
    ##FIXME Does not properly plot tiles that are unoccupied
    ## The problem is that unoccupied tiles do not appear in the model output
    ## data at all. Thus, analyse.R has no idea that they exist, or what temperature
    ## they have...
    m + geom_tile(aes(fill = temp)) + labs(x="Longitude", y="Latitude") +
        scale_fill_continuous(low="lightgrey", high="darkgrey") +
        annotate("rect", xmin=2.5, xmax=3.5, ymin=4.5, ymax=5.5, fill="green", alpha=0.3) +
        annotate("text", x=3, y=5.5, label=paste("t =", timestep)) +
        geom_jitter(data = data, aes(size = abundance, color = lineage, shape = alien)) +
        guides(colour=FALSE, shape=FALSE)
    ggsave(file=paste0(resultdir, "/", simname, "/", simname, "_map_t", timestep, ".jpg"),
           height=4, width=5, units="in", dpi="print")
}

plotTimeSeries = function(rundir, step) {
    data = collateSpeciesTable(rundir)
    ## figure out at which timepoints to create a map
    snapshots = seq(0, worldend, step)
    snapshots[1] = 1
    if (!invasionstart %in% snapshots) snapshots = c(snapshots, invasionstart)
    if (!worldend %in% snapshots) snapshots = c(snapshots, worldend)
    ## plot the maps
    for (s in snapshots) {
        plotMap(subset(data, t == s), s, strsplit(rundir, "/")[[1]][2])
    }
}


### DISPATCH TO THE APPROPRIATE FUNCTIONS
    
visualizeRun = function(outdir, maxt=worldend) {
    plotDiversity(outdir,maxt)
    plotTimeSeries(outdir, round(invasionstart/2))
}
    
analyseAll = function(plotRuns=TRUE,plotAll=TRUE,maxt=worldend,var="invasives") {
    if (plotAll) {
        if (file.exists("experiment_results.dat")) {
            print("Loading analysed results back from file")
            load("experiment_results.dat")
        }
        else {
            results = analyseEstablishment(maxt)
            save(results, file="experiment_results.dat")
        }
        plotEstablishment(results)
        plotFactors(results,var)
    }
    if (plotRuns) {
        for (f in list.files(resultdir)) {
            print(paste("Processing", f))
            simname = f
            outdir = paste0(resultdir, "/", simname)
            if (file.info(outdir)$isdir) visualizeRun(outdir)
        }
    }
}


### CALL THE APPROPRIATE FUNCTIONS
    
# If the simname is given as 'all', do a whole-experiment analysis
if (commandArgs()[length(commandArgs())] == "all") {
    outdir = sub("/all", "", outdir)
    analyseAll(TRUE,TRUE,worldend,"invasives")
    print("Done.")
} else {
    # Otherwise, just look at the specified directory (or the default, if
    # the specified doesn't exist)
    if (!(file.exists(outdir) && file.info(outdir)$isdir)) {
        simname = "tests"
        outdir = paste0(resultdir, "/", simname)
    }
    visualizeRun(outdir)
    print("Done.")
}
