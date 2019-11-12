#!/usr/bin/Rscript
### This is an analysis script for the island invasion model that creates
### graphs of population and diversity development and invasion success

library(ggplot2)
library(ggfortify)
library(reshape2)

## The working directory may be specified via the commandline, otherwise it
## defaults to results/tests
resultdir = "results"
outdir = paste0(resultdir, "/", commandArgs()[length(commandArgs())])

## Preferred output format
outformat = ".eps" ## default: ".jpg", for publication: ".eps"

## The minimum number of cells an alien species must be in to be considered
## invasive (default: 2)
invasiveThreshold = 2

## D-Day and apocalypse
invasionstart = 500
worldend = 1500

### AUXILIARY FUNCTIONS

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

## Extract the information for one lineage and save it separately
extractLineage = function(rundir, species, timestep=worldend)
{
    statFile = grep("stats_", list.files(rundir), value=T)
    stats = read.table(paste0(rundir, "/", statFile[1]), sep="\t", header=T)
    if (timestep > 0) stats = subset(stats, time==timestep)
    lineageData = subset(stats, lineage==species)
    write.csv(lineageData, paste0(rundir, "/", strsplit(rundir,"/")[[1]][2], "_", species, ".csv"))
}

## Convert a lineage string into a (hopefully unique) colour hex
lineageColour = function(lineage) {
    ascii = Reduce(paste0, utf8ToInt(lineage))
    r = substring(ascii, 1, round(nchar(ascii)/3))
    g = substring(ascii, round(nchar(ascii)/3)+1, round(nchar(ascii)/3)*2)
    b = substring(ascii, round(nchar(ascii)/3)*2+1, nchar(ascii)+1)
    rgbVal = as.numeric(sapply(c(r,g,b), function(i) paste0("0.", i)))
    return(rgb(rgbVal[1], rgbVal[2], rgbVal[3]))
}

### ANALYSE THE WHOLE EXPERIMENT

## WARNING: Takes a *long* time to run!
analyseEstablishment = function(timestep=worldend) {
    print("Analysing invasion success")
    ## figure out how many replicates there are
    getRep = function(d) {paste0(substr(strsplit(d, "_")[[1]][1],9,9), strsplit(d, "_")[[1]][2])}
    reps = unique(sapply(list.files(resultdir), getRep))
    nrep = length(reps)
    ## create the results table
    results = array(dim=c(2,2,2,nrep+1,3), dimnames=list(temperature=c("T15", "T35"),
                                                         disturbance=c("1DB", "10DB"),
                                                         propagules=c("1PP", "10PP"),
                                                         replicates=c(reps, "avg", "sum"),
                                                         diversity=c("natives", "aliens", "invasives")))
    for (d in list.files(resultdir)) {
        dir = paste0(resultdir, "/", d)
        if (file.info(dir)$isdir && !grepl("control",d)) {
            print(d)
            ## figure out the scenario
            specs = subset(collateSpeciesTable(dir), t==timestep)
            if (is.null(specs)) next
            repl = getRep(d)
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
                if (length(subset(specs, lineage==a)$lineage) >= invasiveThreshold) {
                    print(paste("Species", a, "is invasive in", d))
                    extractLineage(dir, a)
                    invasives = invasives+1
                }
            }
            nnative = sum(subset(specs, alien==FALSE)$abundance) # needed for ratio analysis
            nalien = sum(subset(specs, alien==TRUE)$abundance)   # i.e. currently unnecessary
            ##FIXME Somehow, NAs end up in this data set -> corrupt/missing data?
            results[temp, dist, prop, repl,] = c(natives, aliens, invasives)
        }
    }
    for (s in dimnames(results)$diversity) {
        ##take the averages
        results["T35","1DB","1PP","avg",s] = mean(results["T35","1DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T35","1DB","10PP","avg",s] = mean(results["T35","1DB","10PP",1:nrep,s], na.rm=TRUE)
        results["T35","10DB","1PP","avg",s] = mean(results["T35","10DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T35","10DB","10PP","avg",s] = mean(results["T35","10DB","10PP",1:nrep,s], na.rm=TRUE)
        results["T15","1DB","1PP","avg",s] = mean(results["T15","1DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T15","1DB","10PP","avg",s] = mean(results["T15","1DB","10PP",1:nrep,s], na.rm=TRUE)
        results["T15","10DB","1PP","avg",s] = mean(results["T15","10DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T15","10DB","10PP","avg",s] = mean(results["T15","10DB","10PP",1:nrep,s], na.rm=TRUE)
        ## take the sums
        results["T35","1DB","1PP","sum",s] = sum(results["T35","1DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T35","1DB","10PP","sum",s] = sum(results["T35","1DB","10PP",1:nrep,s], na.rm=TRUE)
        results["T35","10DB","1PP","sum",s] = sum(results["T35","10DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T35","10DB","10PP","sum",s] = sum(results["T35","10DB","10PP",1:nrep,s], na.rm=TRUE)
        results["T15","1DB","1PP","sum",s] = sum(results["T15","1DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T15","1DB","10PP","sum",s] = sum(results["T15","1DB","10PP",1:nrep,s], na.rm=TRUE)
        results["T15","10DB","1PP","sum",s] = sum(results["T15","10DB","1PP",1:nrep,s], na.rm=TRUE)
        results["T15","10DB","10PP","sum",s] = sum(results["T15","10DB","10PP",1:nrep,s], na.rm=TRUE)
    }
    return(results)
}

plotEstablishment = function(results) {
    print("Plotting establishment matrices")
    for (d in dimnames(results)$diversity) {
        if (outformat == ".eps") {
            setEPS()
            postscript(paste0(d,".eps"), height=4, width=8)
        } else { jpeg(paste0(d,".jpg"), height=480, width=960, quality=100) }
        mppl = results[,,"1PP","avg",d]
        mpph = results[,,"10PP","avg",d]
        par(mfrow=c(1,2), cex=1.2)
        cols = rev(grey.colors(4))
        image(c(15,35),c(1,10),mppl, col=cols,xaxp=c(15,35,1),yaxp=c(1,10,1),
              xlab="Temperature (°C)",ylab="Disturbance (%)", sub="Propagule pressure: 1")
        text(15,1,round(mppl[1,1],2),col="blue",cex=1.5)
        text(15,10,round(mppl[1,2],2),col="blue",cex=1.5)
        text(35,1,round(mppl[2,1],2),col="blue",cex=1.5)
        text(35,10,round(mppl[2,2],2),col="blue",cex=1.5)
        image(c(15,35),c(1,10),mpph, col=cols,xaxp=c(15,35,1),yaxp=c(1,10,1),
              xlab="Temperature (°C)",ylab="Disturbance (%)",
              sub="Propagule pressure: 10")
        text(15,1,round(mpph[1,1],2),col="blue",cex=1.5)
        text(15,10,round(mpph[1,2],2),col="blue",cex=1.5)
        text(35,1,round(mpph[2,1],2),col="blue",cex=1.5)
        text(35,10,round(mpph[2,2],2),col="blue",cex=1.5)
        dev.off()
    }
}

plotInvasives = function(results) {
    print("Plotting invasion graphs")
    nreps = length(dimnames(results)$replicates)-1
    if (outformat == ".eps") {
        setEPS()
        postscript("propagule_pressure.eps")
    } else { jpeg("propagule_pressure.jpg", height=480, width=480, quality=100) }
    barplot(c(sum(results[,,"1PP",1:nreps,"invasives"]),
              sum(results[,,"10PP",1:nreps,"invasives"])),
            col="lightblue", cex.lab=1.4, cex.axis=1.4, cex.names=1.4,
            names.arg=c("1 individual/time step", "10 individuals/time step"),
            ylab="Number of invasive species")
    dev.off()
    if (outformat == ".eps") {
        setEPS()
        postscript("prod_dist.eps")
    } else { jpeg("prod_dist.jpg", height=480, width=480, quality=100) }
    scenarios = array(c(sum(results["T15","1DB",,1:nreps,"invasives"]),
                        sum(results["T15","10DB",,1:nreps,"invasives"]),
                        sum(results["T35","1DB",,1:nreps,"invasives"]),
                        sum(results["T35","10DB",,1:nreps,"invasives"])),
                      dim=c(2,2))
    image(c(1, 10), c(15,35), scenarios,
          col=rev(grey.colors(4)), yaxp=c(15,35,1), xaxp=c(1,10,1),
          ylab="Temperature (°C)", xlab="Disturbance (%)",
          cex.lab=1.4, cex.axis=1.4)
    text(10,35,scenarios[2,2],col="blue",cex=4)
    text(1,35,scenarios[1,2],col="blue",cex=4)
    text(1,15,scenarios[1,1],col="blue",cex=4)
    text(10,15,scenarios[2,1],col="blue",cex=4)
    dev.off()
}

plotFactors = function(results) {
    print("Plotting invasion factor graph")
    res = melt(results[,,,"sum", "invasives"]) ## XXX needs data recalculation, otherwise use "avg"
    ##res$value = round(res$value*60) #only if "avg" is used
    levels(res$temperature) = c("Temperature: 15°C", "Temperature: 35°C")
    levels(res$disturbance) = c("Disturbance: 1%", "Disturbance: 10%")
    ggplot(res, aes(x=temperature, y=value)) +
        facet_grid(disturbance~temperature, as.table=FALSE) +
        geom_col(aes(x=propagules, y=value, fill=propagules),
                 position="dodge", color="black",
                 show.legend=FALSE) +
        scale_y_continuous(name="Number of invasive species",
                           limits=c(0, max(res$value)+1),
                           breaks=seq(0, max(res$value)+1, 2)) +
        scale_x_discrete(name="Propagule pressure per timestep", labels=c("1 individual", "10 individuals")) +
        theme_classic(base_size=14)
    ggsave(file=paste0("factors", outformat), height=4, width=5, units="in", dpi="print")
}

## Superceded by `plotInvasives` (gives a prettier and more informative graph)
plotFactorsOld = function(results, var="invasives") {
    ## Correlate each experimental factor with invasion success
    ## (`var` specifies whether to consider aliens (sensu lato) or invasives (sensu stricto))
    print("Plotting factor boxplots")
    if (outformat == ".eps") {
        setEPS()
        postscript("factors.eps", height=4, width=12)
    } else { jpeg("factors.jpg", height=400, width=1200, quality=100) }
    par(mfrow=c(1,3),cex=1.3)
    nrep = length(dimnames(results)$replicates) - 1
    templ = as.vector(results["T15",,,1:nrep,var])
    temph = as.vector(results["T35",,,1:nrep,var])
    boxplot(templ, temph, names=c("25°C","35°C"), col="lightblue",
            main="Temperature", ylab=paste("Number of", substr(var,1,nchar(var)-1), "species"))
    distl = as.vector(results[,"1DB",,1:nrep,var])
    disth = as.vector(results[,"10DB",,1:nrep,var])
    boxplot(distl, disth, names=c("1% mortality", "10% mortality"), col="lightblue",
            main="Disturbance")
    propl = as.vector(results[,,"1PP",1:nrep,var])
    proph = as.vector(results[,,"10PP",1:nrep,var])
    boxplot(propl, proph, names=c("1 propagule","10 propagules"), col="lightblue",
            main="Propagule pressure")
    dev.off()
}

## WARNING: Takes a *long* time to run!
analyseFitness = function() {
    ## Compare the fitness values of natives, aliens, and invasives across all
    ## runs and populations
    print("Analysing fitness")
    gauss = function(opt, tol, x) {
        ## private function needed to calculate fitness
        ## (copied from `auxfuncts.jl`)
        a = 1 / (tol * sqrt(2 * pi))
        return(a * exp(-(x-opt)^2/(2*tol^2)))
    }
    ## Build a table containing every population's lineage, status, and fitness value
    pops = data.frame()
    for (d in list.files(resultdir)) {
        dir = paste0(resultdir, "/", d)
        if (!file.info(dir)$isdir) next
        print(dir)
        statFile = grep("stats_", list.files(dir), value=T)
        stats = read.table(paste0(dir, "/", statFile[1]), sep="\t", header=T)
        natives = unique(subset(stats, time==invasionstart)$lineage)
        endpops = subset(subset(stats, time==worldend), adults>0)
        for (i in 1:nrow(endpops)) {
            s = endpops[i,]
            if (s$adults == 0) next
            ##FIXME If individuals have not undergone establishment yet, the adaptation values
            ## are set to a default of 1 -> leading to "inverse" graph shapes
            ## Solution: calculate adaptation explicitly
            ##fitness = s$tempadaptionmed + s$precadaptionmed
            fitness = gauss(s$tempoptmed, s$temptolmed, s$temp) +
                gauss(s$precoptmed, s$prectolmed, s$prec)
            if (s$lineage %in% natives) {
                status = "native"
            } else {
                if (dim(subset(endpops, lineage==s$lineage))[1] < invasiveThreshold) {
                    status = "alien"
                } else {
                    status = "invasive"
                }
            }
            pops = rbind(pops, data.frame(run=d, lineage=as.character(s$lineage), status=status, fitness=fitness))
        }
    }
    colnames(pops) = c("run", "lineage", "status", "fitness")
    return(pops)
}

plotFitness = function(pops) {
    ## Visualise this table as a violin plot
    print("Plotting fitness violin graph")
    ggplot(pops, aes(x=status, y=fitness, fill=status)) +
        coord_cartesian(ylim=c(0,2.1)) +
        geom_violin(trim=FALSE, show.legend=FALSE) +
        geom_boxplot(width=0.05, show.legend=FALSE) +
        scale_x_discrete(limits=c("native", "invasive", "alien")) +
        annotate(geom="text", x=1, y=2.125,
                 label=paste(length(unique(subset(pops, status=="native")$lineage)), "species")) +
        annotate(geom="text", x=2, y=2.125,
                 label=paste(length(unique(subset(pops, status=="invasive")$lineage)), "species")) +
        annotate(geom="text", x=3, y=2.125,
                 label=paste(length(unique(subset(pops, status=="alien")$lineage)), "species")) +
        labs(subtitle=paste("n", "=", dim(pops)[1], "populations,", length(unique(pops$lineage)), "species"),
             y="population median fitness") +
        theme_classic()
    ggsave(file=paste0("fitness", outformat), height=4, width=5, units="in", dpi="print")
}

### ANALYSE AN INDIVIDUAL RUN

## Plot population size and diversity indices over time
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
    if (outformat == ".eps") {
        postscript(paste0(outdir, "/", simname, "_population.eps"),
                   height=6, width=round((maxt*(1300/maxt))/120),
                   horizontal=FALSE, onefile=FALSE, paper="special")
    } else {
        jpeg(paste0(outdir, "/", simname, "_population.jpg"), quality=100,
             height=720, width=maxt*(1300/maxt))
    }
    par(cex=1.6)
    plot(data$population[1:maxt], xlab="Time", ylab="Population size",
         ylim=c(0, max(data$population)), col="red", type='l')
    abline(v=invasionstart,lty=2,col="darkgreen")
    dev.off()
    # Plot diversity development
    print("Plotting diversity...")
    ymax = 3 #max(data$alpha, data$beta, data$gamma, data$freespace)
    if (outformat == ".eps") {
        setEPS()
        postscript(paste0(outdir, "/", simname, "_diversity.eps"),
                   height=6, width=round((maxt*(1300/maxt))/120))
    } else {
        jpeg(paste0(outdir, "/", simname, "_diversity.jpg"), quality=100,
             height=720, width=maxt*(1300/maxt))
    }
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
    print(paste0("Plotting map ", simname, " at timestep ", timestep, "..."))
    if (is.null(data)) return()
    data$X = data$X + 1
    data$Y = data$Y + 1
    data$temp = data$temp - 273
    colours = sapply(levels(data$lineage), lineageColour)
    ##XXX Does not properly plot tiles that are unoccupied
    ## The problem is that unoccupied tiles do not appear in the model output
    ## data at all. Thus, analyse.R has no idea that they exist, or what temperature
    ## they have...
    m = ggplot(data, aes(X, Y)) +
        geom_tile(aes(fill = temp)) +
        labs(x="Longitude", y="Latitude", subtitle=paste("t =", timestep, "time steps")) +
        scale_fill_continuous(low="lightgrey", high="darkgrey") +
        ##FIXME The green rectangle doesn't plot under EPS? -> because of the alpha?
        annotate("rect", xmin=2.5, xmax=3.5, ymin=4.5, ymax=5.5, fill="green", alpha=0.3) +
        geom_jitter(data = data, aes(size = abundance, color = lineage, shape = alien)) +
        scale_colour_manual(values=colours) +
        scale_shape_manual(values=c("TRUE"=17, "FALSE"=16)) +
        guides(colour=FALSE, fill=FALSE,
               shape=guide_legend(title="Non-native"),
               size=guide_legend(title="Abundance"))
    ggsave(file=paste0(resultdir, "/", simname, "/", simname, "_map_t", timestep, outformat),
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

plotTraitPCA = function(outdir) {
    ## load data and extract the interesting bits
    ## TODO add  legend
    print("Plotting PCA...")
    statFile = grep("stats_", list.files(outdir), value=T)
    stats = read.table(paste0(outdir, "/", statFile[1]), sep="\t", header=T)
    traits = grep("med", colnames(stats), value=T)
    cnames = c("lineage", traits)
    st = subset(stats, time==worldend)[cnames]
    ## figure out which species are aliens
    natives = unique(subset(stats, time==invasionstart)[,"lineage"])
    alien = !st$lineage %in% natives
    st[,"alien"] = alien
    ## construct the model
    modeldata = st[traits[!traits %in% names(which(apply(st, 2, var)==0))]] # remove zero-variance columns
    model = prcomp(modeldata, scale=TRUE)
    save(model, file=paste0(outdir, "/", "trait_pca_model.dat"))
    ## plot the PCA graph
    shapes = rep(16, length(st$alien))
    shapes[which(st$alien)] = 17
    colours = sapply(as.character(st$lineage), lineageColour)
    autoplot(model, colour=colours, shape=shapes, loadings=TRUE,
             loadings.colour="darkblue", loadings.label=TRUE, loadings.label.size=3)
    if (outformat == ".eps") {
        outname = paste0(outdir, "/", strsplit(outdir, "/")[[1]][2], "_traits.eps")
    } else {
        outname = paste0(outdir, "/", strsplit(outdir, "/")[[1]][2], "_traits.jpg")
    }
    ggsave(outname, height=6, width=9)
}


### DISPATCH TO THE APPROPRIATE FUNCTIONS
    
visualizeRun = function(outdir, maxt=worldend) {
    plotDiversity(outdir,maxt)
    plotTimeSeries(outdir, round(invasionstart/2))
    plotTraitPCA(outdir)
}
    
analyseAll = function(plotRuns=TRUE,plotAll=TRUE,maxt=worldend,var="invasives") {
    if (plotAll) {
        if (file.exists("invasion_results.dat")) {
            print("Loading analysed invasion results back from file")
            load("invasion_results.dat")
        }
        else {
            results = analyseEstablishment(maxt)
            save(results, file="invasion_results.dat")
        }
        plotEstablishment(results)
        plotInvasives(results)
        plotFactors(results,var)
        if (file.exists("fitness_results.dat")) {
            print("Loading analysed fitness results back from file")
            load("fitness_results.dat")
        }
        else {
            pops = analyseFitness()
            save(pops, file="fitness_results.dat")
        }
        plotFitness(pops)
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
    ##analyseAll(FALSE, TRUE,worldend,"invasives")
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
