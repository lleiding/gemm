#!/usr/bin/Rscript
### This is an analysis script for the island invasion model that creates
### graphs of the total and per species populations over time.

library(ggplot2)

## The working directory may be specified via the commandline, otherwise it
## defaults to results/tests
resultdir = "results"
simname = commandArgs()[length(commandArgs())]
outdir = paste0(resultdir, "/", simname)


### ANALYSE THE WHOLE EXPERIMENT

collateSpeciesTable = function(rundir, timestep, compensate=TRUE, showinvaders=TRUE, includeTS=FALSE) {
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
    allspecies$xloc = as.numeric(allspecies$xloc)
    allspecies$yloc = as.numeric(allspecies$yloc)
    allspecies$alien = as.logical(allspecies$alien)
    write.csv2(allspecies, file=paste0(rundir,"/species_",sub(paste0(resultdir,"/"), "", rundir)
                                      ,"_t",timestep,".csv"))
    if (includeTS) {
        return(list("allspecies"=allspecies, "ts"=ts))
    } else {
        return(allspecies)
    }
}

analyseEstablishment = function() {
    print("Analysing invasion success")
    ## create the results table
    results = array(dim=c(2,2,2,6,4), dimnames=list(temperature=c("T25", "T35"),
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
            if (nnative == 0) nnative = 1 #XXX is this acceptable?
            nalien = sum(subset(specs, alien==TRUE)$abundance)
            ratio = nalien / nnative # ratio of total abundances
            results[temp, dist, prop, repl,] = c(natives, aliens, invasives, ratio)
        }
    }
    ##take the averages
    for (s in dimnames(results)$diversity) {
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

plotEstablishment = function(results) {
    print("Plotting establishment matrices")
    for (d in dimnames(results)$diversity) {
        jpeg(paste0(d,".jpg"), width=960)
        mppl = t(results[,,"1PP","avg",d])
        mpph = t(results[,,"10PP","avg",d])
        par(mfrow=c(1,2), cex=1.2)
        maxv = max(mppl,mpph,na.rm=TRUE)
        cols = grey(rev((1:maxv)/maxv))
        image(c(25,35),c(1,10),mppl, col=cols,xaxp=c(25,35,1),yaxp=c(1,10,1),
              xlab="Temperature (°C)",ylab="Disturbance (%)", main="Propagule pressure: 1")
        text(25,1,round(mppl[1,1],2),col="blue",cex=3)
        text(25,10,round(mppl[1,2],2),col="blue",cex=3)
        text(35,1,round(mppl[2,1],2),col="blue",cex=3)
        text(35,10,round(mppl[2,2],2),col="blue",cex=3)
        image(c(25,35),c(1,10),mpph, col=cols,xaxp=c(25,35,1),yaxp=c(1,10,1),
              xlab="Temperature (°C)",ylab="Disturbance (%)",
              main="Propagule pressure: 10")
        text(25,1,round(mpph[1,1],2),col="blue",cex=3)
        text(25,10,round(mpph[1,2],2),col="blue",cex=3)
        text(35,1,round(mpph[2,1],2),col="blue",cex=3)
        text(35,10,round(mpph[2,2],2),col="blue",cex=3)
        dev.off()
    }
}

plotFactors = function(results) {
    print("Plotting factor boxplots")
    jpeg("factors.jpg",height=400, width=1200)
    par(mfrow=c(1,3),cex=1.3)
    templ = as.vector(results["T25",,,1:5,"aliens"])
    temph = as.vector(results["T35",,,1:5,"aliens"])
    boxplot(templ, temph, names=c("25°C","35°C"), col="lightblue",
            main="Temperature", ylab="Number of alien species")
    if (t.test(templ,temph)$p.value < 0.05)
        text(1.5, max(temph,templ,na.rm=TRUE), "*", cex=4)
    distl = as.vector(results[,"1DB",,1:5,"aliens"])
    disth = as.vector(results[,"10DB",,1:5,"aliens"])
    boxplot(distl, disth, names=c("1% mortality", "10% mortality"), col="lightblue",
            main="Disturbance")
    if (t.test(distl,disth)$p.value < 0.05)
        text(1.5, max(disth, distl,na.rm=TRUE), "*", cex=4)
    propl = as.vector(results[,,"1PP",1:5,"aliens"])
    proph = as.vector(results[,,"10PP",1:5,"aliens"])
    boxplot(propl, proph, names=c("1 propagule","10 propagules"), col="lightblue",
            main="Propagule pressure")
    if (t.test(propl,proph)$p.value < 0.05)
        text(1.5, max(proph,propl,na.rm=TRUE), "*", cex=4)
    dev.off()
}

factorAnalysis = function(results) {
    ##TODO ANOVA
    
}

analyseAll = function(plotAll=TRUE,plotRuns=FALSE) {
    if (plotAll) {
        results = analyseEstablishment()
        save(results, file="experiment_results.dat")
        plotEstablishment(results)
        plotFactors(results)
        factorAnalysis(results)
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

### ANALYSE AN INDIVIDUAL RUN

## Plot the distribution of traits in the population at the given timestep (-1 => END)
plotTraits = function(outdir, timestep=-1, toFile=TRUE) {
    print(paste("Plotting traits at timestep", timestep))
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
            ts$dispmean, ts$dispshape,
            names=c("Fitness (x100)", "log(Size)", "log(Seed size)", "log(Reprod. size)",
                    "Chromosome (x0.1)", "Genes (x0.1)", "Temp. tolerance",
                    "Temp. opt. (-293)", "Precip. tolerance",
                    "Precip. optimum", "Compatibility (x10)",
                    "Dispersal mean", "Dispersal shape"))
    legend("top", c(paste("Individuals:", length(ts$counter)),
                    paste("Lineages:", length(unique(ts$lineage)))))
    dev.off()
}

# Plot population size and diversity indices over time
plotDiversity = function(outdir, logfile="diversity.log") {
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

plotMap = function(outdir, timestep=-1, compensate=TRUE, showinvaders=TRUE) {
    ##FIXME Make sure `alien` status is displayed the same each time, adjust legend, create custom mapping
    print(paste0("Plotting map at timestep ", timestep, "..."))
    st = collateSpeciesTable(outdir, timestep, compensate, showinvaders,TRUE)
    ts = st$ts
    allspecies = st$allspecies
    ts$xloc = ts$xloc + 1
    ts$yloc = ts$yloc + 1
    if (is.null(allspecies)) return()
    m = ggplot(ts, aes(xloc, yloc))
    m + geom_tile(aes(fill = temp.C)) + labs(x="Longitude", y="Latitude") +
        scale_fill_continuous(low="lightgrey", high="darkgrey") +
        annotate("rect", xmin=2.5, xmax=3.5, ymin=4.5, ymax=5.5, fill="green", alpha=0.3) +
        #scale_color_gradient(colours = rainbow(5)) +
        #scale_color_discrete(allspecies$lineage,name="Lineages",labels=as.character(allspecies$lineage)) +
        geom_jitter(data = allspecies, aes(size = abundance, color = lineage, shape = alien)) +
        guides(colour=FALSE)
    ggsave(file=paste0(outdir, "/", simname, "_map_t", timestep, ".jpg"), height=8, width=10)
}

plotTimeSeries = function(outdir, step=1) {
    files = grep(".tsv", grep("_t", list.files(outdir), value=T), value=T)
    for (f in files) {
        timestep = as.numeric(substring(strsplit(f, "_")[[1]][2], 2))
        if (timestep %% step == 0 || abs(timestep) == 1)
            plotMap(outdir, timestep, FALSE)
    }
}

visualizeRun = function(outdir) {
    plotDiversity(outdir)
    plotTraits(outdir)
    plotTraits(outdir, 1000)
    plotTimeSeries(outdir, 250)
}


### CALL THE APPROPRIATE FUNCTIONS
    
# If the simname is given as 'all', do a whole-experiment analysis
if (simname == "all") {
    outdir = sub("/all", "", outdir)
    analyseAll(FALSE,TRUE)
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
