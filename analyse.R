#!/usr/bin/Rscript
# This is an analysis script for the island invasion model that creates
# graphs of the total and per species populations over time.

# The working directory may be specified via the commandline, otherwise it
# defaults to results/tests
outdir = paste("results/", commandArgs()[length(commandArgs())], sep="")
if (!(file.exists(outdir) && file.info(outdir)$isdir)) {
    outdir = "results/tests"
}


# Plot the population sizes and diversity indices over time
visualize = function(logfile="diversity.log", toFile=TRUE) {
    logfile = paste(outdir, logfile, sep="/")
    data = read.csv(logfile)
    data$lineages = data$lineages/10 #otherwise the Y axis is too big
    # Plot population sizes
    if (toFile) {
        jpeg(paste(outdir, "population.jpg", sep="/"), height=720, width=length(data$population)*20)
    }
    plot(data$population, xlab="Time", ylab="Population size", ylim=c(0, max(data$population)),
         col="red", type='l')
    if (toFile) dev.off()
    # Plot diversity development
    ymax = max(data$lineages, data$alpha, data$beta, data$gamma)
    if (toFile) {
        jpeg(paste(outdir, "diversity.jpg", sep="/"), height=720,
             width=length(data$population)*20)
    }
    plot(data$lineages, col="orange", type='l', lty=2, ylim=c(0,ymax*1.5),
         xlab="Time", ylab="Diversity")
    lines(data$freespace, col="cyan", type='l', lty=2)
    lines(data$alpha, col="blue", type='l')
    lines(data$beta, col="green", type='l')
    lines(data$gamma, col="red", type='l')
    legend("topright", c("Lineages (x 0.1)", "Free space per tile", "Alpha diversity",
                         "Beta diversity", "Gamma diversity"),
           col=c("orange", "cyan", "blue", "green", "red"), lwd=2)
    if (toFile) dev.off()
}

visualize()
