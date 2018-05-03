#!/usr/bin/Rscript
# This is an analysis script for the island invasion model that creates
# graphs of the total and per species populations over time.

# The working directory may be specified via the commandline, otherwise it
# defaults to results/tests
outdir = paste("results/", commandArgs()[length(commandArgs())], sep="")
if !(file.exists(outdir) && file.info(outdir)$isdir) {
    outdir = "results/tests"
}

readDataSeries = function() {
    # read in all model output
    data = list()
    files = grep(".tsv", list.files(outdir), value=TRUE)
    for (f in 1:length(files)) {
        print(paste("Reading file", files[f]))
        filename = paste(outdir, files[f], sep="/")
        data[[f]] = read.table(filename, header=T)
    }
    names(data) = gsub(".tsv", "", files)
    return(data)
}


