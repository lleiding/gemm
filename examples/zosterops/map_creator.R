#!/usr/bin/Rscript
###
### This script reads in satellite data (as GeoTIFF files) of elevation
### and forest cover, visualises these and converts them into a map file
### for the GeMM model.
###
### (c) Daniel Vedder 2020
### Licensed under the terms of the MIT license.
###

library(rgdal)
library(raster)
library(RStoolbox)
library(ggplot2)

elevation_file = "taita_elevation.tif"
forest_file = "taita_forest_cover.tif"

image_output_file = "taita_map.jpg"
map_output_file = "taita_hills.map"

## Read in the GeoTIFF files and make sure they're workable
loadData = function(el_file=elevation_file, fo_file=forest_file, script=TRUE)
{
    ## Read in the data
    print("Reading GeoTIFF files")
    elevation_data = raster(el_file)
    forest_data = raster(fo_file)
    ## Do some sanity checks
    print("Checking for compatibility")
    if (!all(as.vector(extent(elevation_data)) == as.vector(extent(forest_data)))) {
        print("Mismatching input data extent. Aborting.")
        if (script) quit()
    }
    if (!all(dim(elevation_data) == dim(forest_data))) {
        print("Mismatching input data resolution. Aborting.")
        if (script) quit()
    }
    ## Identify bad data caused by clouds
    print("Fixing bad data")
    elevation_data[which(as.vector(elevation_data) < 0)] = NA
    forest_data[which(as.vector(forest_data) > 100)] = NA
    ## Return both raster objects
    return(c(elevation_data, forest_data))
}

## Take in two rasters (elevation and forest cover) and create a 2D map image
visualiseMap = function(elevation, forest, out=image_output_file)
{
    ggplot(data=forest) +
        geom_raster(aes(fill=taita_forest_cover, x=x, y=y)) +
        geom_contour(data=elevation, binwidth=200, colour="gray30",
                     aes(x=x, y=y, z=taita_elevation)) +
        scale_fill_gradient(high="forestgreen", low="greenyellow", na.value="gray90") +
        guides(fill=guide_legend(title="Forest cover")) +
        labs(caption="Taita Hills, Kenya") +
        coord_fixed() +
        theme_void()
    ggsave(out)
}

## Take in two rasters (elevation and forest cover) and create a 3D map image
visualise3DMap = function(elevation, forest, out=image_output_file)
{
    #TODO create 3D map with rayshader
}

## Take in two rasters (elevation and forest cover) and create a GeMM map file
convertMap = function(elevation, forest, out=map_output_file)
{
    ##TODO
}

run = function()
{
    data = loadData()
    elevation = data[1][[1]]
    forest = data[2][[1]]
    visualiseMap(elevation, forest)
    convertMap(elevation, forest)
}

#run()
