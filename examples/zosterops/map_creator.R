#!/usr/bin/Rscript
###
### This script reads in satellite data (as GeoTIFF files) of elevation
### and forest cover, visualises these and converts them into a map file
### for the GeMM model.
###
### (c) Daniel Vedder 2020
### Licensed under the terms of the MIT license.
###

library(raster)
library(ggplot2)

_elevation_data_file = "taita_elevation.tif"
_forest_data_file = "taita_forest_cover.tif"
_image_output_file = "taita_map.jpg"
_map_output_file = "taita_hills.map"

_elevation_data = NULL
_forest_data = NULL

## Read in a GeoTIFF file and return the main image raster
loadData = function(dataFile)
{
    ##TODO
    
}

## Take in two rasters (elevation and forest cover) and create a 3D map image
visualiseMap = function(elevation=_elevation_data, forest=_forest_data, out=_image_output_file)
{

}

## Take in two rasters (elevation and forest cover) and create a GeMM map file
convertMap = function(elevation=_elevation_data, forest=_forest_data, out=_map_output_file)
{

}

run = function()
{
    _elevation_data = loadData(_elevation_data_file)
    _forest_data = loadData(_forest_data_file)
    visualiseMap()
    convertMap()
}

run()
