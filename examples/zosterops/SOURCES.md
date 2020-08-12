# DATA SOURCES

https://earthexplorer.usgs.gov/

## Forest cover

USGS LP-DAAC GFCC30TC
(https://lpdaac.usgs.gov/products/gfcc30tcv003/)

Landsat 5 & 7 (2015)

```
VALUE 	DESCRIPTION

0-100 	Percent of pixel area covered by tree cover
200 	Water
210 	Cloud
211 	Shadow
220 	Fill Value
```

## Elevation

USGS Digital Elevation SRTM 1 Arc-Second Global
(https://cmr.earthdata.nasa.gov/search/concepts/C1220567890-USGS_LTA.html)

Shuttle Radar Topography Mission (2000)

--

CGIAR-CSI (http://srtm.csi.cgiar.org/)

SRTM 90m DEM Digital Elevation Database, interpolated to fill no-data regions

Jarvis A., H.I. Reuter, A.  Nelson, E. Guevara, 2008, Hole-filled  seamless SRTM
data V4, International  Centre for Tropical  Agriculture (CIAT), available  from
http://srtm.csi.cgiar.org.


# PROCESSING

## Resampling

Cubic convolution resampling (30m to 100m resolution):
	* https://support.esri.com/en/technical-article/000005606
	* https://docs.qgis.org/3.10/en/docs/user_manual/working_with_raster/raster_analysis.html

## Accessing

With R (`raster` package):
	* https://www.earthdatascience.org/courses/earth-analytics/lidar-raster-data-r/introduction-to-spatial-metadata-r/
	* https://www.rdocumentation.org/packages/raster/versions/3.0-12
