# LSS3 - SoSe 2025

## Data

### Sentinel2-export

Gathered data from Earth Engine about the test field

### Sentinel-1 Data: 

for storage reasons data is not supplied. Can be downloaded via Copernicus Dataspace Ecosystem with following attributes:
Orbit 44, Polygon: POLYGON((11.901171 51.141638, 11.924750 51.141638, 11.924750 51.155587, 11.901171 51.155587, 11.901171 51.141638)), Timerange: 01.04 - 24.06.25, Bands VV + VH, Ascending, S1A + S1C

## Processing

### Sentinel-1_SNAP_processing_graph.xml

Processing Chain for SNAP used to prerpocess Sentinel-1 Data to radiometrically corrected linear backscatter VV and VH in sigma nought

### Sentinel-1_processing.R

processing of preprocessed Sentinel-1 Data: calculation of indices, interpolation, plotting

### Sentinel-1_cor_test.R

performes normality test, correlation test and plotting of processed Sentinel-1 Data

