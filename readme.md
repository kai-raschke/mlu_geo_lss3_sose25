# LSS3 - SoSe 2025

## Data

### Sentinel-2-export

Gathered data from Earth Engine about the test field

### Sentinel-1 Data: 

for storage reasons data is not supplied. Can be downloaded via Copernicus Dataspace Ecosystem with following attributes:
Orbit 44, Polygon: POLYGON((11.901171 51.141638, 11.924750 51.141638, 11.924750 51.155587, 11.901171 51.155587, 11.901171 51.141638)), Timerange: 01.04 - 24.06.25, Bands VV + VH, Ascending, S1A + S1C

### Sentinel-2 Copernicus Browser:
Can be downloaded via Copernicus Dataspace Ecosystem with following attributes:
Timerange 3.3. - 24.6.25, 50% Cloudcoverage, Polygon: POLYGON((11.901171 51.141638, 11.924750 51.141638, 11.924750 51.155587, 11.901171 51.155587, 11.901171 51.141638))

## Processing

### Sentinel-1_SNAP_processing_graph.xml

Processing Chain for SNAP used to prerpocess Sentinel-1 Data to radiometrically corrected linear backscatter VV and VH in sigma nought

### Sentinel-1_processing.R

processing of preprocessed Sentinel-1 Data: calculation of indices, interpolation, plotting

### Sentinel-1_cor_test.R

performes normality test, correlation test and plotting of processed Sentinel-1 Data

### Chlor_NDVI_final.R

R script for processing and analysing field data (chlorophyll) and NDVI time series from Sentinel-2 and UAV data. Contains data cleansing, interpolation, correlation tests and visualisations (scatterplots, line diagrams). Results are exported as CSV and PNG.

### NDMI_Plotha_B08_B11.R

NDMI calculation and extraction from Sentinel-2 data

### NDMI_Plotha_B08A_B11.R

NDMI calculation and extraction from Sentinel-2 data

### NDVI_Plotha_B04_B08.R

NDVI calculation and extraction from Sentinel-2 data

### NDVI_Plotha_UAV_B03_B05.R

NDVI calculation and extraction from UAV data

### LSS3_full_14062025.xlsx

Calculation of standard error

### Timeline (interpoliert, mit UAV 28.4.), Analysen & Grafiken_Wasser_rel.xlsx

Time series, descriptive statistical analysis, diagrams

### esu_shapes.shp

Polygon shapefile used for UAV data

###

Point shapefile used for Sentinel-2 data


