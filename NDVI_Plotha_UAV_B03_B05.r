library(terra)

setwd("/Users/sebastianschmidt/Downloads/UAV")

# Dateien einlesen
shp <- vect("/Users/sebastianschmidt/Documents/Uni/Halle/2. Sem/LSS 3 Plotha/LSS3 Gruppenarbeit/Shapes/esu_shape_von_Kai/esu_shapes.shp")
tif_files <- list.files(pattern = "\\.tif$", full.names = TRUE)

# Datenframe für die Ergebnisse anlegen
results <- data.frame(
  filename   = character(length(tif_files)),
  mean_ndvi  = numeric(length(tif_files)),
  stringsAsFactors = FALSE
)

# Schleife über alle Dateien
for (i in seq_along(tif_files)) {
  file <- tif_files[i]
  
  # Raster einlesen
  r <- rast(file)
  
  # NDVI berechnen: (NIR - Red) / (NIR + Red)
  ndvi <- (r[[5]] - r[[3]]) / (r[[5]] + r[[3]])
  
  # Auf das Shapefile-Polygon zuschneiden und maskieren
  ndvi_clip <- crop(ndvi, shp)            # räumliches Zuschneiden :contentReference[oaicite:0]{index=0}
  ndvi_mask <- mask(ndvi_clip, shp)       # außerhalb Polygon → NA
  
  # Mittelwert innerhalb des Polygons extrahieren
  ext <- extract(ndvi_mask, shp, fun = mean, na.rm = TRUE)  # :contentReference[oaicite:1]{index=1}
  
  # ID und Mean: ext hat Spalten ID und layer-Name
  results$filename[i]  <- basename(file)
  results$mean_ndvi[i] <- ext[,2]
}

# Ergebnisse als CSV speichern
write.csv(results,
          file = "ndvi_plotha_uav_b03_b05.csv",
          row.names = FALSE)