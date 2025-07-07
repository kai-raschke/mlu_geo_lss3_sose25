library (terra)
setwd ("E:/Desktop/B04_08")
shp <- vect("E:/Desktop/B04_08/SSU_WGS84.shp")
rasterdateien <- dir(pattern="tiff")
daten <- unique (substr (rasterdateien,1,10)) ## Zeitschnitte aus den Dateinamen ausschneiden
b04 <- rasterdateien[seq(1,40,2)] ## die 1., 3. 5., .... Datei bestimmen
b08 <- rasterdateien[seq(2,40,2)] ## die 2., 4. 6., ... Datei bestimmen

band04 <- rast (b04[1]) ## lies die erste Band 04 Datei ein
band08 <- rast (b08[1]) ## lies die erste Band 08 Datei ein
ndvi <- (band08 - band04)/(band08 + band04) ## berechne daraus den ndvi

for (i in 2:length(b04)){ ## Setze der Reihe nach 2,3,4,5,... fuer i ein
  band04 <- rast (b04[i]) ## lies das i-te Band 4 ein
  band08 <- rast (b08[i]) ## lies das i-te Band 8 ein
  ndvi.neu <- (band08 - band04)/(band08 + band04) ## berechne den ndvi fuer diese beiden Baender
  ndvi.neu <- resample (ndvi.neu, ndvi) ## resample den neuen ndvi auf den ersten ndvi => sind jetzt deckungsgleich
  ndvi <- c(ndvi,ndvi.neu) ## haeng den neuen ndvi an den bestehenden dran => hat jetzt ein weiteres Band
}
names(ndvi)<-daten ## benenne die Baender im ndvi-Stapel mit den Daten

plothandvi <- extract (ndvi, shp) ## Lies die ndvi-Werte fuer die SSUs aus. Punkte ausserhalb bekommen NA als Wert

write.table(plothandvi,"ndvi_plotha_b04_08.csv") ## Schreib die Tabelle auf die Festplatte