library (terra)
setwd ("/Users/sebastianschmidt/Documents/Uni/Halle/2. Sem/LSS 3 - Plotha/LSS3 Gruppenarbeit/Tiff/SentBrowser/B8A_11")
shp <- vect("/Users/sebastianschmidt/Documents/Uni/Halle/2. Sem/LSS 3 - Plotha/LSS3 Gruppenarbeit/Shapes/SSU_WGS84.shp")
rasterdateien <- dir(pattern="tiff")
daten <- unique (substr (rasterdateien,1,10)) ## Zeitschnitte aus den Dateinamen ausschneiden
b11 <- rasterdateien[seq(1,40,2)] ## die 1., 3. 5., .... Datei bestimmen
b8A <- rasterdateien[seq(2,40,2)] ## die 2., 4. 6., ... Datei bestimmen

band8A <- rast (b8A[1]) ## lies die erste Band 8A Datei ein
band11 <- rast (b11[1]) ## lies die erste Band 11 Datei ein
ndmi <- (band8A - band11)/(band8A + band11) ## berechne daraus den ndmi

for (i in 2:length(b8A)){ ## Setze der Reihe nach 2,3,4,5,... fuer i ein
  band8A <- rast (b8A[i]) ## lies das i-te Band 4 ein
  band11 <- rast (b11[i]) ## lies das i-te Band 8 ein
  ndmi.neu <- (band8A - band11)/(band8A + band11) ## berechne den ndmi fuer diese beiden Baender
  ndmi.neu <- resample (ndmi.neu, ndmi) ## resample den neuen ndmi auf den ersten ndmi => sind jetzt deckungsgleich
  ndmi <- c(ndmi,ndmi.neu) ## haeng den neuen ndmi an den bestehenden dran => hat jetzt ein weiteres Band
}
names(ndmi)<-daten ## benenne die Baender im ndmi-Stapel mit den Daten

plothandmi <- extract (ndmi, shp) ## Lies die ndmi-Werte fuer die SSUs aus. Punkte ausserhalb bekommen NA als Wert

write.table(plothandmi,"ndmi_plotha_b8A_11.csv") ## Schreib die Tabelle auf die Festplatte