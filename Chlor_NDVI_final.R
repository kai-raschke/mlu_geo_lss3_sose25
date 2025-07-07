library(terra)      # für Rasterverarbeitung
library(sf)         # für Vektor- und Shapefile-Verarbeitung
library(dplyr)      # für Datenmanipulation
library(ggplot2)    # für Visualisierung
library(openxlsx)

# Felddaten einladen ###############################################################

# Arbeitsverzeichnis anpassen
setwd("~/Studium/MA_Halle/FS2/LSS3/Poster/LSS3_Plotha")

# Pfad zum SSU-Shapefile
ssu_path <- "shapes_SSU/SSU.shp"
esu_path <- "shapes_SSU/ESUs.shp"
#Pfad zum esu_shapes.geojson mit den Polygonen 
esu_polygons_path <- "shapes_SSU/kai_esu_shapes.geojson"

# SSU und ESU laden
ssu <- st_read(ssu_path)
esu <- st_read(esu_path)

# Polygone für die ESUs einladen, mit denen wir Sentinel2 und UAV ndvi berechnen
esu_polygons <- st_read(esu_polygons_path)

# Chlorophyll-Daten einlesen ####################################################

chlor <- read.csv("Chlorophyll_24062025.csv", sep = ";", dec = ",") %>%
  mutate(
    Datum = as.Date(Datum, format = "%d.%m.%Y"),
    ID = paste0("ESU", ESU, "SSU", SSU),
    Quelle = "Chlorophyll"
  )


# Timeline interpoliert von Sebastian mit NDVI-Werten einladen ##################
timeline_NDVI <- read.xlsx("Timeline_interpoliert.xlsx", sheet = 1)


# S2 NDVI aus Timeline extrahieren ##############################################
# nur die Spalten "Datum" und "NDVI" auswählen
timeline_s2 <- timeline_NDVI[, c(1,6)]
timeline_s2$Datum <- as.Date(timeline_s2$Datum, origin = "1899-12-30")

# UAV NDVI aus Timeline extrahieren #############################################
# nur die Spalten "Datum" und "NDVI" auswählen
timeline_UAV <- timeline_NDVI[, c(1,8)]
timeline_UAV$Datum <- as.Date(timeline_UAV$Datum, origin = "1899-12-30")
# Entferne Zeilen, in denen NDVI.(UAV.B03/B05) NA ist
timeline_UAV <- timeline_UAV[!is.na(timeline_UAV$`NDVI.(UAV.B03/B05)`), ]


# SSU-Daten anpassen und säubern ####################################################

# Sicherstellen, dass beide dieselbe Projektion haben
ssu <- st_transform(ssu, crs = st_crs(esu))

# Räumlicher Join
ssu_joined <- st_join(ssu, esu[, c("ESU_No")])

ssu_joined$ESU_No[ssu_joined$SSU_No %in% c("ESU4SSU3", "ESU4SSU6", "ESU4SSU4", "ESU4SSU1", "ESU4SSU2")] <- "4"
ssu_joined$ESU_No[ssu_joined$SSU_No %in% c("ESU2SSU3", "ESU2SSU6", "ESU2SSU4", "ESU2SSU1", "ESU2SSU2")] <- "2"
ssu_joined$ESU_No[ssu_joined$SSU_No %in% c("ESU3SSU2", "ESU3SSU4", "ESU3SSU1", "ESU3SSU6", "ESU3SSU3")] <- "3"

# Jetzt alle ESU1 löschen; für uns nicht relevant 
ssu_joined <- ssu_joined %>% filter(!ESU_No == 1)

# Filtern nach mittigen SSUs
ssu_joined <- ssu_joined[grepl("SSU4$|SSU5$|SSU6$|SSUA$|SSUB$|SSUC$|SSUD$", ssu_joined$SSU_No),]


# Chlorophyll-Daten filtern und Mittelwerte berechnen mit Quelle-Spalte #################


spad_cols <- grep("^SPAD", names(chlor), value = TRUE)

#Konvertiere SPAD-Spalten zu numerisch
chlor[spad_cols] <- lapply(chlor[spad_cols], function(x) as.numeric(as.character(x)))


chlor_filter <- chlor %>%
  filter(grepl("SSU4$|SSU5$|SSU6$|SSUA$|SSUB$|SSUC$|SSUD$", ID))

# Zeilenweise Mittelwert der SPAD-Spalten berechnen
chlor_filter <- chlor_filter %>%
  rowwise() %>%
  mutate(mean_spad = mean(c_across(starts_with("SPAD")), na.rm = TRUE)) %>%
  ungroup()

# Tagesmittelwert berechnen und Quelle hinzufügen
chlor_means <- chlor_filter %>%
  group_by(Datum) %>%
  summarise(mean_chlor = mean(mean_spad, na.rm = TRUE)) %>%
  rename(Datum.chlor = Datum) %>%
  mutate(Quelle = "Chlorophyll") %>%
  filter(!is.na(mean_chlor)) 

# Ausgabe prüfen
print(chlor_means)


# Interpolation beider timelines ################################################ 

# Interpolationsfunktion definieren
interpolate_ndvi <- function(ndvi_df, chlor_dates) {
  approx_result <- approx(
    x = as.numeric(ndvi_df$Datum.ndvi),
    y = ndvi_df$mean_ndvi,
    xout = as.numeric(chlor_dates),
    method = "linear",
    rule = 2  # extrapoliere mit konstantem Wert außerhalb
  )
  
  data.frame(
    Datum.chlor = as.Date(approx_result$x, origin = "1970-01-01"),
    mean_ndvi_interp = approx_result$y
  )
}


# für S2
ndvi_s2_interp_clean <- interpolate_ndvi(
  ndvi_df = timeline_s2 %>% rename(Datum.ndvi = Datum, mean_ndvi = `NDVI.(S2.B04/B08)`),
  chlor_dates = chlor_means$Datum.chlor
)

# für UAV 
ndvi_uav_interp_clean <- interpolate_ndvi(
  ndvi_df = timeline_UAV %>% rename(Datum.ndvi = Datum, mean_ndvi = `NDVI.(UAV.B03/B05)`),
  chlor_dates = chlor_means$Datum.chlor
)


# Interpolierte NDVI-Werte an Chlorophyll-Mittelwerte anhängen ###################

# Sentinel-2 (S2) NDVI an Chlorophyll-Daten anhängen
chlor_interp_s2 <- chlor_means %>%
  left_join(ndvi_s2_interp_clean, by = "Datum.chlor") %>%
  mutate(Quelle = "Sentinel")

# UAV NDVI an Chlorophyll-Daten anhängen
chlor_interp_uav <- chlor_means %>%
  left_join(ndvi_uav_interp_clean, by = "Datum.chlor") %>%
  mutate(Quelle = "UAV")

# Beide Quellen zusammenführen ###################################################
chlor_ndvi_all <- bind_rows(chlor_interp_s2, chlor_interp_uav)

# Ergebnis prüfen
print(chlor_ndvi_all)

# Optional: in CSV exportieren
write.csv(chlor_ndvi_all, file = "Ergebnisse/chlor_ndvi_all.csv", row.names = FALSE)



# auf Normalverteilung prüfen ###################################################


library(dplyr)
library(tidyr)

#Vorbereitung:
# Wide-Format erstellen: mean_ndvi_interp von Sentinel und UAV als separate Spalten, um direkt vergleichen zu können 
chlor_ndvi_wide <- chlor_ndvi_all %>%
  select(Datum.chlor, mean_chlor, Quelle, mean_ndvi_interp) %>%
  pivot_wider(
    names_from = Quelle,
    values_from = mean_ndvi_interp
  )

print(chlor_ndvi_wide)
write.csv(chlor_ndvi_wide, file = "Ergebnisse/chlor_ndvi_wide.csv", row.names = FALSE)


# chlor_ndvi_wide auf Normalverteilung prüfen
shapiro.test(chlor_ndvi_wide$mean_chlor) # nicht normalverteilt
hist(chlor_ndvi_wide$mean_chlor)
qqnorm(chlor_ndvi_wide$mean_chlor); qqline(chlor_ndvi_wide$mean_chlor)

## chlor daten log-transformieren? 
# chlor_ndvi_wide <- chlor_ndvi_wide %>%
#   mutate(mean_chlor_log = log(mean_chlor))
# 
# # Verteilung prüfen
# shapiro.test(chlor_ndvi_wide$mean_chlor_log)
# hist(chlor_ndvi_wide$mean_chlor_log)
# qqnorm(chlor_ndvi_wide$mean_chlor_log); qqline(chlor_ndvi_wide$mean_chlor_log)
## sind weiterhin nicht normalverteilt, daher bei chlor_means bleiben ohne transformation

shapiro.test(chlor_ndvi_wide$Sentinel) # nicht normalverteilt
hist(chlor_ndvi_wide$Sentinel)
qqnorm(chlor_ndvi_wide$Sentinel); qqline(chlor_ndvi_wide$Sentinel)

shapiro.test(chlor_ndvi_wide$UAV) # nicht normalverteilt
hist(chlor_ndvi_wide$UAV)
qqnorm(chlor_ndvi_wide$UAV); qqline(chlor_ndvi_wide$UAV)


# Korrelation Chlorophyll vs. Sentinel NDVI #####################################
# s2 normalverteilt, UAV Daten nicht normalverteilt
cor_sentinel <- cor.test(chlor_ndvi_wide$mean_chlor, chlor_ndvi_wide$Sentinel, method = "spearman")
print(cor_sentinel)

# Korrelation Chlorophyll vs. UAV NDVI
cor_uav <- cor.test(chlor_ndvi_wide$mean_chlor, chlor_ndvi_wide$UAV, method = "spearman")
print(cor_uav)


# Visualisieren #################################################################

library(ggplot2)

# Sentinel und UAV getrennt plotten:

p_sentinel <- ggplot(chlor_ndvi_wide, aes(x = mean_chlor, y = Sentinel)) +
  geom_point(color = "limegreen", size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "chartreuse4") +
  labs(title = "Chlorophyll vs Sentinel2-NDVI",
       x = "Mean Chlorophyll",
       y = "Sentinel2 Mean NDVI") +
  theme_minimal()

# Plot UAV NDVI vs Chlorophyll
p_uav <- ggplot(chlor_ndvi_wide, aes(x = mean_chlor, y = UAV)) +
  geom_point(color = "cornflowerblue", size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "mediumblue") +
  labs(title = "Chlorophyll vs UAV-NDVI",
       x = "Mean Chlorophyll",
       y = "UAV Mean NDVI") +
  theme_minimal()

# Plots anzeigen
print(p_sentinel)
print(p_uav)

# zusammen anzeigen lassen 

#install.packages("gtable")
library(gtable)
#install.packages("patchwork")
library(patchwork)

# z.B. y-Achse von 0 bis 1, in Schritten von 0.1
y_limits <- c(0.55, 1)
y_breaks <- seq(0, 1, by = 0.1)

p_sentinel <- p_sentinel + 
  scale_y_continuous(limits = y_limits, breaks = y_breaks)

p_uav <- p_uav + 
  scale_y_continuous(limits = y_limits, breaks = y_breaks)

(p_sentinel | p_uav)


# p_uav <- p_uav + 
#   geom_smooth(method = "loess", se = TRUE, color = "red", fill = "red", alpha = 0.2)
# 
# # Dann nebeneinander anzeigen (untereinander mit / ) 
# (p_sentinel | p_uav)


library(ggplot2)
library(dplyr)
library(tidyr)


# Boxplot

# Sentinel2, UAV und Chlorophyll zusammen plotten: 
# Daten ins long Format bringen (Sentinel & UAV NDVI)
# chlor_ndvi_long <- chlor_ndvi_wide %>%
#   pivot_longer(cols = c("Sentinel", "UAV"), names_to = "Quelle", values_to = "NDVI")
# 
# ggplot(chlor_ndvi_long, aes(x = Quelle, y = NDVI)) +
#   geom_boxplot(alpha = 0.5, fill = "lightblue") +   # Boxplots NDVI je Quelle
#   geom_jitter(aes(color = mean_chlor), width = 0.2, size = 3) +  # Punkte NDVI gefärbt nach Chlorophyll
#   scale_color_viridis_c(option = "magma") +
#   labs(title = "Verteilung der NDVI-Werte je Quelle mit Chlorophyll-Farbcode",
#        x = "Quelle",
#        y = "NDVI",
#        color = "Chlorophyll") +
#   theme_minimal()


# Scatterplot alles zusammen 

# Sentinel2, UAV und Chlorophyll zusammen plotten: 
# Daten ins long Format bringen (Sentinel & UAV NDVI)
chlor_ndvi_long <- chlor_ndvi_wide %>%
  pivot_longer(cols = c("Sentinel", "UAV"), names_to = "Quelle", values_to = "NDVI")

# Plot mit Chlorophyll auf x-Achse, NDVI auf y-Achse, nach Quelle farblich unterschieden
scatter <- ggplot(chlor_ndvi_long, aes(x = mean_chlor, y = NDVI, color = Quelle)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1) +
  labs(
   # title = "Chlorophyll vs. NDVI (Sentinel-2 & UAV)",
    x = "Mean Chlorophyll",
    y = "Mean NDVI"
  ) +
  scale_color_manual(
    name = "Datenquelle",
    values = c("Sentinel" = "green", "UAV" = "green4")
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey92", size = 0.3),
    panel.grid.minor = element_line(color = "grey97", size = 0.2),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom"
  )

print(scatter)

# Speichern
ggsave("Ergebnisse/chlor_ndvi_scatter.png", plot = scatter, width = 8, height = 3.5, dpi = 300)
# Nehmen diesen Scatterplot fürs Poster 



# Liniendiagramm mit Legende ###################################################

line <- ggplot(chlor_ndvi_wide, aes(x = Datum.chlor)) +
  geom_line(aes(y = mean_chlor / 40, color = "Chlorophyll / 40"), linewidth = 1, lineend = "round") +
  geom_line(aes(y = Sentinel, color = "Sentinel NDVI"), linewidth = 1, lineend = "round") +
  geom_line(aes(y = UAV, color = "UAV NDVI"), linewidth = 1, lineend = "round") +
  scale_color_manual(
    name = "Messgröße",
    values = c("Chlorophyll / 40" = "purple4", 
               "Sentinel NDVI" = "green", 
               "UAV NDVI" = "green4")
  ) +
  scale_x_date(date_labels = "%d.%m.%y") +  # <- Datumsformatierung OHNE **
  labs(x = "Datum", y = "Wert\n(Chlorophyll: SPAD Units / 40, NDVI)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey92", size = 0.3),
    panel.grid.minor = element_line(color = "grey97", size = 0.2),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 8)
  )
print(line)

# Speichern
ggsave("Ergebnisse/chlor_ndvi_line.png", plot = line,  width = 6, height = 3, dpi = 300)

