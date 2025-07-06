packages <- c("terra", "tidyr", "ggplot2", "openxlsx", "dplyr")
lapply(packages, library, character.only = TRUE)


############## read data ###########
NDVI <- read.xlsx("./data/Seb_NDVI.xlsx", sheet = 1)
field_1 <- terra::vect("./data/vect/Field_1_WGS84.gpkg")
ssu_location <- vect(r"(data\vect\shapes\SSU.shp)")
tif_files <- list.files("./data/rast/CDSE_processed", full.names = TRUE)
in_situ_bestand <- read.xlsx(r"(data\LSS3_full_14062025.xlsx)", sheet = 1)
in_situ_bestand_SSUs <- read.xlsx(r"(data\LSS3_full_14062025.xlsx)", sheet = 1)
in_situ_labor <- read.xlsx(r"(data\LSS3_full_14062025.xlsx)", sheet = 5)
in_situ_lai <- read.xlsx(r"(data\LSS3_full_14062025.xlsx)", sheet = 2) # Spaltenname LAI in der Excel Dati setzen!!

if(!crs(ssu_location) == crs(rast(tif_files[1]))) {
  ssu_location <- terra::project(ssu_location, crs(rast(tif_files[1])))
}
### in_situ subset bestand
in_situ_bestand <- in_situ_bestand[, c("Datum", "Height", "SM")]
in_situ_bestand <- in_situ_bestand[!is.na(in_situ_bestand$SM), ]
in_situ_bestand <- in_situ_bestand[ave(in_situ_bestand$Datum, in_situ_bestand$Datum, FUN = seq_along) == 2, ]
in_situ_bestand$Datum <- as.Date(in_situ_bestand$Datum, origin = "1899-12-30")

### in_situ subset SSUs
in_situ_bestand_SSUs <- in_situ_bestand_SSUs[in_situ_bestand_SSUs$SSU %in% c("4","5","6","A","B","C","D"), c("ESU","SSU", "Datum", "SM_final", "Height_final")]
in_situ_bestand_SSUs$Datum <- as.Date(in_situ_bestand_SSUs$Datum, origin = "1899-12-30")
in_situ_bestand_SSUs <- in_situ_bestand_SSUs %>%
  mutate(SSU_No = paste0("ESU", ESU, "SSU", SSU))
names(in_situ_bestand_SSUs) <- c("ESU", "SSU", "Date", "SM_final", "Height_final", "SSU_No")

### in_situ subset labor
in_situ_labor <- in_situ_labor[, c("Datum", "FM", "TM")]
in_situ_labor <- in_situ_labor[!is.na(in_situ_labor$FM), ]
in_situ_labor <- in_situ_labor[ave(in_situ_labor$Datum, in_situ_labor$Datum, FUN = seq_along) == 2, ]
in_situ_labor$Datum <- as.Date(in_situ_labor$Datum, origin = "1899-12-30")
in_situ_labor$moisture <- (in_situ_labor$FM - in_situ_labor$TM) / in_situ_labor$FM
in_situ_labor <- in_situ_labor[, c("Datum", "moisture")]
in_situ_labor[11,2] <- 0.178 # linear interpolieren weil fehlerhaft: trockenmasse höher als feuchtemasse

### in_situ subset LAI
in_situ_lai <- in_situ_lai[, c("Datum", "LAI")] # Spaltenname muss in der Excel Dati manuel gesetzt werden!!
in_situ_lai <- in_situ_lai[!is.na(in_situ_lai$LAI), ]
in_situ_lai <- in_situ_lai[ave(in_situ_lai$Datum, in_situ_lai$Datum, FUN = seq_along) == 2, ]
in_situ_lai$Datum <- as.Date(in_situ_lai$Datum, origin = "1899-12-30")

in_situ <- cbind(in_situ_bestand, PM = in_situ_labor[,-1], lai = in_situ_lai[,-1])

# subset NDVI
NDVI <- NDVI[,c(1,4)]
NDVI$Datum <- as.Date(NDVI$Datum, origin = "1899-12-30")
names(NDVI) <- c("Datum", "NDVI")

S1_dates <- as.Date(character())
reflectance_list <- list()
reflectance_all <- data.frame(SSU_No = ssu_location$SSU_No)
############### process data ###########
### calculate RVI
for (i in seq_along(tif_files)) {
  rast <- rast(tif_files[i])
  names(rast) <- c("VH", "VV")

  # calculate RVI
  rast$RVI <- (4 * rast$VH) / (rast$VH + rast$VV)
  rast$RVI4S1 <- (sqrt(1 - (rast$VV / (rast$VH + rast$VV)))) * (4 * (rast$VH / (rast$VH + rast$VV)))#(1 - (rast$VH / rast$VV)) / (1 + (rast$VH / rast$VV)) #RVI4S1 
  rast$ratio <- rast$VH / rast$VV
  rast <- rast[[c("RVI", "RVI4S1", "ratio")]] # select bands
  satellite <- substr(basename(tif_files[i]), 8, 10)
  date <- substr(basename(tif_files[i]), 25, 32)
  names(rast) <- c(paste0(satellite, "_", date, "_RVI"),
                   paste0(satellite, "_", date, "_RVI4S1"),
                   paste0(satellite, "_", date, "_ratio"))

  #available dates
  S1_dates <- c(av_dates, as.Date(date, format = "%Y%m%d"))
  values <- terra::extract(rast, ssu_location, ID = FALSE)
  #reflectance_list[[i]] <- values
  reflectance_all <- cbind(reflectance_all, values)
}
#safe
write.csv(reflectance_all, "./data/reflectance_all.csv", row.names = FALSE)

# subset of SSUs in the mid of the field
reflectance_all <- reflectance_all[grepl("SSU4$|SSU5$|SSU6$|SSUA$|SSUB$|SSUC$|SSUD$", reflectance_all$SSU_No), ]
reflectance_all <- reflectance_all[!grepl("^ESU1.*", reflectance_all$SSU_No), ] # remove ESU1 SSUs

reflectance_mean <- as.data.frame(t(colMeans(reflectance_all[,-1])))

# longformat
cols_to_pivot <- names(reflectance_all)[2:length(names(reflectance_all))]

reflectance_long <- pivot_longer(
  reflectance_all,
  cols = all_of(cols_to_pivot),
  names_to = c("Satellite", "Date", "Band"),
  names_sep = "_",
  values_to = "Reflectance"
)
reflectance_long$Date <- as.Date(reflectance_long$Date, format = "%Y%m%d")

cols_to_pivot <- names(reflectance_mean)[1:length(names(reflectance_mean))]
reflectance_mean_long <- pivot_longer(
                 reflectance_mean,
                 cols = all_of(cols_to_pivot),
                 names_to = c("Satellite", "Date", "Band"),
                 names_sep = "_",
                 values_to = "Reflectance")

reflectance_mean_long$Date <- as.Date(reflectance_mean_long$Date, format = "%Y%m%d")

################ interpolate radar reflectance ####
RVI <- approx(unique(reflectance_mean_long$Date), reflectance_mean_long$Reflectance[reflectance_mean_long$Band == "RVI"], xout = in_situ_lai$Datum,)
RVI4S1 <- approx(unique(reflectance_mean_long$Date), reflectance_mean_long$Reflectance[reflectance_mean_long$Band == "RVI4S1"], xout = in_situ_lai$Datum,)
ratio <- approx(unique(reflectance_mean_long$Date), reflectance_mean_long$Reflectance[reflectance_mean_long$Band == "ratio"],xout = in_situ_lai$Datum,)

all_mean_data <- cbind(in_situ, RVI = RVI$y, RVI4S1 = RVI4S1$y, ratio = ratio$y)
# save all data
write.csv(all_mean_data, "./data/all_mean_data.csv", row.names = FALSE)

################# interpolate radar reflectance for each SSU ######
SSUs <- unique(reflectance_long$SSU_No)

in_situ_bestand$Datum == in_situ_lai$Datum
reflectance_long_interpol <- data.frame()
for (SSU in SSUs) {
  # subset for each SSU
  reflectance_ssu <- reflectance_long[reflectance_long$SSU_No == SSU, ]
  print(paste("Processing SSU:", SSU))
  print(head(reflectance_ssu))
  # interpolate
  RVI <- approx(reflectance_ssu$Date[reflectance_ssu$Band == "RVI"],
                reflectance_ssu$Reflectance[reflectance_ssu$Band == "RVI"],
                xout = in_situ_bestand$Datum)
  RVI4S1 <- approx(reflectance_ssu$Date[reflectance_ssu$Band == "RVI4S1"],
                   reflectance_ssu$Reflectance[reflectance_ssu$Band == "RVI4S1"],
                   xout = in_situ_bestand$Datum)
  ratio <- approx(reflectance_ssu$Date[reflectance_ssu$Band == "ratio"],
                  reflectance_ssu$Reflectance[reflectance_ssu$Band == "ratio"],
                  xout = in_situ_bestand$Datum)
  # combine all data
  reflectance_long_interpol <- rbind(reflectance_long_interpol,
                                     data.frame(SSU_No = SSU,
                                                Date = rep(in_situ_bestand$Datum, 3),
                                                Band = c(rep("RVI", length(in_situ_bestand$Datum)),
                                                         rep("RVI4S1", length(in_situ_bestand$Datum)),
                                                         rep("ratio", length(in_situ_bestand$Datum))),
                                                Reflectance = c(RVI$y, RVI4S1$y, ratio$y)))
}

# safe
write.csv(reflectance_long_interpol, "./data/reflectance_long_interpol.csv", row.names = FALSE)

# merge data frames
all_data_interpol <- merge.data.frame(reflectance_long_interpol, in_situ_bestand_SSUs, by = c("Date", "SSU_No"), all = TRUE)

# safe
write.csv(all_data_interpol, "./data/all_data_interpol.csv", row.names = FALSE)

# daten prüfen
mean(all_data_interpol$Reflectance[all_data_interpol$Band == "RVI" & all_data_interpol$Date == "2025-04-08"])
mean(all_data_interpol$Reflectance[all_data_interpol$Band == "RVI4S1" & all_data_interpol$Date == "2025-04-08"])

################## plots ####
windows(width = 800, height = 600)
ggplot() +
  geom_line(
    data = in_situ_labor, aes(x = as.Date(Datum, format = "%Y-%m-%d"), y = moisture, color = "Wassergehalt"), size = 1
  ) +
  geom_line(
    data = reflectance_mean_long[reflectance_mean_long$Band == "RVI", ],
    aes(x = as.Date(Date, format = "%Y%m%d"), y = Reflectance, color = "RVI"), size = 1
  ) +
  geom_line(
    data = in_situ_bestand, aes(x = as.Date(Datum, format = "%Y-%m-%d"), y = Height/100, color = "Höhe"), size = 1
  ) +
  geom_line(
    data = in_situ_bestand, aes(x = as.Date(Datum, format = "%Y-%m-%d"), y = SM/100, color = "Bodenfeuchte"), size = 1
  ) +
  labs(x = "Datum",
       y = "RVI / Wassergehalt / Höhe",
       color = "",
  ) +
  scale_x_date(breaks = in_situ$Datum[seq(1, length(in_situ$Datum), by = 2)], date_labels = "%d.%m.%y" ) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),  # 45 Grad, rechtsbündig
    legend.position = "bottom",
    axis.ticks.x = element_line(size = 0.5, color = "grey"),   # x-Achsen-Ticks
    axis.ticks.length = unit(0.15, "cm")
  ) +
  scale_color_manual(values = c("RVI" = "black",
                             #   "ratio" = "yellow",
                                "Bodenfeuchte" = "lightblue",
                                "Höhe" = "red",
                                "Wassergehalt" = "darkblue",
                             #   "RVI4S1" = "#535353",
                                "LAI / 10" = "green"))

data_long <- pivot_longer(
  all_mean_data,
  cols = c("RVI", "RVI4S1", "ratio", "SM", "Height", "PM", "lai"),
  names_to = "Band",
  values_to = "Reflectance"
)
data_long$Datum <- as.Date(data_long$Datum)

# vgl alle Zeitreihen
ggplot() +
  geom_line(data = data_long[data_long$Band == "RVI", ], aes(x = Datum, y = Reflectance, color = "RVI")) + 
  geom_line(data = data_long[data_long$Band == "RVI4S1", ], aes(x = Datum, y = Reflectance, color = "RVI4S1")) +
  geom_line(data = data_long[data_long$Band == "ratio", ], aes(x = Datum, y = Reflectance, color = "ratio")) + 
  geom_line(data = data_long[data_long$Band == "SM", ], aes(x = Datum, y = Reflectance/200, color = "SM")) +
  geom_line(data = data_long[data_long$Band == "Height", ], aes(x = Datum, y = Reflectance/200, color = "Height")) +
  geom_line(data = data_long[data_long$Band == "PM", ], aes(x = Datum, y = Reflectance, color = "PM")) +
  geom_line(data = data_long[data_long$Band == "lai", ], aes(x = Datum, y = Reflectance/4, color = "LAI")) +
  labs(x = "Date",
       y = "Reflectance",
       color = "Data",
       title = "Comparison of Sentinel-1 Radar reflectance VV + VH, Radar Vegetation Index (RVI), LAI, Soil Moisture Plant Moisture and Plant Height") +  
  scale_color_manual(values = c("RVI" = "black",
                                "RVI4S1" = "darkgrey",
                                "ratio" = "yellow",
                                "SM" = "lightblue",
                                "Height" = "red",
                                "PM" = "darkblue",
                                "LAI" = "green")) 

# Scatterplot of RVI and Height
ggplot(all_mean_data, aes(x = RVI, y = Height)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Scatterplot of Reflectance and Height",
       x = "Reflectance",
       y = "Height") +
  theme_minimal() +
  scale_color_discrete() +
  theme(legend.position = "right")
