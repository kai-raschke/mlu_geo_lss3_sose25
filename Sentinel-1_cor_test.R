packages <- c("corrplot", "greybox", "scatterplot3d", "rgl", "dplyr", "plotly", "ggplot2", "htmlwidgets", "patchwork")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)  # Paket installieren
    library(pkg, character.only = TRUE)  # Paket laden
  }
}

########### Lade die Daten ##########
all_mean_data <- read.csv("./data/all_mean_data.csv", stringsAsFactors = FALSE)
# Zusätzlich für kürzeren Zeitraum:
all_mean_data <- all_mean_data[all_mean_data$Datum != "2025-04-08", ] # Filtere den 8. April 2025 aus
all_mean_data <- all_mean_data[all_mean_data$Datum != "2025-06-24", ] # Filtere den 8. April 2025 aus

# auf Normalverteilung prüfen
for (col in c("RVI", "RVI4S1", "ratio", "Height", "lai", "PM", "SM")) {
  test <- shapiro.test(all_mean_data[[col]])
  if (test$p.value > 0.05) {
    cat(paste("Shapiro-Wilk-Test für", col, ": p =", test$p.value, "-> Normalverteilt\n"))
  } else {
    cat(paste("Shapiro-Wilk-Test für", col, ": p =", test$p.value, "-> Nicht normalverteilt, prüfe Log transfromation\n"))
    test <- shapiro.test(log(all_mean_data[[col]]))
    cat(paste("sahpiro-Wilk-Test für Log-Transformation für", col, ": p =", test$p.value, ifelse(
      test$p.value > 0.05,
      "Normalverteilt",
      "Nicht normalverteilt"
    ), "\n"))
    if (test$p.value > 0.05) {
      all_mean_data[[col]] <- log(all_mean_data[[col]])
      cat(paste("Log-Transformation für", col, "erfolgreich durchgeführt.\n"))
    } else {
      cat(paste("Log-Transformation für", col, "nicht erfolgreich, Daten bleiben unverändert. Spearman nutzen.\n"))
    }
  }
}

############# Correlation matrix ########
cor_matrix <- cor(all_mean_data[, c("RVI", "RVI4S1", "ratio", "Height", "lai", "PM", "SM")], use = "complete.obs", method = "spearman") # pearson, wenn normalverteilt
 
cor_matrix <- cor_matrix[4:7, 2:3] # Reduziere die Matrix auf die gewünschten Variablen
windows()

corrplot(cor_matrix, method = "circle",  tl.col = "black", tl.srt = 0, diag = TRUE, addCoef.col = "black", cl.pos = "n",mar = c(1, 1, 3, 1), tl.pos = "d")
        title("Korrelationsmatrix der gemittelten Variablen", line = 4)

cor_results <- cor.test(all_mean_data$RVI, all_mean_data$PM, method = "spearman")
cor_result <- cor.test(all_mean_data$ratio, all_mean_data$PM, method = "spearman")
cor_results <- cor.test(all_mean_data$RVI4S1, all_mean_data$PM, method = "spearman")
cor_result <- cor.test(all_mean_data$ratio, all_mean_data$PM, method = "spearman")
# multiple correlation sprengt den Rahmen
# mcor(all_data[, c("Height", "lai", "PM", "SM")], all_data$RVI)
# mcor(all_data[, c("lai", "PM")], all_data$RVI)
# mcor(all_data[, c("Height", "SM")], all_data$RVI)
# mcor(all_data[, c("Height", "PM", "SM")], all_data$RVI)
# mcor(all_data[, c("PM", "SM")], all_data$RVI)
# mcor(all_data[, c("Height", "lai")], all_data$RVI)
# mcor(all_data[, c("Height", "lai", "PM")], all_data$RVI)
# mcor(all_data[, c("Height", "PM")], all_data$RVI)
# mcor(all_data[1:9, c("Height", "SM")], all_data$RVI[1:9])

############## Linear Regression ##########
lm_model <- lm(Reflectance ~ Height_final + SM_final, data = all_data_interpol[all_data_interpol$Band == "RVI",])
summary(lm_model)
plot(lm_model)
lm_model <- lm(Reflectance ~ Height_final, data = all_data_interpol)

############### Plots ##########
windows()
# 2d scatterplot
ggplot(all_mean_data, aes(x = RVI, y = PM)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Scatterplot of Reflectance and Height",
       x = "RVI4S1",
       y = "Plat Moisture") +
  scale_color_discrete() +
  theme(legend.position = "right")
  
  ggplot(all_mean_data, aes(x = RVI, y = Height)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Scatterplot of RVI and Height",
       x = "RVI",
       y = "Height") 

# Korrelationsmatrix für Sebastians Daten
seb_cor <- read.xlsx("./data/SS_corr.xlsx")
seb_cor <- seb_cor[7:15,1:3]

# Daten umbenennen
seb_cor$X2 = c(
    "NDMI",
    "NDMI",
    "NDVI (S2)",
    "NDMI", 
    "NDVI (S2)",#5
    "NDVI (UAV)",   
    "NDMI", 
    "NDVI (S2)", 
    "NDVI (UAV)", 
    "Poro"    #10
  )
seb_cor$X1 = c(
    "NDVI (S2)",
    "NDVI (UAV)",
    "NDVI (UAV)",
    "Poro",
    "Poro",
    "Poro",
    "PM",
    "PM",
    "PM",
    "PM"
  )

names(seb_cor) <- c("Variable1", "Variable2", "Korrelationswert")
data <- seb_cor
# Leere Matrix mit allen Variablen als Zeilen und Spalten
alle_variablen <- unique(c(data$Variable1, data$Variable2))
mat <- matrix(NA, nrow=length(alle_variablen), ncol=length(alle_variablen),
              dimnames=list(alle_variablen, alle_variablen))

# Werte eintragen
for(i in 1:nrow(data)){
  v1 <- data$Variable1[i]
  v2 <- data$Variable2[i]
  val <- data$Korrelationswert[i]
  mat[v1, v2] <- val
  mat[v2, v1] <- val # Symmetrisch eintragen
}

# Diagonale auf 1 setzen + runden
diag(mat) <- 1
print(round(mat, 3))

# Matrix anzeigen
mat <- mat[2:4,-c(3:4)] 
corrplot(mat, method = "circle", tl.col = "black", tl.srt = 1, 
         addCoef.col = "black", cl.pos = "r", mar = c(1, 1, 3, 1))
