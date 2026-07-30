library(GWmodel)
library(sp)
library(sf)
library(terra)

kernel_type <- "gaussian"
istraffic <- FALSE
temp_col <- "Temp" # Temp_Diff , Temp
tageszeit <- "morgen"
# heatdata_mittag_traffic_wd, heatdata_cleaned_20260630_abend, heatdata_cleaned_20260630_morgen
messdaten <- read.csv("Data/heatdata_cleaned_20260630_morgen.csv", sep=";")
punkte_sf <- st_as_sf(messdaten, coords = c("Lon", "Lat"), crs = 4326)
punkte_sf <- st_transform(punkte_sf, crs = 25832)


ndvi <- rast("Data/gi_ndvi_utm32n.tif")
svf <- rast("Data/svf_gi_200m.tif")
impervious <- rast("Data/impervious_gi.tif")
street_dens <- rast("Data/street_dens_50m.sdat")
traffic <- st_read("Data/avg_traffic.shp")
traffic <- st_transform(traffic, crs = 25832)

traffic_terra <- vect(traffic)

# raster nur auf die stadtgrenze zuschneiden
stadtgrenze <- st_read("Data/stadtgrenze_gi_manuell.shp")
stadtgrenze <- st_transform(stadtgrenze, crs = 25832) 
grenze_terra <- vect(stadtgrenze)

ndvi_crop <- crop(ndvi, grenze_terra)
ndvi_mask <- mask(ndvi_crop, grenze_terra)

# alle Raster auf die gleiche Zellgröße bringen
svf_angepasst <- resample(svf, ndvi_mask, method="bilinear")
impervious_angepasst <- resample(impervious, ndvi_mask, method="bilinear")
street_dens_angepasst <- resample(street_dens, ndvi_mask, method="bilinear")

traffic_raster <- rasterize(traffic_terra, ndvi_mask, field = "relative_s", background = 1)


#in einen Stack packen
if (istraffic) {
  raster_stapel <- c(ndvi_mask, svf_angepasst, impervious_angepasst, street_dens_angepasst, traffic_raster)
  names(raster_stapel) <- c("ndvi", "svf", "impervious", "street_dens", "traffic_speed")
} else {
  raster_stapel <- c(ndvi_mask, svf_angepasst, impervious_angepasst, street_dens_angepasst)
  names(raster_stapel) <- c("ndvi", "svf", "impervious", "street_dens")
}

# an den messpunkten die Rasterdaten extrahieren
raster_werte <- extract(raster_stapel, punkte_sf)

punkte_sf$ndvi <- raster_werte$ndvi
punkte_sf$svf <- raster_werte$svf
punkte_sf$impervious <- raster_werte$impervious
punkte_sf$street_dens <- raster_werte$street_dens

punkte_sf <- punkte_sf[!is.na(punkte_sf$ndvi), ]
rownames(punkte_sf) <- NULL

punkte_sf <- st_cast(punkte_sf, "POINT")

#punkte_sf <- na.omit(punkte_sf)

# in sp type umwandel, da damit gwr arbeiten kann
punkte_sp <- as(punkte_sf, "Spatial")

raster_stapel <- aggregate(raster_stapel, fact= 1)
grid_df <- as.data.frame(raster_stapel, xy = TRUE, na.rm = TRUE)
grid_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = 25832)
grid_sp <- as(grid_sf, "Spatial")

# Die mathematische Gleichung definieren
if (istraffic) {
  gleichung <- get(temp_col) ~ ndvi + svf + impervious + street_dens + traffic_speed
} else {
  gleichung <- get(temp_col) ~ ndvi + svf + impervious + street_dens
}

print("Suche optimale Bandbreite...")
# bw.gwr sucht die optimale Anzahl an Nachbarn
bw_opt <- bw.gwr(gleichung, 
                 data = punkte_sp, 
                 approach = "AICc", 
                 kernel = kernel_type, 
                 adaptive = TRUE)

print(paste("Die optimale Bandbreite (Anzahl Nachbarn) ist:", bw_opt))

# 1. Wir nutzen dein volles, hochauflösendes Grid (grid_df)
# Angenommen, das hat jetzt z.B. 1 Million Zeilen.

# 2. Wie viele Pixel sollen pro Durchgang berechnet werden?
chunk_size <- 3000
anzahl_zeilen <- nrow(grid_df)

# Eine leere Liste vorbereiten, in die wir die fertigen Temperaturen packen
alle_vorhersagen <- numeric(anzahl_zeilen)

print(paste("Starte Vorhersage in Blöcken. Insgesamt", ceiling(anzahl_zeilen/chunk_size), "Blöcke..."))

# 3. Die Schleife: Wir gehen das Grid Block für Block durch
for (i in seq(1, anzahl_zeilen, by = chunk_size)) {
  
  # Definiere Start- und Endpunkt für diesen Block
  start_idx <- i
  end_idx <- min(i + chunk_size - 1, anzahl_zeilen)
  
  print(paste("Berechne Block von", start_idx, "bis", end_idx))
  
  # Schneide genau diesen Block aus dem großen Grid heraus
  grid_chunk_df <- grid_df[start_idx:end_idx, ]
  
  # Wandle nur diesen kleinen Block in ein Spatial-Objekt um
  chunk_sf <- st_as_sf(grid_chunk_df, coords = c("x", "y"), crs = 25832)
  chunk_sp <- as(chunk_sf, "Spatial")
  
  # Führe GWR nur für diesen kleinen Block aus
  chunk_pred <- gwr.predict(gleichung, 
                            data = punkte_sp, 
                            predictdata = chunk_sp, 
                            bw = bw_opt, 
                            kernel = kernel_type, 
                            adaptive = TRUE)
  
  # Speichere die vorhergesagten Temperaturen in unsere lange Liste
  alle_vorhersagen[start_idx:end_idx] <- chunk_pred$SDF$prediction
  
  # Zwinge R, den Arbeitsspeicher für den nächsten Block aufzuräumen (WICHTIG!)
  gc()
}

print("Vorhersage abgeschlossen! Füge Raster zusammen...")

# 4. Die fertigen Ergebnisse wieder an das große Grid anheften
grid_df$Temp_pred <- alle_vorhersagen

# 5. Wie gehabt in ein Raster umwandeln und speichern
heatmap_raster <- rast(grid_df[, c("x", "y", "Temp_pred")], type = "xyz", crs = "EPSG:25832")
plot(heatmap_raster, col = hcl.colors(100, "Inferno"), main = paste("GWR",kernel_type, temp_col, "Heatmap (°C)"))
writeRaster(heatmap_raster, paste0("heatmaps/",trimws(temp_col),"_Gießen_",trimws(tageszeit),"_10m_",trimws(kernel_type),".tif"), overwrite = TRUE)
