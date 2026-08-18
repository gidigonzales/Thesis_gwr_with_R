library(GWmodel)
library(sp)
library(sf)
library(terra)

# ---------------------------------------------------------
# 1. PARAMETER FÜR DIE SCHLEIFEN DEFINIEREN
# ---------------------------------------------------------
kernel_typen <- c("gaussian", "bisquare")
ziel_variablen <- c("Temp", "Temp_Diff")

# WICHTIG: Die Namen müssen exakt zu den Dateinamen passen!
tageszeiten <- c("morgen", "mittag", "abend")

# Ordnerstruktur für Dateien anlegen
# Wir mappen den Tageszeiten-String auf die echten Dateinamen
datei_mapping <- list(
  "morgen" = "Data/heatdata_cleaned_20260630_morgen.csv",
  "mittag" = "Data/heatdata_mittag_traffic_wd.csv",
  "abend"  = "Data/heatdata_cleaned_20260630_abend.csv"
)

# ---------------------------------------------------------
# 2. EINMALIGES VORBEREITEN DER UMWELT-RASTER 
# ---------------------------------------------------------
print("Bereite Hintergrund-Raster vor...")

ndvi <- rast("Data/gi_ndvi_utm32n.tif")
svf <- rast("Data/svf_gi_200m.tif")
impervious <- rast("Data/impervious_gi.tif")
street_dens <- rast("Data/street_dens_50m.sdat")
traffic <- st_read("Data/avg_traffic.shp")
traffic <- st_transform(traffic, crs = 25832)
traffic_terra <- vect(traffic)

stadtgrenze <- st_read("Data/stadtgrenze_gi_manuell.shp")
stadtgrenze <- st_transform(stadtgrenze, crs = 25832) 
grenze_terra <- vect(stadtgrenze)

ndvi_crop <- crop(ndvi, grenze_terra)
ndvi_mask <- mask(ndvi_crop, grenze_terra)

svf_angepasst <- resample(svf, ndvi_mask, method="bilinear")
impervious_angepasst <- resample(impervious, ndvi_mask, method="bilinear")
street_dens_angepasst <- resample(street_dens, ndvi_mask, method="bilinear")

# Traffic rastern (für den Fall, dass es gebraucht wird)
traffic_raster <- rasterize(traffic_terra, ndvi_mask, field = "relative_s", background = 1)

# Einmalig das Prediction-Grid bauen
basis_stapel <- c(ndvi_mask, svf_angepasst, impervious_angepasst, street_dens_angepasst, traffic_raster)
names(basis_stapel) <- c("ndvi", "svf", "impervious", "street_dens", "traffic_speed")
# Aggregation hier einstellen, falls gewünscht (z.B. fact=1 oder fact=2)
basis_stapel <- aggregate(basis_stapel, fact= 1) 

grid_df <- as.data.frame(basis_stapel, xy = TRUE, na.rm = TRUE)
grid_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = 25832)
grid_sp <- as(grid_sf, "Spatial")
anzahl_zeilen <- nrow(grid_df)
chunk_size <- 4000 

# Leere Liste für die fertigen Raster (für den finalen Plot)
alle_fertigen_raster <- list()
# ---------------------------------------------------------
# 3. DIE SCHLEIFEN & PLOTTING
# ---------------------------------------------------------

for (tageszeit in tageszeiten) {
  
  # NEU: Liste für die Raster NUR für die aktuelle Tageszeit leeren
  tageszeit_raster_liste <- list()
  
  # Traffic Logik: Traffic nur mittags nutzen
  istraffic <- ifelse(tageszeit == "mittag", TRUE, FALSE)
  
  # Aktuelle CSV laden
  csv_pfad <- datei_mapping[[tageszeit]]
  print(paste("----- LADE DATEN FÜR:", toupper(tageszeit), "-----"))
  messdaten <- read.csv(csv_pfad, sep=";")
  
  punkte_sf <- st_as_sf(messdaten, coords = c("Lon", "Lat"), crs = 4326)
  punkte_sf <- st_transform(punkte_sf, crs = 25832)
  
  # Rasterwerte extrahieren
  raster_werte <- extract(basis_stapel, punkte_sf)
  punkte_sf$ndvi <- raster_werte$ndvi
  punkte_sf$svf <- raster_werte$svf
  punkte_sf$impervious <- raster_werte$impervious
  punkte_sf$street_dens <- raster_werte$street_dens
  
  punkte_sf <- punkte_sf[!is.na(punkte_sf$ndvi), ]
  rownames(punkte_sf) <- NULL
  punkte_sf <- st_cast(punkte_sf, "POINT")
  punkte_sp <- as(punkte_sf, "Spatial")
  
  for (temp_col in ziel_variablen) {
    
    # Checken, ob die Spalte in der CSV existiert
    if (!(temp_col %in% names(messdaten))) {
      print(paste("Spalte", temp_col, "existiert nicht in", tageszeit, "-> Überspringe."))
      next
    }
    
    # Gleichung bauen
    if (istraffic) {
      gleichung_formel <- paste(temp_col, "~ ndvi + svf + impervious + street_dens + traffic_speed")
    } else {
      gleichung_formel <- paste(temp_col, "~ ndvi + svf + impervious + street_dens")
    }
    gleichung <- as.formula(gleichung_formel)
    
    for (kernel_type in kernel_typen) {
      
      print(paste("--- Rechne:", tageszeit, "|", temp_col, "|", kernel_type, "---"))
      
      bw_opt <- bw.gwr(gleichung, 
                       data = punkte_sp, 
                       approach = "AICc", 
                       kernel = kernel_type, 
                       adaptive = TRUE)
      
      alle_vorhersagen <- numeric(anzahl_zeilen)
      
      for (i in seq(1, anzahl_zeilen, by = chunk_size)) {
        start_idx <- i
        end_idx <- min(i + chunk_size - 1, anzahl_zeilen)
        
        chunk_sf <- st_as_sf(grid_df[start_idx:end_idx, ], coords = c("x", "y"), crs = 25832)
        chunk_sp <- as(chunk_sf, "Spatial")
        
        chunk_pred <- gwr.predict(gleichung, 
                                  data = punkte_sp, 
                                  predictdata = chunk_sp, 
                                  bw = bw_opt, 
                                  kernel = kernel_type, 
                                  adaptive = TRUE)
        
        alle_vorhersagen[start_idx:end_idx] <- chunk_pred$SDF$prediction
        gc()
      }
      
      grid_df_lokal <- grid_df
      grid_df_lokal$Temp_pred <- alle_vorhersagen
      
      aktuelles_raster <- rast(grid_df_lokal[, c("x", "y", "Temp_pred")], type = "xyz", crs = "EPSG:25832")
      
      # Name für den Plot
      plot_titel <- paste(temp_col, kernel_type, sep="_")
      names(aktuelles_raster) <- plot_titel
      
      # In die Liste der AKTUELLEN Tageszeit speichern
      tageszeit_raster_liste[[plot_titel]] <- aktuelles_raster
      
      # TIF abspeichern
      writeRaster(aktuelles_raster, paste0("heatmaps/", tageszeit, "_", plot_titel, ".tif"), overwrite = TRUE)
      
    } # Ende Kernel-Loop
  } # Ende Variablen-Loop
  
  # ---------------------------------------------------------
  # NEU: PLOT FÜR DIE JEWEILIGE TAGESZEIT ERSTELLEN
  # ---------------------------------------------------------
  if (length(tageszeit_raster_liste) > 0) {
    print(paste("Erstelle Plot für", tageszeit, "..."))
    
    tageszeit_stapel <- rast(tageszeit_raster_liste)
    
    # PNG Datei anlegen (z.B. GWR_Uebersicht_morgen.png)
    png(paste0("heatmaps/GWR_Uebersicht_", tageszeit, ".png"), width=1200, height=1000, res=150)
    
    par(oma = c(0, 0, 3, 0))
    # Die Karten für diese Tageszeit im 2x2 Grid plotten
    plot(tageszeit_stapel, 
         
         col = hcl.colors(100, "heat", rev = TRUE),
         main = names(tageszeit_stapel), 
         mar = c(2, 2, 2, 2),
         axes = FALSE)
         
    titel_text <- paste("GWR Vorhersagen - Tageszeit:", toupper(tageszeit))
    mtext(titel_text, side = 3, outer = TRUE, cex = 1.5, font = 2)
    # Plot speichern und schließen
    dev.off()
  }
  
} # Ende Tageszeiten-Loop

print("Alle Durchläufe und Plots sind fertig!")