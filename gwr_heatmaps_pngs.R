library(GWmodel)
library(sp)
library(sf)
library(terra)

# Schleifenparameter
kernel_typen <- c("bisquare") 
ziel_variablen <- c("Temp", "Temp_Diff")

# Die Namen müssen exakt zu den Dateinamen passen!
tageszeiten <- c("morgen", "mittag", "abend")


datei_mapping <- list(
  "morgen" = "Data/heatdata_cleaned_20260630_morgen.csv",
  "mittag" = "Data/heatdata_mittag_traffic_wd.csv",
  "abend"  = "Data/heatdata_cleaned_20260630_abend.csv"
)

if(!dir.exists("heatmaps")) dir.create("heatmaps")
# vorbereiten der Raster
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

# Traffic rastern
traffic_raster <- rasterize(traffic_terra, ndvi_mask, field = "relative_s", background = 1)

# prediction grid erstellen
basis_stapel <- c(ndvi_mask, svf_angepasst, impervious_angepasst, street_dens_angepasst, traffic_raster)
names(basis_stapel) <- c("ndvi", "svf", "impervious", "street_dens", "traffic_speed")

basis_stapel <- aggregate(basis_stapel, fact= 1) 

grid_df <- as.data.frame(basis_stapel, xy = TRUE, na.rm = TRUE)
grid_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = 25832)
grid_sp <- as(grid_sf, "Spatial")
anzahl_zeilen <- nrow(grid_df)
chunk_size <- 4000 

alle_fertigen_raster <- list()

# berechnung gwr

for (tageszeit in tageszeiten) {
  
  istraffic <- ifelse(tageszeit == "mittag", TRUE, FALSE)
  
  csv_pfad <- datei_mapping[[tageszeit]]
  print(paste("----- LADE DATEN FÜR:", toupper(tageszeit), "-----"))
  messdaten <- read.csv(csv_pfad, sep=";")
  
  punkte_sf <- st_as_sf(messdaten, coords = c("Lon", "Lat"), crs = 4326)
  punkte_sf <- st_transform(punkte_sf, crs = 25832)
  
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
      
      # plot name
      plot_titel <- paste(tools::toTitleCase(tageszeit), temp_col, sep=" - ")
      names(aktuelles_raster) <- plot_titel
      
      # liste für finalen plot
      listen_key <- paste(tageszeit, temp_col, sep="_")
      alle_fertigen_raster[[listen_key]] <- aktuelles_raster
      
      writeRaster(aktuelles_raster, paste0("heatmaps/", tageszeit, "_", temp_col, "_", kernel_type, ".tif"), overwrite = TRUE)
      
    } # Ende Kernel-Loop
  } # Ende Variablen-Loop
} 

# plotting
if (length(alle_fertigen_raster) > 0) {
  print("Erstelle finalen Gesamt-Plot...")
  
  finaler_stapel <- rast(alle_fertigen_raster)
  
  png("heatmaps/GWR_Tagesverlauf_Uebersicht.png", width=1600, height=2200, res=200)
  
  par(oma = c(0, 0, 4, 0)) # Platz für Haupttitel oben 
  
  # nr = 3 (Zeilen) und nc = 2 (Spalten) 
  plot(finaler_stapel, 
       nr = 3, nc = 2,
       col = hcl.colors(100, "heat", rev = TRUE),
       main = names(finaler_stapel), 
       mar = c(2, 2, 3, 4), #  Platz rechts für die Legende lassen
       axes = FALSE)
  
  titel_text <- paste("GWR Vorhersagen ( Kernel:", kernel_typen[1], ")")
  mtext(titel_text, side = 3, outer = TRUE, line = 2, cex = 1.8, font = 2)
  
  dev.off()
  print("Plot gespeichert unter: heatmaps/GWR_Tagesverlauf_Uebersicht.png")
}

print("Alle Durchläufe und Plots sind fertig!")