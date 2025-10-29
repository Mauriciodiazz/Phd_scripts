
# Projecting time-specific niche models -----------------------------------

library(tidyverse)
library(tidyterra)
library(terra)
library(tenm)

# Lista de nombres de la especie
spp.list <- list.files('M:/Phd_PC/Script_Phd_PC/inputs/', full.names = T) 
spp.name <- basename(spp.list) |> 
  str_remove('_ENMm')

# Cargar todas las capas del mes
month.code <- c(paste0('-0', 1:9), '-10', '-11', '-12')

y <- 2
for(y in 2:3){
  
  temp_env <- new.env()
  load(paste0(spp.list[y],'/ENM_',spp.name[y],'.Rdata'), envir = temp_env)
  
  mod_sel.list <- temp_env$mod_sel.list
  
  # m <- 1
  for(m in 1:12){
    m.files <- list.files('D:/CHELSA_ym_5k/', pattern = month.code[m], full.names = T)
    
    # Cargar puntos
    pts <- mod_sel.list[[m]]$temporal_df
    
    pts.vec <- vect(pts, geom = c('lon', 'lat'), crs='EPSG:4326')
    
    # Convertir grados a metros primero si el CRS está en lat/lon
    
    # Definición LCC recomendada para cubrir América
    lambert_america <- "+proj=lcc +lat_1=-10 +lat_2=65 +lat_0=0 +lon_0=-100 +datum=WGS84 +units=m +no_defs"
    
    pts.vec2 <- project(pts.vec, lambert_america)  # CRS UTM adecuado para México central, cambia según tu zona
    
    # Crear buffer de 250,000 metros (250 km - 50 ceeldas de 5km)
    buffer_250km <- buffer(pts.vec2, width = 250000) |> 
      aggregate() |> 
      project('EPSG:4326')
    
    
    rast.list <- list()
    
    # Abrir raster por cada tiempo (month-year)
    for (x in 1:length(m.files)) { # 
      r <- rast(list.files(m.files[x], full.names=T))
      
      
      r.crop <- crop(r, buffer_250km, mask = T)
      # plot(r.crop)
      names(r.crop) <- c('prec', "tmax", "tmin")
      
      #Hay un problema con CHELSAcruts y hay algunos pocos pixeles de precipitacion que no estan correspondidos por pixeles de temperatura (hay mas pixeles de prec que de temp), esto genera que en la extraccion existan NAs. Para solucionar esto, toca hacer un mask de las variables con alguna de temperatura (min o max) para homogeneizarlo
      
      r.crop2 <- crop(r.crop, r.crop[[2]], mask = T)
      
      # Esto de abajo es una prueba de que funciono lo anterior
      # # a ver, si yo obtengo los centroides de esos pixeles y extraigo los valores, todos los pixeles deberian tener 3 valores. Veamos
      #
      # pre.pts <- as.points(r.crop[[1]])
      # a <- extract(r.crop, pre.pts, ID=F, xy=T)
      # a[which(is.na(a$tmin)),]
      
      
      suit <- tenm::predict(
        mod_sel.list[[m]],
        model_variables = c("prec", "tmax", "tmin"),
        layers = r.crop2,
        layers_ext = ".tif$",
        output = "suitability")
      
      rast.list[[x]] <- ntbox::bin_model(
        model = raster::raster(suit), 
        pts[, c(2, 3)], 
        percent = 10)
    }
    
    rast.list |> 
      raster::stack() |> 
      rast() |> 
      sum() |> 
      terra::writeRaster(paste0("./borrar/datos_phd/", spp.name[y], "_proj_", m, ".tif"), overwrite=TRUE)
    
  }
  
}
