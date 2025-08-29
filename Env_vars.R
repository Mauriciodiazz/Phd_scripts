## Cutting, saving and resampling variables

library(terra)
library(tidyverse)

setwd("D:/") # Setting disk work directory

#1. Cut variables to America extent
chelsa.files <- list.files("./CHELSA_1k/", full.names = F)

extent<-ext(-180,-10,-60,85)

for (x in 1:length(chelsa.files)) {
  print(chelsa.files[x])
  
  env.files <- 
    list.files(paste0("./CHELSA_1k/", chelsa.files[x]), full.names = T)
  env.names <- 
    list.files(paste0("./CHELSA_1k/", chelsa.files[x]), full.names = F) |> 
    str_remove("CHELSAcruts_") |> 
    str_remove("_V\\.1\\.0\\.tif$")
  
  # barra de progreso
  pb <- txtProgressBar(min = 0,      # Valor mínimo de la barra de progreso
                       max = length(env.files), # Valor máximo de la barra de progreso
                       style = 3,    # Estilo de la barra (también style = 1 y style = 2)
                       width = 50,   # Ancho de la barra. Por defecto: getOption("width")
                       char = "=")   # Carácter usado para crear la barra
  
  for (y in 1:length(env.files)) {
    env.rast<-rast(env.files[y])
    
    out<-crop(env.rast, 
              extent, 
              filename=paste0("./CHELSA_ame_1k/", chelsa.files[x], "_ame/", env.names[y], "_ame.tif"))
    
    setTxtProgressBar(pb, y)
  }
  close(pb) # Cerramos la conexión	
}


# Mean monthly layers -----------------------------------------------------

# Variable counter
chelsa.files.m <- list.files("./CHELSA_ame_5k/", full.names = F)[-1]
# Month counter
env.month <- 
  paste0("_", seq(1,12), "_")

#	list.files("./CHELSA_ame_1k/prec_ame/")
for (x in 1:length(chelsa.files.m)) {
  for (y in 1:length(env.month)) {
    # Open raster stack files
    env.files.m <- 
      list.files(paste0("./CHELSA_ame_5k/", chelsa.files.m[x]), pattern = env.month[y], full.names = T)
    env.rast.m<-rast(env.files.m)
    # Mean calculation
    out.m<-mean(env.rast.m)
    # Save output
    writeRaster(out.m, paste0("./CHELSA_ame_5k/mean_ame_5k/", chelsa.files.m[x], env.month[y], "mean.tif"))
    print(paste0(chelsa.files.m[x], env.month[y]," done"))
  }
}


# Remuestrear a 5km -------------------------------------------------------
# 1k raster folder direction
chelsa.files.m <- list.files("./CHELSA_ame_1k/", full.names = F)[-1]

# 5k raster folder direction
chelsa.files.m5 <- list.files("./CHELSA_ame_5k/", full.names = F)[-1]

# resample 5k resolution template
template<-rast("./CHELSA_ame_1k/tmin_ame_1k/tmin_1_1901_ame.tif")
res(template) <- 0.045


for (x in 1:length(chelsa.files.m)) {
  env.files.m <- 
    list.files(paste0("./CHELSA_ame_1k/", chelsa.files.m[x]), full.names = T)
  env.names.m <- 
    list.files(paste0("./CHELSA_ame_1k/", chelsa.files.m[x]), full.names = F) |> str_remove(".tif")
  
  # Progress bar (optional)
  pb <- txtProgressBar(min = 0,      # Valor mínimo de la barra de progreso
                       max = length(env.files.m), # Valor máximo de la barra de progreso
                       style = 3,    # Estilo de la barra (también style = 1 y style = 2)
                       width = 50,   # Ancho de la barra. Por defecto: getOption("width")
                       char = "=")   # Carácter usado para crear la barra
  
  for (y in 1:length(env.files.m)) {
    # Open raster stack files
    env.rast.m<-rast(env.files.m[y])
    # Resampling calculation
    out.m<-resample(env.rast.m, template, 
                    method="bilinear", 
                    filename=paste0("./CHELSA_ame_5k/", 
                                    chelsa.files.m5[x], "/", env.names.m[y], "_5k.tif"))
    # Show progress bar (optional)
    setTxtProgressBar(pb, y)
  }
  # Close progress bar conection
  close(pb) 
}


# To tenm niche time-specific models I need to structure the variables per year and month
# year -> year-month -> variables (prec.tif, tmax.tif, tmin.tif)
setwd("D:/") # Setting disk work directory

# Rutas de las carpetas originales
prec_folder <- "./CHELSA_ame_5k/prec_ame_5k"
tmax_folder <- "./CHELSA_ame_5k/tmax_ame_5k"
tmin_folder <- "./CHELSA_ame_5k/tmin_ame_5k"

# Ruta base para la nueva estructura
base_path <- "./CHELSA_ym_5k"

#1. creating year-month file

years <- seq(1901,2016)

month.orig <- 1:12
months <- sprintf("%02d", 1:12)


for (y in 1:length(years)) {
 for (m in 1:length(months)) {
    # Nombre de la carpeta con formato "YYYY-MM"
    period_folder <- file.path(base_path, paste0(years[y], "-", months[m]))
    dir_create(period_folder)
    
    # Definir rutas de los archivos originales
    prec_file <- file.path(prec_folder, paste0("prec_", month.orig[m], "_", years[y], "_ame_5k.tif"))
    tmax_file <- file.path(tmax_folder, paste0("tmax_", month.orig[m], "_", years[y], "_ame_5k.tif"))
    tmin_file <- file.path(tmin_folder, paste0("tmin_", month.orig[m], "_", years[y], "_ame_5k.tif"))
    
    # Copiar y renombrar los archivos en la nueva carpeta
    if (file.exists(prec_file)) {
      file_copy(prec_file, file.path(period_folder, "prec.tif"))
    }
    if (file.exists(tmax_file)) {
      file_copy(tmax_file, file.path(period_folder, "tmax.tif"))
    }
    if (file.exists(tmin_file)) {
      file_copy(tmin_file, file.path(period_folder, "tmin.tif"))
    }

  }
}


#######
library(terra)

f<-1 #files
dir.list<-list.files("D:/prueba/", full.names = F)

path<- abex$temporal_df$layers_path |> str_remove("D:/CHELSA_ym_5k/")

for (f in 1:length(path)) {
 pb <- txtProgressBar(min = 0,      # Valor mínimo de la barra de progreso
                       max = length(path), # Valor máximo de la barra de progreso
                       style = 3,    # Estilo de la barra (también style = 1 y style = 2)
                       width = 50,   # Ancho de la barra. Por defecto: getOption("width")
                       char = "=")   # Carácter usado para crear la barra
  for (r in 1:3) {
    raster <- rast(list.files(abex$temporal_df$layers_path[f], full.names =T, pattern=".tif$")[r])
    
    # NAflag(raster)
    # is.na(raster)
    raster[is.nan(raster)] <- -9999
    
    writeRaster(raster, list.files(paste0("D:/prueba/", path[f]), full.names =T, pattern=".tif$")[r], overwrite = T)
    
  }
    setTxtProgressBar(pb, f)
}

r1<-rast(list.files(paste0("D:/prueba/",dir.list[f]), full.names=T, pattern=".tif$")[1])
r2<-rast(list.files(paste0("D:/prueba/",dir.list[f]), full.names=T, pattern=".tif$")[2])
r1
r2
NAflag(r1)
NAflag(r2)

abex$temporal_df$layers_path |> unique()

plot(r1)
