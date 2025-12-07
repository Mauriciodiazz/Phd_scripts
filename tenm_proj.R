
# Projecting time-specific niche models -----------------------------------

library(tidyverse)
library(tidyterra)
library(terra)
library(tenm)

# Month -------------------------------------------------------------------

spp.list.m <- list.files('./species/mig_ENM/monthly/', full.names = T)
spp.name.m <- basename(spp.list.m) |> 
  str_remove('_ENMm')

# Cargar todas las capas del mes
month.code <- c(paste0('-0', 1:9), '-10', '-11', '-12')

 # y <- 78
for(y in c(78)){ # set chrys 72 - set palma 78
  
temp_env <- new.env()
load(paste0(spp.list.m[y],'/ENM_',spp.name.m[y],'.Rdata'), envir = temp_env)

mod_sel.list <- temp_env$mod_sel.list
r.stack.list <- list()

spp.ext <- data.frame(spp=NA, month=1:12, occ.out=NA, occ.tot=NA)

# m <- 1
for(m in 1:12){
  print(paste(spp.name.m[y],"month", m))
m.files <- list.files('D:/CHELSA_ym_5k/', pattern = month.code[m], full.names = T)

# Cargar puntos
pts <- mod_sel.list[[m]]$temporal_df

pts.vec <- vect(pts, geom = c('lon', 'lat'), crs='EPSG:4326')

# Convertir grados a metros primero si el CRS está en lat/lon

# Definición LCC recomendada para cubrir América
lambert_america <- "+proj=lcc +lat_1=-10 +lat_2=65 +lat_0=0 +lon_0=-100 +datum=WGS84 +units=m +no_defs"

pts.vec2 <- project(pts.vec, lambert_america)  # CRS UTM adecuado para México central, cambia según tu zona

# Crear buffer de 250,000 metros (250 km - 50 celdas de 5km)
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
  # # ------
  # pre.pts <- as.points(r.crop2[[1]]) # Pixeles de precipitacion (hay mas en crop)
  # a <- extract(r.crop2, pre.pts, ID=F, xy=T)
  # a[which(is.na(a$tmin)),]
  # 
  # # ------
  
  suit <- tenm::predict(
    object = mod_sel.list[[m]],
    model_variables = c("prec", "tmax", "tmin"),
    layers = r.crop2,
    layers_ext = ".tif$",
    output = "suitability")
  
  rast.list[[x]] <- ntbox::bin_model(
    model = raster::raster(suit), 
    pts[, c(2, 3)], 
    percent = 10)
}

r.stack <- rast.list |> 
  raster::stack() |> 
  rast() 


r.stack.list[[m]] <- r.stack
writeRaster(r.stack, paste0("./borrar/bin_mods/", spp.name.m[y], "_", m, ".tif"), overwrite=T)

}

# save(list = c("r.stack.list"),
#            file = paste0("./borrar/bin_mods/", spp.name.m[y], "_proj.Rdata"))

}


# Leer los raster y calcular el mejor mapa
r.proj.path <- list.files("./borrar/bin_mods/", pattern=".tif$", full.names=T)
spp.list.m
spp.name.m

spp.ext.list <- list()

y <- 78
# Cargar los puntos de presencia
r.proj.list <- r.proj.path[grep(spp.name.m[y], r.proj.path)]

temp_env <- new.env()
load(paste0(spp.list.m[y],'/ENM_',spp.name.m[y],'.Rdata'), envir = temp_env)

spp.ext2 <- data.frame(spp=NA, month=1:12, occ.in.tr=NA, occ.in.te=NA, occ.tot=NA, thr = NA)

abex.list <- temp_env$abex.list

# Umbrales
# 0.75*terra::nlyr(r.proj) #87
# 0.90*terra::nlyr(r.proj) #104
# 0.95*terra::nlyr(r.proj) #110

# thr <- 87
for (thr in c(104,110)) {
# m <- 1
for (m in 1:12) {
  r.proj <- rast(r.proj.list[m])

# Cargar los puntos de presencia
pts <- abex.list[[m]]$temporal_df

pts.vec <- vect(pts, geom = c('lon', 'lat'), crs='EPSG:4326')

mod_sel.list <- temp_env$mod_sel.list

r_final <- sum(r.proj) > thr
r_final <- classify(r_final, cbind(0, 0, 0))  # para asegurarte que sea 0/1

 writeRaster(r_final, paste0("./borrar/bin_mods/bin_sums/", spp.name.m[y], "_", thr, "_", m, ".tif"), overwrite=TRUE)

spp.ext2[m,1] <-  basename(r.proj.list[m]) |> str_sub(1,9)
spp.ext2[m,3] <-
  terra::extract(r_final, pts.vec, ID=F, bind=T) |>
  select(sum, trian_test) |>
  as_tibble() |> 
  filter(sum != 0) |>
  filter(trian_test=="Train") |> 
  nrow()
spp.ext2[m,4] <-
  terra::extract(r_final, pts.vec, ID=F, bind=T) |>
  select(sum, trian_test) |>
  as_tibble() |> 
  filter(sum != 0) |>
  filter(trian_test=="Test") |> 
  nrow()

spp.ext2[m,5] <- length(pts.vec)
spp.ext2[m,6] <- thr
}

spp.ext.list[[thr]] <- spp.ext2
}


do.call(rbind, spp.ext.list) |> 
  mutate(prop = occ.in.tr/occ.tot*100,
         thr = as.factor(thr),
         month = as.factor(month)) |> 
  summarise(mean=mean(prop),.by=thr)
  ggplot(aes(y=prop, x=month, fill=thr)) +
  geom_bar(stat='identity', position=position_dodge())

# plot(r_final)
# plot(pts.vec, add=T)

maps.path <- list.files("./borrar/bin_mods/bin_sums/", full.names = T)

x_104 <- maps.path[grepl("104", maps.path)]
x_110 <- maps.path[grepl("110", maps.path)]

col_104 <- c('#4a72b0','#4a72b0','#4a72b0',
             '#6ea96e','#6ea96e','#6ea96e',
             "#da9500","#da9500","#da9500",
             "#294029","#294029","#294029")

col_110 <- c('#2f486f','#2f486f','#2f486f',
             '#426642','#426642','#426642',
             "#664600","#664600","#664600",
             "#0e160e","#0e160e","#0e160e")

### Mapas

x <- 1

r_104 <- rast(x_104[x])
r_110 <- rast(x_110[x])

ggplot() +
  geom_spatraster(data = r_104, aes(fill = sum)) +
  scale_fill_gradient(
    low = NA,        # no pinta los ceros
    high = col_104[x],   # color para los unos
    limits = c(0, 1),
    na.value = NA) 
  
