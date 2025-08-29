

# Comparación de nichos ---------------------------------------------------

# Comparación de modelos de elipsoides con Elipsenm

# Paquetes
library(ellipsenm)
library(tidyverse)
library(rgl)
#en caso de que no instale ellipsenm, correr lo siguiente:

#leemos los tres archivos .csv con las siguientes columnas: pca1, pca2, pca3, Long y Lat



Reg_pcadom<-read.csv("G:/Mi unidad/Bioclimatologia/Proyectos/Protonicho/analisis/pcadom.csv")
head(Reg_pcadom)

Reg_pcahis<-read.csv("G:/Mi unidad/Bioclimatologia/Proyectos/Protonicho/analisis/pcahis.csv")
head(Reg_pcahis)

Reg_pcaita<-read.csv("G:/Mi unidad/Bioclimatologia/Proyectos/Protonicho/analisis/pcaita.csv")
head(Reg_pcaita)

#Creamos los elipsoides

elip_dom<-overlap_object(data=Reg_pcadom, species="dom",longitude="Long", latitude = "Lat", method = "mve1", level = 99)
elip_his<-overlap_object(data=Reg_pcahis, species="his",longitude="Long", latitude = "Lat", method = "mve1", level = 99)
elip_ita<-overlap_object(data=Reg_pcaita, species="ita",longitude="Long", latitude = "Lat", method = "mve1", level = 99)

# sobreposición de las elipses de las tres especies 
overlap_t<-ellipsoid_overlap(elip_dom,elip_his, elip_ita)

# graficamos las tres elipses de las especies
plot_overlap(object=overlap_t, niches = c(1,3), data = T,  data_col = c("black","blue","red"),  background = T, change_labels = FALSE, xlab = "PCA1", ylab = "PCA2", zlab = "PCA3", legend = F, niche_col = c("black","blue","red"))

setwd("G:/Mi unidad/Bioclimatologia/Proyectos/Protonicho/analisis")
svg(filename ="elipses.svg")

rgl.postscript('elipses.pdf', fmt = 'pdf')

rgl.postscript("elipse.pdf",fmt="pdf")


##### promedios

bg_swd.list[[m]] |> 
  write.table("./borrar/Sel_call_back.txt", sep="\t", dec=".", row.names = F)

# 1. Leer shapefile de puntos
pts <- vect(bg_swd.list[[m]], geom=c("lon", "lat"))  # capa de puntos

# 2. Crear raster vacío con extensión y resolución deseada
# (ajusta res() al tamaño de celda que quieras)

r <- rast("./borrar/clima/capas_mean/1/prec_ame_5k_1_mean.tif")

r <- rast(ext(pts), res = 0.01, crs = crs(pts))

# 3. Rasterizar usando una columna de atributos (ej: "valor")
prec.bg<-rasterize(pts, r, field = "prec")
tmax.bg<-rasterize(pts, r, field = "tmax")
tmin.bg<-rasterize(pts, r, field = "tmin")
# 
#   writeRaster("./borrar/clima/rasterize_prec.tif", overwrite = TRUE)
# rasterize(pts, r, field = "tmax") |> 
#   writeRaster("./borrar/clima/rasterize_tmin.tif", overwrite = TRUE)
# rasterize(pts, r, field = "tmin") |> 
#   writeRaster("./borrar/clima/rasterize_tmax.tif", overwrite = TRUE)


# 4. Guardar el resultado (opcional)

# 5. Ver
plot(r_val)
points(pts, col = "red")


a<-rast("D:/CHELSA_ym_5k/1901-01/prec.tif")
b<-rast("./outputs/borrar/background.tif")
a

global(a, fun = "notNA", na.rm = TRUE)
global(b, fun = "notNA", na.rm = TRUE)
