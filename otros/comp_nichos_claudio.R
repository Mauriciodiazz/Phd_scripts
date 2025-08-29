#comparación de modelos de elipsoides entre Passer domesticus, Passer hispanicus y Passer Italiae

# Paquetes
library(ellipsenm)
library(dplyr)
library(rgl)
#en caso de que no instale ellipsenm, correr lo siguiente:

#if(!require(devtools)){
#  install.packages("devtools")
#}
#if(!require(ellipsenm)){
# devtools::install_github("marlonecobos/ellipsenm")
#}

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
