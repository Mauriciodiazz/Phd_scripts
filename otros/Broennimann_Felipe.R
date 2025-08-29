#Curso MNE 2023
#Comparación de nichos con Broennimann et al. 2012, Di Cola et al. 2017 
#comparaci?n entre Pseudomyrmex ferrugineus (hormiga) y Acacia cornigera


setwd("D:/Proyectos/kuenm_MNE2022_INECOL_UDEA/Comp_nichos/")

#instalar las sigueinte paquetería
# install.packages("ecospat")
# install.packages("raster")
# install.packages("dismo")

library(ecospat)
library(raster)
library(dismo)

#Crear una carpeta que contenga una archivo .csv con los registros de la especie 1 (tres columnas 1 especie, 2 longitud 3 altitud), registros con la sp 2, subcarpeta con las variables de las sp1, subcarpeta con las variables de la sp 2 (las subcarpetas deberán contener el mismo número de variables y se deber?n llamar igual).

#registros sp 1 (sp, x, y) 
rega<- read.csv("Comp/reg/P_ferrugineus.csv", header=T)
head(rega)
#registros de la sp 2
regb<- read.csv("Comp/reg/V_cornigera.csv", header=T)


#perfil clim?tico de la entidad 1

#enlistar las variables clim?ticas del tama?o de la m en formato .asc o .tif
clima<- list.files(path ="Comp/varfer/", pattern = ".tif",full.names = TRUE)
#compiliar todas las variables
varclima<- stack(clima)
#convertir el raster (pixeles) a puntos
climapts <- rasterToPoints(varclima[[1]], fun=NULL, spatial=TRUE)
#extraer los valores de cada una de las variables 
climada<- extract (varclima,climapts)
#formateamos el archivo anterior para crear una tabla con columna x y y
climafin<- data.frame(coordinates(climapts),climada)
head(climafin)


#perfil clim?tico de la entidad 2

climab<- list.files(path="Comp/varcor/",pattern = ".tif",full.names = TRUE)
varclimab<- stack(climab)
climabpts <- rasterToPoints(varclimab[[1]], fun=NULL, spatial=TRUE)
climadb<- extract (varclimab,climabpts)
climabfin<- data.frame(coordinates(climabpts),climadb)



#Crear un perfil clim?tico de ambas especies 
ambosclimas<- na.omit(rbind(climafin,climabfin))
row.names(ambosclimas) <- 1:nrow(ambosclimas)
ambosclimas <- ambosclimas[]
head(ambosclimas)
tail(ambosclimas)

#A?adir los valores de cada una de las variables a cada uno de los puntos de registro de las especies (del principio)

sp1<-na.exclude(ecospat.sample.envar(dfsp =rega, colspxy = 2:3, colspkept = 2:3, dfvar = ambosclimas,colvarxy = 1:2,colvar = "all", resolution = 0.1))
sp2<-na.exclude(ecospat.sample.envar(dfsp =regb, colspxy = 2:3, colspkept = 2:3, dfvar = ambosclimas,colvarxy = 1:2,colvar = "all", resolution = 0.1))

####################PCA-DEL_AMBIENTE###################

#Generar una base de datos que incluya los valores de las variables de clima (que incluye las m de ambas especies generado previamente) mas los valores de la sp1 y luego la sp2

data <-rbind.data.frame(ambosclimas[,3:7],sp1[,3:7],sp2[,3:7])
w <-c(rep(1,nrow(ambosclimas)), rep(0,nrow(sp1)), rep(0,nrow(sp2)))

#si marca un eror en este punto, es posible que se deba instalar la sig. paqueter?a
#install.packages("ade4")

library(ade4)
pca.cal<-dudi.pca(data, row.w = w, center = T, scale=T, scannf = F, nf = 2)

#Filas en las que est?n los datos de clima  de cada una de las especies

row.clima<-1:nrow(climafin)
row.climab<-(nrow(climafin)+1):(nrow(climafin)+nrow(climabfin))
row.climaab <-1:nrow(ambosclimas)
row.sp1a<-(1+nrow(climafin)+nrow(climabfin)):(nrow(climafin)+nrow(climabfin)+nrow(sp1))
row.sp2b<-(1+nrow(climafin)+nrow(climabfin)+nrow(sp1)):(nrow(climafin)+nrow(climabfin)+nrow(sp1)+nrow (sp2))

#coordenadas en cada uno de los ejes del PCA de todos los datos(emes y especies)

scores.clima<-pca.cal$li[row.clima,]
scores.climab<-pca.cal$li[row.climab,]
scores.climaab<-pca.cal$li[row.climaab,]
scores.sp1a<-pca.cal$li[row.sp1a,]
scores.sp2b<-pca.cal$li[row.sp2b,]

#contribuci?n de las varaibles a cada componente

contribucion<- ecospat.plot.contrib(contrib=pca.cal$co, eigen = pca.cal$eig)


################################################################################
######Superficie de densidad de registros y ?rea de cada especie################

z1<-ecospat.grid.clim.dyn(scores.climaab, scores.clima, scores.sp1a, R=100)
z2<-ecospat.grid.clim.dyn(scores.climaab, scores.climab,scores.sp2b, R=100)

#gr?ficos de PCA 
ecospat.plot.niche (z1, title="fer", name.axis1="PCA1", name.axis2="PCA2", cor=F)
ecospat.plot.niche (z2, title="cor", name.axis1="PCA1", name.axis2="PCA2", cor=F)

#Obtenci?n del ?ndice de similitud de Schoener 1970 Valor observado "D":
ecospat.niche.overlap (z1=z1,z2=z2,cor=TRUE)

ecospat.plot.niche.dyn(z1, z2, quant=0.25, interest=1, title= "Niche overlap fer vs cor", name.axis1="PC1", name.axis2="PC2")


########################################################################
#################### TEST DE EQUIVALENCIA DE NICHO ######################
########################################################################

#a<-ecospat.niche.equivalency.test(z1,z2,rep=100)

#ecospat.plot.overlap.test(a,"D","Equivalencia barVSmac")

########################################################################
#################### TEST DE SIMILARIDAD DE NICHO Di Cola et al 2017####
########################################################################

#corremos el Similarity test, en este caso higher asumimos conservadurismo de nicho

b<-ecospat.niche.similarity.test(z1,z2,rep=100, overlap.alternative = "higher")

ecospat.plot.overlap.test(b,"D","fer a cor")


#################### An?lisis con base en di Cola et al 2017 lower ######################
#en este caso estamos probando divergencia de nicho 
#sim.test<-ecospat.niche.similarity.test(z1,z2, rep=100, alternative = "lower")
# Plot Similarity test



