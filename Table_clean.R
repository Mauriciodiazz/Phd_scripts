#Cleaning species
library(tidyverse)
library(ntbox)
library(sf)
library(terra)

# Species data
data<-read_tsv("./species/GBIF/0039607-240906103802322.csv", quote = "") |> 
  select(1, 4:10, 13,16:18, 22:25,30:33,36:38)
data |> 
  head()

spp<-levels(factor(data$species))
spp

# Rasters
bio1<-rast("inputs/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif")
bio12<-rast("inputs/CHELSA/CHELSA_bio12_1981-2010_V.2.1.tif")
ocean<-rast("inputs/CHELSA/CHELSAcruts_prec_12_2016_V.1.0.tif")
# Variables terra stack
chelsa<-c(bio1, bio12, ocean)
names(chelsa)<-c("CHELSA_B1","CHELSA_B12", "CHELSA_ocean")


#Esto es para evitar que la notación científica (8e-5) cree conflicto con la función decimalPlaces
options(scipen=10) 

decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

for(i in 1:length(spp)){
print(spp[i])
# 1. Seleccionar cada especie
spp.name<-spp[i] |> 
  str_split_1(pattern=" ") |> 
  substr(c(1,1), c(3,5)) |> 
  paste(collapse = '_')

data.sp<-data |> 
  filter(species==spp[i])
# 1.1 guardar en una tabla los registros crudos
print(paste(spp[i],"writing_raw"))
data.sp |> 
  write.table(paste0("./species/tablas_crudas/",spp.name,".txt"), 
              dec=".", sep="\t", row.names = F, quote = F)
print(paste(spp[i],"writed_raw"))
#2. Eliminar filas
#2.1 Filas Sin mes y sin cordenadas menor a dos decimales
mon.data<-data.sp[which(!is.na(data.sp$month)),]
# keep only the ones with more than 2 decimals
mon.data2<-mon.data[sapply(mon.data$decimalLatitude, decimalplaces) > 2 & 
            sapply(mon.data$decimalLongitude, decimalplaces) > 2, ]


#2.2 Eliminar filas con coordenadas duplicadas por mes
print(paste(spp[i],"month filtering"))
mon.list<-list()
for(y in 1:12){
    mon.list[[y]]<- mon.data2 |> 
      filter(month==y) |> 
      clean_dup(longitude= "decimalLongitude" , 
            latitude = "decimalLatitude", 
            threshold = 0.008333) #0.008333=1km
}
print(paste(spp[i],"month filtered"))

mon.list.sd<-as.data.frame(do.call(rbind, mon.list))
mon.list.sd |> 
  write.table(paste0("./species/tablas_sindups/",spp.name,"_sd.txt"), 
              dec=".", sep="\t", row.names = F, quote = F)

#3. extraer variables ambientales
print(paste(spp[i],"variables extraction"))
spp.shp<-mon.list.sd |> 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),crs=st_crs(4326))

#3.1 Extraer datos ambientales: bio1 - bio12 - ocean (no terrestres)
spp.chel<-extract(chelsa, spp.shp, ID=F, bind=T) |> 
  st_as_sf()

spp.chel[which(spp.chel$CHELSA_ocean!=-32768),] |> 
  write_sf(paste0("./species/mig_shapes/",spp.name, ".shp"))
print(paste(spp[i],"variables extracted"))

# 4. Plotting length of monthly data
spp.chel |> 
  as.data.frame() |> 
  summarise(n=n(), .by=month) |> 
  ggplot(aes(x=factor(month, level = c(1:12)), y=n))+
    geom_bar(stat='identity')+
  labs(title=spp[i], y="numero de registros")+
  theme(axis.title.x = element_blank())

ggsave(filename = paste0("./species/monthly_graphs/",spp.name,".png"),
       width = 10,
       height = 10, #alto
       scale=2,
       units ="cm",
       dpi = 200)
}

