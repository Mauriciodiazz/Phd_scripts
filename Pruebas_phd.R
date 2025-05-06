## Analisis de completitud 

# Probando iNEXT
#install.packages("iNEXT")
library(iNEXT)
library(tidyverse)

# Intentemoslo para 20 especies
spp.path<-list.files("./species/tablas_sindups/", full.names = T)

rows<-sample(1:length(spp.path), 50, replace = FALSE)

spp.list<-list()

for(x in 1:length(rows)){
  spp.list[[x]]<- read_tsv(paste0(spp.path[rows[x]])) |> 
    select(species, decimalLatitude, decimalLongitude, month, year)
}

spp_df <- bind_rows(spp.list)

spp_df |> 
  na.omit() |> 
  filter(month==1 & decimalLatitude<30) |> 
  mutate(rangos = cut(year, breaks = seq(1800, 2030, by = 10), include.lowest = TRUE),
         franjas= case_when(
           decimalLatitude > 25 &  decimalLatitude < 30 ~ "F1",
           decimalLatitude > 20 &  decimalLatitude < 25 ~ "F2",
           decimalLatitude > 15 &  decimalLatitude < 20 ~ "F3",
           decimalLatitude > 10 &  decimalLatitude < 15 ~ "F4",
           decimalLatitude > 5 &  decimalLatitude < 10 ~ "F5",
           decimalLatitude > 0 &  decimalLatitude < 5 ~ "F6",
           decimalLatitude > -5 &  decimalLatitude < 0 ~ "F7",
           decimalLatitude > -10 &  decimalLatitude < -5 ~ "F8",
           decimalLatitude > -15 &  decimalLatitude < -10 ~ "F9",
           decimalLatitude > -20 &  decimalLatitude < -15 ~ "F10",
           decimalLatitude > -25 &  decimalLatitude < -20 ~ "F11",
           decimalLatitude > -30 &  decimalLatitude < -25 ~ "F12",
           decimalLatitude > -35 &  decimalLatitude < -30 ~ "F13",
           decimalLatitude > -40 &  decimalLatitude < -35 ~ "F14",
           decimalLatitude > -45 &  decimalLatitude < -40 ~ "F15",
           decimalLatitude > -50 &  decimalLatitude < -45 ~ "F16",
           decimalLatitude > -55 &  decimalLatitude < -50 ~ "F17",
           decimalLatitude > -60 &  decimalLatitude < -55 ~ "F18",
           decimalLatitude > -65 &  decimalLatitude < -60 ~ "F19")) |>
  summarise(n = n(), .by = c(rangos, franjas)) |>
  pivot_wider(names_from = franjas, values_from = n) |>
  mutate(across(everything(), ~ replace_na(.x, 0))) |> 
  column_to_rownames(var = "rangos") |> 
  DataInfo()


data[is.na(data$franjas),]

# ahora el mapa
library(sf)
library(tidyverse)
america<-read_sf("./inputs/america_shapes/america_cont/America.shp")

lats<- seq(-65, 30, by= 5)

# Crear una lista de polígonos rectangulares
franjas <- lapply(lats, function(lat) {
  st_polygon(list(matrix(c(-170, lat,   # Punto inferior izquierdo
                           -10, lat,   # Punto inferior derecho
                           -10, lat + 5, # Punto superior derecho + incremento latitud
                           -170, lat + 5, # Punto superior izquierdo + incremento latitud
                           -170, lat), # Cierre del polígono
                         ncol = 2, byrow = TRUE)))
})


# Convertir la lista en un objeto sf
franjas_sf <- st_sf(geometry = st_sfc(franjas, crs = 4326)) |> 
  mutate(lat_band = paste(lats, lats + 5, sep = " to "),
         franja=paste0("F", 20:1))

franjas_sf2<-franjas_sf[-20,]

franjas_sf2 |> 
  select(lat_band) |> 
  as_tibble()
#

franjas_sf2 |> 
  ggplot() +
  geom_sf() +
  geom_sf(data=america, aes(geometry=geometry, alpha=0))


# análisis de completitud usando KnowBR package 
library(KnowBR)
library(sp)
library(sf)


data(States)

States |> class()
data<-read.table("./outputs/borrar/10grados/centroids10.txt", header=T, dec=".", sep="\t")[,-1]
shape<-read_sf("./outputs/borrar/10grados/10_grid.shp")
a<-as(shape, "Spatial")
data(adworld)

data2<-
  data |> 
  rename("Longitude"=x, "Latitude"=y) |> 
  pivot_longer(cols = !c(Longitude, Latitude), names_to = "Species", values_to = "count") |> 
  filter(count!=0) |>
  sample_n(size=500) |> 
  relocate(Species) |> 
  as.data.frame()

KnowBPolygon(data=data2, shape=States, admAreas=T, shapenames="NAME", minLon=-170,
             maxLon=-30, minLat=-100, maxLat=100, colscale=rev(heat.colors(100)), jpg=FALSE,
             curve="Rational", save="CSV", 
             file1="Species per site",
             file2="Estimators", 
             file3="Standard error of the estimators", 
             na="NA",
             dec=".", row.names=FALSE)

KnowB(data=data2, format="A", cell=60, minLon=-170,
             maxLon=-30, minLat=-100, maxLat=100, colscale=rev(heat.colors(100)), jpg=FALSE,
             curve="Rational", save="CSV", 
             file1="Species per site",
             file2="Estimators", 
             file3="Standard error of the estimators", 
             na="NA",
             dec=".", row.names=FALSE)


#Download records from GBIF of the flowering plants of the family Polygonaceae

library(rgbif)
r.gbif<-occ_search(scientificName = "Polygonaceae", limit=500, return='data',
                    hasCoordinate=TRUE)

#Data frame with the format A required by the function KnowBPolygon

records<-data.frame(r.gbif$data$species, r.gbif$data$decimalLongitude, r.gbif$data$decimalLatitude)
names(records)<-c("Species","Longitude","Latitude")

#A column is added to the records with the number of counts
#(format A), assuming 1 count per record

dim<-dim(records)
Counts<-rep(2,dim[1])
records<-cbind(records,Counts)
records$Counts<-1

#Running the function
setwd("./temps/knowbr_borrar/")
data(States) #State Boundaries of the United States
data(adworld)
KnowBPolygon(data=records, format="A", shape=States, admAreas=T, shapenames="NAME", minLon=-130,
             maxLon=-70, minLat=25, maxLat=50, colscale=rev(heat.colors(100)), jpg=FALSE,
             curve="Rational", save="csv", 
             file1="Species per site",
             file2="Estimators", 
             file3="Standard error of the estimators", 
             na="NA",
             dec=".", row.names=FALSE)

## End(Not run)
plot(States, col="lightgrey", border="black", lwd=0.5, main="States of the United States")
points(records$Longitude, records$Latitude, pch=20, col="red", cex=0.5)

States |> 
  st_as_sf() |>
  ggplot() +
  geom_sf() +
  geom_sf_text(aes(label = NAME))+
  geom_point(data=records, aes(x=Longitude, y=Latitude, color=Species), size=2) +
  theme_minimal() +
  xlim(-130, -65) +
  ylim(25, 50) + 
  theme(legend.position = "none")

States |> 
  st_as_sf() |> 
  as.data.frame() |> 
  select(!geometry) |> 
  arrange(STUSPS) |> 
  head()

estimators<-read.csv("./Estimators.CSV")
estimators
spxsite<-read.csv("./Species per site.CSV")
SEestim<-read.csv("./Standard error of the estimators.CSV")

unique(records$Species)
States@plotOrder


#### grafico porcentaje
library(terra)


grids<-list.files("./outputs/borrar/phd_riqobs/", pattern="_2.tif", full.names = T)

df.grid<-data.frame(grid=NA, below=NA, above=NA)
for (i in 1:length(grids)) {
  grid<-rast(grids[i])
  df.grid[i,1]<-grids[i] |> str_remove("./outputs/borrar/phd_riqobs/dif_") |>  str_remove("_2.tif")
  df.grid[i,2]<-length(grid[grid<100])
  df.grid[i,3]<-length(grid[grid>100])
}


df.grid |> 
  mutate(total=(below + above),
         bel.por=below*100/total,
         abo.por=above*100/total) |> 
  pivot_longer(cols=c(bel.por, abo.por), names_to = "type", values_to = "percentage") |> 
  ggplot(aes(x=factor(grid, levels=c("10","5","1","0.5","0.083")), y=percentage, fill=type)) +
  geom_bar(stat='identity') +
  labs(x="", y="%") +
  #  scale_fill_manual(values=c("darkred", "orange")) +
  scale_fill_discrete(name="% Obs.", labels=c(">100 (Esp < Obs)", "<100 (Esp > Obs)"), type=c("darkred", "lightblue"))+
  theme_classic()

ggsave(filename = "./outputs/borrar/phd_riqobs/dif_perc.png",
       width = 13,
       height = 10, #alto
       scale=1,
       units ="cm",
       dpi = 200)

### Curva de acumulación de individuos por especie
library(terra)
library(sf)
library(tidyverse)

grid.shp<-read_sf("./outputs/borrar/10grados/10_grid.shp")
plot(grid.shp)

spp.list<-list.files("./species/mig_shapes/ok/", full.names=T, pattern = ".shp")
spp.names<-list.files("./species/mig_shapes/ok/", full.names=F, pattern = ".shp") |> str_remove(".shp")


data.list<-list()

for (x in 1:length(spp.list)) { 
  # Abrir shape de la especie
  spp <- read_sf(spp.list[x])
data<-data.frame(matrix(ncol=12, nrow=nrow(grid.shp)))
names(data)<-month.name
  
  # Filtro mensual
  for(m in 1:12){
    # Vector mensual
    spp.month <-
      spp |>
      filter(month == m) 
    # Conteo de valore spor cuadrante
    qxspp <- st_contains_properly(st_make_valid(grid.shp), spp.month)# arroja el número de la fila que coincide en cada polígono de grid.shp
    
    data[,m]<-
    qxspp |>
      lapply(FUN = length) |>
      unlist()
  }
data.list[[x]]<-data
  print(x)
}

library(iNEXT)


a<- spp |> 
  as_tibble() |> 
  summarise(n=n(), .by=month) |> 
  as.data.frame()

out<- iNEXT(a$n, q = 0, datatype = "abundance")
ggiNEXT(out, type = 2)

DataInfo(a)
data(bird)

spp
datos<- data.frame(matrix(nrow=spp$year |> unique() |> length(), ncol=13))
names(datos) <- c("year", month.name)

year.sp<-spp$year |> unique()
for (y in 1:length(year.sp)) {
  spp.y<-
    spp |> 
    filter(year==year.sp[y])
  datos[y,1]<-year.sp[y]
for(m in 2:13){
  # Vector mensual
  datos[y,m] <-
    spp.y |>
    filter(month == m) |> 
    nrow()
}
}

datos2<-datos[,2:13]
row.names(datos2)<-datos[,1]
datos3<-datos2+1

DataInfo(datos3, datatype="abundance")

out<-iNEXT(datos3, q=0, datatype="abundance", endpoint = 500)
ggiNEXT(out, type=2, facet.var="Assemblage", color.var = "Assemblage")
