## Analisis de completitud 

# Probando iNEXT
#install.packages("iNEXT")
library(iNEXT)
library(tidyverse)

spp.path<-list.files("./Phd/species/species/tablas_sindups/", full.names = T, pattern=".txt")

# rows<-sample(1:length(spp.path), 50, replace = FALSE)

spp.list<-list()

for(x in 1:length(spp.path)){
  spp.list[[x]]<- read_tsv(paste0(spp.path[x])) |> 
    select(species, decimalLatitude, decimalLongitude, month, year)
}

spp_df <- bind_rows(spp.list)

spp_df2<- spp_df |> 
  na.omit() |> 
  filter(decimalLatitude<30) |> #month==1 & 
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
           decimalLatitude > -65 &  decimalLatitude < -60 ~ "F19")) 

franjas<-paste0("F",1:19)

data.comp<-as.data.frame(matrix(ncol=12, nrow = length(franjas)))
row.names(data.comp)<-franjas
names(data.comp)<-month.name

# F1
for (f in 1:nrow(data.comp)) {
  for (m in 1:ncol(data.comp)) {
    complet<-
      spp_df2 |> 
      filter(franjas==franjas[f]) |> # franja
      summarise(n = n(), .by = c(species, month, rangos)) |>  
      filter(month==m) |> # mes
      pivot_wider(names_from = rangos, values_from = n) |> 
      mutate(across(everything(), ~ replace_na(.x, 0))) |> 
      column_to_rownames(var = "species") |> 
      select(!month) |>
      DataInfo(datatype = "abundance") |> 
      mutate(rango_limpio = str_replace_all(Assemblage, "[\\(\\)\\[\\]]", ""),
             inicio = as.numeric(str_extract(rango_limpio, "^[^,]+")),  # Extrae el primer número
             fin = as.numeric(str_extract(rango_limpio, "[^,]+$")),     # Extrae el segundo número
             rango_texto = paste(as.integer(inicio), "to", as.integer(fin))) |> 
      select(!c(rango_limpio, inicio, fin))
    
    # complet |> 
    # 	summarise(median=median(SC), mean=median(SC))
    
    data.comp[f,m]<- median(complet$SC)
    
    complet |> 
      ggplot(aes(x=rango_texto, y=SC)) +
      geom_bar(stat = 'identity') +
      labs(title=paste0("Franja 1 (30-25) - Enero", " SC_median: ", median(complet$SC))) +
      theme(axis.text.x = element_text(angle = 90),
            axis.title.x = element_blank()) 
    
    ggsave(filename = paste0("./borrar/graphs_completeness/",franjas[f],"_",m,".png"),
           width = 10,
           height = 10,
           scale=2,
           units = "cm",
           dpi=200)
  }
}

str(a)

summarise(n = n(), .by = c(rangos, franjas)) |>
  column_to_rownames(var = "rangos") |> 
  
  
  data[is.na(data$franjas),]

library(dplyr)
library(stringr)

# Ejemplo de datos
df <- tibble(rango = c("(2.01e+03,2.02e+03]", "(2.02e+03,2.03e+03]", "(1.95e+03,2.00e+03]"))

df <- df %>%
  mutate(
    # Eliminar paréntesis y corchetes
    rango_limpio = str_replace_all(rango, "[\\(\\)\\[\\]]", ""),
    
    # Separar los límites del rango
    inicio = as.numeric(str_extract(rango_limpio, "^[^,]+")),  # Extrae el primer número
    fin = as.numeric(str_extract(rango_limpio, "[^,]+$")),     # Extrae el segundo número
    
    # Convertir a formato "YYYY to YYYY"
    rango_texto = paste(as.integer(inicio), "to", as.integer(fin))
  ) %>%
  select(rango, rango_texto)  # Mantiene solo columnas relevantes

print(df)

## Ahora con la idea del profe
load("./scripts/Rdatas/spp.collapsed.Rdata")

# 1. Hacer un mapa de riqueza esperada con los polígonos de la IUCN de mis 243 spp 
library(terra)
library(letsR)
library(sf)
library(tidyverse)

iucn.spp2<-vect("./Phd/species/species/mig_shapes/birdlife_mig2.shp")
iucn.spp<-iucn.spp2[-13,] #G. granadensis

#emp.iucn<-vect("./Phd/borrar_phd/spp_borrar/Emp_birdlife.shp")

PAM<-lets.presab(iucn.spp, 
                 xmn = -170,  
                 xmx = -10,
                 ymn = -100,
                 ymx = 100,
                 reso=0.083,
                 crs = "+proj=longlat +datum=WGS84")

PAM$Richness_Raster |> plot()

PAM$Richness_Raster |> 
  writeRaster("./Phd/species/species/mig_shapes/0.5grados/birdlife_rich0.083.tif",overwrite=TRUE)#

# 2. Hacer le mapa de riqueza observada
# iucn.shp<-read_sf("./Phd/species/species/mig_shapes/5grados/birdlife_rich5.shp")
iucn.shp<-read_sf("./Phd/species/species/mig_shapes/0.083grados/0.083_grid.shp")
plot(iucn.shp)

# Contar cuantas especies hay en cada cuadrante

# Cargar listas de especies

# Total -------------------------------------------------------------------
spp.list<-list.files("./Phd/species/species/mig_shapes/mig_shapes/", pattern=".shp",full.names = T)[-c(91,32,93,95,147,113,210)] # Retiro especies con inconsistencias taxonómicas
spp.names<-list.files("./Phd/species/species/mig_shapes/mig_shapes/", pattern=".shp",full.names = F)[-c(91,32,93,95,147,113,210)] |> 
  str_remove(".shp")


ab.table<-as.data.frame(matrix(ncol=length(spp.list), nrow = nrow(iucn.shp)))
names(ab.table) <- spp.names
dim(ab.table)

sf_use_s2(FALSE)
for (x in 1:length(spp.list)) { 
  spp <- #read_sf(spp.list[x])
    
    # Intercepto los puntos dentro de cada centroide
    qxspp <- st_intersects(st_make_valid(iucn.shp), spp) # arroja el número de la fila
  
  # Abundancia
  ab.table[, x] <-
    qxspp |>
    lapply(FUN = length) |>
    unlist()
  
  print(x)
}

write.table(ab.table, "./Phd/species/species/mig_shapes/ab.table0.083.txt", sep="\t", dec=".", row.names = F)
read.table("./Phd/species/species/mig_shapes/ab.table0.5.txt", sep="\t", dec=".", header=T)
# load("./scripts/Rdatas/Obs.richness.Rdata")
# rm(list=ls()[which(ls()!=C("ab.table10","ab.table5"))])
# save.image("./scripts/Rdatas/Obs.richness.Rdata")


# Obs para Breeding ----------------------------------------------------------------
# 1. Necesito saber cuántas y cuáles especies tienen polígono de breeding
bl.bred<-vect("./Phd/species/species/mig_shapes/temps/birdlife_mig2_bred.shp")
bl.bred$seasonal |> unique() #2= breeding

# 2. Cargo la lista de especies observadas
spp.list<-list.files("./Phd/species/species/mig_shapes/mig_shapes/", pattern=".shp",full.names = T)[-c(91,32,93,95,147,113,210)] # Retiro especies con inconsistencias taxonómicas
spp.names<-list.files("./Phd/species/species/mig_shapes/mig_shapes/", pattern=".shp",full.names = F)[-c(91,32,93,95,147,113,210)] |> 
  str_remove(".shp")

a<-bl.bred |> 
  as.data.frame() |> 
  select(sci_name) |> 
  mutate(codigo = str_split(sci_name, " ", simplify = TRUE)) # separar en palabras

paste(substr(a[,1], 1, 3), substr(a[,2], 1, 5), sep="_") |>      # extraer y unir partes
  unique() |> 
  sort()


# ---------------------------
gbif.10<-read.table("./Phd/species/species/mig_shapes/ab.table10.txt", sep="\t", dec=".", header=T)
gbif.5<-read.table("./Phd/species/species/mig_shapes/ab.table5.txt", sep="\t", dec=".", header=T)
gbif.1<-read.table("./Phd/species/species/mig_shapes/ab.table1.txt", sep="\t", dec=".", header=T)
gbif.05<-read.table("./Phd/species/species/mig_shapes/ab.table0.5.txt", sep="\t", dec=".", header=T)
#gbif.008<-read.table("./Phd/species/species/mig_shapes/ab.table0.083.txt", sep="\t", dec=".", header=T)


# Para 10 grados ----------------------------------------------------------
iucn.shp10<-read_sf("./Phd/species/species/mig_shapes/10grados/10_grid.shp")
bl.rich10<-rast("./Phd/species/species/mig_shapes/10grados/birdlife_rich10.tif")
bl.rich10[bl.rich10$lyr.1==0]<-NA

pam.gbif10<-
  gbif.10 |> 
  mutate(across(
    .cols = everything(),
    .fns = function(x)
      ifelse(x > 0, 1, 0)))

iucn.shp10$abun.obs10<-
  apply(gbif.10, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp10$rich.obs10<-
  apply(pam.gbif10, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp10$rich.obs10[which(iucn.shp10$rich.obs10==0)]<-NA

# Rasterizar
plot(iucn.shp10[,4])

obs.rast10<-rasterize(vect(iucn.shp10), bl.rich10, field="rich.obs10")

dif10<-obs.rast10-bl.rich10
dif10<-obs.rast10*100/bl.rich10

writeRaster(dif10,"./Phd/species/species/mig_shapes/10grados/dif_10_2.tif",overwrite=T)

# Para 5 grados -----------------------------------------------------------
iucn.shp5<-read_sf("./Phd/species/species/mig_shapes/5grados/5_grid.shp")
bl.rich5<-rast("./Phd/species/species/mig_shapes/5grados/birdlife_rich5.tif")
bl.rich5[bl.rich5$lyr.1==0]<-NA

pam.gbif5<-
  gbif.5 |> 
  mutate(across(.cols = everything(), .fns = function(x) ifelse(x > 0, 1, 0)))

iucn.shp5$abun.obs5<-
  apply(gbif.5, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp5$rich.obs5<-
  apply(pam.gbif5, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp5$rich.obs5[which(iucn.shp5$rich.obs5==0)]<-NA

# Rasterizar
plot(iucn.shp5[,4])

obs.rast5<-rasterize(vect(iucn.shp5), bl.rich5, field="rich.obs5")

dif5<-obs.rast5-bl.rich5
dif5<-obs.rast5*100/bl.rich5
plot(dif5)

writeRaster(dif5,"./Phd/species/species/mig_shapes/5grados/dif_5_2.tif",overwrite=T)

# Para 1 grado -----------------------------------------------------------
iucn.shp1<-read_sf("./Phd/species/species/mig_shapes/1grado/1_grid.shp")
bl.rich1<-rast("./Phd/species/species/mig_shapes/1grado/birdlife_rich1.tif")

bl.rich1[bl.rich1$lyr.1==0]<-NA

pam.gbif1<-
  gbif.1 |> 
  mutate(across(.cols = everything(), .fns = function(x) ifelse(x > 0, 1, 0)))

iucn.shp1$abun.obs1<-
  apply(gbif.1, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp1$rich.obs1<-
  apply(pam.gbif1, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp1$rich.obs1[which(iucn.shp1$rich.obs1==0)]<-NA

# Rasterizar
plot(iucn.shp1[,4])

obs.rast1<-rasterize(vect(iucn.shp1), bl.rich1, field="rich.obs1")

dif1<-obs.rast1-bl.rich1
dif1<-obs.rast1*100/bl.rich1
plot(dif1)

writeRaster(dif1,"./Phd/species/species/mig_shapes/1grado/dif_1_2.tif",overwrite=T)

# Para 0.5 grado -----------------------------------------------------------
iucn.shp0.5<-read_sf("./Phd/species/species/mig_shapes/0.5grados/0.5_grid.shp")
bl.rich0.5<-rast("./Phd/species/species/mig_shapes/0.5grados/birdlife_rich0.5.tif")

bl.rich0.5[bl.rich0.5$lyr.1==0]<-NA

pam.gbif0.5<-
  gbif.05 |> 
  mutate(across(.cols = everything(), .fns = function(x) ifelse(x > 0, 1, 0)))

iucn.shp0.5$abun.obs0.5<-
  apply(gbif.05, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp0.5$rich.obs0.5<-
  apply(pam.gbif0.5, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp0.5$rich.obs0.5[which(iucn.shp0.5$rich.obs0.5==0)]<-NA

# Rasterizar
plot(iucn.shp0.5[,4])

obs.rast0.5<-rasterize(vect(iucn.shp0.5), bl.rich0.5, field="rich.obs0.5")

dif0.5<-obs.rast0.5-bl.rich0.5
dif0.5<-obs.rast0.5*100/bl.rich0.5

plot(dif0.5)

writeRaster(dif0.5,"./Phd/species/species/mig_shapes/0.5grados/dif_0.5_2.tif",overwrite=T)

# Para 0.083 grado -----------------------------------------------------------
iucn.shp0.08<-read_sf("./Phd/species/species/mig_shapes/0.083grados/0.083_grid.shp")
bl.rich0.08<-rast("./Phd/species/species/mig_shapes/0.083grados/birdlife_rich0.083.tif")

bl.rich0.08[bl.rich0.08$lyr.1==0]<-NA

pam.gbif0.083<-
  gbif.008 |> 
  mutate(across(.cols = everything(), .fns = function(x) ifelse(x > 0, 1, 0)))

iucn.shp0.08$abun.obs0.083<-
  apply(gbif.008, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp0.08$rich.obs0.083<-
  apply(pam.gbif0.083, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp0.08$rich.obs0.083[which(iucn.shp0.08$rich.obs0.083==0)]<-NA

# Rasterizar
# plot(iucn.shp0.083[,4])

obs.rast0.083<-rasterize(vect(iucn.shp0.08), bl.rich0.08, field="rich.obs0.083")
writeRaster(obs.rast0.083,"./Phd/species/species/mig_shapes/0.083grados/obs_0.083.tif",overwrite=T)

dif0.083<-((obs.rast0.083+1)-(bl.rich0.08+1))/(bl.rich0.08+1)
dif0.083<-obs.rast0.083*100/bl.rich0.08
plot(dif0.083)

dif0.083 |> 
  as.data.frame() |> 
  filter(rich.obs0.083<100) |> 
  ggplot(aes(x=rich.obs0.083))+
  geom_histogram()+
  labs(x="porcentaje", x="frecuencia")

writeRaster(dif0.083,"./Phd/species/species/mig_shapes/0.083grados/dif_0.083_2.tif",overwrite=T)


# Ahora voy a hacer lo que propone Backtrom et al. 2024

# 1. hay que calcular el "inventory completeness"

# Para 1 grado -----------------------------------------------------------
iucn.shp1<-read_sf("./Phd/species/species/mig_shapes/1grado/1_grid.shp")
bl.rich1<-rast("./Phd/species/species/mig_shapes/1grado/birdlife_rich1.tif") #Esperado

bl.rich1[bl.rich1$lyr.1==0]<-NA

pam.gbif1<-
  gbif.1 |> 
  mutate(across(.cols = everything(), .fns = function(x) ifelse(x > 0, 1, 0)))

iucn.shp1$abun.obs1<-
  apply(gbif.1, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp1$rich.obs1<-
  apply(pam.gbif1, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp1$rich.obs1[which(iucn.shp1$rich.obs1==0)]<-NA

# Rasterizar
plot(iucn.shp1[,4])

obs.rast1<-rasterize(vect(iucn.shp1), bl.rich1, field="rich.obs1")

IC_1<-obs.rast1/bl.rich1
plot(IC_1)

writeRaster(IC_1,"./Phd/species/species/mig_shapes/inventory_completeness/IC_1.tif",overwrite=T)


# Para 0.5 grado -----------------------------------------------------------
iucn.shp0.5<-read_sf("./Phd/species/species/mig_shapes/0.5grados/0.5_grid.shp")
bl.rich0.5<-rast("./Phd/species/species/mig_shapes/0.5grados/birdlife_rich0.5.tif")

bl.rich0.5[bl.rich0.5$lyr.1==0]<-NA

pam.gbif0.5<-
  gbif.05 |> 
  mutate(across(.cols = everything(), .fns = function(x) ifelse(x > 0, 1, 0)))

iucn.shp0.5$abun.obs0.5<-
  apply(gbif.05, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp0.5$rich.obs0.5<-
  apply(pam.gbif0.5, MARGIN = 1, FUN = sum, na.rm=T)

iucn.shp0.5$rich.obs0.5[which(iucn.shp0.5$rich.obs0.5==0)]<-NA

# Rasterizar
plot(iucn.shp0.5[,4])

obs.rast0.5<-rasterize(vect(iucn.shp0.5), bl.rich0.5, field="rich.obs0.5")

IC_05<-obs.rast0.5/bl.rich0.5


plot(IC_05)

writeRaster(IC_05,"./Phd/species/species/mig_shapes/inventory_completeness/IC_05.tif",overwrite=T)


# Mean Inventory Completeness ---------------------------------------------

# MIC - 1 grado
spp.list
spp.names

ic_rast1<- rast("./Phd/species/species/mig_shapes/inventory_completeness/IC_1.tif")
ic_rast05<- rast("./Phd/species/species/mig_shapes/inventory_completeness/IC_05.tif")
mic_table<- data.frame(species=spp.names, MIC_1=NA,  sd_1=NA, MIC_05=NA, sd_05=NA)

pb<- txtProgressBar(min=0, max= nrow(spp.names), style=3)

for (x in 1:length(spp.names)) {
  # 1. Open species shapefiles
  spp<-vect(spp.list[x])
  ext1<- terra::extract(ic_rast1, spp, ID=F)
  ext05<- terra::extract(ic_rast05, spp, ID=F)
  ic_table[x,2]	<- mean(ext1$rich.obs1, na.rm=T)
  ic_table[x,3]	<- sd(ext1$rich.obs1, na.rm=T)
  ic_table[x,4]	<- mean(ext05$rich.obs0.5, na.rm=T)
  ic_table[x,5]	<- sd(ext05$rich.obs0.5, na.rm=T)
  
  setTxtProgressBar(pb,x)
}


# Total Inventory completeness --------------------------------------------
iucn.spp
ic_rast1
ic_rast05

tic_table<- data.frame(species=NA, TIC_1=NA,  sd_1=NA, TIC_05=NA, sd_05=NA)

pb<- txtProgressBar(min=0, max= nrow(iucn.spp), style=3)

for (x in 1:nrow(iucn.spp)) {
  sp<-iucn.spp[x,3] |> 
    as.data.frame() |> 
    pull() |> 
    strsplit(" ") |> 
    unlist()
  
  tic_table[x,1] <-	paste0(substr(sp[1], 1, 3), "_", substr(sp[2], 1, 5))
  
  tic_table[x,2]<- crop(ic_rast1, iucn.spp[x], mask=T) |> global("mean", na.rm=T) |> pull()
  tic_table[x,3]<- crop(ic_rast1, iucn.spp[x], mask=T) |> global("sd", na.rm=T) |> pull()
  tic_table[x,4]<- crop(ic_rast05, iucn.spp[x], mask=T) |> global("mean", na.rm=T) |> pull()
  tic_table[x,5]<- crop(ic_rast05, iucn.spp[x], mask=T) |> global("sd", na.rm=T) |> pull()
  
  setTxtProgressBar(pb,x)
}

# Unir ambos en una sola tabla
tic_table
mic_table

ic_table<-inner_join(tic_table, mic_table, by=join_by(species))

# save(list=c("ic_table","tic_table", "mic_table"), file="D:/Mauro/scripts/Rdatas/completenes.Rdata")
# load("D:/Mauro/scripts/Rdatas/completenes.Rdata")

min(ic_table$TIC_1)

TIC1_g<- 
  ic_table |> 
  ggplot(aes(x=species, y=TIC_1)) +
  geom_bar(stat='identity', fill="darkcyan") +
  # geom_errorbar(aes(ymin=TIC_1-sd_1.x, ymax=TIC_1+sd_1.x), width=.2,
  # 							position=position_dodge(.9)) +
  geom_hline(yintercept = c(1,0.5), col="red", linetype='dashed') +
  geom_hline(yintercept = min(ic_table$TIC_1), col="darkred", linetype='dashed') +
  labs(x="Especies", y="Total Inventory Completeness", title="1°")+
  theme_classic()+
  theme(axis.text.x = element_blank()) 

TIC05_g<- 
  ic_table |> 
  ggplot(aes(x=species, y=TIC_05)) +
  geom_bar(stat='identity', fill="darkcyan") +
  geom_hline(yintercept = c(1,0.5), col="red", linetype='dashed') +
  geom_hline(yintercept = min(ic_table$TIC_05), col="darkred", linetype='dashed') +
  labs(x="Especies", y="", title="0.5°")+
  theme_classic()+
  theme(axis.text.x = element_blank()) 

MIC1_g<- 
  ic_table |> 
  ggplot(aes(x=species, y=MIC_1)) +
  geom_bar(stat='identity', fill="darkolivegreen") +
  geom_hline(yintercept = c(1,0.5), col="red", linetype='dashed') +
  geom_hline(yintercept = min(ic_table$MIC_1), col="darkred", linetype='dashed') +
  labs(x="Especies", y="Mean Inventory Completeness", title="1°")+
  theme_classic()+
  theme(axis.text.x = element_blank()) 

MIC05_g<- 
  ic_table |> 
  ggplot(aes(x=species, y=MIC_05)) +
  geom_bar(stat='identity', fill="darkolivegreen") +
  geom_hline(yintercept = c(1,0.5), col="red", linetype='dashed') +
  geom_hline(yintercept = min(ic_table$MIC_05), col="darkred", linetype='dashed') +
  labs(x="Especies", y="", title="0.5°")+
  theme_classic()+
  theme(axis.text.x = element_blank()) 

#library(patchwork)
(TIC1_g | TIC05_g )/(MIC1_g | MIC05_g)


#### Idea del Ox, contar por cuadrante cuantos meses ha sido muestreada una especie, esto lo voy a hacer para 4 especies que tienen poblaciones residentes chiquitas
library(sf)
library(terra)
library(tidyverse)
library(tidyterra)

# Cargar especie

spp <- read_sf("./species/mig_shapes/ok/Set_cerul.shp")
grid <- read_sf("./inputs/america_shapes/amrica_1dgr/1_grid.shp") |> #"./outputs/borrar/10grados/10_grid.shp"
  mutate(month=NA)

#sf_use_s2(FALSE)

# Intercepto los puntos dentro de cada centroide
qxspp <- st_intersects(st_make_valid(grid), spp) # arroja el número de la fila

pb<- txtProgressBar(min=0, max= nrow(grid), style=3) # Process bar

for (i in 1:nrow(grid)) {
  value<-  spp[qxspp[[i]],] |> 
    as.data.frame() |> 
    select(month) |> 
    unique() |> 
    nrow()
  if(value!=0){
    grid[i,3]<- value
  }
  else{grid[i,3]<- NA}
  
  setTxtProgressBar(pb,i)
}

# rasterizar
grid1.template<-rast("./outputs/borrar/10grados/grid10.tif")

grid1.r<-rasterize(vect(grid), grid1.template, field="month")
# residentes
# plot(grid1.r, xlim=c(-115,-90), ylim=c(15,35)) # Calothorax lucifer
# plot(grid1.r, xlim=c(-128,-84), ylim=c(8,47)) # Cardellina rubrifrons
# plot(grid1.r, xlim=c(-180,-40), ylim=c(10,75)) # Catharus guttatus
# plot(grid1.r, xlim=c(-140,-50), ylim=c(5,60)) # Icterus bullockii

# ejemplos migratorios full
# ver setophaga cerulea
plot(grid1.r, xlim=c(-140,-45), ylim=c(-20,70)) # Cardellina canadensis
# plot(grid1.r, xlim=c(-130,-45), ylim=c(-5,55)) # Vermivora chrysoptera
# plot(grid1.r, xlim=c(-140,-45), ylim=c(-5,70)) # Setophaga castanea
plot(grid1.r)


# Plotting
amer <- vect("./inputs/america_shapes/america_cont/America.shp")

colores <- c(
  "#e4e4e4",
  "#cecece",
  "#b8b8b8",
  "#a3a3a3",
  "#8e8e8e",
  "#7a7a7a",
  "#666666",
  "#535353",
  "#414141",
  "#303030",
  "#1f1f1f",
  "#B32929")


ggplot(amer)+
  geom_spatraster(data=grid1.r, aes(fill=month)) +
  geom_spatvector(fill=NA)+
  coord_sf(# st_bbox calcula el extent de la capa mas pequeña
    xlim = c(st_bbox(spp)$xmin, st_bbox(spp)$xmax), 
    ylim = c(st_bbox(spp)$ymin, st_bbox(spp)$ymax), 
    expand = T) + # expand = T genera un plot mas grande que los puntos y ya no queda tan preciso
  scale_fill_gradientn(colours = colores[1:global(grid1.r, fun="max", na.rm=T)$max], 
                       na.value = "transparent") +
  labs(title="Setophaga cerulea**")+
  theme_minimal()

ggsave(filename = "./outputs/images/Set_cerul_month.jpg",
       width = 10,
       height = 10, #alto
       scale=1.5,
       units ="cm",
       dpi = 200)


### Mensualidad por franja latitudinal ----------------
spp.list<-list.files("./species/mig_shapes/ok/", pattern=".shp$", full.names = T)
spp.names<-list.files("./species/mig_shapes/ok/", pattern=".shp$") |> str_remove(".shp")
amer


i<-6

spp1 <- read_sf(spp.list[i])

coords<-
  spp1 |> 
  st_coordinates()

spp <- spp1 |> 
  mutate(lon=coords[,1],
         lat=coords[,2]) |> 
  mutate(franjas= case_when(
    lat > 60 &  lat < 75   ~ "F1",
    lat > 45 &  lat < 60   ~ "F2",
    lat > 30 &  lat < 45   ~ "F3",
    lat > 15 &  lat < 30   ~ "F4",
    lat >0 &  lat < 15     ~ "F5",
    lat > -15 &  lat < 0   ~ "F6",
    lat > -30 &  lat < -15 ~ "F7",
    lat > -45 &  lat < -30 ~ "F8",
    lat > -60 &  lat < -45 ~ "F9"))

spp.graf<-
  spp |> 
  as.data.frame() |> 
  summarise(n=n(), .by=c(franjas, month)) |> 
  ggplot(aes(x=as.factor(month), y=n)) +
  geom_bar(stat='identity') +
  labs(y="Registros", x="", title=spp.names[i]) +
  facet_wrap(~franjas, scale="free_y")

mapa<-ggplot() +
  geom_sf(data=st_as_sf(amer)) + 
  geom_sf(data=spp, size=0.5, alpha=0.4, color="gray20") +
  scale_color_manual(values=)+
  coord_sf(# st_bbox calcula el extent de la capa mas pequeña
    xlim = c(st_bbox(spp)$xmin, st_bbox(spp)$xmax), 
    ylim = c(st_bbox(spp)$ymin, st_bbox(spp)$ymax), 
    expand = T) + # expand = T genera un plot mas grande que los puntos y ya no queda tan preciso
  geom_hline(yintercept = c(75,60,45,30,15,0,-15,-30,-45,-60), color="red", linetype='dotted')

#library(patchwork)
spp.graf+mapa

ggsave(filename = paste0("./outputs/images/", spp.names[i],"_franj.jpg"),
       width = 10,
       height = 5, #alto
       scale=3,
       units ="cm",
       dpi = 200)



