

# comparando monhts y season ----------------------------------------------
# Broennimann
list.m<-list()
list.s<-list()
df.m <- list.files("./species/mig_ENM/overlaps_tables/month/", full.names=T)[1:8]
df.s <- list.files("./species/mig_ENM/overlaps_tables/season/", full.names=T)[1:8]

for (x in 1:length(df.m)) {
  a <- read.table(df.m[x], sep="\t", dec = ".", header=T)
  list.m[[x]] <- a
  b <- read.table(df.s[x], sep="\t", dec = ".", header=T)
  list.s[[x]] <- b
}

data.ovrlp.m <- do.call(rbind, list.m)
data.ovrlp.s <- do.call(rbind, list.s)

# write.table(data.ovrlp.m, "./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt",
#             sep="\t", dec = ".", row.names=F)
# write.table(data.ovrlp.s, "./species/mig_ENM/overlaps_tables/data.ovrlp_S_df.txt",
#             sep="\t", dec = ".", row.names=F)

# data.ovrlp.m<-read.table("./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt", 
#                           sep="\t", dec = ".", header=T)
# 
# data.ovrlp.s<-read.table("./species/mig_ENM/overlaps_tables/data.ovrlp_S_df.txt", 
#                           sep="\t", dec = ".", header=T)

## Cargar valores de jaccard
# Month
spp.list.m <- list.files('./species/mig_ENM/monthly/', full.names=T)[1:8]
uni_ovl.m <- list()

# x <- 1
for(x in 1:length(spp.list.m)) {
  uni_ovl.m[[x]] <- 
  read.table(paste0(spp.list.m[x], '/m_overlap/', basename(spp.list.m[x]), '_m_uni_ovl.txt'), sep='\t', dec='.', header=T) |> 
  # select(spp, niche.comp, overlap, p_value) |> 
  mutate(spp = basename(spp.list.m[x]) |> str_remove('_ENMm'),
         niche.comp = str_replace_all(niche.comp, "Niche_", "m")) |> 
    mutate(niche.comp = str_replace_all(niche.comp, "vs_", "m")) |> 
  rename('overlap'= niche.comp,
         'jaccard' = overlap)}

uni_ovl.df.m <- do.call(rbind, uni_ovl.m) |>
  mutate(overlap=factor(overlap, levels = unique(overlap)),
         color_flag=ifelse(p_value > 0.05, ">0.05", "<0.05"))

 m<-
  data.ovrlp.m |> 
  mutate(overlap=factor(overlap, levels = unique(overlap)),
         color_flag=ifelse(p.D > 0.05, ">0.05", "<0.05")) |> 
  #  summarise(n=n(), .by=c(color_flag, spp))
  ggplot(aes(x=obs.D, y=fct_rev(overlap), color=color_flag)) +
  geom_point(size=2) +
  geom_point(data=uni_ovl.df.m, aes(x=jaccard, y=overlap), size=2, shape = 2) +
  geom_vline(xintercept = 0.5, linetype ="dotted") +
  geom_vline(xintercept = c(0,1), linetype ="dashed") +
  scale_color_manual(values = c("<0.05" = "black",
                                ">0.05" = "red")) +
  xlim(0,1)+
  labs(x="", y="Monthly overlap", color="p-value", size="D observed")+
  facet_wrap(~spp, ncol=19) +
    theme_test() 

# Ahora hay que repetirlo pero para las estaciones 
  spp.list.s <- list.files('./species/mig_ENM/seasonaly/', full.names=T)[1:8]
  uni_ovl.s <- list()
  
  # x <- 1
  for(x in 1:length(spp.list.s)) {
     uni_ovl.s[[x]] <- 
      read.table(paste0(spp.list.s[x], '/S_overlap/', basename(spp.list.s[x]) |> str_remove('_ENMs'), '_s_uni_ovl.txt'), sep='\t', dec='.', header=T) |> 
      # select(spp, niche.comp, overlap, p_value) |> 
      mutate(spp = basename(spp.list.s[x]) |> str_remove('_ENMs'),
             niche.comp = c('S1_S2', 'S1_S3', 'S1_S4', 'S2_S3', 'S2_S4', 'S3_S4')) |> 
      rename('overlap'= niche.comp,
             'jaccard' = overlap)
     }
  
  uni_ovl.df.s <- do.call(rbind, uni_ovl.s) |>
    mutate(overlap=factor(overlap, levels = unique(overlap)),
           color_flag=ifelse(p_value > 0.05, ">0.05", "<0.05"))
  
  
s<-
  data.ovrlp.s |> 
  mutate(overlap=factor(overlap, levels = unique(overlap)),
         color_flag=ifelse(p.D > 0.05, ">0.05", "<0.05")) |> 
  #  summarise(n=n(), .by=c(color_flag, spp))
  ggplot(aes(x=obs.D, y=fct_rev(overlap), color=color_flag)) +
  geom_point(size=2) +
  geom_point(data=uni_ovl.df.s, aes(x=jaccard, y=overlap), size=2, shape = 2) +
  scale_color_manual(values = c("<0.05" = "black",
                                ">0.05" = "red")) +
  xlim(0,1)+
  geom_vline(xintercept = 0.5, linetype ="dotted") +
  geom_vline(xintercept = c(0,1), linetype ="dashed") +
  #  geom_vline(xintercept = 0.05, color="red", linetype ="dotted") +
  labs(x="Overlap", y="Season overlap", color="p-value")+
  scale_y_discrete(label= c("Ve-Ot","Pr-Ot","Pr-Ve", "In-Ot", "In-Ve", "In-Pr"))+ # va de abajo hacia arriba
  #  scale_size(range = c(0, 4)) +  # ajusta rango de tamaños
  facet_wrap(~spp, ncol=19)  +
  theme_test() 

#library(patchwork)
m/s
# cuadrar: 
# 1. decimales del eje X
# 2. nombres del eje Y y quizá el tamaño

ggsave(file = "./outputs/s_m_overlaps.svg",
       width = 1706,
       height = 1271,
       scale=3,
       units ="px")



#### Diferencias con los jaccard de 100k y 50k

temp100k<- new.env() # Esto me crea un ambiente en un conjunto temporal
temp50k<- new.env() # Esto me crea un ambiente en un conjunto temporal
load("./species/mig_ENM/monthly/Sel_calli_ENMm/m_overlap/Sel_calli_ENMm_m_jaccard100k.Rdata", envir = temp100k)
load("./species/mig_ENM/monthly/Sel_calli_ENMm/m_overlap/Sel_calli_ENMm_m_jaccard50k.Rdata", envir = temp50k)

full.overlap100K<- temp100k$full.overlap[[1]]
full.overlap100K$type<-"100k"

full.overlap50K <- temp50k$full.overlap[[1]]
full.overlap50K$type<-"50k"

union.overlap100K<- temp100k$union.overlap[[1]]
union.overlap100K$type<-"100k"

union.overlap50K <- temp50k$union.overlap[[1]]
union.overlap50K$type<-"50k"


rbind(full.overlap50K,full.overlap100K) |> 
#  select(spp, overlap, type) |> 
  ggplot(aes(x=type, y=p_value)) +
  geom_boxplot()

rbind(union.overlap50K,union.overlap100K) |> 
 # select(spp, overlap, type) |> 
  ggplot(aes(x=p_value, y=niche.comp, color=type))+
  geom_point() +
  geom_vline(xintercept = 0.05, color="red", linetype="dashed")

full.overlap50K$overlap - full.overlap100K$overlap

temp100k$union.overlap

wilcox.test(union.overlap50K$p_value,union.overlap50K$p_value)



### MAPAS + ELIPSOIDES
library(tidyverse)
library(tidyterra)
library(terra)

maps.list <- list.files('./outputs/mapas_spp/', full.names = T)
maps.nums <- maps.list |> 
  str_extract("\\d+") |> 
  as.numeric()

america <- 
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |> 
  subset(region_un == "Americas") |> 
  vect()

month <- month.name
month[maps.nums[x]]

x <- 1
for (x in 25:length(maps.list)) {
  
map.r <- rast(maps.list[x])

map.r.rec <- classify(map.r, rbind(c(-Inf, 115.9999, NA)))
names(map.r.rec) <- 'layer.1'

 a <- 
  ggplot() +
  geom_spatraster(data = map.r.rec, aes(fill = layer.1)) +
  scale_fill_viridis_c(na.value = 'transparent') +
  geom_spatvector(data = america, fill = NA, color = 'black') +
  xlim(c(-140,-55)) +
  ylim(c(5, 80)) +
  labs(x = 'Longitude', y = 'Latitude', title = month[maps.nums[x]]) +
  theme_classic() + 
  theme(legend.position = 'none')

ggsave(a,file = paste0("./outputs/mapas_spp/", basename(maps.list[x]) |> str_remove('.tif'),".png"),
       width = 1706,
       height = 1271,
       scale=1,
       units ="px")
}

##### Mapas MOP

rast.list <- list.files('./outputs/MOP_R2/', full.names = T)

# Vector con los nombres originales
 rast.names <-  # la posicion 26 corresponde a la comparacion del las estaciones, la pongo a mano
  basename(rast.list[-26]) |> 
    str_remove('MOP_realms_') |> 
    str_remove('.tif') 
 
# Vector con los meses
# month <- month.name
month

# Extraer los números antes y después del "_"
num1 <- as.numeric(sub("m(\\d+)_.*", "\\1", rast.names))
num2 <- as.numeric(sub(".*_(\\d+)", "\\1", rast.names))

# Crear nuevo vector con nombres de los meses
nuevo <- paste(month[num1], month[num2], sep = " (Neart) vs ")
nuevo2 <- c(nuevo, 'Winter (Neart) vs Summer')

x <- 1
for (x in 1:length(rast.list)) {
  
mop <- rast(rast.list[x])
names(mop) <- 'layer.1'
plot(mop)  

b <- ggplot() +
  geom_spatraster(data = mop, aes(fill=layer.1)) +
  scale_fill_viridis_c(na.value = 'transparent') + 
  labs(title = paste(nuevo2[x], '(Neotr)'), fill = 'MOP') + 
  theme_classic()

ggsave(b, file = paste0("./outputs/MOP_R2/", basename(rast.list[x]) |> str_remove('.tif'),".png"),
       width = 1706,
       height = 1271,
       scale=1,
       units ="px")
}




### guardar los valores de uni_olv y full_ovl por especie para las estaciones
spp.list.s <- list.files('./species/mig_ENM/seasonaly/', full.names=T)

x <- 1
for (x in 1:length(spp.list.s)) {
load(paste0(spp.list.s[x], '/S_overlap/', basename(spp.list.s[x]) |> str_remove('_ENMs'), '_jaccard.Rdata'))  

  Sup_all@full_overlap |> 
    mutate(spp=basename(spp.list.s[x]) |> str_remove('_ENMs') 
           |> str_remove("_ENMs")) |> 
    rownames_to_column(var = "niche.comp") |> 
    relocate(spp) |> 
    write.table(paste0(spp.list.s[x], '/S_overlap/', basename(spp.list.s[x]) |> str_remove('_ENMs'), "_s_full_ovl.txt"), 
              sep="\t", dec = ".", row.names=F)
  
  Sup_all@union_overlap |> 
    mutate(spp=basename(spp.list.s[x]) |> str_remove('_ENMs') 
           |> str_remove("_ENMs")) |> 
    rownames_to_column(var = "niche.comp") |> 
    relocate(spp) |> 
    write.table(paste0(spp.list.s[x], '/S_overlap/', basename(spp.list.s[x]) |> str_remove('_ENMs'), "_s_uni_ovl.txt"), 
              sep="\t", dec = ".", row.names=F)
  
  }
