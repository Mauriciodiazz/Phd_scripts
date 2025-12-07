### Figuras


# comparando monhts y season ----------------------------------------------
library(tidyverse)

# Broennimann
list.m<-list()
list.s<-list()
df.m <- list.files("./species/mig_ENM/overlaps_tables/month/", full.names=T)
df.s <- list.files("./species/mig_ENM/overlaps_tables/season/", full.names=T)

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

# ## Cargar valores de jaccard
# # Month
# spp.list.m <- list.files('./species/mig_ENM/monthly/', full.names=T)[c(65,68:85)] # Especies con jaccard mensual
# uni_ovl.m <- list()
# 
# # x <- 1
# for(x in 1:length(spp.list.m)) {
#   uni_ovl.m[[x]] <- 
#   read.table(paste0(spp.list.m[x], '/m_overlap/', basename(spp.list.m[x]), '_m_uni_ovl.txt'), sep='\t', dec='.', header=T) |> 
#   # select(spp, niche.comp, overlap, p_value) |> 
#   mutate(spp = basename(spp.list.m[x]) |> str_remove('_ENMm'),
#          niche.comp = str_replace_all(niche.comp, "Niche_", "m")) |> 
#     mutate(niche.comp = str_replace_all(niche.comp, "vs_", "m")) |> 
#   rename('overlap'= niche.comp,
#          'jaccard' = overlap)}
# 
# uni_ovl.df.m <- do.call(rbind, uni_ovl.m) |>
#   mutate(overlap=factor(overlap, levels = unique(overlap)),
#          color_flag=ifelse(p_value > 0.05, ">0.05", "<0.05"))

# Graph congreso ---------------------------------------
# Cargo base de datos donde estan todos los nombres (visual)
spp.fl <- read.table('./inputs/final-spp-list.txt', header=T, sep='\t', dec='.') |> 
  separate(Scientific.Name, into = c('gen', 'sp'), sep = ' ', remove=F) |> 
  mutate(spp=paste(str_sub(gen, 1,3), str_sub(sp, 1,5), sep = '_')) |> 
  select(!c(Taxonomic.Order, gen, sp))

# Orden de especies
sp.ord <- 
  spp.fl$Scientific.Name[c(73,83,76,79,81,82)]
# "Set_chrys" "Set_stria" "Set_magno" "Set_palma" "Set_petec" "Set_rutic"

eca.fct <- data.ovrlp.m |> 
  # Solo para parulidos
  filter(str_starts(spp, 'Set_')) |> 
  summarise(mean=mean(obs.D), .by=c(spp)) |> 
  arrange(mean) |> 
  select(spp) |> 
  pull()

m<-
  data.ovrlp.m |> 
  left_join(spp.fl, by='spp') |> 
  relocate(Order, Family, Scientific.Name) |> 
  mutate(overlap = factor(overlap, levels = unique(overlap)),
         color_flag = ifelse(p.D > 0.05, ">0.05", "<0.05"),
         spp = factor(spp, levels=eca.fct),
         overlap2 = str_replace_all(overlap, c('m1_' = 'Ene_', 
                                               'm2' = 'Feb', 
                                               'm3' = 'Mar',
                                               'm4' = 'Abr',
                                               'm5' = 'May',
                                               'm6' = 'Jun',
                                               'm7' = 'Jul',
                                               'm8' = 'Ago',
                                               'm9' = 'Sep',
                                               'm10' = 'Oct',
                                               'm11' = 'Nov',
                                               'm12' = 'Dic',
                                               '_' = ' vs '))) |>  
  mutate(overlap2 = factor(overlap2, levels = unique(overlap2))) |> 
  # Solo graficar 6 especies, de manera ilustrativa y ajustarles ese orden
  filter(spp==eca.fct[1] | spp==eca.fct[2] | spp==eca.fct[9] | 
           spp==eca.fct[16] | spp==eca.fct[17] | spp==eca.fct[18]) |> 
  mutate(Scientific.Name = factor(Scientific.Name, levels=sp.ord)) |> 
  ggplot(aes(x=obs.D, y=fct_rev(overlap2), color=color_flag)) +
  geom_point(size=2) +
  # geom_point(data=uni_ovl.df.m, aes(x=jaccard, y=overlap), size=2, shape = 2) + # Esta linea plotea el resultado de jaccard
  geom_vline(xintercept = 0.5, linetype ="dotted") +
  geom_vline(xintercept = c(0,1), linetype ="dashed") +
  scale_color_manual(values = c("<0.05" = "black",
                                ">0.05" = "red")) +
  xlim(0,1)+
  labs(x="D de Schoener", y="Superposición mensual", color="P-value")+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", '0.75', "1")) +
  facet_wrap(~Scientific.Name, ncol=19) +
  theme_test()

# # Jaccard para las estaciones
#   spp.list.s <- list.files('./species/mig_ENM/seasonaly/', full.names=T)[c(65,68:85)]
#   uni_ovl.s <- list()
#   
#   # x <- 1
#   for(x in 1:length(spp.list.s)) {
#      uni_ovl.s[[x]] <- 
#       read.table(paste0(spp.list.s[x], '/S_overlap/', basename(spp.list.s[x]) |> str_remove('_ENMs'), '_s_uni_ovl.txt'), sep='\t', dec='.', header=T) |> 
#       # select(spp, niche.comp, overlap, p_value) |> 
#       mutate(spp = basename(spp.list.s[x]) |> str_remove('_ENMs'),
#              niche.comp = c('S1_S2', 'S1_S3', 'S1_S4', 'S2_S3', 'S2_S4', 'S3_S4')) |> 
#       rename('overlap'= niche.comp,
#              'jaccard' = overlap)
#      }
#   
#   uni_ovl.df.s <- do.call(rbind, uni_ovl.s) |>
#     mutate(overlap=factor(overlap, levels = unique(overlap)),
#            color_flag=ifelse(p_value > 0.05, ">0.05", "<0.05"))


s<-
  data.ovrlp.s |>
  left_join(spp.fl, by='spp') |> 
  relocate(Order, Family, Scientific.Name) |> 
  mutate(overlap=factor(overlap, levels = unique(overlap)),
         color_flag=ifelse(p.D > 0.05, ">0.05", "<0.05"),
         spp=factor(spp, levels=eca.fct)) |> 
  # Solo graficar 6 especies, de manera ilustrativa y ajustarles ese orden
  filter(spp==eca.fct[1] | spp==eca.fct[2] | spp==eca.fct[9] | 
           spp==eca.fct[16] | spp==eca.fct[17] | spp==eca.fct[18]) |> 
  mutate(Scientific.Name = factor(Scientific.Name, levels=sp.ord)) |> 
  #  summarise(n=n(), .by=c(color_flag, spp))
  ggplot(aes(x=obs.D, y=fct_rev(overlap), color=color_flag)) +
  geom_point(size=2) +
  # geom_point(data=uni_ovl.df.s, aes(x=jaccard, y=overlap), size=2, shape = 2) +
  scale_color_manual(values = c("<0.05" = "black",
                                ">0.05" = "red")) +
  xlim(0,1)+
  geom_vline(xintercept = 0.5, linetype ="dotted") +
  geom_vline(xintercept = c(0,1), linetype ="dashed") +
  #  geom_vline(xintercept = 0.05, color="red", linetype ="dotted") +
  labs(x="D de Schoener", y="Superposición estacional", color="P-value")+
  scale_y_discrete(label= c("Ve-Ot","Pr-Ot","Pr-Ve", "In-Ot", "In-Ve", "In-Pr"))+ # va de abajo hacia arriba
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", '0.75', "1")) +
  #  scale_size(range = c(0, 4)) +  # ajusta rango de tamaños
  facet_wrap(~Scientific.Name, ncol=19)  +
  theme_test() +
  theme(legend.position = 'none')

#library(patchwork)
# m/s
s+m
# cuadrar: 

ggsave(file = "./outputs/images/s_m_ovl.png",
       width = 20,
       height = 10,
       scale=2.3,
       units ="cm")


# Environmental Coincidence Average ---------------------------------------

data.ovrlp.m |> 
  summarise(eca=mean(obs.D), sd.eca=sd(obs.D), .by=c(spp)) |> 
  arrange(eca) |> 
  mutate(spp = str_replace_all(spp, '_', '. ')) |> 
  ggplot(aes(y=eca, x=fct_reorder(spp, eca))) +
  geom_point(size=2) +
  # geom_bar(stat='identity',  fill = 'lightblue4')+
  geom_errorbar(aes(ymin = eca - sd.eca, ymax = eca + sd.eca), width = 0.5) +
  # geom_line(aes(x=orden, group=spp)) +
  # geom_hline(yintercept = 0.5, linetype ="dotted", col = 'red') +
  labs(x="Especies", y="Promedio de coincidencia ambiental intermensual") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

ggsave(file = "./outputs/images/ieca.png",
       width = 10,
       height = 7,
       scale=3,
       units ="cm")


# Pruebas con las familias de los bichos
spp.fl <- read.table('./inputs/final-spp-list.txt', header=T, sep='\t', dec='.') |> 
  separate(Scientific.Name, into = c('gen', 'sp'), sep = ' ', remove=F) |> 
  mutate(spp=paste(str_sub(gen, 1,3), str_sub(sp, 1,5), sep = '_')) |> 
  select(!c(Taxonomic.Order, gen, sp))

data.ovrlp.m |> 
  summarise(eca=mean(obs.D), sd.eca=sd(obs.D), .by=c(spp)) |> 
  arrange(eca) |> 
  left_join(spp.fl, by='spp') |> 
  relocate(Order, Family, Scientific.Name) |> 
  # summarise(n=n(), .by = Order)
  mutate(spp = str_replace_all(spp, '_', '. ')) |> 
  ggplot(aes(y=eca, x=fct_reorder(spp, eca), fill=Family)) +
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin = eca - sd.eca, ymax = eca + sd.eca), width = 0.5) +
  geom_hline(yintercept = 0.5, linetype ="dotted", col = 'red') +
  labs(x="", y="Coincidencia Ambiental Promedio (ECA)") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Correlaciones entre Jaccard y Shoenner ----------------------------------

uni_ovl.df.m
data.ovrlp.m

j.d.df <- data.ovrlp.m |> 
  as_tibble() |> 
  select(spp, overlap, obs.D, p.D) |> 
  rename('compar' = overlap,
         'overlap' = obs.D,
         'p_value' = p.D) |> 
  mutate(type = 'shoenner') |> 
  rbind(uni_ovl.df.m |> 
          as_tibble() |> 
          select(spp, overlap, jaccard, p_value) |> 
          rename('compar' = overlap,
                 'overlap' = jaccard) |> 
          mutate(type = 'jaccard'))

# Como se ve la relacion por especie
j.d.df |> 
  # filter(spp == 'Sel_calli') |> 
  select(!p_value) |> 
  pivot_wider(names_from = 'type', values_from = 'overlap') |> 
  # select(shoenner, jaccard) |> 
  ggplot(aes(x = shoenner, y = jaccard, color = spp)) +
  geom_point() +
  geom_smooth(method = "lm", se=F) +
  theme_classic() +
  labs(title = 'Month overlap (species)', color = 'Species')


j.d.df |> 
  # filter(spp == 'Sel_calli') |> 
  select(!p_value) |> 
  pivot_wider(names_from = 'type', values_from = 'overlap') |> 
  ggplot(aes(x = shoenner, y = jaccard, color = compar)) +
  geom_point() +
  geom_smooth(method = "lm", se=F) +
  labs(title = 'Month overlap (mx vs my)') +
  theme_classic() +
  theme(legend.position = 'none') 


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

rast.list <- list.files('./outputs/MOP_R2/Nea_Nea/', full.names = T)

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

