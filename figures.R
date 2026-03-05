### Figures and formal analysis

library(tidyverse)

# Broenniman results ------------------------------------------------------

list.m<-list()
list.s<-list()
df.m <- list.files("./species/mig_ENM/overlaps_tables/B_ovl_month/", full.names=T)
df.s <- list.files("./species/mig_ENM/overlaps_tables/B_ovl_season/", full.names=T)

for (x in 1:length(df.m)) {
  a <- read.table(df.m[x], sep="\t", dec = ".", header=T)
  list.m[[x]] <- a
  b <- read.table(df.s[x], sep="\t", dec = ".", header=T)
  list.s[[x]] <- b
}

data.ovrlp.m <- do.call(rbind, list.m) |> 
  tibble()
data.ovrlp.s <- do.call(rbind, list.s) |> 
  tibble()

# write.table(data.ovrlp.m, "./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt",
#             sep="\t", dec = ".", row.names=F)
# write.table(data.ovrlp.s, "./species/mig_ENM/overlaps_tables/data.ovrlp_S_df.txt",
#             sep="\t", dec = ".", row.names=F)

data.ovrlp.m<-read.table("./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt",
                         sep="\t", dec = ".", header=T)

data.ovrlp.s<-read.table("./species/mig_ENM/overlaps_tables/data.ovrlp_S_df.txt",
                         sep="\t", dec = ".", header=T)

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

esa.fct <- data.ovrlp.m |> 
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
         color_flag = ifelse(p.D > 0.05, "p>0.05", "p<0.05"),
         spp = factor(spp, levels=esa.fct),
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
  filter(spp==esa.fct[1] | spp==esa.fct[2] | spp==esa.fct[9] | 
           spp==esa.fct[16] | spp==esa.fct[17] | spp==esa.fct[18]) |> 
  mutate(Scientific.Name = factor(Scientific.Name, levels=sp.ord)) |> 
  ggplot(aes(x=obs.D, y=fct_rev(overlap2), color=color_flag)) +
  geom_point(size=2) +
  # geom_point(data=uni_ovl.df.m, aes(x=jaccard, y=overlap), size=2, shape = 2) + # Esta linea plotea el resultado de jaccard
  geom_vline(xintercept = 0.5, linetype ="dotted") +
  geom_vline(xintercept = c(0,1), linetype ="dashed") +
  scale_color_manual(values = c("p<0.05" = "black",
                                "p>0.05" = "red")) +
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
         color_flag=ifelse(p.D > 0.05, "p>0.05", "p<0.05"),
         spp=factor(spp, levels=esa.fct)) |> 
  # Solo graficar 6 especies, de manera ilustrativa y ajustarles ese orden
  filter(spp==esa.fct[1] | spp==esa.fct[2] | spp==esa.fct[9] | 
           spp==esa.fct[16] | spp==esa.fct[17] | spp==esa.fct[18]) |> 
  mutate(Scientific.Name = factor(Scientific.Name, levels=sp.ord)) |> 
  #  summarise(n=n(), .by=c(color_flag, spp))
  ggplot(aes(x=obs.D, y=fct_rev(overlap), color=color_flag)) +
  geom_point(size=2) +
  # geom_point(data=uni_ovl.df.s, aes(x=jaccard, y=overlap), size=2, shape = 2) +
  scale_color_manual(values = c("p<0.05" = "black",
                                "p>0.05" = "red")) +
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

# library(patchwork)
s/m
s+m

ggsave(file = "./outputs/images/s_m_ovl.svg",
       width = 20,
       height = 10,
       scale=2.3,
       units ="cm")


# ¿Para cuántas spp p<0.05?
data.ovrlp.s |> 
  filter(p.D<0.05) |>
  # filter(p.D<0.05) |>
  # filter(spp=="Car_canad")# |>
  select(spp) |> 
  unique() 

spp.sum <-
  data.ovrlp.s |> 
  mutate(type="season",
         spp= spp |> str_replace("_", ". ")) |> 
  rbind(data.ovrlp.m |> 
          mutate(type="month", 
                 spp= spp |> str_replace("_", ". "))) |> 
  select(spp, overlap, obs.D, p.D, type) |> 
  summarise(n=sum(p.D < 0.05, na.rm = TRUE), sd=sd(obs.D), .by=c(spp, type)) |> 
  mutate(len = case_when(type == "season" ~ 6,
                         type == "month" ~ 66),
         prop = n/len) |> tibble()

# ord.spp <- 
  spp.sum |> 
  filter(type=="season") |> 
  filter(n==6) |>
    print(n=30)
  arrange(n)# |> 
  pull(spp)

spp.sum |> 
  mutate(spp =  factor(spp, levels=rev(ord.spp))) |> 
  ggplot(aes(y=spp, x=prop)) +
  geom_bar(stat="identity") +
  facet_wrap(~factor(type, levels=c("season", "month")), 
             labeller = as_labeller(c("season"="Estación","month"="Mes"))) +
  labs(y="", x="Proporción de pares de nichos diferentes") +
  theme_classic() 

ggsave(file = "./outputs/images/ovals_num.png",
       width = 10,
       height = 15,
       scale=2,
       units ="cm")


# Anual Environmental Coincidence Average ---------------------------------------
library(harmonicmeanp)

## P-value harmonic mean  -------------------------------------------------------

l <- list()
nam <- data.ovrlp.m$spp |> unique()
sp.df.m <- data.frame(spp=nam, p.m=NA)
sp.df.s <- data.frame(spp=nam, p.s=NA)

for(x in 1:103){
  p.m <- data.ovrlp.m |> 
    filter(spp==nam[x]) |> 
    pull(p.D)
  sp.df.m[x,2] <- harmonicmeanp::p.hmp(p.m, L=length(p.m))[1]
  
  p.s <- data.ovrlp.s |> 
    filter(spp==nam[x]) |> 
    pull(p.D)
  sp.df.s[x,2] <- harmonicmeanp::p.hmp(p.s, L=length(p.s))[1]
}

s.ps <- data.ovrlp.s |> 
  summarise(esa=mean(obs.D), 
            sd.esa=sd(obs.D), .by=c(spp)) |> 
  left_join(sp.df.s, by="spp") |> 
  mutate(spp = str_replace_all(spp, '_', '. '),
         p.cod = case_when(
           p.s < 0.05 ~ "p<0.05",
           p.s > 0.05 ~ "p>0.05")) 

esa.df <- data.ovrlp.m |> 
  summarise(esa=mean(obs.D), 
            sd.esa=sd(obs.D), .by=c(spp)) |> 
  left_join(sp.df.m, by="spp") |> 
  arrange(esa) |> 
  mutate(p.cod = case_when(
    p.m < 0.05 ~ "p<0.05",
    p.m > 0.05 ~ "p>0.05")) |> 
  tibble()

# esa.df |>  write.table('./species/mig_ENM/overlaps_tables/esa.df.txt', sep="\t", dec = ".", row.names=F)

esa.g <-
  esa.df |> 
  mutate(spp = str_replace_all(spp, '_', '. ')) |> 
  ggplot(aes(y=esa, x=fct_reorder(spp, esa), color=p.cod)) +
  geom_point(size=2) +
  # geom_point(data=s.ps, aes(y=esa, x=fct_reorder(spp, esa), color=p.cod), shape=5) +
  # geom_bar(stat='identity',  fill = 'lightblue4')+
  geom_errorbar(aes(ymin = esa - sd.esa, ymax = esa + sd.esa), width = 0.5) +
  scale_color_manual(values = c("p<0.05" = "black",
                                "p>0.05" = "red")) +
  # geom_line(aes(x=orden, group=spp)) +
  # geom_hline(yintercept = 0.5, linetype ="dotted", col = 'red') +
  labs(x="", y="Environmental Similarity Average", color="P-value(c)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file = "./outputs/images/iesa.svg",
       width = 10,
       height = 7,
       scale=3,
       units ="cm")


# Jackard and intersection volume -----------------------------------------

# This analysis was the result of ellipsoid comparison using the overllip package (Osorio-Olvera, 2020). This analysis calculates a Jackard index and the volume of the ellipsoid's intersection between all ellipsoids. Seasonal analysis was removed from this part and only the monthly results have been used.

jack <- read.table("./species/mig_ENM/overlaps_tables/jack_res_Tm.txt", sep="\t", dec=".", header=T)
int <- read.table("./species/mig_ENM/overlaps_tables/int_vols_Tm.txt", sep="\t", dec=".", header=T)

spp.order <- esa.df |> 
  tibble() |> 
  arrange(esa) |> 
  select(spp) |> 
  pull()

# int |> 
jack |> 
  select(1:13) |> 
  tibble() |> 
  pivot_longer(cols = !species, names_to = 'months', values_to = 'values') |> 
  mutate(months = factor(months, levels = c("total_vs_Dec", "total_vs_Jan","total_vs_Feb", "total_vs_Mar", "total_vs_Apr", "total_vs_May", "total_vs_Jun", "total_vs_Jul", "total_vs_Aug", "total_vs_Sep", "total_vs_Oct", "total_vs_Nov")),
         species = factor(species, levels=spp.order |> str_replace('. ', '_'))) |> 
  ggplot(aes(x=species, y=values, color=months)) +
  geom_point() +
  scale_color_manual(values=c("#4a72b0", "#4a72b0", "#4a72b0", "#6ea96e", "#6ea96e", "#6ea96e", "#da9500", "#da9500", "#da9500", "#294029","#294029","#294029")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

jack |> 
  select(1:13) |> 
  tibble() |> 
  pivot_longer(cols = !species, names_to = 'months', values_to = 'values') |> 
  summarise(mean.j=mean(values), sd.j=sd(values), .by=species) |> 
  rename('spp' = species) |> 
  left_join(esa.spp, by='spp') |> 
  mutate(spp = factor(spp, levels=spp.order |> str_replace('. ', '_'))) |> 
  ggplot(aes(x=spp, y=mean.j)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean.j - sd.j, ymax = mean.j + sd.j), width = 0.5) +
  geom_point(aes(x=spp, y=esa, color=p.cod), shape=17) +
  geom_errorbar(aes(ymin = esa - sd.esa, ymax = esa + sd.esa), width = 0.5, linetype='dashed') +
  scale_color_manual(values = c("p<0.05" = "black",
                                "p>0.05" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Mensual proportion of the total environment -----------------------------------

# This is the result of ellipsoid overlap using the ellipsenm package. The Overlap column represents the environmental shared in proportion with on environment of reference Seasonal analysis was removed from this part and only the monthly results have been used.

ovl.m <- list.files('./species/mig_ENM/overlaps_tables/prop_ovl_month/', full.names = T)
# ovl.s <- list.files('./species/mig_ENM/overlaps_tables/prop_ovl_season/', full.names = T)
spp.order

olv.m.list <- list()
# olv.s.list <- list()

# x <- 1
for(x in 1:length(ovl.m)){
  
  olv.m.list[[x]] <- read.table(ovl.m[x],
                                sep = "\t",
                                dec = ".",
                                header = T)
  # olv.s.list[[x]] <- read.table(ovl.s[x],
  #                               sep = "\t",
  #                               dec = ".",
  #                               header = T)
}


olv.m.df <- do.call(rbind, olv.m.list)
# olv.s.df <- do.call(rbind, olv.s.list)

olv.m.g <-
  olv.m.df |> 
  tibble() |> 
  filter(str_starts(ovl.type, 'Niche_1_')) |> # Niche 1 es el nicho completo (total) 
  select(spp, ovl.type, overlap, prop_size_niche_2_vs_1) |> 
  mutate(spp = spp |> str_replace('_', '. ')) |> 
  mutate(ovl.type = factor(ovl.type, levels = c("Niche_1_vs_13", "Niche_1_vs_2", "Niche_1_vs_3", "Niche_1_vs_4", "Niche_1_vs_5", "Niche_1_vs_6", "Niche_1_vs_7", "Niche_1_vs_8", "Niche_1_vs_9", "Niche_1_vs_10", "Niche_1_vs_11", "Niche_1_vs_12")),
         spp = factor(spp, levels=spp.order |> str_replace('_', '. '))) |>
  summarise(mean=mean(overlap), sd=sd(overlap), .by=spp) |> 
  ggplot(aes(x=spp, y=mean)) +
  # ggplot(aes(x=spp, y=overlap, color=ovl.type)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5) +
  theme_minimal() +
  labs(x='Species',  y='Mensual average proportion of the Grinnellian niche') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# y= Mensual average proportion of the total Grinnellian niche

# library(patchwork)

esa.g /olv.m.g

ggsave(file = "./outputs/images/iesa_prom.jpg",
       width = 10,
       height = 7,
       scale=3,
       units ="cm")

# Linear model (ESA ~ Overlap) --------------------------------------------

# esa.df <- read_tsv('./species/mig_ENM/overlaps_tables/esa.df.txt')
esa.df

ovl.mean <- olv.m.df |> 
  tibble() |> 
  filter(str_starts(ovl.type, 'Niche_1_')) |> # Niche 1 es el nicho completo (total) 
  select(spp, ovl.type, overlap, prop_size_niche_2_vs_1) |> 
  summarise(mean.ovl=mean(overlap), sd.ovl=sd(overlap), .by=spp)

olv.esa.df <- ovl.mean |> 
  left_join(esa.df, by='spp')

m <- lm(olv.esa.df$esa~olv.esa.df$mean.ovl)
summary(m)
plot(m)

# Spearman correlation
cor(olv.esa.df[,c(2,4)], method='spearman')

# Normality
ks.test(m$residuals, rnorm(length(m$residuals), 0, sd(m$residuals)))

# Homocedasticity
# puede calcularse a través de elevar al cuadrado los residuales y ajustarlos con los valores del modelo
lm(m$residuals^2~m$fitted.values) |> summary()
# elvalor de la pendiente es p-value: 0.3057 por lo tanto no hay evidencia suficiente para rechazar la hipotesis nula en donde asume que no hay diferencias entre las varianzas. Por lo tanto hay homocedasticidad



# Broenniman graphs -------------------------------------------------------

# months
pca.path <- list.files("./species/mig_ENM/overlaps_tables/PCA_m/", full.names = T)
pca.names <- basename(pca.path) |> str_remove_all("_m_PCA.txt")

# seasons
# pca.path <- list.files("./species/mig_ENM/overlaps_tables/PCA_s/", full.names = T)
# pca.names <- basename(pca.path) |> str_remove_all("_s_PCA.txt")
# basename(pca.path[1])

# x <- 39
for (x in 1:length(pca.path)) {
  
spp <- read.table(pca.path[x],header=T, sep='\t', dec='.') |> 
  filter(str_starts(elip, "abex_"))


# #---
# # Estos graficos de densidad calculan la probabilidad y el 0.1 es la mas externa. Agrupa esos valores y realiza un CH. por eso queda mocho
# spp |> 
#   filter(elip =="abex_Lei_criss_3") |> 
#   ggplot(aes(x=Axis1, y=Axis2)) + 
#   geom_point() +
#   stat_density2d(geom="polygon", alpha=0.10, na.rm = T, contour_var = "ndensity", breaks = c(seq(0.1,0.9, by=0.05)))
# 
# #---

pca.fig <-
  ggplot()+
  stat_density2d(data=spp,
                 # filter(elip == paste0("abex_", pca.names[x], "_1") | elip == paste0("abex_", pca.names[x], "_2")),
                 aes(x = Axis1,y = Axis2, color=elip, fill=elip), 
                 geom="polygon", alpha=0.10, na.rm = T, contour_var = "ndensity")+
    theme_bw()+
    labs(title=pca.names[x] |> str_replace("_", ". "),x="PC1", y="PC2")+
    geom_hline(yintercept=0, linetype="dashed", color = "gray30", linewidth=0.5,alpha=0.70)+
    geom_vline(xintercept=0, linetype="dashed", color = "gray30", linewidth=0.5, alpha=0.70)+
  
  # Months
  scale_fill_manual(values=c("#4a72b0", "#4a72b0", "#6ea96e", "#6ea96e", "#6ea96e", "#da9500", "#da9500", "#da9500", "#294029","#294029","#294029", "#4a72b0")) +
    scale_color_manual(values=c("#4a72b0", "#4a72b0", "#6ea96e", "#6ea96e", "#6ea96e", "#da9500", "#da9500", "#da9500", "#294029","#294029","#294029", "#4a72b0")) +
  
# # Seasons
#     scale_fill_manual(values=c("#4a72b0", "#6ea96e", "#da9500", "#294029")) +
#     scale_color_manual(values=c("#4a72b0", "#6ea96e", "#da9500", "#294029")) +
    theme(legend.position="none")

ggsave(pca.fig,
       file = paste0("./outputs/images/PCA/M_", pca.names[x] ,"_PCA.png"),
       width = 5,
       height = 5,
       scale=3,
       units ="cm")

print(paste(basename(pca.path[x]) |> str_remove("_s_PCA.txt"), x, "done"))

}
  

  #....
dens <- 
  ggplot( data = spp |> filter(str_starts(elip, "abex_")),
  aes(x = Axis1, y = Axis2, color=elip, fill=elip)) +
  stat_density_2d( aes(fill = after_stat(level)), geom = "polygon")

dens_data <- ggplot_build(dens)$data[[1]]

outer_poly <- dens_data %>%
  filter(level == min(level))

ggplot() +
  geom_polygon(
    data = outer_poly,
    aes(x, y, group = group),
    fill = "steelblue",
    alpha = 0.15)

  #....

 # Pruebas con las familias de los bichos
spp.fl <- read.table('./inputs/final-spp-list.txt', header=T, sep='\t', dec='.') |> 
  separate(Scientific.Name, into = c('gen', 'sp'), sep = ' ', remove=F) |> 
  mutate(spp=paste(str_sub(gen, 1,3), str_sub(sp, 1,5), sep = '_')) |> 
  select(!c(Taxonomic.Order, gen, sp))

data.ovrlp.m |> 
  summarise(esa=mean(obs.D), sd.esa=sd(obs.D), .by=c(spp)) |> 
  arrange(esa) |> 
  left_join(spp.fl, by='spp') |> 
  relocate(Order, Family, Scientific.Name) |> 
  # summarise(n=n(), .by = Order)
  mutate(spp = str_replace_all(spp, '_', '. ')) |> 
  ggplot(aes(y=esa, x=fct_reorder(spp, esa), fill=Family)) +
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin = esa - sd.esa, ymax = esa + sd.esa), width = 0.5) +
  geom_hline(yintercept = 0.5, linetype ="dotted", col = 'red') +
  labs(x="", y="Coincidencia Ambiental Promedio (esa)") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Broenniman PCA graphs 2 -------------------------------------------------

# Totales vs Mensual
library(tidyverse)
library(geometry) # Para calcular áreas de polígonos convexos
library(patchwork)

pca.list.s <- list.files("./species/mig_ENM/overlaps_tables/PCA_s/", full.names = T)
pca.list.m <- list.files("./species/mig_ENM/overlaps_tables/PCA_m/", full.names = T)

x <- 1
for (x in 1:length(pca.list.s)) {
spp.pca.s <- read_tsv(pca.list.s[x]) |> 
  filter(str_starts(elip, 'abex')) |> 
  mutate(elip2 = elip |> str_remove(paste0("abex_", basename(pca.list.s[x]) |> 
                                             str_remove("s_PCA.txt")))) |> 
  mutate(elip2 = case_when(
    elip2 == 1 ~ "W",
    elip2 == 2 ~ "T1",
    elip2 == 3 ~ "S",
    elip2 == 4 ~ "T2")) |> 
  mutate(elip2=factor(elip2, levels=c("W", "T1", "S", "T2")))

spp.pca.m <- read_tsv(pca.list.m[x]) |> 
  filter(str_starts(elip, 'abex')) |> 
  mutate(elip2 = elip |> str_remove(paste0("abex_", basename(pca.list.m[x]) |> 
                                             str_remove("m_PCA.txt")))) |> 
  mutate(elip2 = case_when(
    elip2 == 1 ~ "Jan",
    elip2 == 2 ~ "Feb",
    elip2 == 3 ~ "Mar",
    elip2 == 4 ~ "Apr",
    elip2 == 5 ~ "May",
    elip2 == 6 ~ "Jun",
    elip2 == 7 ~ "Jul",
    elip2 == 8 ~ "Aug",
    elip2 == 9 ~ "Sep",
    elip2 == 10 ~ "Oct",
    elip2 == 11 ~ "Nov",
    elip2 == 12 ~ "Dec")) 

spp.g.s <-
  ggplot(spp.pca.s, aes(x = Axis1, y = Axis2)) +
  # El fondo con todos los puntos en gris
  geom_point(alpha = 0.2,  aes(shape = elip2)) +
  # Los polígonos por mes (puedes filtrar por uno específico o usar facets)
  geom_polygon(data = spp.pca.s |>  
                 group_by(elip2) |>  
                 slice(chull(Axis1, Axis2)),
               aes(color = factor(elip2), fill= factor(elip2)), alpha = 0.1) +
  geom_polygon(data = spp.pca.s |> 
                 slice(chull(Axis1, Axis2)), 
               alpha = 0.4, fill=NA, color="red", linetype="dashed") +
  scale_color_manual(values=c("#4a72b0", "#6ea96e", "#da9500", "#294029")) +
  scale_fill_manual(values=c("#4a72b0", "#6ea96e", "#da9500", "#294029")) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray30", linewidth=0.5,alpha=0.70)+
  geom_vline(xintercept=0, linetype="dashed", color = "gray30", linewidth=0.5, alpha=0.70)+
  theme_minimal() +
  labs(title = paste(basename(pca.list.s[x]) |> str_remove("_s_PCA.txt")),
       fill = "Season", color="Season", shape="Season")

spp.g.m <- 
  ggplot(spp.pca.m, aes(x = Axis1, y = Axis2)) +
  # El fondo con todos los puntos en gris
  geom_point(alpha = 0.2,  aes(shape = elip2)) +
  # Los polígonos por mes (puedes filtrar por uno específico o usar facets)
  geom_polygon(data = spp.pca.m |>  
                 group_by(elip2) |>  
                 slice(chull(Axis1, Axis2)),
               aes(color = factor(elip2), fill= factor(elip2)), alpha = 0.1) +
  geom_polygon(data = spp.pca.m |> 
                 slice(chull(Axis1, Axis2)), 
               alpha = 0.4, fill=NA, color="red", linetype="dashed") +
  scale_color_manual(values=c("#4a72b0","#4a72b0","#4a72b0", "#6ea96e","#6ea96e","#6ea96e", "#da9500","#da9500","#da9500", "#294029", "#294029", "#294029")) +
  scale_fill_manual(values=c("#4a72b0","#4a72b0","#4a72b0", "#6ea96e","#6ea96e","#6ea96e", "#da9500","#da9500","#da9500", "#294029", "#294029", "#294029")) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray30", linewidth=0.5,alpha=0.70)+
  geom_vline(xintercept=0, linetype="dashed", color = "gray30", linewidth=0.5, alpha=0.70)+
  scale_shape_manual(values=c(19,19,19, 17,17,17, 15,15,15, 3,3,3)) +
  theme_minimal() +
  labs(fill = "Month", color="Month", shape="Month")

g <- spp.g.s + spp.g.m

ggsave(g,
       file = paste0("./outputs/images/PCA_2/", basename(pca.list.s[x]) |> str_remove("_s_PCA.txt"),"_PCA.png"),
       width = 15,
       height = 10,
       scale=1.8,
       units ="cm")

print(paste(basename(pca.list.s[x]) |> str_remove("_s_PCA.txt"), "saved"))
}

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
# Install the tmap package
install.packages("tmap")

# Load the package
library(tmap)


rast.list <- list.files('./outputs/MOP_R2/Nea_Nea/', full.names = T, pattern=".tif")
rast.list2 <- list.files('./outputs/MOP_R2/Nea_Neo/', full.names = T, pattern=".tif")

# Vector con los nombres originales
rast.names <-  # la posicion 26 corresponde a la comparacion del las estaciones, la pongo a mano
  basename(rast.list[26]) |> 
  str_remove('MOP_realms_') |> 
  str_remove('.tif') 

# Vector con los meses
# month <- month.name
month

# Extraer los números antes y después del "_"
num1 <- as.numeric(sub("m(\\d+)_.*", "\\1", rast.names) |> str_remove_all("Nea_Nea_"))
num2 <- as.numeric(sub(".*_(\\d+)", "\\1", rast.names))

# Crear nuevo vector con nombres de los meses
nuevo <- paste(month[num1], month[num2], sep = " (Neart) vs ")
nuevo2 <- c(nuevo, 'Winter (Neart) vs Summer')

x <- 1
for (x in 1:length(rast.list)) {
  
  mop <- rast(rast.list[x])
  mop2 <- rast(rast.list2[x])
  names(mop2) <- 'layer.1'
  names(mop) <- 'layer.1'
  plot(mop)  


  b <-
    ggplot() +
      geom_raster(data = mop, aes(x = x, y = y, fill=layer.1)) + 
      coord_equal() +
    # geom_spatraster(data = mop, aes(x = x, y = y, fill=layer.1))# +
    # coord_sf(datum = st_crs(4326)) 
    scale_fill_viridis_c(na.value = 'transparent') + 
    # labs(title = paste(nuevo2[x], '(Neotr)'), fill = 'MOP', x="Longitude", y="Latitude") + 
    labs(fill = 'MOP', x="Longitude", y="Latitude") + 
    theme_classic()
  
  # c <-
  #   ggplot() +
  #     geom_raster(data = mop2, aes(x = x, y = y, fill=layer.1)) + 
  #     coord_equal() +
  #   # geom_spatraster(data = mop, aes(x = x, y = y, fill=layer.1))# +
  #   # coord_sf(datum = st_crs(4326)) 
  #   scale_fill_viridis_c(na.value = 'transparent') + 
  #   labs(fill = 'MOP', x="Longitude", y="Latitude") + 
  #   theme_classic()
  
  
  # d <- b/c + plot_annotation(tag_levels = 'A')
  
  ggsave(b, file = paste0("./outputs/MOP_R2/", basename(rast.list[x]) |> str_remove('.tif'),".png"),
         width = 1706,
         height = 1271,
         scale=1,
         units ="px")
  
  # ggsave(d, file = "./outputs/MOP_R2/MOP_graphs2.png",
  #        width = 10,
  #        height = 15,
  #        scale=1,
  #        units ="cm")
}



# Amplitudes --------------------------------------------------------------
library(tidyverse)

# Rangos del PCA
pca.paths.m <- list.files("./species/mig_ENM/overlaps_tables/PCA_m/", full.names=T)
pca.paths.s <- list.files("./species/mig_ENM/overlaps_tables/PCA_s/", full.names=T)

pca.m <- list()
pca.s <- list()

for (x in 1:length(pca.paths.m)) {
  pca.m[[x]] <- 
    read_tsv(pca.paths.m[x]) |> 
    filter(str_starts(elip,"abex_"))
  pca.s[[x]] <- 
    read_tsv(pca.paths.s[x]) |> 
    filter(str_starts(elip,"abex_"))
  
  print(paste(basename(pca.paths.s[x] |> str_remove("_s_PCA.txt")), x, "ready"))
}

pca.df.m <- 
  do.call(rbind, pca.m)

pca.df.s <- 
  do.call(rbind, pca.s)

amp.m<- pca.df.m |> 
  summarise(min.a1=min(Axis1), 
            max.a1=max(Axis1), 
            min.a2=min(Axis2), 
            max.a2=max(Axis2), 
            .by=elip) |> 
  extract(elip, 
          into = c("spp", "month"),
          regex = "(.*)_([0-9]+$)") |> 
  mutate(spp = spp |> str_remove_all("abex_"),
         cod = case_when(
           month == 1 ~ "In", 
           month == 2 ~ "In", 
           month == 3 ~ "In", 
           month == 4 ~ "T1", 
           month == 5 ~ "T1", 
           month == 6 ~ "T1", 
           month == 7 ~ "Ve", 
           month == 8 ~ "Ve", 
           month == 9 ~ "Ve", 
           month == 10 ~ "T2", 
           month == 11 ~ "T2", 
           month == 12 ~ "T2")) |> 
  summarise(min.a1.tot=min(min.a1),
            max.a1.tot=max(max.a1),
            min.a2.tot=min(min.a2),
            max.a2.tot=max(max.a2),
            .by=c(spp, cod)) |>
  # Calculo de la amplitud para cada temporada por especie
  mutate(
    # Primero detectamos los extremos de cada fila
    min_extremo = pmin(min.a1.tot, min.a2.tot),
    max_extremo = pmax(max.a1.tot, max.a2.tot),
    # min_extremo = pmin(min.a1, min.a2),
    # max_extremo = pmax(max.a1, max.a2),
    # Calculamos la diferencia
    amplitud = max_extremo - min_extremo) |> 
  select(spp, cod, min_extremo, max_extremo, amplitud)


amp.m |> 
  group_by(spp) |> 
  # filter(amplitud == max(amplitud)) |>
  filter(amplitud == min(amplitud)) |>
  ungroup() |> 
  summarise(n=n()/103*100, n2=n(), .by=cod) |> 
  arrange(-n)

amp.m |> 
  group_by(spp) |> 
  filter(amplitud == max(amplitud)) |> 
  select(spp, cod, amplitud) |> 
  print(n=110)

amp.m |> 
  ggplot(aes(x=spp, y=amplitud, fill=cod)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#4a72b0", "#6ea96e", "#294029", "#da9500")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="", y="Amplitud", fill="Estación") 

ggsave(file = "./outputs/images/amplitud.svg",
       width = 20,
       height = 15,
       scale=1.5,
       units ="cm")


# Revision ----------------------------------------------------------------
library(tidyverse)

rev <- read_tsv('./outputs/tablas/review_fig.txt')
rev.comp <- read_tsv('./outputs/tablas/review_fig_complete.txt', locale = locale(encoding = "Latin1"))

rev |> 
  pivot_longer(cols = !Especies, names_to = 'papers', values_to = 'type') |> 
  na.omit() |> 
  ggplot(aes(x=Especies, y=papers, color=type, shape=type)) +
  geom_point(size=3) +
  scale_color_manual(values = c('aquamarine3', 'coral3', 'orange'))+
  theme_classic() + 
  labs(y='', color='Tipo de migracion', shape='Tipo de migracion') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='italic'),
        axis.text.y = element_text(face='bold'))
  
ggsave(file = "./outputs/images/review.png",
       width = 20,
       height = 15,
       scale=1.4,
       units ="cm")

rev.comp |> 
  rename('type.sum'=`Migration type summary`,
         'papers'=paper,
         'type'=`Migration type`) |>
  mutate(type.sum=case_when(
    type.sum == 'Other' ~ 'Otro',
    type.sum == 'Switcher' ~ 'Permutadora',
    type.sum == 'Tracker' ~ 'Seguidora',
    type.sum == 'No clasificada' ~ 'No clasificada')) |> 
  ggplot(aes(x=Species, y=papers, color=type.sum, shape=type.sum)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(4, 19, 17, 15)) + # Asigna forma a cada grupo
  scale_color_manual(values = c('gray8','aquamarine3', 'coral3', 'orange'))+
  theme_classic() + 
  labs(y='', color='Tipo de migracion', shape='Tipo de migracion') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face='italic'),
        axis.text.y = element_text(face='bold'))


ggsave(file = "./outputs/images/review_com.png",
       width = 20,
       height = 10,
       scale=1.8,
       units ="cm")



