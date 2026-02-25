


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

## Cualcular tipo de migracion segun Cohen et al. 2023
data <- read.table('./borrar/cohenetal2023_appendix_1.txt', sep='\t', dec='.', header=T)
head(data)


data |> 
  ggplot(aes(x=TD_Adj.Maha.dist, y=TD_Adj.Det.ratio)) +
  geom_point() +
  geom_hline(yintercept = median(data$TD_Adj.Det.ratio)) +
  geom_vline(xintercept = median(data$TD_Adj.Maha.dist))


data |> 
  rename('N.pos.simi' = TD_Adj.Maha.dist,
         'N.bre.simi' = TD_Adj.Det.ratio,) |> 
  select(Species, N.pos.simi, N.bre.simi) |> 
  mutate(mig.clas = case_when(
    N.pos.simi < median(N.pos.simi) & N.bre.simi < median(N.bre.simi) ~ 'Switcher',
    N.pos.simi < median(N.pos.simi) & N.bre.simi > median(N.bre.simi) ~ 'Bre.Tracker',
    N.pos.simi > median(N.pos.simi) & N.bre.simi < median(N.bre.simi) ~ 'Pos.Tracker',
    N.pos.simi > median(N.pos.simi) & N.bre.simi > median(N.bre.simi) ~ 'Tracker')) |> 
  arrange(Species) |> 
  filter(str_detect(Species, 'Vireo'))

# Comparación de nichos con Broennimann et al. 2012, Di Cola et al. 2017 

# Esta es una nueva version, en donde se calcula el PCA para los registros de toda la especie y se optienen los valores correspondientes a cada mes. Previamente, había calculado un valor de PCA por elipsoide y al calcular el ECI los valores de D no eran comparables.

# Librerias
library(ecospat)
library(terra)
library(tidyverse)

# Path de las especies
month.names<-list.files("./species/mig_ENM/monthly/", full.names = F) |> str_remove("_ENMs")
month.folder<-list.files("./species/mig_ENM/monthly/", full.names = T)

# Tabla que organiza el bucle (combinación 1 entre 12)
combs.m<-data.frame(a= rep(1:11, 11:1), b= unlist(lapply(2:12, function(i) i:12)))

data.ovrlp.list<-list()
# f<-65
for (f in 1:length(month.names)) { #1:length(month.names)
  
  # Cargo el Rdata que contiene los ENM y el background ---------------------
  
  temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
  load(list.files(month.folder[f], pattern="*.Rdata$", full.names = T), envir = temp_env)
  
  # Objetos de temp_env que se van a usar
  bg_swd.list<-temp_env$bg_swd.list
  abex.list<-temp_env$abex.list
  
  dir.create(paste0(month.folder[f], "/M_overlap"))
  
  data.overlap <- data.frame(
    spp = NA,
    overlap=NA,
    obs.D = NA,
    obs.I = NA,
    obs.expan = NA,
    obs.stabi = NA,
    obs.unfil = NA,
    p.D = NA,
    p.I = NA,
    p.expan = NA,
    p.stabi = NA,
    p.unfil = NA)
  
  # PCA total (toda la especie) ---------------------------------------------

    # Background
  bgd <- 
    do.call(rbind, bg_swd.list) |> 
    select(lon, lat, prec, tmax, tmin) |> 
    rownames_to_column(var = "elip") |> 
    mutate(elip = elip |> str_sub(1,19) |> str_remove_all("\\."))
    
    # Ambiente
  env2 <- list()
  
   # Extracting environmental values from temporal data frame
  for(p in 1:12){ # 1 to 12 (months / length of analysis (12))
    env2[[p]] <- abex.list[[p]]$temporal_df |> 
      mutate(elip = names(abex.list)[p])
  }
  
  env <-
    do.call(rbind, env2) |> 
    select(lon, lat, prec, tmax, tmin, elip)
    
#  Bind background and environmental dataframes
  data <- rbind.data.frame(bgd[,c(1,4:6)], env[,c(6,3:5)])
  
  # Weight vector. Occurrences= 0 and background (survey sites)=1
  w <- c(rep(1, nrow(bgd)), rep(0, nrow(env)))
  
  pca.cal <- ade4::dudi.pca(data[,2:4],
                            row.w = w,
                            center = T,
                            scale = T,
                            scannf = F,
                            nf = 2) # Produce solo dos PC
  
  # Ellipsoid ID  
  pca.cal2 <- pca.cal$li |> 
    mutate(elip = data$elip)
  
  # Extraction ellipsoid names
  
  bg_name <-  unique(pca.cal2$elip)[1:12] # It will be different for months 
  abex_name <-  unique(pca.cal2$elip)[13:24] # It will be different for months
 
  # Select PCA values for each ellipsoid
  
  # x<-1
  for (x in 1:nrow(combs.m)) {

# PCA climatic scores  ----------------------------------------------------
    
  # Backgrounds
  scores.clima1 <- pca.cal2 |> 
    filter(elip == bg_name[combs.m[x,1]]) |> 
    select(!elip)
  
  scores.clima2 <- pca.cal2 |> 
    filter(elip == bg_name[combs.m[x,2]]) |> 
    select(!elip)
  
  scores.clima12 <- rbind(scores.clima1, scores.clima2) |> 
    add_row(Axis1=min(pca.cal2$Axis1), Axis2=min(pca.cal2$Axis2)) |> 
    add_row(Axis1=max(pca.cal2$Axis1), Axis2=max(pca.cal2$Axis2))

  # Environments
  scores.sp1a <- pca.cal2 |> 
    filter(elip == abex_name[combs.m[x,1]]) |> 
    select(!elip)
  
  scores.sp2b <- pca.cal2 |> 
    filter(elip == abex_name[combs.m[x,2]]) |> 
    select(!elip)
  
# 
#  # Variable contribution to PCA analysis
#  contribucion<- ecospat.plot.contrib(contrib=pca.cal$co,
#                                         eigen = pca.cal$eig)
    
    
    # Superficie de densidad de registros -------------------------------------
    
    z1<-ecospat.grid.clim.dyn(glob=scores.clima12, 
                              glob1=scores.clima1, 
                              sp=scores.sp1a, 
                              R = 100)
    # head(z1)
    
    z2<-ecospat.grid.clim.dyn(glob=scores.clima12,
                              glob1=scores.clima2, 
                              sp=scores.sp2b,
                              R = 100)
    # head(z2)
    
    # PCA individual Graphs ---------------------------------------------------------
    
    # First ellipsoid
    
    png(paste0(month.folder[f], "/M_overlap/", 
               month.names[f], "_PCA_", "m", combs.m[x,1] ,".png"),
        width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche (z1,
                        title = paste0(month.names[f] |> str_replace("_", " "), 
                                       "-m", combs.m[x,1]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
     dev.off()
    
    # Second ellipsoid
    png(paste0(month.folder[f], "/M_overlap/", 
               month.names[f], "_PCA_", "m", combs.m[x,2] ,".png"), 
        width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche (z2,
                        title = paste0(month.names[f] |> str_replace("_", " "), 
                                       "-m", combs.m[x,2]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
    dev.off()
    

    # Overlap graph -----------------------------------------------------------

    png(paste0(month.folder[f], "/M_overlap/", month.names[f], 
               "_B.ovl_", "m", combs.m[x,1], "m", combs.m[x,2]  ,".png"), width = 800, height = 800, res = 100)
    
    par(mar=c(5,3,2,2))
    ecospat.plot.niche.dyn(z1,
                           z2,
                           quant = 0.25,
                           interest = 1,
                           title = paste0("Niche overlap ", month.names[f] |> str_replace("_", " "), 
                                          "-m", combs.m[x,1], "_", "m", combs.m[x,2]),
                           name.axis1 = "PC1",
                           name.axis2 = "PC2")
    dev.off()
    
    # D and I values:
    # ecospat.niche.overlap(z1=z1,z2=z2,cor=TRUE)
    
    # Similarity test (di Cola et al 2017 - higher) ----------------------
    
    sim.test <- ecospat.niche.similarity.test(
      z1,
      z2,
      rep = 100,
      overlap.alternative = "higher",
      rand.type = 1)
    
    data.overlap[x,1] <- month.names[f]
    data.overlap[x,2] <- paste0("M", combs.m[x, 1], "_M", combs.m[x, 2])
    data.overlap[x, 3:7] <- unlist(sim.test$obs)
    data.overlap[x,8] <- sim.test$p.D
    data.overlap[x,9] <- sim.test$p.I
    data.overlap[x,10] <- sim.test$p.expansion
    data.overlap[x,11] <- sim.test$p.stability
    data.overlap[x,12] <- sim.test$p.unfilling
    
    
    ## Plot Similarity test -  evalua si la similaridad es mas similar de la esperada por azar
    
    png(paste0(month.folder[f], "/M_overlap/", month.names[f], "_ST_", "m", combs.m[x,1], "m", combs.m[x,2]  ,".png"), width = 800, height = 800, res = 100)
    
    ecospat.plot.overlap.test(sim.test, "D", paste0("Similarity greater rand.type 1", "\n", 
                                                    month.names[f] |> str_replace("_", " "), 
                                                    "-m",combs.m[x,1],"_m", combs.m[x,2]))
    
    dev.off()
    
    print(paste(basename(month.folder[f]), combs.m[x,1], "v", combs.m[x,2], "ready"))
  }
  
  print(paste(basename(month.folder[f]), f, "done"))
  
  write.table(data.overlap,
              paste0("./species/mig_ENM/overlaps_tables/month/", month.names[f], "_m_ovl.txt"),
              sep="\t", dec = ".", row.names=F)
  
  
  data.ovrlp.list[[f]]<-data.overlap
  
  rm(list = setdiff(ls(), c("month.names", "month.folder", "combs.s", "data.ovrlp.list")))
}
# hay que hacer los graficos de todos los elipsoides juntos en esta forma de kernells


# En caso de haber borrado el objeto data.ovrlp.list
df.path <- list.files("./species/mig_ENM/overlaps_tables/season/", full.names = T)

for (x in 1:length(df.path)) {
  data.ovrlp.list[[x]] <- read.table(df.path[x], sep="\t", dec=".", header=T)
}

data.ovrlp.df<-do.call(rbind, data.ovrlp.list)

data.ovrlp.df|> head()

# write.table(data.ovrlp.df, "./species/mig_ENM/overlaps_tables/data.ovrlp_s_df.txt",
#             sep="\t", dec = ".", row.names=F)

data.ovrlp.df |> 
  mutate(overlap=factor(overlap, levels = unique(overlap)),
         color_flag=ifelse(p.D > 0.05, ">0.05", "<0.05")) |> 
  #  summarise(n=n(), .by=c(color_flag, spp))
  ggplot(aes(x=obs.D, y=fct_rev(overlap), color=color_flag)) +
  geom_point(size=2) +
  scale_color_manual(values = c("<0.05" = "black",
                                ">0.05" = "red")) +
  xlim(0,1)+
  #  geom_vline(xintercept = 0.05, color="red", linetype ="dotted") +
  labs(x="Shoener's D", y="Monthly overlap", color="p-value", size="D observed")+
  #  scale_size(range = c(0, 4)) +  # ajusta rango de tamaños
  facet_wrap(~spp, ncol=10) 


# Comparar los resultados
old.paths <- list.files("./species/mig_ENM/overlaps_tables/season/", full.names=T)[c(1,65,71,72,78, 79)]
new.paths <- list.files('./borrar/elpsborrar/', full.names=T)

old.l <- list()
new.l <- list()

for (x in 1:length(old.paths)) {
  old.l[[x]] <- read.table(old.paths[x], dec='.', sep='\t', header=T)
  new.l[[x]] <- read.table(new.paths[x], dec='.', sep='\t', header=T)
}


old <- do.call(rbind, old.l) |> 
  mutate(type='old')

new <- do.call(rbind, new.l) |> 
  mutate(type='new')

rbind(old,new) |> 
  mutate(color_flag = ifelse(p.D > 0.05, ">0.05", "<0.05")) |> 
  ggplot(aes(x=obs.D, y=overlap, color = color_flag, shape=type)) +
  geom_point() +
  facet_wrap(~spp)

# Como que no hay muchas diferencias... Habria que probarlo mensual y ver el ECA como es



# Path de las especies
season.names<-list.files("./species/mig_ENM/seasonaly/", full.names = F) |> str_remove("_ENMs")
season.folder<-list.files("./species/mig_ENM/seasonaly/", full.names = T)

### Esto hay que borrarlo cuando corra

# f <- 1
for (f in 1:length(season.names)) { #1:length(season.names)
  
  ## Cargo el Rdata que contiene los ENM y el background ---------------------
  
  temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
  load(list.files(season.folder[f], pattern="*.Rdata$", full.names = T), envir = temp_env)
  
  # Objetos de temp_env que se van a usar
  bg_swd.list<-temp_env$bg_swd.list
  abex.list<-temp_env$abex.list
  
  # PCA total (toda la especie) ---------------------------------------------
  
  # Background
  bgd <- 
    do.call(rbind, bg_swd.list) |> 
    select(lon, lat, prec, tmax, tmin) |> 
    rownames_to_column(var = "elip") |> 
    mutate(elip = elip |> str_sub(1,19) |> str_remove_all("\\."))
  
  # Ambiente
  env2 <- list()
  
  # Extracting environmental values from temporal data frame
  for(p in 1:4){ # 1 to 12 (months / length of analysis (12))
    env2[[p]] <- abex.list[[p]]$temporal_df |> 
      mutate(elip = names(abex.list)[p])
  }
  
  env <-
    do.call(rbind, env2) |> 
    select(lon, lat, prec, tmax, tmin, elip)
  
  #  Bind background and environmental dataframes
  data <- rbind.data.frame(bgd[,c(1,4:6)], env[,c(6,3:5)])
  
  # Weight vector. Occurrences= 0 and background (survey sites)=1
  w <- c(rep(1, nrow(bgd)), rep(0, nrow(env)))
  
  pca.cal <- ade4::dudi.pca(data[,2:4],
                            row.w = w,
                            center = T,
                            scale = T,
                            scannf = F,
                            nf = 2) # Produce solo dos PC
  
  # Ellipsoid ID  
  pca.cal2 <- pca.cal$li |> 
    mutate(elip = data$elip)
  
  write.table(pca.cal2,
              paste0("./species/mig_ENM/overlaps_tables/PCA_s/", season.names[f] |> str_remove_all("_ENMm"),
                     "_s_PCA.txt"),
              sep="\t", dec = ".", row.names=F)
  
  
  print(paste(season.names[f], f, "done"))
  }


# comparcion cohen vs yo
library(tidyverse)

data.ovrlp.m<-read.table("./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt",sep="\t", dec = ".", header=T)

mio <- data.ovrlp.m |> 
  summarise(eca=mean(obs.D), sd.eca=sd(obs.D), .by=c(spp)) |> 
  mutate(spp = spp |> str_remove_all("_ENMm"))


cohen <- read.table("./borrar/cohenetal2023_appendix_1.txt", sep="\t", dec = ".", header=T)
head(cohen)


cohen2 <- cohen |> 
  tibble() |> 
  select(Species, TD_Adj.Maha.dist) |> 
  separate(col=Species,
           into = c("genero", "especie"), 
           sep = " ", 
           remove = FALSE) |> 
  mutate(let_gen = str_sub(genero, start = 1, end = 3),
    let_esp = str_sub(especie, start = 1, end = 5),
    spp = paste(let_gen, let_esp, sep="_"),
    count = 1) 


mio |> 
  left_join(cohen2, by = "spp") |> 
  tibble() |>
  na.omit() |>
  select(spp, eca, TD_Adj.Maha.dist) |> 
  mutate(eca = scale(eca), TD_Adj.Maha.dist = scale(TD_Adj.Maha.dist)) |> 
  ggplot(aes(x=TD_Adj.Maha.dist, y=eca)) +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  labs(x="Cohen index", y="Environmental Coincidence Average")

a <- mio |> 
  left_join(cohen2, by = "spp") |> 
  tibble() |>
  na.omit() |>
  select(spp, eca, TD_Adj.Maha.dist) |> 
  mutate(eca = scale(eca), TD_Adj.Maha.dist = scale(TD_Adj.Maha.dist)) |> 
  select(!spp)

a.lm <- lm(a$eca~a$TD_Adj.Maha.dist)
plot(a.lm)
summary(a.lm)
  
  cor(method="spearman")


pca.cal2 |>
# pca.cal3 |>
  tibble() |> 
  summarise(
    min.a1=min(Axis1), min.a2=min(Axis2),
    max.a1=max(Axis1), max.a2=max(Axis2), .by = elip)

pca.cal3$elip |> unique()


#### Valores de los modelos
library(tidyverse)
spp.list.s <- list.files('./species/mig_ENM/seasonaly/', full.names=T)
spp.list.a <- list.files('./species/mig_ENM/annual/', full.names=T)
spp.list.m <- list.files('./species/mig_ENM/monthly/', full.names=T)

list.enm.s <- list()
list.enm.a <- list()
list.enm.m <- list()

for (x in 1:length(spp.list.s)) {
  list.enm.s[[x]] <- 
    read.table(paste0(spp.list.s[x], '/ENM_mods_s_', 
                      basename(spp.list.s[x] |> 
                                 str_remove('_ENMs')), '.txt'), 
               sep='\t', dec='.', header=T)
  list.enm.a[[x]] <- 
    read.table(paste0(spp.list.a[x], '/ENM_mods_a_', 
                      basename(spp.list.s[x] |> 
                                 str_remove('_ENMs')), '.txt'), 
               sep='\t', dec='.', header=T)
  list.enm.m[[x]] <- 
    read.table(paste0(spp.list.m[x], '/ENM_mods_', 
                      basename(spp.list.s[x] |> 
                                 str_remove('_ENMs')), '.txt'), 
               sep='\t', dec='.', header=T)
}

spp.df.s <- do.call(rbind, list.enm.s)
spp.df.a <- do.call(rbind, list.enm.a)
spp.df.m <- do.call(rbind, list.enm.m)

spp.df.s |> 
  filter(om_rate_train > 0.06) |>
  nrow()
  # head()
spp.df.a |> 
  filter(om_rate_train > 0.06) |> 
  nrow()
  
s.g <- spp.df.s |> 
  ggplot(aes(x=om_rate_test)) +
  geom_histogram() +
  geom_vline(xintercept = c(0.05, 0.1, mean(spp.df.s$om_rate_test)), 
             color=c("red", "blue", "black"), linetype="dashed") +
  labs(title = "Omision rate test (Season)") +
  theme_classic()
a.g <- spp.df.a |> 
  ggplot(aes(x=om_rate_test)) +
  geom_histogram() +
  geom_vline(xintercept = c(0.05, 0.1, mean(spp.df.a$om_rate_test)), 
             color=c("red", "blue", "black"), linetype="dashed") +
  labs(title = "Omision rate test (Annual)") +
  theme_classic()
m.g <- spp.df.m |> 
  ggplot(aes(x=om_rate_test)) +
  geom_histogram() +
  geom_vline(xintercept = c(0.05, 0.1, mean(spp.df.m$om_rate_test)), 
             color=c("red", "blue", "black"), linetype="dashed") +
  labs(title = "Omision rate test (Month)") +
  theme_classic()

# library(patchwork)
s.g / a.g / m.g + plot_annotation(tag_levels = "A")


spp.df.m[which(spp.df.m$om_rate_test==max(spp.df.m$om_rate_test)),]

x <- 39
occ.list.vect <- list()
occ.list.thin <- list()
occ.df.thin <- data.frame(spp=NA, En=NA, Fe=NA, Ma=NA, Ab=NA, May=NA, Jun=NA, Jul=NA, Ag=NA, Se=NA, Oc=NA, No=NA, Di=NA)
for (x in 1:103) {
  
temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
load(list.files(spp.list.m[x], pattern="*.Rdata$", full.names = T), 
     envir = temp_env)

occ.df.thin[x,1] <- spp.names[x]
for (y in 1:12) {
 occ.df.thin[x,y+1] <- 
  temp_env$abex.list[[y]]$temporal_df |> 
    nrow()
}

# Solo tiene 9 registros?
spp<-vect(spp.list[x])

occ.list.vect[[x]] <- 
  spp |> 
  as.data.frame() |> 
  tibble() |> 
  summarise(n=n(), .by = month) |> 
  mutate(species=spp.names[x]) |> 
  pivot_wider(names_from = month, values_from = n)

print(paste(spp.names[x], x, "ready"))
}

do.call(rbind, occ.list.vect)|> write.table("./borrar/occ.df.vect.txt",
                                            sep="\t",
                                            dec=".",
                                            row.names=F)

occ.df.thin |> write.table("./borrar/occ.df.thin.txt",
                           sep="\t",
                           dec=".",
                           row.names=F)
