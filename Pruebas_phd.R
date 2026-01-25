


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

#### minimos y maximos ####
  
  library(ecospat)
  library(terra)
  library(sf)
  library(tidyverse)
  
  # Path de las especies
  season.names<-list.files("./species/mig_ENM/seasonaly/", full.names = F) |> str_remove("_ENMs")
  season.folder<-list.files("./species/mig_ENM/seasonaly/", full.names = T)
  
  # Tabla que organiza el bucle (combinación 1 entre 4)
  combs.s<-data.frame(a= rep(1:3, 3:1), b= unlist(lapply(2:4, function(i) i:4)))
  
  # # ---- borrar
  # 
  # a <- data.frame(spp=NA, dif=NA)
  # a.list <- list()
  # #---
  
  data.ovrlp.list<-list()
  
  # f<-103
  for (f in c(1:21,25:103)) { #1:length(season.names)
  
    # 1. Cargo el Rdata que contiene los ENM y el background ---------------------
    
    temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
    load(list.files(season.folder[f], pattern="*.Rdata$", full.names = T), envir = temp_env)
    
    # Objetos de temp_env que se van a usar
    
    abex.list<-temp_env$abex.list
    
    abbg.list <- list()
    for(x in 1:length(temp_env$abbg.list)) { 
      abbg.list[[x]] <-  temp_env$abbg.list[[x]]$env_bg
      names(abbg.list)[x]<-paste0("abbg_", season.names[f], "_", x)
    }
    
    # Se crea el directorio donde se guardan los datos
    dir.create(paste0(season.folder[f], "/S_overlap"))
    
    # Recipiente de los datos
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
    
    # 2. PCA total (toda la especie) ---------------------------------------------
    
    ##### 2.1 Background ----------------------------------------------------------

    bgd <- 
      do.call(rbind, abbg.list) |> 
      select(lon, lat, prec, tmax, tmin, ID_YEAR) |> 
      rownames_to_column(var = "elip") |> 
      mutate(elip = elip |> str_remove("\\..*"))


    ##### 2.2 Ambiente ------------------------------------------------------------
    env2 <- list()
    
    # Extracting environmental values from temporal data frame
    for(p in 1:length(abex.list)){ # 1 to 4 (seasons / length of analysis (4))
      env2[[p]] <- abex.list[[p]]$temporal_df |> 
        mutate(elip = names(abex.list)[p])
    }
    
    env <-
      do.call(rbind, env2) |> 
      select(lon, lat, prec, tmax, tmin, elip, layers_path) |> 
      rename("ID_YEAR"=layers_path)
    
    ##### 3. Construcción del PCA ----------------------------------------------------

    #  Bind background and environmental dataframes
    data <- rbind.data.frame(bgd[,c(1,4:7)], env[,c(6,3:5,7)])
    
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
      mutate(elip = data$elip,
             ID_YEAR = data$ID_YEAR)
    
    # Extraction ellipsoid names
    
    bg_name <-  unique(pca.cal2$elip)[1:4] # It will be different for months 
    abex_name <-  unique(pca.cal2$elip)[5:8] # It will be different for months
    
    # Determinación de si el rango del bg es mas que el rango de los puntos
    
    min.max.bg <- pca.cal2 |> 
      filter(str_starts(elip, "abbg_")) |> 
      dplyr::summarise(min.a1=min(Axis1), 
                       max.a1=max(Axis1), 
                       min.a2=min(Axis2), 
                       max.a2=max(Axis2), .by = elip) |> 
      arrange(elip)
    
    min.max.sp <- pca.cal2 |> 
      filter(str_starts(elip, "abex_")) |> 
      dplyr::summarise(min.a1=min(Axis1), 
                       max.a1=max(Axis1), 
                       min.a2=min(Axis2), 
                       max.a2=max(Axis2), .by = elip) |> 
      arrange(elip)
    
    min.max.logic <- data.frame( 
    min.a1 = min.max.bg$min.a1 > min.max.sp$min.a1,
    max.a1 = min.max.bg$max.a1 < min.max.sp$max.a1,
    min.a2 = min.max.bg$min.a2 > min.max.sp$min.a2,
    max.a2 = min.max.bg$max.a2 < min.max.sp$max.a2) |> as.matrix()

    if (any(min.max.bg$min.a1 > min.max.sp$min.a1) ||
        any(min.max.bg$min.a2 > min.max.sp$min.a2) ||
        any(min.max.bg$max.a1 < min.max.sp$max.a1) ||
        any(min.max.bg$max.a2 < min.max.sp$max.a2)) {
      
      elips.prob <- 
       which(min.max.logic, arr.ind = TRUE) |> 
        as_tibble() |> 
        mutate(elip = bg_name[row]) |>  # Mapeamos el ID original
        select(elip) |> 
        pull() |> 
        unique()
      
      print(data.frame(which=paste(length(elips.prob), "correct needed")))
      
  #    ## --- borrar
  #    # Esto es para crear una tabla que me diga donde hay estos conflictos
  #    a.list[[f]] <- data.frame(spp=elips.prob, ID=rep(f, length(elips.prob)))
  #   }
  #   print(paste(season.names[f], f, "done"))
  # }
  # do.call(rbind, a.list) |> 
  #   write.table("./borrar/bg_conflict.S.txt", sep="\t", dec=".", row.names=F)
  #    elips.prob <- c("abbg_Vir_solit_3", "abbg_Vir_solit_2")
  #    
  #    ## ---
     
      # Este bucle identifica los espacios donde hay conflictos (elips.prob) y hace la corrección para cada conjunto de       datos. Corta los raster correspondiente al tiempo de cada registro, calcula el min y el max y obtiene un min-max        final que se añade a bgd.  
     # e <- 1
     for (e in 1:length(elips.prob)) {
       
      bgd.filter <-  bgd |> 
         filter(elip == elips.prob[e])

      # 1. Construcción de los Cuadrados (Polígonos)
      # La lógica es: xmin, ymin -> xmax, ymin -> xmax, ymax -> xmin, ymax -> cerrar
      eme.t <-
        env |>
        filter(elip==elips.prob[e] |> str_replace("abbg","abex")) |> 
        tibble() |> 
        st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
        mutate(
          # Construimos la cadena de texto que define el polígono (WKT). 
          # El cuadrado debe ser de 4.5° de lado, pero se crea como si fuera un buffer de 2.25
          wkt_geometry = glue::glue(
            "POLYGON (({lon - 2.25} {lat - 2.25}, 
                 {lon + 2.25} {lat - 2.25}, 
                 {lon + 2.25} {lat + 2.25}, 
                 {lon - 2.25} {lat + 2.25}, 
                 {lon - 2.25} {lat - 2.25}))")) |> # El último punto debe ser igual al primero para cerrar
        # Convertimos esa columna de texto en geometría real sf
        st_drop_geometry() %>% # Soltamos la geometría de puntos original
        st_as_sf(wkt = "wkt_geometry", crs = 4326) |> 
        summarise() |> 
        vect()
      
      # 2. Obtener los paths correspondientes a las fechas de los registros
      
      dat.path <- 
        pca.cal2 |> 
        filter(elip == elips.prob[e]) |> 
        select(ID_YEAR) |> 
        unique() |> pull()
      
      print(paste(elips.prob[e], length(dat.path), "rast paths"))
      
      # 2.1 Abrir raster de cada path
      
      min.max.list <- list()
      
      # r <- 1L
      for (r in 1:length(dat.path)) {# length(dat.path)
        r.stack <- rast(list.files(dat.path[r] |> str_replace("E:", "D:"), full.names = T))
        
        # Extracción de los valores dentro de M
        min.max.list[[r]] <- 
          r.stack |> 
          terra::extract(eme.t, ID=F) |> 
          rename("prec" = 1, "tmax" = 2, "tmin" = 3) |> 
          apply(MARGIN=2, range, na.rm=T) # [1,] min / [2,] max 
      }
      
      # Esto es lo que se adicionaría al bg filtrado (bgd.filter)
      min.max.df <- 
        do.call(rbind, min.max.list) |> 
        apply(MARGIN=2, range, na.rm=T) |> 
        data.frame()
       
      bgd.filter2 <- 
        bgd.filter |> 
        tibble() |> 
        # min
          add_row(elip=elips.prob[e], 
                  lon=bgd.filter$lon[1],
                  lat=bgd.filter$lat[1],
                  prec=min.max.df$prec[1], 
                  tmax=min.max.df$tmax[1],
                  tmin=min.max.df$tmin[1],
                  ID_YEAR=bgd.filter$ID_YEAR[1])  |> 
        # max
        add_row(elip=elips.prob[e], 
                lon=bgd.filter$lon[1],
                lat=bgd.filter$lat[1],
                prec=min.max.df$prec[2], 
                tmax=min.max.df$tmax[2],
                tmin=min.max.df$tmin[2],
                ID_YEAR=bgd.filter$ID_YEAR[1])
      
      bgd <- bgd |> 
        filter(elip != elips.prob[e]) |>
        bind_rows(bgd.filter2) |>
        tibble()
      
      print(paste(elips.prob[e], "correction done"))
       }
     

      # 4. Calculo de PCA corregido (bgd) ------------------------------------------
      
    
     #  Bind background and environmental dataframes
     data <- rbind.data.frame(bgd[,c(1,4:7)], env[,c(6,3:5,7)])

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
       mutate(elip = data$elip,
              ID_YEAR = data$ID_YEAR)
     
     # Prueba de que funcionó

     min.max.bg <- pca.cal2 |>
       filter(str_starts(elip, "abbg_")) |>
       dplyr::summarise(min.a1=min(Axis1),
                        max.a1=max(Axis1),
                        min.a2=min(Axis2),
                        max.a2=max(Axis2), .by = elip) |> 
       arrange(elip)

     min.max.sp <- pca.cal2 |>
       filter(str_starts(elip, "abex_")) |>
       dplyr::summarise(min.a1=min(Axis1),
                        max.a1=max(Axis1),
                        min.a2=min(Axis2),
                        max.a2=max(Axis2), .by = elip) |> 
       arrange(elip)

     min.max.logic <- data.frame(
       min.a1 = min.max.bg$min.a1 > min.max.sp$min.a1,
       max.a1 = min.max.bg$max.a1 < min.max.sp$max.a1,
       min.a2 = min.max.bg$min.a2 > min.max.sp$min.a2,
       max.a2 = min.max.bg$max.a2 < min.max.sp$max.a2) |> as.matrix()
     
     # Si no funciona, entonces esos puntos extremos se eliminan de los registros
     
     if (any(min.max.bg$min.a1 > min.max.sp$min.a1) ||
         any(min.max.bg$max.a1 < min.max.sp$max.a1) ||
         any(min.max.bg$min.a2 > min.max.sp$min.a2) ||
         any(min.max.bg$max.a2 < min.max.sp$max.a2)) {

       print(paste("records need correction"))
       
       limites_bg.pca <- 
         pca.cal2 |> 
         group_by(elip) |> 
         summarise(across(
           c(Axis1, Axis2), 
           list(min = \(x) min(x, na.rm = TRUE), 
                max = \(x) max(x, na.rm = TRUE)),
           .names = "{.col}_{.fn}" # Esto genera nombres como tmin_min, tmin_max, etc.
         )) |> 
         filter(str_starts(elip, "abbg_")) |> 
         mutate(elip2 = elip |> str_remove("abbg_")) |> 
         select(!elip)
       
       pca.cal2 <- 
       pca.cal2 |> 
         tibble() |> 
         filter(str_starts(elip, "abex_")) |> 
         mutate(elip2 = elip |> str_remove("abex_")) |> 
         left_join(limites_bg.pca, by = "elip2") |> 
           # Para cada variable, encerramos el valor entre su min y max de categoría
         mutate(Axis1 = pmax(pmin(Axis1, Axis1_max), Axis1_min),
                Axis2 = pmax(pmin(Axis2, Axis2_max), Axis2_min)) |> 
         # 3. Limpieza: eliminamos las columnas de límites que ya no necesitamos
         select(Axis1, Axis2, elip, ID_YEAR) |> 
         rbind(pca.cal2 |> filter(str_starts(elip, "abbg_")))
       

       # # Prueba de que funcionó x2
       # 
       # min.max.bg <- pca.cal2 |>
       #   filter(str_starts(elip, "abbg_")) |>
       #   dplyr::summarise(min.a1=min(Axis1),
       #                    max.a1=max(Axis1),
       #                    min.a2=min(Axis2),
       #                    max.a2=max(Axis2), .by = elip) |> 
       #   arrange(elip)
       # 
       # min.max.sp <- pca.cal2 |>
       #   filter(str_starts(elip, "abex_")) |>
       #   dplyr::summarise(min.a1=min(Axis1),
       #                    max.a1=max(Axis1),
       #                    min.a2=min(Axis2),
       #                    max.a2=max(Axis2), .by = elip) |> 
       #   arrange(elip)
       # 
       # min.max.logic <- data.frame(
       #   min.a1 = min.max.bg$min.a1 > min.max.sp$min.a1,
       #   max.a1 = min.max.bg$max.a1 < min.max.sp$max.a1,
       #   min.a2 = min.max.bg$min.a2 > min.max.sp$min.a2,
       #   max.a2 = min.max.bg$max.a2 < min.max.sp$max.a2) |> as.matrix()
       
     }

     
    } 
    
    
    # 5. Select PCA values for each ellipsoid ------------------------------------

    # x<-1
    for (x in 1:nrow(combs.s)) {
      
      ##### 5.1 PCA climatic scores  ----------------------------------------------------
      
      # Backgrounds
      scores.clima1 <- pca.cal2 |> 
        filter(elip == bg_name[combs.s[x,1]]) |> 
        select(!c(elip, ID_YEAR))
      
      scores.clima2 <- pca.cal2 |> 
        filter(elip == bg_name[combs.s[x,2]]) |> 
        select(!c(elip, ID_YEAR))
      
      scores.clima12 <- rbind(scores.clima1, scores.clima2)
      
      # Environments
      scores.sp1a <- pca.cal2 |> 
        filter(elip == abex_name[combs.s[x,1]]) |> 
        select(!c(elip, ID_YEAR))
      
      scores.sp2b <- pca.cal2 |> 
        filter(elip == abex_name[combs.s[x,2]]) |> 
        select(!c(elip, ID_YEAR))
      
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
      
      png(paste0(season.folder[f], "/S_overlap/", 
                 season.names[f], "_PCA_", "s", combs.s[x,1] ,".png"),
          width = 800, height = 800, res = 100)
      par(mar=c(5,3,2,2))
      ecospat.plot.niche (z1,
                          title = paste0(season.names[f] |> str_replace("_", " "), 
                                         "-s", combs.s[x,1]),
                          name.axis1 = "PC1",
                          name.axis2 = "PC2",
                          cor = F)
      dev.off()
      
      # Second ellipsoid
      png(paste0(season.folder[f], "/S_overlap/", 
                 season.names[f], "_PCA_", "s", combs.s[x,2] ,".png"), 
          width = 800, height = 800, res = 100)
      par(mar=c(5,3,2,2))
      ecospat.plot.niche (z2,
                          title = paste0(season.names[f] |> str_replace("_", " "), 
                                         "-s", combs.s[x,2]),
                          name.axis1 = "PC1",
                          name.axis2 = "PC2",
                          cor = F)
      dev.off()
      
      
      # Overlap graph -----------------------------------------------------------
      
      png(paste0(season.folder[f], "/S_overlap/", season.names[f], 
                 "_B.ovl_", "s", combs.s[x,1], "s", combs.s[x,2]  ,".png"), width = 800, height = 800, res = 100)
      
      par(mar=c(5,3,2,2))
      ecospat.plot.niche.dyn(z1,
                             z2,
                             quant = 0.25,
                             interest = 1,
                             title = paste0("Niche overlap ", season.names[f] |> str_replace("_", " "), 
                                            "-s", combs.s[x,1], "_", "s", combs.s[x,2]),
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
      
      data.overlap[x,1] <- season.names[f]
      data.overlap[x,2] <- paste0("S", combs.s[x, 1], "_S", combs.s[x, 2])
      data.overlap[x, 3:7] <- unlist(sim.test$obs)
      data.overlap[x,8] <- sim.test$p.D
      data.overlap[x,9] <- sim.test$p.I
      data.overlap[x,10] <- sim.test$p.expansion
      data.overlap[x,11] <- sim.test$p.stability
      data.overlap[x,12] <- sim.test$p.unfilling
      
      
      ## Plot Similarity test -  evalua si la similaridad es mas similar de la esperada por azar
      
      png(paste0(season.folder[f], "/S_overlap/", season.names[f], "_ST_", "s", combs.s[x,1], "s", combs.s[x,2]  ,".png"), width = 800, height = 800, res = 100)
      
      ecospat.plot.overlap.test(sim.test, "D", paste0("Similarity greater rand.type 1", "\n", 
                                                      season.names[f] |> str_replace("_", " "), 
                                                      "-s",combs.s[x,1],"_s", combs.s[x,2]))
      
      dev.off()
      
      print(paste(basename(season.folder[f]), combs.s[x,1], "v", combs.s[x,2], "ready"))
    }
    
    write.table(data.overlap, 
                 paste0("./species/mig_ENM/overlaps_tables/B_ovl_season/",
                        season.names[f],
                        "_s_ovl.txt"), sep = "\t", dec = ".", row.names = F)
    data.ovrlp.list[[f]] <- data.overlap
    
    print(paste(basename(season.folder[f]), f, "done"))
    rm(list = setdiff(ls(), c("season.names", "season.folder", "combs.s", "data.ovrlp.list")))
    
  }

