


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

# Librerias
library(ecospat)
library(terra)
library(tidyverse)

# Path de las especies
season.names<-list.files("./species/mig_ENM/seasonaly/", full.names = F) |> str_remove("_ENMs")
season.folder<-list.files("./species/mig_ENM/seasonaly/", full.names = T)

# Tabla que organiza el bucle (combinación 1 entre 4)
combs.s<-data.frame(a= rep(1:3, 3:1), b= unlist(lapply(2:4, function(i) i:4)))

data.ovrlp.list<-list()
# f<-65
for (f in c(79)) { #1:length(season.names)
  
  # Cargo el Rdata que contiene los ENM y el background ---------------------
  
  temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
  load(list.files(season.folder[f], pattern="*.Rdata$", full.names = T), envir = temp_env)
  
  # Objetos de temp_env que se van a usar
  bg_swd.list<-temp_env$bg_swd.list
  abex.list<-temp_env$abex.list
  
  # dir.create(paste0(season.folder[f], "/S_overlap"))
  
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
    mutate(elip = elip |> str_sub(1,18))
  
  # Ambiente
  # p <- 1 # 1 to 4 (seasons / length of analysis)
  env2 <- list()
  
  for(p in 1:4){
    env2[[p]] <- abex.list[[p]]$temporal_df |> 
      mutate(elip = names(abex.list)[p])
  }
  
  env <-
    do.call(rbind, env2) |> 
    select(lon, lat, prec, tmax, tmin, elip)
    
#  unir todo en un solo df
  data <- rbind.data.frame(bgd[,c(1,4:6)], env[,c(6,3:5)])
  
  # Vector de peso 0 para las ocurrencias y 1 para los sitios del ?rea de estudio 
  w <- c(rep(1, nrow(bgd)), rep(0, nrow(env)))
  
  pca.cal <- ade4::dudi.pca(data[,2:4],
                            row.w = w,
                            center = T,
                            scale = T,
                            scannf = F,
                            nf = 2) # Produce solo dos PC
  pca.cal2 <- pca.cal$li |> 
    mutate(elip = data$elip)
  
  pca.cal2$elip |> unique()
  
  bg_name <-  unique(pca.cal2$elip)[1:4] # ESTO CAMBIARIA PARA LOS MESES
  abex_name <-  unique(pca.cal2$elip)[5:8] # ESTO CAMBIARIA PARA LOS MESES
 
  # Select PCA values for each ellipsoid
  
  # x<-1
  for (x in 1:nrow(combs.s)) {
  
  # PCA climatic scores 
  # Backgrounds
  scores.clima1 <- pca.cal2 |> 
    filter(elip == bg_name[combs.s[x,1]]) |> 
    select(!elip)
  
  scores.clima2 <- pca.cal2 |> 
    filter(elip == bg_name[combs.s[x,2]]) |> 
    select(!elip)
  
  scores.clima12 <- rbind(scores.clima1, scores.clima2) |> 
    add_row(Axis1=min(pca.cal2$Axis1), Axis2=min(pca.cal2$Axis2)) |> 
    add_row(Axis1=max(pca.cal2$Axis1), Axis2=max(pca.cal2$Axis2))

  # Environments
  scores.sp1a <- pca.cal2 |> 
    filter(elip == abex_name[combs.s[x,1]]) |> 
    select(!elip)
  
  scores.sp2b <- pca.cal2 |> 
    filter(elip == abex_name[combs.s[x,2]]) |> 
    select(!elip)
  
# 
#     # Contribuci?n de las varaibles a cada componente
#     contribucion<- ecospat.plot.contrib(contrib=pca.cal$co,
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
    
    # Graficos de PCA ---------------------------------------------------------
    
    # Primer elipsoide
    # png(paste0(season.folder[f], "/S_overlap/", season.names[f], "_PCA_", "s", combs.s[x,1] ,".png"),
    #     width = 800, height = 800, res = 100)
    # par(mar=c(5,3,2,2))
    ecospat.plot.niche (z1,
                        title = paste0(season.names[f], " s", combs.s[x,1]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
    # dev.off()
    
    # Segundo elipsoide
    # png(paste0(season.folder[f], "/S_overlap/", season.names[f], "_PCA_", "s", combs.s[x,2] ,".png"), width = 800, height = 800, res = 100)
    # par(mar=c(5,3,2,2))
    ecospat.plot.niche (z2,
                        title = paste0(season.names[f], " s", combs.s[x,2]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
    # dev.off()
    
    # Overlap
    # png(paste0(season.folder[f], "/S_overlap/", season.names[f], "_ovlp_", "s", combs.s[x,1], "s", combs.s[x,2]  ,".png"), width = 800, height = 800, res = 100)
    # par(mar=c(5,3,2,2))
    ecospat.plot.niche.dyn(z1,
                           z2,
                           quant = 0.25,
                           interest = 1,
                           title = paste0("Niche overlap ", season.names[f], 
                                          " s", combs.s[x,1], " vs ", "s", combs.s[x,2]),
                           name.axis1 = "PC1",
                           name.axis2 = "PC2")
    # dev.off()
    
    # Para obtener los valores de D e I se corre lo siguiente:
    # ecospat.niche.overlap(z1=z1,z2=z2,cor=TRUE)
    
    # Test de similaridad (di Cola et al 2017 - higher) ----------------------
    
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
    
    # png(paste0(season.folder[f], "/S_overlap/", season.names[f], "_ST_", "s", combs.s[x,1], "s", combs.s[x,2]  ,".png"), width = 800, height = 800, res = 100)
    
    ecospat.plot.overlap.test(sim.test, "D", paste0("Similarity greater rand.type 1", "\n", 
                                                    season.names[f], "_s",combs.s[x,1],"_s", combs.s[x,2]))
    
    # dev.off()
    
    # print(paste(basename(season.folder[f]), combs.s[x,1], "v", combs.s[x,2], "ready"))
  }
  
  print(paste(basename(season.folder[f]), f, "done"))
  
  # write.table(data.overlap, 
  #             paste0("./species/mig_ENM/overlaps_tables/season/", season.names[f], "_s_ovl.txt"),
  #             sep="\t", dec = ".", row.names=F)
  write.table(data.overlap, 
              paste0("./borrar/elpsborrar/", season.names[f], "_s_ovl.txt"),
              sep="\t", dec = ".", row.names=F)
  
  # data.ovrlp.list[[f]]<-data.overlap
  
  rm(list = setdiff(ls(), c("season.names", "season.folder", "combs.s", "data.ovrlp.list")))
}

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

