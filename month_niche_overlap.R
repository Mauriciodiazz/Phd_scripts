
# Comparación de nichos con Broennimann et al. 2012, Di Cola et al. 2017 

# Librerias
library(ecospat)
library(terra)
library(tidyverse)

# Path de las especies
month.names<-list.files("./species/mig_ENM/monthly/", full.names = F) |> str_remove("_ENMm")
month.folder<-list.files("./species/mig_ENM/monthly/", full.names = T)

# Tabla que organiza el bucle (combinación 1 entre 12)
combs.m<-data.frame(a= rep(1:11, 11:1), b= unlist(lapply(2:12, function(i) i:12)))

data.ovrlp.list<-list()

# f<-1
for (f in 8:length(month.names)) { #1:length(month.names)
  
  # Cargo el Rdata que contiene los ENM y el background ---------------------
  temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
  load(list.files(month.folder[f], pattern="*.Rdata$", full.names = T), envir = temp_env)
  
  # Objetos de temp_env que se van a usar
  bg_swd.list<-temp_env$bg_swd.list
  abex.list<-temp_env$abex.list
  
  dir.create(paste0(month.folder[f], "/m_overlap"))
  
  data.overlap <- data.frame(
    spp = NA,
    overlap = NA,
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
  
  # x<-1
  for (x in 1:nrow(combs.m)) {
    
    # Perfiles climaticos -----------------------------------------------------
    
    ### Background --------------------------------------------------------------
  
    climafin1<-bg_swd.list[[combs.m[x,1]]] |> select(lon, lat, prec, tmax, tmin) # enero
    #head(climafin1)
    climafin2<-bg_swd.list[[combs.m[x,2]]] |> select(lon, lat, prec, tmax, tmin) # febrero
    #head(climafin2)
    # 
    # 1. Crear un perfil climatico de ambas especies 
    ambosclimas<- na.omit(rbind(climafin1,climafin2))
    
    
    # tENM spp ----------------------------------------------------------------

    # Esto es el SDW que genera tenm para cada especie
    sp1 <- abex.list[[combs.m[x, 1]]]$temporal_df |>
      mutate(el.name = paste0(month.names[f], "_", x)) |>
      select(el.name, lon, lat, prec, tmax, tmin) |> #, prec, tmax, tmin
      as.data.frame()
    #head(sp1)

    sp2 <- abex.list[[combs.m[x, 2]]]$temporal_df |>
      mutate(el.name = paste0(month.names[f], "_", x + 1)) |>
      select(el.name, lon, lat, prec, tmax, tmin) |> #
      as.data.frame()
    #head(sp2)
    
    # PCA ---------------------------------------------------------------------
    
    # Generar una base de datos que incluya los valores de las variables de clima (que incluye las m de ambas especies generado previamente o sus backgrounds) mas los valores de las variables para la sp1 y la sp2
    
    data <- rbind.data.frame(ambosclimas[, 3:5], sp1[, 4:6], sp2[, 4:6])
    
    # Vector de peso 0 para las ocurrencias y 1 para los sitios del ?rea de estudio 
    w <- c(rep(1, nrow(ambosclimas)), rep(0, nrow(sp1)), rep(0, nrow(sp2)))
    
    # head(w)
    # tail(w)
    
    # El PCA se realiza con todos los datos del ?rea de estudio (es decir, ambas emes). Las presencias no son usadas para la calibraci?n del PCA pero sus coordenadas en los PCA son calculadas
    
    # library(ade4)
    pca.cal <- ade4::dudi.pca(data,
                              row.w = w,
                              center = T,
                              scale = T,
                              scannf = F,
                              nf = 2)
    #head(pca.cal)
    
    # Filas en las que estan los datos de clima  de cada una de las especies
    
    row.clima1 <- 1:nrow(climafin1) # Background elipsoide 1
    row.clima2 <- (nrow(climafin1) + 1):(nrow(climafin1) + nrow(climafin2)) # Background elipsoide 2
    row.clima12 <- 1:nrow(ambosclimas) # Background total
    row.sp1a <- (1 + nrow(climafin1) + nrow(climafin2)):(nrow(climafin1) + 
                                                           nrow(climafin2) +
                                                           nrow(sp1))
    row.sp2b <- (1 + nrow(climafin1) + nrow(climafin2) + nrow(sp1)):(nrow(climafin1) +
                                                                       nrow(climafin2) + 
                                                                       nrow(sp1) + nrow(sp2))
    
    # Filtro de las coordenadas en cada uno de los ejes del PCA de todos los datos (emes y especies)
    scores.clima1 <- pca.cal$li[row.clima1, ]
    scores.clima2 <- pca.cal$li[row.clima2, ]
    # scores.clima12 <- pca.cal$li[row.clima12, ]
    scores.sp1a <- pca.cal$li[row.sp1a, ]
    scores.sp2b <- pca.cal$li[row.sp2b, ]
    
    # Con un perfil de background normal esto funcionaría bien. Sin embargo, cuando se crea el background con el paquete tenm, a veces los rangos de los PCA del backgorund son menores que los ejes de la especie, lo que genera problemas en la generación de los kernel. Es decir, min(scores.clima12)>min(scores.sp) y debe ser al contrario: el clima global debe incluir el clima de la especie
    
    scores.clima12<- pca.cal$li[row.clima12, ] |> 
      add_row(Axis1=min(pca.cal$li$Axis1), Axis2=min(pca.cal$li$Axis2)) |> 
      add_row(Axis1=max(pca.cal$li$Axis1), Axis2=max(pca.cal$li$Axis2))
    
    # Contribuci?n de las varaibles a cada componente
    contribucion<- ecospat.plot.contrib(contrib=pca.cal$co,
                                        eigen = pca.cal$eig)
    
    
    # Superficie de densidad de registros -------------------------------------
    
    z1<-ecospat.grid.clim.dyn(glob=scores.clima12, 
                              glob1=scores.clima1, 
                              sp=scores.sp1a, 
                              R = 100)
    #head(z1)
    
    z2<-ecospat.grid.clim.dyn(glob=scores.clima12,
                              glob1=scores.clima2, 
                              sp=scores.sp2b,
                              R = 100)
    #head(z2)
    
    # Graficos de PCA ---------------------------------------------------------
    
    # Primer elipsoide
    png(paste0(month.folder[f], "/m_overlap/", month.names[f], "_PCA_", "m", combs.m[x,1] ,".png"),
        width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche (z1,
                        title = paste0(month.names[f], " m", combs.m[x,1]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
    dev.off()
    
    # Segundo elipsoide
    png(paste0(month.folder[f], "/m_overlap/", month.names[f], "_PCA_", "m", combs.m[x,2] ,".png"), 
        width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche (z2,
                        title = paste0(month.names[f], " m", combs.m[x,2]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
    dev.off()
    
    # Overlap
    png(paste0(month.folder[f], "/m_overlap/", month.names[f], "_ovlp_", "m_", combs.m[x,1], "_",combs.m[x,2],".png"), width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche.dyn(z1,
                           z2,
                           quant = 0.25,
                           interest = 1,
                           title = paste0("Niche overlap ", month.names[f], 
                                          " m", combs.m[x,1], " vs ", "m", combs.m[x,2]),
                           name.axis1 = "PC1",
                           name.axis2 = "PC2")
    dev.off()
    
    # Para obtener los valores de D e I se corre lo siguiente:
    # ecospat.niche.overlap(z1=z1,z2=z2,cor=TRUE)
    
    # Test de similaridad (di Cola et al 2017 - higher) ----------------------
    
    sim.test <- ecospat.niche.similarity.test(
      z1,
      z2,
      rep = 100,
      overlap.alternative = "higher",
      rand.type = 1)

    data.overlap[x,1] <- month.names[f] |> str_remove("_ENMm")
    data.overlap[x,2] <- paste0("m", combs.m[x, 1], "_m", combs.m[x, 2])
    data.overlap[x, 3:7] <- unlist(sim.test$obs)
    data.overlap[x,8] <- sim.test$p.D
    data.overlap[x,9] <- sim.test$p.I
    data.overlap[x,10] <- sim.test$p.expansion
    data.overlap[x,11] <- sim.test$p.stability
    data.overlap[x,12] <- sim.test$p.unfilling
    
    ## Plot Similarity test -  evalua si la similaridad es mas similar de la esperada por azar
    
    png(paste0(month.folder[f], "/m_overlap/", 
               month.names[f], "_ST_", "m_", combs.m[x,1], "_", combs.m[x,2]  ,".png"), 
        width = 800, height = 800, res = 100)
    
    ecospat.plot.overlap.test(sim.test, "D", paste0("Similarity greater rand.type 1", "\n", 
                                                    month.names[f], "_m",combs.m[x,1],"_m", combs.m[x,2]))
    
    dev.off()
    
    print(paste(basename(month.folder[f]), combs.m[x,1], "v", combs.m[x,2], "ready"))
  }
  
  print(paste(basename(month.folder[f]), f, "done"))
  
  write.table(data.overlap, paste0("./species/mig_ENM/overlaps_tables/month/", month.names[f], "_m_ovl.txt"),
              sep="\t", dec = ".", row.names=F)
  
  data.ovrlp.list[[f]]<-data.overlap
  
  rm(list = setdiff(ls(), c("month.names", "month.folder", "combs.m", "data.ovrlp.list")))
}

# En caso de haber borrado el objeto data.ovrlp.list
df.path <- list.files("./species/mig_ENM/overlaps_tables/month/", full.names = T)

for (x in 1:length(df.path)) {
  data.ovrlp.list[[x]] <- read.table(df.path[x], sep="\t", dec=".", header=T)
}

data.ovrlp.df<-do.call(rbind, data.ovrlp.list)
# write.table(data.ovrlp.df, "./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt", sep="\t", dec = ".", row.names=F)

data.ovrlp.df<-read.table("./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt", sep="\t", dec = ".", header=T)

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
  facet_wrap(~spp, ncol=4) 


