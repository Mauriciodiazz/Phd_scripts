
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
# f<-1
for (f in 1:length(season.names)) {

# Cargo el Rdata que contiene los ENM y el background ---------------------

  load(list.files(season.folder[f], pattern="*.Rdata$", full.names = T))
  
  data.overlap<-data.frame(spp.s=NA,obs.D=NA, obs.I=NA, obs.expan=NA, obs.stabi=NA, obs.unfil=NA, p.D=NA, p.I=NA, p.expan=NA, p.stabi=NA, p.unfil=NA)
  
  # x<-1
  for (x in 1:nrow(combs.s)) {

# Perfiles climaticos -----------------------------------------------------

### Background --------------------------------------------------------------
    
    climafin1<-bg_swd.list[[combs.s[x,1]]] |> select(lon, lat, prec, tmax, tmin) # enero
    #head(climafin1)
    climafin2<-bg_swd.list[[combs.s[x,2]]] |> select(lon, lat, prec, tmax, tmin) # febrero
    head(climafin2)

    # 1. Crear un perfil climatico de ambas especies 
    ambosclimas<- na.omit(rbind(climafin1,climafin2))


# tENM spp ----------------------------------------------------------------
    # Esto es el SDW que genera tenm para cada especie
    sp1 <- abex.list[[combs.s[x, 1]]]$temporal_df |>
      mutate(el.name = paste0(season.names[f], "_", x)) |>
      select(el.name, lon, lat, prec, tmax, tmin) |> #, prec, tmax, tmin
      as.data.frame()
    #head(sp1)

    sp2 <- abex.list[[combs.s[x, 2]]]$temporal_df |>
      mutate(el.name = paste0(season.names[f], "_", x + 1)) |>
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
    scores.clima12 <- pca.cal$li[row.clima12, ]
    scores.sp1a <- pca.cal$li[row.sp1a, ]
    scores.sp2b <- pca.cal$li[row.sp2b, ]


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
    
# si de esto se genera un NA es porque los climas de scores.clima no contienen valores de la especie
# 
# scores.clima12[scores.clima12$Axis2 %in% scores.sp1a$Axis2,] |> nrow()
# scores.clima12[scores.clima12$Axis2 %in% scores.sp2b$Axis2,] |> nrow()


# Graficos de PCA ---------------------------------------------------------

# Primer elipsoide
    png(paste0(season.folder[f], "/", season.names[f], "_PCA_", "s", combs.s[x,1] ,".png"),
        width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche (z1,
                        title = paste0(season.names[f], " s", combs.s[x,1]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
    dev.off()

# Segundo elipsoide
    png(paste0(season.folder[f], "/", season.names[f], "_PCA_", "s", combs.s[x,2] ,".png"), 
        width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche (z2,
                        title = paste0(season.names[f], " s", combs.s[x,2]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
    dev.off()
    
# Overlap
    png(paste0(season.folder[f], "/", season.names[f], "_ovlp_", "s", combs.s[x,1], "s", combs.s[x,2]  ,".png"), width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche.dyn(z1,
                           z2,
                           quant = 0.25,
                           interest = 1,
                           title = paste0("Niche overlap ", season.names[f], 
                                          " s", combs.s[x,1], " vs ", "s", combs.s[x,2]),
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
    
    data.overlap[x, 1] <- paste0(season.names[f], "_s", combs.s[x, 1], "_s", combs.s[x, 2])
    data.overlap[x, 2:6] <- unlist(sim.test$obs)
    data.overlap$p.D[x] <- sim.test$p.D
    data.overlap$p.I[x] <- sim.test$p.I
    data.overlap$p.expan[x] <- sim.test$p.expansion
    data.overlap$p.stabi[x] <- sim.test$p.stability
    data.overlap$p.unfil[x] <- sim.test$p.unfilling


## Plot Similarity test

png(paste0(season.folder[f], "/", season.names[f], "_ST_", "s", combs.s[x,1], "s", combs.s[x,2]  ,".png"), width = 800, height = 800, res = 100)

ecospat.plot.overlap.test(sim.test, "D", paste0("Similarity greater rand.type 1", "\n", 
                                                season.names[f], "_s",combs.s[x,1],"_s", combs.s[x,2]))

dev.off()

  }

  data.ovrlp.list[[f]]<-data.overlap
}
#################### An?lisis con base en di Cola et al 2017 lower ######################




#eq.test<-ecospat.niche.equivalency.test(z1,z2, rep=100, alternative = "lower")

#sim.test<-ecospat.niche.similarity.test(z1,z2, rep=100, alternative = "lower", rand.type = 2)

## Plot equivalencity test
#ecospat.plot.overlap.test(eq.test, "D", "Equivalency lower")
## Plot Similarity test
#ecospat.plot.overlap.test(sim.test, "D", "Similarity lower rand.type 2")


