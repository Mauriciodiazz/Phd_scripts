
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
for (f in 1:length(season.names)) { #1:length(season.names)

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
  
  # x<-1
  for (x in 1:nrow(combs.s)) {

# Perfiles climaticos -----------------------------------------------------

### Background --------------------------------------------------------------
    
    climafin1<-bg_swd.list[[combs.s[x,1]]] |> select(lon, lat, prec, tmax, tmin) # enero
    #head(climafin1)
    climafin2<-bg_swd.list[[combs.s[x,2]]] |> select(lon, lat, prec, tmax, tmin) # febrero
    #head(climafin2)

    # 1. Crear un perfil climatico de ambos climas 
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
    png(paste0(season.folder[f], "/S_overlap/", season.names[f], "_PCA_", "s", combs.s[x,1] ,".png"),
        width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche (z1,
                        title = paste0(season.names[f], " s", combs.s[x,1]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
    dev.off()

# Segundo elipsoide
    png(paste0(season.folder[f], "/S_overlap/", season.names[f], "_PCA_", "s", combs.s[x,2] ,".png"), 
        width = 800, height = 800, res = 100)
    par(mar=c(5,3,2,2))
    ecospat.plot.niche (z2,
                        title = paste0(season.names[f], " s", combs.s[x,2]),
                        name.axis1 = "PC1",
                        name.axis2 = "PC2",
                        cor = F)
    dev.off()
    
# Overlap
    png(paste0(season.folder[f], "/S_overlap/", season.names[f], "_ovlp_", "s", combs.s[x,1], "s", combs.s[x,2]  ,".png"), width = 800, height = 800, res = 100)
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
                                                season.names[f], "_s",combs.s[x,1],"_s", combs.s[x,2]))

dev.off()

print(paste(basename(season.folder[f]), combs.s[x,1], "v", combs.s[x,2], "ready"))
  }
  
print(paste(basename(season.folder[f]), f, "done"))

write.table(data.overlap, 
            paste0("./species/mig_ENM/overlaps_tables/season/", season.names[f], "_s_ovl.txt"),
            sep="\t", dec = ".", row.names=F)

  data.ovrlp.list[[f]]<-data.overlap
  
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



# Jaccard overlap ---------------------------------------------------------

# Superposicion de nichos. Elipses de volumen minimo
# Curso: Modelos de nichos ecologicos y de distribucion de especies
# Inecol, 2025

# Instalar paquetes que no esten instalados ---------------------------------------------------------
paquetes <- c("remotes","terra","tidyverse","patchwork")
instalar <- paquetes[!(paquetes %in% rownames(installed.packages()))]

if (length(instalar) > 0) {
  install.packages(instalar)
}

# Instalar ellipsenm2 si no esta instalado ---------------------------------------------------------
if(!require(ellipsenm)){
  remotes::install_github("marlonecobos/ellipsenm2", force = TRUE)
}

# ================================= Similitud de nichos =============================================

library(terra)
library(ellipsenm)
library(tidyverse)
library(patchwork)
library(rgl)

# Preparar los objetos de superposicion de cada especie para realizar los analsis ---------------------
# Hay una vaina y es que el background en algunos casos es menor al tamaño de las especies. Eso genera un error, cuando se corran los modelos de todas las especies hay que tener cuidado con eso y asegurarse de que el tamaño del background sea mayor o igual a la cantidad de puntos. Por esto Set_ameri, Set_petec y Set_rutic no irian en estos resultados preliminares

# season.names<-list.files("./species/mig_ENM/seasonaly/", full.names = F) |> str_remove("_ENMs")
# season.folder<-list.files("./species/mig_ENM/seasonaly/", full.names = T)
season.names
season.folder

# Tables
vol.df <- data.frame(spp=NA, wint=NA, spri=NA, summ=NA, autu=NA)
# ---
a <- data.frame(spp=NA, seas=NA, nrow.sp=NA, nrow.bg=NA, 
                min.prec.sp=NA,
                min.prec.bg=NA,
                max.prec.sp=NA,
                max.prec.bg=NA,
                min.tmin.sp=NA,
                min.tmin.bg=NA,
                max.tmin.sp=NA,
                max.tmin.bg=NA,
                min.tmax.sp=NA,
                min.tmax.bg=NA,
                max.tmax.sp=NA,
                max.tmax.bg=NA)
c <- list()
# ---

for (f in 1:103) { #1 Y 12 length(season.folder)
  temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
  load(list.files(season.folder[f], pattern="*.Rdata$", full.names = T), envir = temp_env)
  
  abex.list<-temp_env$abex.list
  bg_swd.list<-temp_env$bg_swd.list
  
  # Lists
  Nich_ovobj <- list()
  full.overlap <- list()
  union.overlap <- list()
  
  
  for (x in 1:length(abex.list)) {
    # Niche tenm models -------------------------------------------------------
    
    niche.swd <- abex.list[[x]]$temporal_df |>
      mutate(el.name = paste0(season.names[f], "_", x)) |>
      select(el.name, lon, lat, prec, tmax, tmin)
    
    # Backgorund --------------------------------------------------------------
    
    climafin<-bg_swd.list[[x]] |> select(prec, tmax, tmin)
    
    # Overlap objects list----------------------------------------------------------
    
    Nich_ovobj[[x]] <- ellipsenm::overlap_object(data = niche.swd,
                                                 species =  "el.name",
                                                 longitude = "lon",
                                                 latitude = "lat",
                                                 method = "mve1",
                                                 level = 95,
                                                 variables = climafin)
    
#     # ---
#     # OJO que esto solo son pruebas para ver lo del background
#     a[x,1] <- season.names[f]
#     a[x,2] <- x
#     a[x,3] <- nrow(niche.swd)
#     a[x,4] <- nrow(climafin)
#     a[x,5] <- min(niche.swd$prec) #spp
#     a[x,6] <- min(climafin$prec) #bg
#     a[x,7] <- max(niche.swd$prec) #spp
#     a[x,8] <- max(climafin$prec) #bg
#     a[x,9] <- min(niche.swd$tmin) #spp
#     a[x,10] <- min(climafin$tmin) #bg
#     a[x,11] <- max(niche.swd$tmin) #spp
#     a[x,12] <- max(climafin$tmin) #bg
#     a[x,13] <- min(niche.swd$tmax) #spp
#     a[x,14] <- min(climafin$tmax) #bg
#     a[x,15] <- max(niche.swd$tmax) #spp
#     a[x,16] <- max(climafin$tmax) #bg
#     # ---
   } ##### Este siempre va
#   c[[f]] <- a
#   }
# 
# # El bg no debe ser mayor que el vector de las especies. Esto debe de dar cero
#   do.call(rbind, c) |>
#     mutate(res=nrow.sp-nrow.bg) |>
#     filter(res>0)
# 
# # ¿En alguna ocasión los rangos del bg son menores que los rangos de la especie?
#   do.call(rbind, c) |> 
#     mutate(res.min.prec= min.prec.bg - min.prec.sp,
#            res.max.prec= max.prec.sp - max.prec.bg,
#            res.min.tmax= min.tmax.bg - min.tmax.sp,
#            res.max.tmax= max.tmax.sp - max.tmax.bg,
#            res.min.tmin= min.tmin.bg - min.tmin.sp,
#            res.max.tmin= max.tmin.sp - max.tmin.bg) |> 
#     select(1,17:22) |> 
#     pivot_longer(!spp, names_to = "type", values_to = "values") |> 
#     filter(values>0) |> 
#     # summarise(n=n(), .by=spp)
#     ggplot(aes(x=values, y=type)) + 
#     geom_point() +
#     geom_vline(xintercept = 0, color="red", linetype="dashed")
#   
  # Esto debe ir en una lista
  # Estimación de las superposiciones de nicho entre pares de especies con prueba de significancia ------
  Sup_all <- ellipsoid_overlap(Nich_ovobj[[1]], Nich_ovobj[[2]], Nich_ovobj[[3]], Nich_ovobj[[4]],
                               n_points = 100000,
                               significance_test = TRUE, iterations = 1000,
                               overlap_type = "all", confidence_limit = 0.05)
  
  full.overlap[[f]] <- 
    Sup_all@full_overlap |> 
    mutate(spp=basename(season.folder[f]) 
           |> str_remove("_ENMs")) |> 
    rownames_to_column(var = "niche.comp") |> 
    relocate(spp)
  
  write.table(vol.df, paste0(season.folder[f], "/S_overlap/", season.names[f], "_full_ovl.txt"), 
              sep="\t", dec = ".", row.names=F)
  
  union.overlap[[f]] <- 
    Sup_all@union_overlap |> 
    mutate(spp=basename(season.folder[f]) 
           |> str_remove("_ENMs")) |> 
    rownames_to_column(var = "niche.comp") |> 
    relocate(spp)
  
  write.table(vol.df, paste0(season.folder[f], "/S_overlap/", season.names[f], "_uni_ovl.txt"), 
              sep="\t", dec = ".", row.names=F)
  
  # Volumen
  vol.df[f,1]<-season.names[f]
  vol.df[f,2]<-Sup_all@ellipsoids$Niche_1@niche_volume
  vol.df[f,3]<-Sup_all@ellipsoids$Niche_2@niche_volume
  vol.df[f,4]<-Sup_all@ellipsoids$Niche_3@niche_volume
  vol.df[f,5]<-Sup_all@ellipsoids$Niche_4@niche_volume
  
  # Graficos prueba de superposicion de nichos -----------------------------------------------------------
  # Titulos
  mains<-
    paste(rep(c("Union","Full"),each=6), c("winter vs spring", 
                                           "winter vs summer", 
                                           "winter vs autum",
                                           "spring vs summer",
                                           "spring vs autum",
                                           "summer vs autum")) |> 
    matrix(byrow = F, ncol=2) |> 
    as.data.frame() |> 
    rename("union"=V1,"full"=V2)
  
  # Layout
  png(paste0(season.folder[f], "/S_overlap/", season.names[f], "_jaccard.png"), width = 800, height = 800, res = 100)
  #x11(width=10, height=10)
  par(mfrow = c(6,2))
  
  # Grafico
  for(z in 1:6) {
    hist(Sup_all@significance_results$union_random[[z]]$overlap,
         main = mains[z,1], xlab = "Overlap", xlim = c(0, 1),
         col = "azure3")
    abline(v = Sup_all@union_overlap$overlap[[z]], col = "darkred", lwd = 2)
    
    # Full winter vs spring
    hist(Sup_all@significance_results$full_random[[z]]$overlap, 
         main = mains[z,2], xlab = "Overlap", xlim = c(0, 1),
         col="azure3")
    abline(v = Sup_all@full_overlap$overlap[[z]], col = "darkred", lwd = 2)
  }
  par(mfrow = c(1,1))
  dev.off() 
  
  # Guardar .Rdata
  save(list = c("full.overlap", "union.overlap", "Sup_all", "vol.df"), 
       file = paste0(season.folder[f], "/S_overlap/", season.names[f],"_jaccard.Rdata"))
  
  print(paste(season.names[f], f, "done"))
}

# Grafico de los elipsoides ----------------------------------------------------------------------------
# ellipsenm::plot_overlap(object=Sup_all, niches = c(1,4), data = F, 
#                                        niche_col = c("#4a72b0", "#6ea96e", "#da9500", "#294029"), 
#                                        data_col = c("#548235", "#2a678c", "#e6960e", "red"), 
#                                        background = F, background_type = "full", 
#                                        background_col = viridis::viridis, proportion = 1,
#                                        change_labels = F, legend = T)






