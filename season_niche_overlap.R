
# Comparación de nichos con Broennimann et al. 2012, Di Cola et al. 2017 

# Esta es una nueva version, en donde se calcula el PCA para los registros de toda la especie y se optienen los valores correspondientes a cada mes. Previamente, había calculado un valor de PCA por elipsoide y al calcular el ECI los valores de D no eran comparables. 
# Hablar con el profe y borrar la primera versión. Esto permitirá que los PCA sean comparables para cacular el ECA (Environmaental Coincidence Average)

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

# registros borrados en los filtros
spp.fil <- data.frame(spp=NA, fil1=NA, fill2=NA)

# f <- 7
for (f in 1:length(season.names)) { #1:length(season.names)
  
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
  # dir.create(paste0(season.folder[f], "/S_overlap"))
  
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
  
  # ----
  spp.fil[f,1] <- season.names[f]
    
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
        
        # Mi pc
        # r.stack <- rast(list.files(dat.path[r] |> 
        #                              str_replace("E:", "D:"), 
        #                            full.names = T))
        
        # PC del lab
        r.stack <- rast(list.files(dat.path[r], 
                                   full.names = T))
        
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
    
    #---
    spp.fil[f,2] <- pca.cal2 |> 
      filter(str_starts(elip, "abex_")) |> 
      nrow()
    #---
    
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
      
      print(paste("records for need correction"))
      
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
      
      #---
      spp.fil[f,3] <- pca.cal2 |> 
        filter(str_starts(elip, "abex_")) |> 
        nrow()
      #---
      
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
  
  write.table(pca.cal2,
              paste0("./species/mig_ENM/overlaps_tables/PCA_s/", season.names[f] |> str_remove_all("_ENMm"),
                     "_s_PCA.txt"),
              sep="\t", dec = ".", row.names=F)
  # # -----
  # } # final "prematuro" del bucle inicial
  # # -----
  
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

#-----------------


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


# Jackard overlap using overllip package (Osorio-Olvera, 2020) ------------

# devtools::install_github("luismurao/overllip")

library(overllip)
library(tidyverse)

pca.list.s <- list.files("./species/mig_ENM/overlaps_tables/PCA_s/", full.names = T)
pca.list.m <- list.files("./species/mig_ENM/overlaps_tables/PCA_m/", full.names = T)

jack_res <- as.data.frame(matrix(ncol=11, nrow=103))
int_vols <- as.data.frame(matrix(ncol=11, nrow=103))
names(jack_res) <- c("species", "total_vs_W", "total_vs_T1", "total_vs_S",
                     "total_vs_T2", "W_vs_T1", "W_vs_S", "W_vs_T2", "T1_vs_S",
                     "T1_vs_T2", "S_vs_T2")
names(int_vols) <- names(jack_res)

list.ovlps <- list()

# x <- 67
for (x in 1:103) {
  
  jack_res[x,1] <- basename(pca.list.s[x]) |> 
    str_remove("_s_PCA.txt")
  int_vols[x,1] <- basename(pca.list.s[x]) |> 
    str_remove("_s_PCA.txt")
  
  
  # Load species data -------------------------------------------------------
  pca.s <- read_tsv(pca.list.s[x])
  
  # Filter background --------------------------------------------------------
  
  env.pca.s <- 
    pca.s |> 
    filter(str_starts(elip, 'abbg')) |> 
    mutate(elip2 = elip |> str_remove(paste0("abbg_", basename(pca.list.s[x]) |> 
                                               str_remove("s_PCA.txt")))) |> 
    mutate(elip2 = case_when(
      elip2 == 1 ~ "W",
      elip2 == 2 ~ "T1",
      elip2 == 3 ~ "S",
      elip2 == 4 ~ "T2")) |> 
    mutate(elip2=factor(elip2, levels=c("W", "T1", "S", "T2")))
  
  # Filter total species ellipsoid ------------------------------------------
  
  spp.pca.s <- pca.s |> 
    filter(str_starts(elip, 'abex')) |> 
    mutate(elip2 = elip |> str_remove(paste0("abex_", basename(pca.list.s[x]) |> 
                                               str_remove("s_PCA.txt")))) |> 
    mutate(elip2 = case_when(
      elip2 == 1 ~ "W",
      elip2 == 2 ~ "T1",
      elip2 == 3 ~ "S",
      elip2 == 4 ~ "T2")) |> 
    mutate(elip2=factor(elip2, levels=c("W", "T1", "S", "T2")))
  
  tot <- spp.pca.s[,1:2]
  cen.tot <- colMeans(tot) 
  sig.tot <- cov(tot)
  
  # Overllip total object ---------------------------------------------------
  
  ellip_total <- overllip::ellipsoid_data(centroid=cen.tot,
                                          sigma = sig.tot, 
                                          cf=0.95)
  
  # Filter period species ellipsoids ----------------------------------------
  
  periods <- spp.pca.s$elip2 |> unique() |> as.character()
  ellip_list.s <- list()
  
  for(p in 1:length(periods)){
    est <- spp.pca.s |> 
      filter(elip2==periods[p]) |> 
      select(Axis1,Axis2)
    cen.est <- colMeans(est) 
    sig.est <- cov(est) 
    
    ellip_list.s[[p]] <- overllip::ellipsoid_data(centroid=cen.est,
                                                  sigma = sig.est, 
                                                  cf=0.95)
    names(ellip_list.s)[p] <- periods[p]
  }
  
  
  # Ellipsoids stack --------------------------------------------------------
  
  ellipsoid_stack <- overllip::stack(ellip_total,
                                     ellip_list.s[[1]],
                                     ellip_list.s[[2]],
                                     ellip_list.s[[3]],
                                     ellip_list.s[[4]],
                                     ellipsoid_names = c("total",periods))
  
  # Unioning data (Hypercube) -----------------------------------------------
  
  env_data <- overllip::hypercube(ellipsoids = ellipsoid_stack,
                                  n = 10000000)
  
  
  # Ellipsoid comparison ----------------------------------------------------
  
  rmat <- overllip::stack_overlap(ellipsoid_stack = ellipsoid_stack,
                                  env_data =  env_data,
                                  parallel = T,
                                  ncores = 6)
  
  # Save data ---------------------------------------------------------------
  
  jack_res[x,2:11] <- rmat@Jackard_indices
  int_vols[x,2:11] <- rmat@intersection_volumes
  
  list.ovlps[[x]] <- rmat
  names(list.ovlps)[x] <-  basename(pca.list.s[x]) |> str_remove("_s_PCA.txt")
  
  print(paste(basename(pca.list.s[x]) |> str_remove("_s_PCA.txt"), x, "overlaped done"))
  
  jack_res |>
    write.table(paste0("./species/mig_ENM/overlaps_tables/jack_res_Ts.txt"), 
                sep="\t",
                dec=".",
                row.names=F)
  int_vols |>
    write.table(paste0("./species/mig_ENM/overlaps_tables/int_vols_Ts.txt"), 
                sep="\t",
                dec=".",
                row.names=F)
  
  rm(list = setdiff(ls(), c("jack_res", "int_vols", "list.ovlps", "pca.list.s")))
}

save.image("./species/mig_ENM/overlaps_tables/jack_int_s.Rdata")



# Proporción de Nichos estacionales (ellipsenm) ---------------------------

library(ellipsenm)
library(tidyverse)
library(rgl)


### Cargar lista de especies ------------------------------------------------

season.names<-list.files("./species/mig_ENM/seasonaly/", full.names = F) |> str_remove("_ENMs")
season.folder<-list.files("./species/mig_ENM/seasonaly/", full.names = T)


# f <- 1
for(f in 5:103){#length(season.folder)

### Cargar Rdata de la especie ---------------------------------------------

  temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
  load(list.files(season.folder[f], pattern="*.Rdata$", full.names = T), envir = temp_env)
  
  
  abex.list <- temp_env$abex.list 
  
  list.temp <- list()
  
  for(x in 1:length(abex.list)){
    list.temp[[x]] <- abex.list[[x]]$temporal_df |> 
      mutate(spp = season.names[f],
             period = x) |> 
      relocate(c(spp, period))
  }
  
  abex.tot <- do.call(rbind, list.temp)
  
  s.nam <- abex.tot$period |> unique()
  
  abex.tot.ovl <- abex.tot |> 
    select(spp, lon, lat, prec, tmax, tmin)
  

### Elipsoide total ---------------------------------------------------------

  elip_tot<-overlap_object(data=abex.tot.ovl,
                           species="spp",
                           longitude="lon",
                           latitude = "lat", 
                           method = "mve1", 
                           level = 99)
  

### Elipsoides por periodo --------------------------------------------------
  
  elip_s <- list()
  # y <- 2
  for(y in 1:length(s.nam)){
    abex.s <- abex.tot |> 
      filter(period == y) |> 
      select(spp, lon, lat, prec, tmax, tmin)
    
    elip_s[[y]]<-overlap_object(data=abex.s,
                                species="spp",
                                longitude="lon",
                                latitude = "lat", 
                                method = "mve1", 
                                level = 99)
  }
  

### Sobreposición de los ambientes ------------------------------------------

  overlap_t <- ellipsoid_overlap(elip_tot,
                                 elip_s[[1]], # 1 siempre es el total
                                 elip_s[[2]],
                                 elip_s[[3]],
                                 elip_s[[4]])

### Almacenar información ---------------------------------------------------

  olv.data <- overlap_t@full_overlap |> 
    mutate(spp=season.names[f]) |> 
    rownames_to_column(var="ovl.type") |> 
    relocate(spp, ovl.type) 
  
  write.table(olv.data,
              paste0("./species/mig_ENM/overlaps_tables/prop_ovl_season/", season.names[f] |> str_remove_all("_ENMs"),
                     "_s_prop.txt"),
              sep="\t", dec = ".", row.names=F)
  print(paste(season.names[f], f, "done"))
  
}