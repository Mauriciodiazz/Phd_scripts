
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
  
  rm(list = setdiff(ls(), c("month.names", "month.folder", "combs.m")))
}

# Segunda version (Final?) ------------------------------------------------

# Esta es una nueva version, en donde se calcula el PCA para los registros de toda la especie y se optienen los valores correspondientes a cada mes. Previamente, había calculado un valor de PCA por elipsoide y al calcular el ECI los valores de D no eran comparables. 
# Hablar con el profe y borrar la primera versión. Esto permitirá que los PCA sean comparables para cacular el ECA (Environmaental Coincidence Average)

library(ecospat)
library(terra)
library(sf)
library(tidyverse)

# Path de las especies
month.names<-list.files("./species/mig_ENM/monthly/", full.names = F) |> str_remove("_ENMm")
month.folder<-list.files("./species/mig_ENM/monthly/", full.names = T)

# Tabla que organiza el bucle (combinación 1 entre 4)
combs.m<-data.frame(a= rep(1:11, 11:1), b= unlist(lapply(2:12, function(i) i:12)))

# # ---- borrar
# 
# a <- data.frame(spp=NA, dif=NA)
# a.list <- list()
# #---

data.ovrlp.list<-list()

# f<-97
for (f in 1:length(month.names)) { #1:length(month.names)
  
  # 1. Cargo el Rdata que contiene los ENM y el background ---------------------
  
  temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
  load(list.files(month.folder[f], pattern="*.Rdata$", full.names = T), envir = temp_env)
  
  # Objetos de temp_env que se van a usar
  
  abex.list<-temp_env$abex.list
  
  abbg.list <- list()
  for(x in 1:length(temp_env$abbg.list)) { 
    abbg.list[[x]] <-  temp_env$abbg.list[[x]]$env_bg
    names(abbg.list)[x]<-paste0("abbg_", month.names[f], "_", x)
  }
  
  # Se crea el directorio donde se guardan los datos
  # dir.create(paste0(month.folder[f], "/M_overlap"))
  
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
  
  bg_name <-  unique(pca.cal2$elip)[1:12] # It will be different for months 
  abex_name <-  unique(pca.cal2$elip)[13:24] # It will be different for months
  
  # Determinación de si el rango del bg es mas que el rango de los puntos
  
  min.max.bg <- pca.cal2 |> 
    filter(str_starts(elip, "abbg_")) |> 
    dplyr::summarise(min.a1=min(Axis1), 
                     max.a1=max(Axis1), 
                     min.a2=min(Axis2), 
                     max.a2=max(Axis2), .by = elip) |> 
    mutate(ID= str_extract_all(elip, "[0-9]+\\.?[0-9]*") |> as.numeric()) |> 
    arrange(ID) |> 
    select(!ID)
  
  min.max.sp <- pca.cal2 |> 
    filter(str_starts(elip, "abex_")) |> 
    dplyr::summarise(min.a1=min(Axis1), 
                     max.a1=max(Axis1), 
                     min.a2=min(Axis2), 
                     max.a2=max(Axis2), .by = elip) |> 
    mutate(ID= str_extract_all(elip, "[0-9]+\\.?[0-9]*") |> as.numeric()) |> 
    arrange(ID) |> 
    select(!ID)
  
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
    #   print(paste(month.names[f], f, "done"))
    # }
    # do.call(rbind, a.list) |>
    #   write.table("./borrar/bg_conflict.M.txt", sep="\t", dec=".", row.names=F)
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
    
    # Prueba de que funcionó
    
    min.max.bg <- pca.cal2 |>
      filter(str_starts(elip, "abbg_")) |>
      dplyr::summarise(min.a1=min(Axis1),
                       max.a1=max(Axis1),
                       min.a2=min(Axis2),
                       max.a2=max(Axis2), .by = elip) |> 
      mutate(ID= str_extract_all(elip, "[0-9]+\\.?[0-9]*") |> as.numeric()) |> 
      arrange(ID) |> 
      select(!ID)
    
    min.max.sp <- pca.cal2 |>
      filter(str_starts(elip, "abex_")) |>
      dplyr::summarise(min.a1=min(Axis1),
                       max.a1=max(Axis1),
                       min.a2=min(Axis2),
                       max.a2=max(Axis2), .by = elip) |> 
      mutate(ID= str_extract_all(elip, "[0-9]+\\.?[0-9]*") |> as.numeric()) |> 
      arrange(ID) |> 
      select(!ID)
    
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
      
      print(paste("records needs for correction"))
      
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
        mutate(ID= str_extract_all(elip2, "[0-9]+\\.?[0-9]*") |> as.numeric()) |> 
        arrange(ID) |> 
        select(!c(ID, elip))

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
  
  write.table(pca.cal2,
              paste0("./species/mig_ENM/overlaps_tables/PCA_m/", month.names[f],
                     "_m_PCA.txt"),
              sep="\t", dec = ".", row.names=F)
#   # -----
# } # final "prematuro" del bucle inicial
# # -----

  # 5. Select PCA values for each ellipsoid ------------------------------------
  
  # x<-1
  for (x in 1:nrow(combs.m)) {
    
    ##### 5.1 PCA climatic scores  ----------------------------------------------------
    
    # Backgrounds
    scores.clima1 <- pca.cal2 |> 
      filter(elip == bg_name[combs.m[x,1]]) |> 
      select(!c(elip, ID_YEAR))
    
    scores.clima2 <- pca.cal2 |> 
      filter(elip == bg_name[combs.m[x,2]]) |> 
      select(!c(elip, ID_YEAR))
    
    scores.clima12 <- rbind(scores.clima1, scores.clima2)
    
    # Environments
    scores.sp1a <- pca.cal2 |> 
      filter(elip == abex_name[combs.m[x,1]]) |> 
      select(!c(elip, ID_YEAR))
    
    scores.sp2b <- pca.cal2 |> 
      filter(elip == abex_name[combs.m[x,2]]) |> 
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
    data.overlap[x,2] <- paste0("m", combs.m[x, 1], "_m", combs.m[x, 2])
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
  
  write.table(data.overlap, 
              paste0("./species/mig_ENM/overlaps_tables/B_ovl_month/",
                     month.names[f],
                     "_m_ovl.txt"), sep = "\t", dec = ".", row.names = F)
  data.ovrlp.list[[f]] <- data.overlap
  
  print(paste(basename(month.folder[f]), f, "done"))
  rm(list = setdiff(ls(), c("month.names", "month.folder", "combs.m", "data.ovrlp.list")))
  
}
### ----



# BD completa de superposicion -----------------------------------------------
df.path <- list.files("./species/mig_ENM/overlaps_tables/month/", full.names = T)

data.ovrlp.list<-list()

for (x in 1:length(df.path)) {
  data.ovrlp.list[[x]] <- read.table(df.path[x], sep="\t", dec=".", header=T)
}

data.ovrlp.df<-do.call(rbind, data.ovrlp.list)
# write.table(data.ovrlp.df, "./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt", sep="\t", dec = ".", row.names=F)

# data.ovrlp.df<-read.table("./species/mig_ENM/overlaps_tables/data.ovrlp_m_df.txt", sep="\t", dec = ".", header=T)

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


# Comparar los resultados -------------------------------------------------

old.paths <- list.files("D:/borrar/mig_ENM_prev_bkp/overlaps_tables/month/", full.names=T)
new.paths <- list.files('F:/Phd_DD/data_Phd_DD/species/mig_ENM/overlaps_tables/month/', full.names=T)

old.l <- list()
new.l <- list()

for (x in 1:length(old.paths)) {
  old.l[[x]] <- read.table(old.paths[x], dec='.', sep='\t', header=T)
  new.l[[x]] <- read.table(new.paths[x], dec='.', sep='\t', header=T)
}


old <- do.call(rbind, old.l) |> 
  mutate(type='old',
         overlap = overlap |> str_replace_all("m", "M")) |> 
  tibble()

new <- do.call(rbind, new.l) |> 
  mutate(type='new',
         spp = spp |> str_remove_all("_ENMm")) |> 
  tibble()


# Differences graph -------------------------------------------------------

tibble(spp = new$spp, overlap = new$overlap, dif=new$obs.D-old$obs.D) |> 
  ggplot(aes(x=dif, y=spp)) +
  geom_boxplot() +
  geom_vline(xintercept = 0, color="red", linetype="dashed") +
  labs(x="Dif. (New-Old)", y="Species") +
  theme_minimal()

wilcox.test(new$obs.D, old$obs.D)


# coincidencia ambiental (prueba) -----------------------------------------

new |>
# old |> 
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

ggsave(file = "F:/Phd_DD/data_Phd_DD/outputs/images/ieca_new.png",
       width = 10,
       height = 7,
       scale=3,
       units ="cm")





