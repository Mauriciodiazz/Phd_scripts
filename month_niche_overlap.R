
# Comparación de nichos con Broennimann et al. 2012, Di Cola et al. 2017 

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
  
  ## 2.1 Background ----------------------------------------------------------
  
  bgd <- 
    do.call(rbind, abbg.list) |> 
    select(lon, lat, prec, tmax, tmin, ID_YEAR) |> 
    rownames_to_column(var = "elip") |> 
    mutate(elip = elip |> str_remove("\\..*"))
  
  
  ## 2.2 Ambiente ------------------------------------------------------------
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
  
  # 3. Construcción del PCA ----------------------------------------------------
  
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
    
    #    # --- borrar
    #    # Esto es para crear una tabla que me diga donde hay estos conflictos
    #    a.list[[f]] <- data.frame(spp=elips.prob, ID=rep(f, length(elips.prob)))
    #   }
    #   print(paste(month.names[f], f, "done"))
    # }
    # do.call(rbind, a.list) |>
    #   write.table("./borrar/bg_conflict.M.txt", sep="\t", dec=".", row.names=F)
    #    elips.prob <- c("abbg_Vir_solit_3", "abbg_Vir_solit_2")
    # 
    #    # ---
    
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

  # 5. Select PCA values for each ellipsoid ------------------------------------
  
  # x<-1
  for (x in 1:nrow(combs.m)) {
    
    ## 5.1 PCA climatic scores  ----------------------------------------------------
    
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
    
    
    # 5.2 Superficie de densidad de registros -------------------------------------
    
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
    
    # 5.3 PCA individual Graphs ---------------------------------------------------------
    
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
    
    
    # 5.4 Overlap graph -----------------------------------------------------------
    
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
    
    # 5.5 Similarity test (di Cola et al 2017 - higher) ----------------------
    
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



# Jackard overlap using overllip package (Osorio-Olvera et al., 2020) ------------

library(overllip)
library(tidyverse)

pca.list.m <- list.files("./species/mig_ENM/overlaps_tables/PCA_m/", full.names = T)

jack_res <- as.data.frame(matrix(ncol=79, nrow=103))
int_vols <- as.data.frame(matrix(ncol=79, nrow=103))
names(jack_res) <- c("species", "total_vs_Jan", "total_vs_Feb", "total_vs_Mar", "total_vs_Apr", "total_vs_May", "total_vs_Jun", "total_vs_Jul", "total_vs_Aug", "total_vs_Sep", "total_vs_Oct", "total_vs_Nov", "total_vs_Dec", "Jan_vs_Feb", "Jan_vs_Mar", "Jan_vs_Apr", "Jan_vs_May", "Jan_vs_Jun", "Jan_vs_Jul", "Jan_vs_Aug", "Jan_vs_Sep", "Jan_vs_Oct", "Jan_vs_Nov", "Jan_vs_Dec", "Feb_vs_Mar", "Feb_vs_Apr",
"Feb_vs_May", "Feb_vs_Jun", "Feb_vs_Jul", "Feb_vs_Aug", "Feb_vs_Sep", "Feb_vs_Oct", "Feb_vs_Nov", "Feb_vs_Dec", "Mar_vs_Apr", "Mar_vs_May", "Mar_vs_Jun", "Mar_vs_Jul", "Mar_vs_Aug", "Mar_vs_Sep", "Mar_vs_Oct", "Mar_vs_Nov", "Mar_vs_Dec", "Apr_vs_May", "Apr_vs_Jun", "Apr_vs_Jul", "Apr_vs_Aug", "Apr_vs_Sep", "Apr_vs_Oct", "Apr_vs_Nov", "Apr_vs_Dec", "May_vs_Jun", "May_vs_Jul", "May_vs_Aug", "May_vs_Sep", "May_vs_Oct", "May_vs_Nov", "May_vs_Dec", "Jun_vs_Jul", "Jun_vs_Aug", "Jun_vs_Sep", "Jun_vs_Oct", "Jun_vs_Nov", "Jun_vs_Dec", "Jul_vs_Aug", "Jul_vs_Sep", "Jul_vs_Oct", "Jul_vs_Nov", "Jul_vs_Dec", "Aug_vs_Sep", "Aug_vs_Oct", "Aug_vs_Nov", "Aug_vs_Dec", "Sep_vs_Oct", "Sep_vs_Nov", "Sep_vs_Dec", "Oct_vs_Nov", "Oct_vs_Dec", "Nov_vs_Dec")
names(int_vols) <- names(jack_res)

list.ovlps <- list()

# x <- 1
for (x in 82:103) {
  
  jack_res[x,1] <- basename(pca.list.m[x]) |> 
    str_remove("_m_PCA.txt")
  int_vols[x,1] <- basename(pca.list.m[x]) |> 
    str_remove("_m_PCA.txt")
  
  ## Load species data -------------------------------------------------------
  
  pca.s <- read_tsv(pca.list.m[x])
  
  ## Filter background --------------------------------------------------------
  
  env.pca.s <- 
    pca.s |> 
    filter(str_starts(elip, 'abbg')) |> 
    mutate(elip2 = elip |> str_remove(paste0("abbg_", basename(pca.list.m[x]) |> 
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
      elip2 == 12 ~ "Dec")) |> 
    mutate(elip2=factor(elip2, 
                        levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
                                 "Sep","Oct","Nov","Dec")))
  
  ## Filter total species ellipsoid ------------------------------------------
  
  spp.pca.s <- pca.s |> 
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
      elip2 == 12 ~ "Dec")) |> 
    mutate(elip2=factor(elip2, 
                        levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",
                                 "Sep","Oct","Nov","Dec")))
  
  tot <- spp.pca.s[,1:2]
  cen.tot <- colMeans(tot) 
  sig.tot <- cov(tot)
  
  ## Overllip total object ---------------------------------------------------
  
  ellip_total <- overllip::ellipsoid_data(centroid=cen.tot,
                                          sigma = sig.tot, 
                                          cf=0.95)
  
  ## Filter period species ellipsoids ----------------------------------------
  
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
  
  
  ## Ellipsoids stack --------------------------------------------------------
  
  ellipsoid_stack <- overllip::stack(ellip_total,
                                     ellip_list.s[[1]],
                                     ellip_list.s[[2]],
                                     ellip_list.s[[3]],
                                     ellip_list.s[[4]],
                                     ellip_list.s[[5]],
                                     ellip_list.s[[6]],
                                     ellip_list.s[[7]],
                                     ellip_list.s[[8]],
                                     ellip_list.s[[9]],
                                     ellip_list.s[[10]],
                                     ellip_list.s[[11]],
                                     ellip_list.s[[12]],
                                     ellipsoid_names = c("total",periods))
  
  ## Unioning data (Hypercube) -----------------------------------------------
  
  env_data <- overllip::hypercube(ellipsoids = ellipsoid_stack,
                                  n = 10000000)
  
  ## Ellipsoid comparison ----------------------------------------------------
  
  rmat <- overllip::stack_overlap(ellipsoid_stack = ellipsoid_stack,
                                  env_data =  env_data,
                                  parallel = T,
                                  ncores = 6)
  
  ## Save data ---------------------------------------------------------------
  
  jack_res[x,2:79] <- rmat@Jackard_indices
  int_vols[x,2:79] <- rmat@intersection_volumes
  
  list.ovlps[[x]] <- rmat
  names(list.ovlps)[x] <-  basename(pca.list.m[x]) |> str_remove("_m_PCA.txt")
  
  print(paste(basename(pca.list.m[x]) |> 
                str_remove("_m_PCA.txt"), x, "overlaped done"))
  
  jack_res |>
    write.table("./species/mig_ENM/overlaps_tables/jack_res_Tm2.txt", 
                sep="\t",
                dec=".",
                row.names=F)
  int_vols |>
    write.table("./species/mig_ENM/overlaps_tables/int_vols_Tm2.txt", 
                sep="\t",
                dec=".",
                row.names=F)
  
  rm(list = setdiff(ls(), c("jack_res", "int_vols", "list.ovlps", "pca.list.m")))
}

# save.image("./species/mig_ENM/overlaps_tables/jack_int_m.Rdata")



# Proporción de Nichos mensual (ellipsenm) -----------------------------------------

month.names<-list.files("./species/mig_ENM/monthly/", full.names = F) |> str_remove("_ENMm")
month.folder<-list.files("./species/mig_ENM/monthly/", full.names = T)


# f <- 1
for(f in 72:103){#length(month.folder)
  
  ## Cargar Rdata de la especie ---------------------------------------------
  
  temp_env <- new.env() # Esto me crea un ambiente en un conjunto temporal
  load(list.files(month.folder[f], pattern="*.Rdata$", full.names = T), envir = temp_env)
  
  abex.list <- temp_env$abex.list 
  
  list.temp <- list()
  
  for(x in 1:length(abex.list)){
    list.temp[[x]] <- abex.list[[x]]$temporal_df |> 
      mutate(spp = month.names[f],
             period = x) |> 
      relocate(c(spp, period))
  }
  
  abex.tot <- do.call(rbind, list.temp)
  
  s.nam <- abex.tot$period |> unique() #ANTES DE CORRERLO REVISAR ESTE, TIENE QUE TENER 12 OBJETOS
  
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
    abex.m <- abex.tot |> 
      filter(period == y) |> 
      select(spp, lon, lat, prec, tmax, tmin)
    
    elip_s[[y]]<-overlap_object(data=abex.m,
                                species="spp",
                                longitude="lon",
                                latitude = "lat", 
                                method = "mve1", 
                                level = 99)
  }
  
  ### Sobreposición de los ambientes ------------------------------------------
  
  overlap_t <- ellipsoid_overlap(elip_tot,
                                 elip_s[[1]],
                                 elip_s[[2]],
                                 elip_s[[3]],
                                 elip_s[[4]],
                                 elip_s[[5]],
                                 elip_s[[6]],
                                 elip_s[[7]],
                                 elip_s[[8]],
                                 elip_s[[9]],
                                 elip_s[[10]],
                                 elip_s[[11]],
                                 elip_s[[12]])
  
  ### Almacenar información ---------------------------------------------------
  
  olv.data <- overlap_t@full_overlap |> 
    mutate(spp=month.names[f]) |> 
    rownames_to_column(var="ovl.type") |> 
    relocate(spp, ovl.type) 
  
  write.table(olv.data,
              paste0("./species/mig_ENM/overlaps_tables/prop_ovl_month/", month.names[f] |> str_remove_all("_ENMm"),
                     "_m_prop.txt"),
              sep="\t", dec = ".", row.names=F)
  print(paste(month.names[f], f, "done"))
  
}
