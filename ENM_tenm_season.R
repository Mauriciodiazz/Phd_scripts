# ----------------- TENM bucle -----------------------

library(tenm)
library(terra)
library(tidyverse)
library(sf)

# Open species file
tenm_mask <- terra::rast(list.files("E:/CHELSA_ym_5k/1901-01/", full.names=T, pattern=".tif$"))
names(tenm_mask)<-names(tenm_mask) |> 
  str_remove("CHELSAcruts_") |> 
  str_remove("_V.1.0")

# Funcion para generar el background --------------------------------------
# las celdas vacías de CHELSAcruts estan indexadas como valores nulos o NA (NaN) y eso genera conflicto en la conformación del backgorund de este paquete. La función establece pixeles al rededor de los puntos de presencia para obtener el backgorund, lo que hice fue modificarlo de manera que cuando detectara valores ausentes en al menos una de las capas en esa coordenada eliminara el id de esas celdas (que era el mismo basicamente), porque estaba cogiendo el valor de los pixeles que tuvieran ambiente así solo un raster tuviera valor ambiental.

# Modificación funciones
# Aquella que genera los puntos de backgorund
backg_CHcr<-function (this_species, buffer_ngbs = NULL, buffer_distance = 1000, 
                      n_bg = 50000, process_ngbs_by = 100) 
{
  stopifnot(inherits(this_species, "sp.temporal.env"))
  tdf <- this_species$temporal_df
  samp_prop <- layer_dates <- var_name <- layer_val <- NULL
  df_samps <- dplyr::summarise(dplyr::group_by(tdf, layer_dates), 
                               samp_prop = dplyr::n()/nrow(tdf), n_samples = ceiling(samp_prop * 
                                                                                       n_bg))
  paths_layers <- purrr::map_df(split(tdf, tdf$layers_path), 
                                function(x) {
                                  layer_path <- list.files(x$layers_path[1], pattern = this_species$layers_ext, 
                                                           full.names = TRUE)[1]
                                  layer_path_df <- data.frame(layers_path = x$layers_path, 
                                                              layer_path)
                                  return(layer_path_df)
                                })
  tdf$layer_path <- paths_layers$layer_path
  ddL <- split(tdf, tdf$layer_path)
  if (!is.null(buffer_ngbs)) {
    cells_to_samp <- purrr::map(seq_along(ddL), function(z) {
      cell_ids <- tenm::cells2samp(data = ddL[[z]], longitude = this_species$lon_lat_vars[1], 
                                   latitude = this_species$lon_lat_vars[2], cell_ids = ddL[[z]]$cell_ids_year, 
                                   buffer_ngbs = buffer_ngbs, raster_mask = terra::rast(ddL[[z]]$layer_path[1]), 
                                   n_bg = df_samps$n_samples[z], process_ngbs_by = process_ngbs_by, 
                                   progress = FALSE)
      return(cell_ids)
    })
  }
  else if (is.null(buffer_ngbs) && is.numeric(buffer_distance)) {
    cells_to_samp <- purrr::map(seq_along(ddL), function(z) {
      occ_va <- ddL[[z]][, c(this_species$lon_lat_vars[1], 
                             this_species$lon_lat_vars[2])]
      r <- terra::rast(ddL[[z]]$layer_path)
      crsL <- terra::crs(r)
      vec_occ <- terra::vect(as.matrix(occ_va), crs = crsL)
      buff_te <- terra::buffer(vec_occ, width = buffer_distance)
      cell2sa <- terra::cells(r, buff_te)[, 2]
      nsamples <- ifelse(df_samps$n_samples[z] > length(cell2sa), 
                         length(cell2sa), df_samps$n_samples[z])
      cell_ids <- sample(cell2sa, nsamples)
      return(cell_ids)
    })
  }
  gc()
  dir_paths <- unique(tdf$layers_path)
  all_layers <- purrr::map_df(seq_along(dir_paths), function(x) {
    cp <- list.files(dir_paths[x], pattern = this_species$layers_ext, 
                     full.names = TRUE, recursive = FALSE)
    data.frame(dir_paths = dir_paths[x], layers_path = cp)
  })
  names(cells_to_samp) <- dir_paths
  ex_date <- furrr::future_map_dfr(seq_len(nrow(all_layers)), 
                                   function(x) {
                                     cellids <- cells_to_samp[[all_layers$dir_paths[x]]]
                                     # env_layers <- terra::rast(all_layers$layers_path[x])
                                     # En vez de muestrear una capa, las extraigo todas y así evito que un pixel                                        tenga solo un valor de las variables y tenga climas "incompletos"
                                     env_layers <- terra::rast(list.files(all_layers$dir_paths[x],
                                                                          full.names=T, 
                                                                          pattern=this_species$layers_ext))
                                     layer_val <- stats::na.omit(env_layers[cellids])
                                     rm_ids <- stats::na.action(layer_val)
                                     if (length(rm_ids) > 0) {
                                       cellids <- cellids[-rm_ids]
                                     }
                                     if (nrow(layer_val) == 0L) 
                                       return()
                                     snam <- paste0("env_layers@", names(attributes(env_layers))[1])
                                     snam <- eval(parse(text = paste0(snam, "$get_sourcenames()")))
                                     # df1 <- data.frame(ID_YEAR = all_layers$dir_paths[x],
                                     #                   cellids, layer_val = layer_val[[1]], var_name = snam)
                                     # Crear un df que contiene las rutas y los ID de las celdas
                                     df3 <- data.frame(ID_YEAR = all_layers$dir_paths[x],
                                                       cellids)
                                     # unirlo con los valores extraidos de cada variable
                                     df2<-cbind(df3, layer_val)
                                     # cambiar el nombre y darle el asignado a snam
                                     names(df2)[-c(1,2)]<-snam
                                     # dejarlo en el mismo formato
                                     df1<- df2 |> 
                                       tidyr::pivot_longer(!c(ID_YEAR, cellids), 
                                                           values_to = "layer_val",
                                                           names_to = "var_name")
                                     return(df1)
                                   }, .progress = TRUE, .options = furrr::furrr_options(seed = NULL))
  gc() 
  # ex_date genera valores duplicados, con esto tengo valores unicos por variable
  ex_date<- ex_date |> 
    dplyr::group_by(cellids) |> 
    dplyr::distinct(ID_YEAR, cellids, var_name, layer_val) |> 
    as.data.frame()
  bg_env <- tidyr::pivot_wider(ex_date, values_fn = list, 
                               names_from = var_name, values_from = layer_val)
  r1 <- terra::rast(all_layers$layers_path[1])
  xys <- terra::xyFromCell(r1, bg_env$cellids)
  colnames(xys) <- this_species$lon_lat_vars
  bg_env <- furrr::future_map_dfr(seq_along(bg_env$ID_YEAR), 
                                  function(x) {
                                    ID_YEAR <- rep(bg_env$ID_YEAR[[x]], length(bg_env[[3]][[x]]))
                                    df_year <- purrr::map_dfc(seq_along(bg_env[-(1:2)]), 
                                                              function(y) {
                                                                data <- bg_env[-(1:2)]
                                                                value <- data[[y]][[x]]
                                                                df1 <- data.frame(value)
                                                                names(df1) <- names(data[y])
                                                                return(df1)
                                                              })
                                    df_res <- data.frame(ID_YEAR, df_year)
                                    return(df_res)
                                  }, .progress = TRUE, .options = furrr::furrr_options(seed = NULL, 
                                                                                       globals = c("bg_env")))
  bg_env <- data.frame(ID_YEAR = bg_env[, 1], xys, bg_env[, 
                                                          -1])
  sp.temp.data.env <- list(temporal_df = tdf, sp_date_var = this_species$sp_date_var, 
                           lon_lat_vars = this_species$lon_lat_vars, layers_ext = this_species$layers_ext, 
                           env_bg = bg_env)
  class(sp.temp.data.env) <- c("sp.temporal.env", "sp.temporal.bg")
  return(sp.temp.data.env)
}
# Aquella que genera el SDW de los puntos de backgorund
tdf2swd2<-function (this_species, sp_name = "sp") 
{
  cls <- class(this_species)
  if (length(cls) == 1L && methods::is(this_species, "sp.temporal.env")) {
    df_base <- this_species$temporal_df
    swd_df <- data.frame(sp_name, df_base[, -c(4, 5, 6, 
                                               ncol(df_base))])
  }
  else if ("sp.temporal.bg" %in% cls) {
    df_base <- this_species$env_bg
    years <- basename(df_base$ID_YEAR) #as.numeric(basename(df_base$ID_YEAR))  = as.numeric("1995-11") generan NA
    swd_df <- data.frame(sp_name = "background", df_base[,c(2, 3)], year = years, df_base[, -c(1, 2, 3)])
  }
  else {
    stop("this_species should be of class sp.temporal.env or sp.temporal.bg")
  }
  return(swd_df)
}

spp.list<-list.files("./species/mig_shapes/Selected/", full.names = T, pattern=".shp$")
spp.names<-list.files("./species/mig_shapes/Selected/", full.names = F, pattern=".shp$") |> str_remove(".shp")

#s<-68 # ok :65, 68 - Setophagas: 68-85
#spp.list[c(65,68:71)]

for (s in c(68:85)) {
  print(spp.names[s])
  
  spp<-vect(spp.list[s])
  
  spp2<-terra::extract(tenm_mask, spp, ID=F, bind=T)
  
  # lists
  #m<-1
  abt.list<-list()
  abex.list<-list()
  abbg.list<-list()
  bg_swd.list<-list()
  
  mod_sel.list<-list()
  mods_table.list<-list()
  
  # spp folder
  # OJO!!!!!!!!!!!!!!
  spp.folder<-paste0(spp.names[s],"_ENMs")
  # load(paste0("./species/mig_ENM/seasonaly/", spp.folder,"/ENM_s_", spp.names[s],".Rdata"))
  #dir.create(paste0("./species/mig_ENM/seasonaly/", spp.folder)) # Seasonal ENM 
  
  
  for (m in 1:4) {
    # Transformin it as data.frame
    spp.df<-
      spp2 |>
      st_as_sf() |> 
      st_coordinates() |> 
      cbind(st_as_sf(spp2)) |>
      drop_na("prec_1_1901", "tmax_1_1901", "tmin_1_1901") |> 
      rename("lon" = X, "lat" = Y) |>
      as.data.frame() |>
      mutate(season=case_when(
        month %in% c(12,1,2) ~ 1,
        month %in% 3:5 ~ 2,
        month %in% 6:8 ~ 3,
        month %in% 9:11 ~ 4)) |> 
      select(species, lon, lat, year, month, season) |>
      filter(season==m) |>  ### SEASON FILTERING
      filter(year>=1901 & year<=2016) |> 
      unite("date", year, month, sep="-", remove=F) |> 
      mutate(date = str_replace(date, "-(\\d)$", "-0\\1")) |> # date format YYYY-MM
      mutate(date=paste(date, "01", sep="-")) 
    
    # We indicate the path where our time-specific modeling layers are located (variables).
    tempora_layers_dir <- "E:/CHELSA_ym_5k/"
    
    # sp_temporal_data will allow us to work with time-specific data
    abt <- tenm::sp_temporal_data(occs = spp.df,
                                  longitude = "lon",
                                  latitude = "lat",
                                  sp_date_var = "date",
                                  occ_date_format="ymd",
                                  layers_date_format= "ym",
                                  layers_by_date_dir = tempora_layers_dir,
                                  layers_ext="*.tif$")
    
    # Save in a list object per specie, which will containing the monthly objects
    abt.list[[m]]<-abt
    names(abt.list)[m]<-paste0("abt_", spp.names[s], "_", m)
    
    # Time-specific spatial data thinning -------------------------------------
    
    # Clean duplicates using a raster mask
    abtc <- tenm::clean_dup_by_date(this_species = abt,
                                    by_mask = TRUE,
                                    threshold = terra::res(tenm_mask)[1],
                                    raster_mask = tenm_mask[[1]], ##
                                    n_ngbs = 0)
    
    # # Check number of records
    # head(tidyr::as_tibble(abtc$temporal_df))
    # nrow(tidyr::as_tibble(abtc$temporal_df))
    # 
    # nrow(abtc$temporal_df)
    # nrow(abt$temporal_df)
    
    # Time-specific environmental data extraction -----------------------------
    print(paste("abex", m))
    future::plan("multisession",workers=5) # Allow that ex_by_date could be run in parallel
    abex <- tenm::ex_by_date(this_species = abtc,
                             train_prop=0.7) # train=0.7 and test=0.3
    future::plan("sequential")
    
    abex.list[[m]]<-abex
    names(abex.list)[m]<-paste0("abex_", spp.names[s], "_", m)
    
    # Time-specific background generation -------------------------------------
    
    print(paste("abbg", m))
    future::plan("multisession",workers=5)
    abbg <- backg_CHcr(this_species = abex,
                       buffer_ngbs=50,
                       n_bg=50000, ### 
                       buffer_distance=5000)
    future::plan("sequential")
    
    abbg.list[[m]]<-abbg
    names(abbg.list)[m]<-paste0("abbg_", spp.names[s], "_", m)

    # Exporting time-specific information as Samples With Data (SWD) format ---------
    
    # SWD table for occurrence records
    occ_swd <- tdf2swd2(this_species = abex, 
                       sp_name = "Sel_calli")
    
    occ_swd |> 
      mutate(season=m) |> 
      relocate(season, .after=sp_name) |> 
      write.table(paste0("./species/mig_ENM/seasonaly/", spp.folder, "/", 
                         spp.names[s], "_ENMs_", m, ".txt"), 
                  sep="\t",
                  dec=".",
                  row.names=F)
    
    # SWD table for background data
    bg_swd <- tdf2swd2(this_species = abbg) #missing year values
    # head(tidyr::as_tibble(occ_swd))
    # head(tidyr::as_tibble(bg_swd))
    
    bg_swd.list[[m]]<-bg_swd
    names(bg_swd.list)[m]<-paste0("bg_swd_", spp.names[s], "_", m)
    
    # Time-specific model calibration and selection ---------------------------
    vars2fit<- c("prec", "tmin", "tmax")
    
    mod_sel <- tenm::tenm_selection(this_species = abbg,
                                    omr_criteria =1,
                                    ellipsoid_level= 0.975, #by default
                                    vars2fit = vars2fit,
                                    nvars_to_fit=3,
                                    proc = T,
                                    RandomPercent = 50,
                                    NoOfIteration=1000,
                                    parallel=TRUE,
                                    n_cores=4)
    
    #mod_sel$mods_table #resultados de la evalaución
    
    mod_sel.list[[m]]<-mod_sel
    mods_table.list[[m]]<-
      mod_sel$mods_table |>
      mutate(season=m,
             species=spp.names[s]) |> 
      relocate(species, season)
    
  }
 
  # drop list
  mods_table.df<-do.call(rbind, mods_table.list) |> 
    select(!c(non_pred_test_ids,non_pred_train_ids))
  
  write.table(mods_table.df, 
              paste0("./species/mig_ENM/seasonaly/", spp.folder,"/ENM_mods_s_", spp.names[s],".txt"),
              sep="\t",
              dec=".",
              row.names=F)

  # Plotting ellipsoid ------------------------------------------------------
  
  x1 <- abex.list[[1]][[1]]$tmin
  y1 <- abex.list[[1]][[1]]$tmax
  z1 <- abex.list[[1]][[1]]$prec
  
  x2 <- abex.list[[2]][[1]]$tmin
  y2 <- abex.list[[2]][[1]]$tmax
  z2 <- abex.list[[2]][[1]]$prec
  
  x3 <- abex.list[[3]][[1]]$tmin
  y3 <- abex.list[[3]][[1]]$tmax
  z3 <- abex.list[[3]][[1]]$prec
  
  x4 <- abex.list[[4]][[1]]$tmin
  y4 <- abex.list[[4]][[1]]$tmax
  z4 <- abex.list[[4]][[1]]$prec
  
  
  # 3D ellipsoid
  rgl::par3d(windowRect = c(150, 150, 1000, 1000))  # c(xmin, ymin, xmax, ymax)
  tenm::plot_ellipsoid(x = x1, y=y1, z=z1, semiaxes= FALSE, xlab="tmin", ylab="tmax", zlab="prec", col="#4a72b0")
  # rgl::plot3d(x=x1,y=y1,z=z1, col="blue",add=T)
  # rgl::rgl.snapshot(paste0("./species/mig_ENM/seasonaly/", spp.folder,"/ENM_1_", spp.names[s],".png"))
  tenm::plot_ellipsoid(x = x2, y=y2, z=z2 ,semiaxes= FALSE, col="#6ea96e", mve=T, add=TRUE)
  # rgl::plot3d(x=x2,y=y2,z=z2, col="#6ea96e",add=T)
  # rgl::rgl.snapshot(paste0("./species/mig_ENM/seasonaly/", spp.folder,"/ENM_2_", spp.names[s],".png"))
  tenm::plot_ellipsoid(x = x3, y=y3, z=z3 ,semiaxes= FALSE, col="#da9500", mve=T, add=TRUE)
  # rgl::plot3d(x=x3,y=y3,z=z3, col="#da9500",add=T)
  # rgl::rgl.snapshot(paste0("./species/mig_ENM/seasonaly/", spp.folder,"/ENM_3_", spp.names[s],".png"))
  tenm::plot_ellipsoid(x = x4, y=y4, z=z4 ,semiaxes= FALSE, col="red", mve=T, add=TRUE)
  # rgl::plot3d(x=x4,y=y4,z=z4, col="red",add=T)
  # rgl::rgl.snapshot(paste0("./species/mig_ENM/seasonaly/", spp.folder,"/ENM_4_", spp.names[s],".png"))
  
  save.image(paste0("./species/mig_ENM/seasonaly/", spp.folder,"/ENM_s_", spp.names[s],".Rdata"))
  
  #delete all except tenm_mask and the background function from environment
  print(paste(spp.names[s], "done"))
  rm(list = setdiff(ls(), c("tenm_mask", "backg_CHcr", "spp.list", "spp.names", "tdf2swd2")))
}

# Projecting time-specific niche models -----------------------------------
# # This section does not run... projection in geography
# 
# Note: The tenm's predict function reads the CHELSAcruts layer with no data as infinite values. To run this section, we need to open them manually and assign the NaN values of -9999 - Due to an infinite values error. But the environmental space of suit_proj will be very large (because its maximum value is -9999). However, whether the data used in this section comes from Worldclim or CHELSA timeseries, this issue surely does not occur. 

# # For one time (e.g., 1901-01)
# env_layers_proj <- list.dirs(tempora_layers_dir, recursive = F)[1]
# 
# vars.stack<- rast(list.files(env_layers_proj, full.names=T))
# names(vars.stack)<-c("prec", "tmax", "tmin")
# NAflag(vars.stack)<- -9999
# var_list<-list()
# var_list[[1]] <- vars.stack
# 
# suit_proj <-
#   tenm::predict(mod_sel,
#                 model_variables = c("prec", "tmax", "tmin"),
#                 layers=var_list,
#                 # layers_path = env_layers_proj,
#                 # layers_ext = ".tif$",
#                 output = "suitability")
# 
# terra::plot(suit_proj, main="Prediction")

### Reviewing model tables
season.folder<-list.files("./species/mig_ENM/seasonaly/", full.names = T)
models.list<-list()

for (x in 1:length(season.folder)) {
models.list[[x]]<-read.table(list.files(season.folder[x], pattern = "ENM_mods_s_", full.names = T),
           sep="\t",
           dec=".",
           header=T)
}

do.call(rbind, models.list) |> 
   filter(om_rate_train>0.05) |> 
  arrange(om_rate_train)


# Elipsenm - Comparing niches ------------------

# Installing and loading packages
library(ellipsenm)

season.folder

x<-1
mig_enm.s<-
  list.files(season.folder[x], pattern=".txt$",full.names = T)
mig_enm.s<-mig_enm.s[-grep("ENM_mods_s.", mig_enm.s)]

elips.list<-list()
for (m in 1:length(mig_enm.s)) {
  occurrences <- read.table(mig_enm.s[m], 
                            sep="\t",
                            dec=".",
                            header=T) |>  select(sp_name, lon, lat, prec, tmax, tmin)
  ####
  
  bg_swd.list[[m]] |> 
    write.table("./borrar/Sel_call_back.txt", sep="\t", dec=".", row.names = F)
  
  # 1. Leer shapefile de puntos
  pts <- vect(bg_swd.list[[m]], geom=c("lon", "lat"))  # capa de puntos
  
  # 2. Crear raster vacío con extensión y resolución deseada
  # (ajusta res() al tamaño de celda que quieras)
  
  r <- rast("./borrar/clima/capas_mean/1/prec_ame_5k_1_mean.tif")
  
  #r <- rast(ext(pts), res = 0.01, crs = crs(pts))
  
  # 3. Rasterizar usando una columna de atributos (ej: "valor")
  prec.bg<-rasterize(pts, r, field = "prec")
  tmax.bg<-rasterize(pts, r, field = "tmax")
  tmin.bg<-rasterize(pts, r, field = "tmin")
  
  vars<- c(prec.bg, tmax.bg, tmin.bg)
  names()
  
  ####
  elips.list[[m]]<-overlap_object(data=occurrences, 
                                  species="sp_name",
                                  longitude="lon", 
                                  latitude = "lat", 
                                  method = "mve1",
                                  level = 99)
}

#Creamos los elipsoides

overlap_t<-
  ellipsoid_overlap(elips.list[[1]], 
                    elips.list[[2]], 
                    elips.list[[3]],
                    elips.list[[4]],
                    overlap_type = "all",
                    significance_test =F,
                    iterations=2)
overlap_t@full_overlap

plot_overlap(overlap_t,
            niches = c(1,2),
            data = F,
            background = F,
            proportion = 1,
            background_type = "full",
            niche_col = rainbow(4),
            data_col = rainbow(4))



# Shoener's D and Hellinger's I -------------------------------------------

library(ecospat)

# -------------------------------------------------------------------
# 1) PCA del espacio ambiental (PCA-env)
#    'glob' será el "fondo" ambiental. Si no tienes background,
#    una opción práctica es usar el pool de ambas especies.
# -------------------------------------------------------------------
x<-1

occ1<-read.table(mig_enm.s[1], sep="\t", dec=".", header=T) |> select(prec, tmax, tmin)
occ2<-read.table(mig_enm.s[2], sep="\t", dec=".", header=T) |> select(prec, tmax, tmin)
occ3<-read.table(mig_enm.s[3], sep="\t", dec=".", header=T) |> select(prec, tmax, tmin)
occ4<-read.table(mig_enm.s[4], sep="\t", dec=".", header=T) |> select(prec, tmax, tmin)


glob_env <- bind_rows(occ1, occ2) %>% as.data.frame()

# ecospat.pca.env calcula PCA usando 'glob' y proyecta sp1 y sp2

z1 <- ecospat.grid.clim.dyn(
  glob = abbg$env_bg |> select(prec, tmax),           # scores del fondo en PCA
  glob1 = abbg$env_bg |> select(prec, tmax),          # mismo fondo
  sp = abex.list[[1]]$temporal_df |> select(prec, tmax),          # scores de sp1 en PCA
  R = 100)                   # resolución

z2 <- ecospat.grid.clim.dyn(
  glob = abbg$env_bg |> select(prec, tmax),           # scores del fondo en PCA
  glob1 = abbg$env_bg |> select(prec, tmax),          # mismo fondo
  sp = abex.list[[2]]$temporal_df |> select(prec, tmax),     # scores de sp1 en PCA
  R = 100)                   # resolución


# -------------------------------------------------------------------
# 3) Solapamiento de nicho: Schoener's D
# -------------------------------------------------------------------
ov <- ecospat.niche.overlap(z1, z2, cor = TRUE)  # retorna D (y Hellinger I si aplica)
ov["D"]


ecospat.plot.niche.dyn(z1, z2, quant=0.25, interest=1,
                       title= "Niche overlap S1 vs S2",
                       name.axis1="Prec", name.axis2="Tmax")


