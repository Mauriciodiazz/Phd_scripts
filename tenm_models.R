### tempo-specific niche models -  tenm package

# if (!require('devtools')) install.packages('devtools')
# devtools::install_github("luismurao/tenm")
# # If you want to build vignette, install pandoc before and then
# devtools::install_github('luismurao/tenm',build_vignettes=TRUE)

# Libraries
library(tenm)
library(terra)
library(tidyverse)
library(sf)

# Open species file
tenm_mask <- terra::rast(list.files("D:/CHELSA_ym_5k/1901-01/", full.names=T, pattern=".tif$"))
names(tenm_mask)<-names(tenm_mask) |> 
  str_remove("CHELSAcruts_") |> 
  str_remove("_V.1.0")

spp<-vect("./species/mig_shapes/Selected/Sel_calli.shp")

spp2<-terra::extract(tenm_mask, spp, ID=F, bind=T)


# Transformin it as data.frame
spp.df<-
  spp2 |>
  st_as_sf() |> 
  st_coordinates() |> 
  cbind(st_as_sf(spp2)) |>
  drop_na("prec_1_1901", "tmax_1_1901", "tmin_1_1901") |>
  rename("lon" = X, "lat" = Y) |>
  as.data.frame() |>
  select(species, lon, lat, year, month) |>
  filter(year>=1901 & year<=2016) |> 
  unite("date", year, month, sep="-", remove=F) |> 
  mutate(date = str_replace(date, "-(\\d)$", "-0\\1")) # date format YYYY-MM

#write.csv(spp.df, "./outputs/borrar/Sel_calli.csv")
head(spp.df)
range(spp.df$year) # Inspecting temporal range (1901-2016)

# We indicate the path where our time-specific modeling layers are located (variables).
tempora_layers_dir <- "D:/CHELSA_ym_5k/"

# Exploring the directory structure
# list.dirs(tempora_layers_dir, recursive = F)

# sp_temporal_data will allow us to work with time-specific data
abt <- tenm::sp_temporal_data(occs = spp.df,
                              longitude = "lon",
                              latitude = "lat",
                              sp_date_var = "date",
                              occ_date_format="ym",
                              layers_date_format= "ym",
                              layers_by_date_dir = tempora_layers_dir,
                              layers_ext="*.tif$")

# Time-specific spatial data thinning -------------------------------------

# Clean duplicates using a raster mask
tenm_mask

abtc <- tenm::clean_dup_by_date(this_species = abt,
                                by_mask = TRUE,
                                threshold = terra::res(tenm_mask)[1],
                                raster_mask = tenm_mask[1],
                                n_ngbs = 0)

# Check number of records
head(tidyr::as_tibble(abtc$temporal_df))
nrow(tidyr::as_tibble(abtc$temporal_df))

nrow(abtc$temporal_df)
nrow(abt$temporal_df)

# Time-specific environmental data extraction -----------------------------

future::plan("multisession",workers=2) # Allow that ex_by_date could be run in parallel
abex <- tenm::ex_by_date(this_species = abtc,
                         train_prop=0.7) # train=0.7 and test=0.3
future::plan("sequential")

# Time-specific background generation -------------------------------------

future::plan("multisession",workers=2)
abbg <- tenm::bg_by_date(this_species = abex,
                         buffer_ngbs=10, 
                         n_bg=10000,
                         buffer_distance=5000)
future::plan("sequential")

# Exporting time-specific information as Samples With Data (SWD) format ---------

# SWD table for occurrence records
occ_swd <- tdf2swd(this_species = abex, 
                   sp_name = "Sel_calli")

# SWD table for background data
bg_swd <- tdf2swd(this_species = abbg) #missing year values
head(tidyr::as_tibble(occ_swd))
head(tidyr::as_tibble(bg_swd))


# Time-specific model calibration and selection ---------------------------
vars2fit<- c("prec", "tmin", "tmax")

mod_sel <- tenm::tenm_selection(this_species = abbg,
                                omr_criteria =0.1,
                                ellipsoid_level=0.999, #0.975 by default
                                vars2fit = vars2fit,
                                nvars_to_fit=3,
                                proc = T,
                                RandomPercent = 50,
                                NoOfIteration=1000,
                                parallel=TRUE,
                                n_cores=4)

names(mod_sel)
mod_sel$mods_table #resultados de la evalaución

# Centroid ----------------------------------------------------------------
mod <- tenm::cov_center(data = abex$env_data,
                        mve = TRUE,
                        level = 0.975,
                        vars = c("prec","tmax","tmin"))

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

# Plotting ellipsoid

x <- abex$temporal_df$tmin
y <- abex$temporal_df$tmax
z <- abex$temporal_df$prec

# 3D ellipsoid
tenm::plot_ellipsoid(x = x, y=y, z=z ,semiaxes= FALSE, xlab="tmin", ylab="tmax", zlab="prec")
# # Adding a plot
# tenm::plot_ellipsoid(x = x+100, y=y, z=z ,semiaxes= FALSE,add=TRUE)


