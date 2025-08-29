# Amplitud de nicho, latitud y elevacion
# Felipe Toro, Daniel Valencia, Fabricio Villalobos, Mauricio Diaz, Juliana Herrera, Gabriel Moulatlet ....


library(tidyverse)
library(sf)
library(rgbif)

sf_use_s2(use_s2 = F)

# Obtener el listado de especies a partir de los poligonos IUCN -----------
mams<- read_sf("Rodentia/MDD_Rodentia.gpkg")
names(mams)

# Separo rondentia de todo el set de datos. En caso de que se use la iucn como fuente de datos, seria necesario hacer un filtro por la columna order_
rodentia<- mams

# grupo<- filter(poligs, order_=="Orden o grupo a filtrar")

# Ahora filtro por los poligonos que estan en America y saco el listado de especie en un vector para luego buscar los registros de presencia
America<- read_sf("América.shp")
rods.america<- st_intersection(rodentia, America)

head(rods.america)
spp_multipolygon <- rods.america %>%
  group_by(sciname) %>%
  summarize(geometry = st_union(geom)) 

spp_multipolygon<- arrange(spp_multipolygon)
head(spp_multipolygon)
spp.rodentia<- unique(spp_multipolygon$sciname) # 1119 especies


# Ahora puedo buscar los registros de las especies del grupo, esto puede tardar algunos dias dependiendo del numero de registros y especies. Se dejan algunos campos para tener una identidad de cada registro descargado de gbif que pueda servir de alguna manera
# Instala y carga los paquetes necesarios
install.packages("rgbif")
install.packages("httr")
install.packages("retry")
library(rgbif)
library(httr)
library(retry)


# Define la función d.gbif
d.gbif <- function(spp, output_file, max_retries = 3, 
                   timeout_seconds = 60, retry_delay = 5) {
  rgbif_data <- list()  # Lista para almacenar los datos
  
  # Configura el tiempo de espera global para las solicitudes
  set_config(timeout(timeout_seconds))
  
  for (i in 1:length(spp)) {
    species_name <- spp[i]
    retry_count <- 0
    success <- FALSE
    
    while (retry_count < max_retries & !success) {
      retry_count <- retry_count + 1
      
      tryCatch({
        # Realiza la solicitud
        occs1 <- occ_search(
          scientificName = species_name, 
          hasCoordinate = TRUE, 
          limit = 30000, 
          fields = c("institutionCode", "catalogNumber", "collectionCode",
                     "country", "stateProvince", "decimalLatitude", 
                     "decimalLongitude", "basisOfRecord", "species")
        )
        
        # Si la solicitud fue exitosa, agrega los datos a la lista y marca como exitoso
        rgbif_data[[species_name]] <- occs1
        success <- TRUE
        
        # Guarda los datos en un archivo después de cada descarga exitosa
        save(rgbif_data, file = output_file)
        
      }, error = function(e) {
        message(paste("Error al descargar datos para", species_name, "en intento", retry_count, ":", e$message))
        if (retry_count >= max_retries) {
          message(paste("No se pudo descargar datos para", species_name, "después de", max_retries, "intentos."))
        }
      })
      
      # Espera un tiempo antes de reintentar
      if (!success) {
        Sys.sleep(retry_delay)  # Espera antes del siguiente intento
      }
    }
  }
  
  # Extrae solo los dataframes de cada elemento en rgbif_data
  data_frames_list <- lapply(rgbif_data, function(x) x$data)
  
  return(data_frames_list)
}

rods1<- d.gbif(spp=spp.rodentia, output_file = "Rodentia_raw_rgbif.RData") # ajustar el nombre del archivo según el grupo
save(list=c("rods1"), file="Rodentia_raw_dfs.RData")

# Convertir cada dataframe en la lista a un objeto sf y hacer la intersección con su polígono respectivo

# Función para procesar cada especie
process_species <- function(species_name) {
  
  coords_columns <- c("decimalLongitude", "decimalLatitude") 
  
  # Obtener los datos para la especie
  species_data <- rods1[[species_name]]
  
  # Mensajes de depuración
  cat("Processing species:", species_name, "\n")
  cat("Type of species_data:", class(species_data), "\n")
  
  # Verificar si el objeto es un tibble y tiene registros
  if (inherits(species_data, "tbl") && nrow(species_data) > 0) {
    
    # Convertir el tibble a objeto sf
    species_sf <- tryCatch(
      st_as_sf(species_data, coords = coords_columns, crs = 4326),
      error = function(e) {
        cat("Error converting to sf:", e$message, "\n")
        NULL
      }
    )
    
    if (is.null(species_sf)) {
      return(NULL)
    }
    
    # Filtrar el polígono de la especie correspondiente
    species_polygon <- spp_multipolygon %>%
      filter(sciname == species_name)
    
    if (nrow(species_polygon) == 0) {
      return(NULL)
    }
    
    # Intersección del sf de puntos con el polígono
    intersection <- st_intersection(species_sf, species_polygon)
    
    # Devolver la intersección solo si tiene registros
    if (nrow(intersection) > 0) {
      return(intersection)
    }
  }
  
  return(NULL)
}

# Aplicar la función a cada especie
intersected_list <- lapply(names(rods1), process_species)
names(intersected_list)<- spp.rodentia

# Filtrar elementos NULL de la lista resultante y mantener solo aquellos que no son NULL
valid_intersections <- !sapply(intersected_list, is.null)

# Aplicar el filtro para eliminar los elementos con 0 registros después de la intersección
intersected_list2 <- intersected_list[valid_intersections]

# Eliminar los elementos que queden con 1 o ningun registro después de la intersección
intersected_list3 <- intersected_list2[sapply(intersected_list2, nrow) > 1]

# Por ultimo hago un subset de los polígonos de distribución e las especies que quedaron al final

spp_polyg<- spp_multipolygon %>% filter(sciname %in% names(intersected_list3))

# Guardo el rdata (guardo filtrado original y final, así como los polígonos)
save(list = c("intersected_list", "intersected_list3", "spp_polyg"), file = "Rodentia_filtoccs_ranges.RData")
