# Paquetes a verificar e instalar si es necesario
#paquetes <- c("shiny", "leaflet", "leaflet.extras", "sf", "tidyverse", "terra", "DT", "plotly", #"spThin")
#
## Comprobar e instalar los paquetes
#for (paquete in paquetes) {
#  if (!require(paquete, character.only = TRUE)) {
#    install.packages(paquete)
#  }
#  library(paquete, character.only = TRUE)
#}

# Cargar las Librerías necesarias
library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)
library(tidyverse)
library(terra)
library(DT)
library(plotly)
library(spThin)

# ------------------------- CARGAR Y FILTRAR REGISTROS ----------------------------

registros <- read_csv("files/Occ_limpio_2/Leptodactylidae/registros_crudos_leptodactylidae.csv")
#registros <- read_csv("files/especies/H_frenatus/observations-585349.csv") #iNaturalist

sp <- unique(registros$species)[14]

# Filtrar los registros de la especie deseada (GBIF)
occ <- registros %>%
  dplyr::filter(species == sp) %>% 
  dplyr::filter(basisOfRecord != "FOSSIL_SPECIMEN", # Elimina fósiles 
                institutionCode != "PBDB") %>% 
  dplyr::select(species, decimalLatitude, decimalLongitude) %>%
  drop_na() %>% # Elimina NAs
  dplyr::distinct() %>% # Elimina duplicados
  dplyr::rename(latitude = decimalLatitude, longitude = decimalLongitude)


## Juntar con registros de inaturalist
#reg2 <- read_csv("files/Occ_limpio_2/Atelopus fronterizo/observations-586211.csv")
#reg_intrlst <- reg2 %>% 
#  rename(species = scientific_name) %>% 
#  select(species, latitude, longitude) %>% 
#  na.omit() %>% 
#  distinct()
##occ <- rbind(occ,reg_intrlst)
#occ <- reg_intrlst
#sp <- unique(occ$species)


#occ <- read_csv("files/Occ_limpio_2/Atelopus spurrelli/Atelopus spurrelli.csv")

# Filtrar registros iNaturalist
#occ <- registros %>% 
#  dplyr::select(scientific_name, longitude, latitude) %>% 
#  rename(species = scientific_name) %>% 
#  drop_na() %>% # Elimina NAs
#  dplyr::distinct()


# Aplicar el thinning espacial para eliminar puntos dentro del radio especificado
set.seed(123)  # Para replicabilidad (opcional)
thinned <- spThin::thin(
  loc.data = occ %>% select(species, latitude, longitude),
  lat.col = "latitude",
  long.col = "longitude",
  spec.col = "species",
  thin.par = 1,  # Distancia mínima en km (ajústalo según necesites)
  reps = 1,       
  locs.thinned.list.return = TRUE,
  write.files = FALSE
)

occ_filt <- thinned[[1]] %>% 
  tibble::as_tibble() %>%
  dplyr::rename(latitude = Latitude, longitude = Longitude) %>% 
  dplyr::mutate(species = unique(occ$species), .before = longitude,
                IDU = row_number())


# Convertir los puntos afinados de nuevo a un objeto sf
thin_points_sf <- sf::st_as_sf(occ_filt, coords = c("longitude", "latitude"), crs = 4326)



#Crear carpetas para los registros y la M
ruta_sp <- paste0("files/Occ_limpio_2/Leptodactylidae/", sp)
if (!dir.exists(ruta_sp)) {
  dir.create(ruta_sp, recursive = TRUE)
}

ruta_M <- paste0("files/Occ_limpio_2/Leptodactylidae/", sp, "/M")
if (!dir.exists(ruta_M)) {
  dir.create(ruta_M, recursive = TRUE)
}


# ------------------------- INTERFAZ 1: MAPA ----------------------------

ui1 <- fluidPage(
  titlePanel("Seleccionar puntos en mapa"),
  leafletOutput("mapa", height = 600),
  fluidRow(
    column(6, actionButton("reset", "Restablecer selección")),
    column(6, downloadButton("download", "Descargar selección"))
  ),
  DTOutput("selected_table")
)

server1 <- function(input, output, session) {
  selected_data <- reactiveVal(data.frame())
  
  output$mapa <- renderLeaflet({
    leaflet(thin_points_sf) %>%
      addProviderTiles(providers$Esri.WorldShadedRelief, group = "Relieve") %>%
      addProviderTiles(providers$OpenTopoMap, group = "Topográfico") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Satelital") %>%
      addProviderTiles(providers$OpenStreetMap, group = "OSM") %>%
      addCircleMarkers(radius = 5, color = "red", 
                       popup = ~paste("Especie:", species, "<br>IDU:", IDU),
                       label = ~paste(species, "- IDU:", IDU)) %>%
      addDrawToolbar(
        targetGroup = "drawn_polygons",
        polygonOptions = drawPolygonOptions(
          shapeOptions = drawShapeOptions(fillOpacity = 0.4, color = "red")
        ),
        editOptions = editToolbarOptions(edit = TRUE, remove = TRUE)
      ) %>%
      addLayersControl(
        baseGroups = c("Relieve", "Topográfico", "Satelital", "OSM"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  observeEvent(input$mapa_draw_new_feature, {
    feature <- input$mapa_draw_new_feature
    if (!is.null(feature)) {
      coords <- feature$geometry$coordinates[[1]] %>% 
        map(~ c(.x[[1]], .x[[2]])) %>% 
        do.call(rbind, .) %>% 
        as.data.frame()
      names(coords) <- c("longitude", "latitude")
      
      poly_sf <- st_as_sf(coords, coords = c("longitude", "latitude"), crs = 4326) %>%
        summarise(geometry = st_combine(geometry)) %>%
        st_cast("POLYGON")
      
      selected <- thin_points_sf[st_intersects(thin_points_sf, poly_sf, sparse = FALSE), ]
      selected_df <- selected %>%
        mutate(longitude = st_coordinates(.)[, 1],
               latitude = st_coordinates(.)[, 2]) %>%
        st_drop_geometry()
      
      selected_data(selected_df)
      
      leafletProxy("mapa") %>%
        clearShapes() %>%
        addPolygons(data = poly_sf, fillOpacity = 0.4, color = "red") %>%
        clearMarkers() %>%
        addCircleMarkers(
          data = selected,
          radius = 5,
          color = "red",
          popup = ~paste("Especie:", species, "<br>IDU:", IDU),
          label = ~paste(species, "- IDU:", IDU)
        )
    }
  })
  
  output$selected_table <- renderDT({
    datatable(selected_data(), options = list(pageLength = 10))
  })
  
  output$download <- downloadHandler(
    filename = function() { "datos_seleccionados.csv" },
    content = function(file) {
      write_csv(selected_data(), file)
    }
  )
  
  observeEvent(input$reset, {
    selected_data(data.frame())
    leafletProxy("mapa") %>%
      clearShapes() %>%
      clearMarkers() %>%
      addCircleMarkers(
        data = thin_points_sf,
        radius = 5,
        color = "red",
        popup = ~paste("Especie:", species, "<br>IDU:", IDU),
        label = ~paste(species, "- IDU:", IDU)
      )
  })
}

shinyApp(ui1, server1)

# ------------------- INTERFAZ 2: SCATTERPLOT AMBIENTAL ----------------------

registros <- read_csv("files/Occ_limpio_2/Leptodactylidae/registros_crudos_leptodactylidae.csv")
#registros <- read_csv("files/especies/H_frenatus/observations-585349.csv") #iNaturalist

sp <- unique(registros$species)[14]

# Leer datos filtrados
sp_occ <- read_csv(paste0("files/Occ_limpio_2/Leptodactylidae/", sp, "/occ_spat.csv"))
vars <- terra::rast(list.files("C:/Users/ASUS/Desktop/Bios_CHELSA/", pattern = "\\.tif$", full.names = TRUE))
elev <- rast("D:/Variables_WC_2.1/wc2.1_30s_elev/wc2.1_30s_elev.tif")
#vars <- c(vars, elev) # si provienen de la misma fuente (mismo extent)
sp_occ_extract1 <- terra::extract(vars, sp_occ[, c("longitude", "latitude")])
sp_occ_extract2 <- terra::extract(elev, sp_occ[, c("longitude", "latitude")])
df <- bind_cols(sp_occ, sp_occ_extract1[-1], sp_occ_extract2[-1]) %>% 
  drop_na()

#--- Filtrado ambiental v 2.0 con qqplot
# Estoy intentando hacer una función que yo le de el numero de la especie, el numero del mes y la variable y me haga el shiny de cada especie por mes para poder hacer el filtro ambiental

spp.list<-list.files("./species/mig_shapes/ok/", pattern=".shp", full.names = T)
spp.names<-list.files("./species/mig_shapes/ok/", pattern=".shp") |> str_remove(".shp")

env.filt <- function(spp, n.month, variable){

spp <- 1 #spp # Número de la especie
n.month<- 1 #month # Mes (número)
variable <- "CHELSA_B12" #variable # CHELSA_B1 o CHELSA_B12

df <- vect(spp.list[spp]) |> 
  as.data.frame() |> 
  select(month, CHELSA_B1, CHELSA_B12, gbifID) |> 
  filter(month==n.month)




# Definir las variables que quieres visualizar

if(variable=="CHELSA_B12"){
  
df.12<-
  df |> 
  dplyr::arrange(CHELSA_B12)

var_y <- "CHELSA_B12" # Variable real

# Para ahcer el QQ-plot hay que construir un eje x teórico que permita contrastar los datos reales vs los esperados
# 1. Ordenar los valores de y
n <- 
  df.12$CHELSA_B12 |> 
  sort() |> 
  length()

p <- (1:n - 0.5) / n

df.12$B12.teor <- qnorm(p)

var_x <- "B12.teor" 

#plot(df.12$B12.teor, df.12$CHELSA_B12) # probar que si este bien

# var_x <- "CHELSA_B1" # Precipitación anual (Bio12)

# Ajustes panel interactivo
ui2 <- fluidPage(
  titlePanel("Selección en scatterplot ambiental"),
  plotlyOutput("scatter", height = 600),
  br(),
  downloadButton("download_selected", "Descargar selección (.csv)"),
  br(), br(),
  DTOutput("selected_table")
)

# Ajustar la información a visualizar
server2 <- function(input, output, session) {
  output$scatter <- renderPlotly({
    plot_ly(
      data = df.12,
      x = as.formula(paste0("~", var_x)),  # Usamos la variable definida
      y = as.formula(paste0("~", var_y)),  # Usamos la variable definida
      type = "scatter",
      mode = "markers",
      text = ~paste("Especie:", spp.names[x], "<br>ID:", gbifID),
      marker = list(color = 'blue', size = 6)
    ) %>%
      layout(dragmode = "lasso")
  })
  
  selected_points <- reactive({
    eventdata <- event_data("plotly_selected")
    if (is.null(eventdata)) return(NULL)
    df[eventdata$pointNumber + 1, ]
  })
  
  output$selected_table <- renderDT({
    req(selected_points())
    datatable(selected_points(), options = list(pageLength = 10))
  })
  
  output$download_selected <- downloadHandler(
    filename = function() { paste0("seleccion_scatterplot_", Sys.Date(), ".csv") },
    content = function(file) {
      write_csv(selected_points(), file)
    }
  )
}

shinyApp(ui2, server2)





} else( df.1<-
          df |> 
          dplyr::arrange(CHELSA_B1))


}
