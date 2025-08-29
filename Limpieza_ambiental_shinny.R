# # Paquetes a verificar e instalar si es necesario
# paquetes <- c("shiny", "leaflet", "leaflet.extras", "sf", "tidyverse", "terra", "DT", "plotly", "spThin")
# 
# # Comprobar e instalar los paquetes
# for (paquete in paquetes) {
#  if (!require(paquete, character.only = TRUE)) {
#    install.packages(paquete)
#  }
#  library(paquete, character.only = TRUE)
# }

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

# ------------------- INTERFAZ: SCATTERPLOT AMBIENTAL ----------------------
# Este script genera una función en donde se le da el número de la especie, el número del mes y el nombre de la variable para generar un Q-Q plot ambiental de cada variable

spp.list<-list.files("./species/mig_shapes/ok/", pattern=".shp$", full.names = T)
spp.names<-list.files("./species/mig_shapes/ok/", pattern=".shp$") |> str_remove(".shp")

# spp: numeric
# n-month: numeric (1-12)
# variable: character ("CHELSA_B1"  "CHELSA_B12")

env.filt <- function(spp, n.month, variable){

spp <- spp 
n.month<- n.month 
variable <- variable 

# 2. Definir las variables que quieres visualizar

if(variable=="CHELSA_B12"){
  
df.12<-
  vect(spp.list[spp]) |> 
  as.data.frame() |> 
  select(month, CHELSA_B1, CHELSA_B12, gbifID) |> 
  filter(month==n.month) |> 
  dplyr::arrange(CHELSA_B12)

var_y <- "CHELSA_B12" # Variable real

# 3. construir variable teórica del QQ-plot
# 3.1 Ordenar los valores del eje y
n <- 
  df.12$CHELSA_B12 |> 
  length()
#3.2 Cuantiles
p <- (1:n - 0.5) / n

#3.3 Valores de probabilidad según una distribución normal
df.12$B12.teor <- qnorm(p)

#3.4 Definición de la variable teórica
var_x <- "B12.teor"

#plot(df.12$B12.teor, df.12$CHELSA_B12) # probar que si este bien

#5. Ajustes panel interactivo
ui_B12 <- fluidPage(
  titlePanel(paste("Selección en scatterplot ambiental", spp.names[spp], "(",spp,")", "month=",n.month, "var=",variable)),
  plotlyOutput("scatter", height = 600),
  br(),
  downloadButton("download_selected", "Descargar selección (.csv)"),
  br(), br(),
  DTOutput("selected_table")
)

# Ajustar la información a visualizar
server_B12 <- function(input, output, session) {
  output$scatter <- renderPlotly({
    plot_ly(
      data = df.12,
      x = as.formula(paste0("~", var_x)),  # Usamos la variable definida
      y = as.formula(paste0("~", var_y)),  # Usamos la variable definida
      type = "scatter",
      mode = "markers",
      text = ~paste("Especie:", spp.names[spp], "(", spp ,")", "<br>ID:", gbifID),
      marker = list(color = 'blue', size = 6)
    ) %>%
      layout(dragmode = "select")
  })
  
  selected_points <- reactive({
    eventdata <- event_data("plotly_selected")
    if (is.null(eventdata)) return(NULL)
    df.12[eventdata$pointNumber + 1, ]
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

return(shinyApp(ui_B12, server_B12))

} else {
  df.1 <-
    vect(spp.list[spp]) |> 
    as.data.frame() |> 
    select(month, CHELSA_B1, CHELSA_B12, gbifID) |> 
    filter(month==n.month) |>
    dplyr::arrange(CHELSA_B1)
  
  var_y <- "CHELSA_B1" # Variable real
  
  # 3. construir variable teórica del QQ-plot
  # 3.1 Ordenar los valores del eje y
  n <- 
    df.1 |> 
    nrow()
  
  #3.2 Cuantiles
  p <- (1:n - 0.5) / n
  
  #3.3 Valores de probabilidad según una distribución normal
  df.1$B1.teor <- qnorm(p)
  
  #3.4 Definición de la variable teórica
  var_x <- "B1.teor"
  
  #plot(df.12$B12.teor, df.12$CHELSA_B12) # probar que si este bien
  
  #5. Ajustes panel interactivo
  ui_B1 <- fluidPage(
    titlePanel(paste("Selección en scatterplot ambiental", spp.names[spp], "(",spp,")", "month=",n.month, "var=",variable)),
    plotlyOutput("scatter", height = 600),
    br(),
    downloadButton("download_selected", "Descargar selección (.csv)"),
    br(), br(),
    DTOutput("selected_table")
  )
  
  # Ajustar la información a visualizar
  server_B1 <- function(input, output, session) {
    output$scatter <- renderPlotly({
      plot_ly(
        data = df.1,
        x = as.formula(paste0("~", var_x)),  # Usamos la variable definida
        y = as.formula(paste0("~", var_y)),  # Usamos la variable definida
        type = "scatter",
        mode = "markers",
        text = ~paste("Especie:", spp.names[spp], "(", spp ,")", "<br>ID:", gbifID),
        marker = list(color = 'blue', size = 6)
      ) %>%
        layout(dragmode = "select")
    })
    
    selected_points <- reactive({
      eventdata <- event_data("plotly_selected")
      if (is.null(eventdata)) return(NULL)
      df.1[eventdata$pointNumber + 1, ]
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
  
  return(shinyApp(ui_B1, server_B1))
}
}

# Este código genera un gráfico que contiene los qqplots mensuales por especie (x)
x <- 91
# x11()
vect(spp.list[x]) |> 
  as.data.frame() |> 
  select(month, CHELSA_B1, CHELSA_B12, gbifID) |> 
  mutate(month = case_when(
    month == 1 ~ "Enero",
    month == 2 ~ "Febrero",
    month == 3 ~ "Marzo",
    month == 4 ~ "Abril",
    month == 5 ~ "Mayo",
    month == 6 ~ "Junio",
    month == 7 ~ "Julio",
    month == 8 ~ "Agosto",
    month == 9 ~ "Septiembre",
    month == 10 ~ "Octubre",
    month == 11 ~ "Noviembre",
    month == 12 ~ "Diciembre")) |> 
  mutate(month= factor(month, 
                       levels = c("Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio", "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre"))) |> 
  pivot_longer(!c(month,gbifID), names_to = "variable", values_to = "value") |> 
  #  sample_n(size=100) |> 
  ggplot(aes(sample=value)) +
  stat_qq_line() +
  stat_qq(col="#2b75a6") +
  facet_wrap(~month+variable, scales="free") +
  labs(x="", y="", title=paste(spp.names[x],"(", x, ")")) +
  theme_classic()

env.filt(spp=x, n.month = 7, variable = "CHELSA_B12")


