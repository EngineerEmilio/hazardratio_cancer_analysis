# Shiny
library(shiny)
library(bslib)
library(shinydashboard)
# Modelaje
library(modeldata)
# Widgets
library(plotly)
library(tidyverse)

#Manipulacion
library('dplyr')
library(ggplot2)
library(readxl)
library(survival)
library(shinythemes)
# install.packages("gtsummary")
# install.packages("shinythemes")
library(gtsummary)

# --------  RUTAS Emilio ----------
setwd("C:/Users/1097513297/Documents/BEDU/R_works/Proyecto_Cancer")
getwd()
archivo_excel <- "Genes para superviencia.xlsx"

df.list <- lapply(excel_sheets(archivo_excel), function(sheet) {
  read_excel(archivo_excel, sheet = sheet)
})

# --------  RUTAS Misael ----------


# --------  RUTAS Kalaumari ----------

# ----------------------- Parte UI

ui <- dashboardPage(
  dashboardHeader(title = "Proyecto Estadistico de Survival Analysis "),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Histogram", tabName = "distPlot", icon = icon("bar-chart")),
      menuItem("Kaplan-Meier Curves", tabName = "kmPlot", icon = icon("line-chart")),
      menuItem("Cox Proportional Hazard Model", tabName = "coxPlot", icon = icon("area-chart")),
      menuItem("Time Series Plot", tabName = "timeSeriesPlot", icon = icon("line-chart")),
      selectInput("dataframe", "Select Dataframe:",
                  choices = c("CDHR1" = "df.list[[2]]",
                              "DSC3" = "df.list[[3]]",
                              "GJB5" = "df.list[[4]]",
                              "S1PR5 " = "df.list[[5]]",
                              "GPC1" = "df.list[[6]]"),
                  selected = "df.list[[2]]")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "distPlot", plotOutput("distPlot")),
      tabItem(tabName = "kmPlot", plotOutput("kmPlot")),
      tabItem(tabName = "coxPlot", tableOutput("coxPlot")),
      tabItem(tabName = "timeSeriesPlot", plotOutput("timeSeriesPlot"))
    ),
    mainPanel(textOutput("title")),
    
    fluidRow(
      column(12, box(title = "Créditos",
               "Dashboard de Análisis de Supervivencia es un proyecto para certificación de Bedu",
               background = "light-blue",  
               solidHeader = FALSE)
             )
           )
    
  )
)

#------------------  Corresponde al server

server <- function(input, output, session) {

  # Valores tipo React
  selected_df <- reactive({
    eval(parse(text = input$dataframe))
  })
  
  
  # Histograma
  output$distPlot <- renderPlot({
    
    color <- 
      if(input$dataframe == "df.list[[2]]"){
        "yellow"
      } else{
      if(input$dataframe == "df.list[[3]]"){
        "green"
      }else{
      if (input$dataframe == "df.list[[4]]"){
        "blue"
      } else {
      if (input$dataframe == "df.list[[5]]"){
          "orange"
      } else {
        if(input$dataframe == "df.list[[6]]"){
            "red"
              }
            }
          }
        }
      }
    
    death.DS <- selected_df()[selected_df()$Event == 1, ]
    x <- death.DS$`Time (months)`
    hist_data <- hist(x, breaks = 60, plot = FALSE)
    

    plot(hist_data, main = "Histograma de muertes", xlab = "Tiempo (meses)",
         ylab = "Frecuencia muertes", col = color, border = "black")
  })
  
  # Kaplan-Meier Curves 
  
  output$kmPlot <- renderPlot({
    hi.exp <- selected_df()[selected_df()$`Expression (1=high)` < 1, ]
    low.exp <- selected_df()[selected_df()$`Expression (1=high)` > 0, ]
    hi.survdata <- Surv(hi.exp$`Time (months)`, hi.exp$Event)
    low.survdata <- Surv(low.exp$`Time (months)`, low.exp$Event)
    K_M.hi <- survfit(hi.survdata ~ 1)
    K_M.low <- survfit(low.survdata ~ 1)
    plot(K_M.hi, main = "Curva de Kaplan-Meier", xlab = "Tiempo (meses)", ylab = "Probabilidad de supervivencia", col = 2, lwd = 2)
    lines(K_M.low, lwd = 2, col = "blue")
    legend('bottomright', inset = 0.05, c("Expresión Alta", "Expresión Baja"), lty = 1, col = c("red", "blue"), title = "Gráficas")
  })
  
  # Cox Proportional Hazard Model 
  
  output$coxPlot <- renderTable({
    coxph(Surv(selected_df()$`Time (months)` , selected_df()$Event) ~ selected_df()$`Expression (1=high)`, 
          data = selected_df()) %>% tbl_regression(exp = TRUE) 
  })
  
  # Series de tiempo
  
  output$timeSeriesPlot <- renderPlot({
    death.DS <- selected_df()[selected_df()$Event == 1, ]
    x <- death.DS$`Time (months)`
    hist_data <- hist(x, breaks = 60, plot = FALSE)
    tiempo <- hist_data$breaks[-1]
    frecuencia <- hist_data$counts
    ggplot(data.frame(tiempo = tiempo, frecuencia = frecuencia), aes(x = tiempo, y = frecuencia)) +
      geom_line() +
      labs(title = "Serie de Tiempo de Muertes", x = "Tiempo (meses)", y = "Número de muertes")
  })
  
  output$title <- renderText({
    color <- 
      if(input$dataframe == "df.list[[2]]"){
        paste("CDHR1: Codifica para la proteína cadherina-1, relacionada con la adhesión celular y la diferenciación tisular.")
      } else{
        if(input$dataframe == "df.list[[3]]"){
          paste("DSC3: Codifica para la desmocolina-3, una proteína de unión celular involucrada en la cohesión y mantenimiento de la integridad de la epidermis.")
        }else{
          if (input$dataframe == "df.list[[4]]"){
            paste("GJB5: Codifica para la proteína conexina 31, que forma parte de las uniones comunicantes y está involucrada en la comunicación intercelular.")
          } else {
            if (input$dataframe == "df.list[[5]]"){
              paste("S1PR5: Codifica para el receptor 5 de esfingosina-1-fosfato, que participa en la regulación de procesos celulares como la migración y la supervivencia.")
            } else {
              if(input$dataframe == "df.list[[6]]"){
                paste("GPC1: Codifica para la glicoproteína de células germinales 1, que se ha asociado con la progresión del cáncer y la angiogénesis.")
              }
            }
          }
        }
      }
    
    
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)