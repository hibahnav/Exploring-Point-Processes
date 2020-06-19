library(shiny)
library(tidyverse)
library(vroom)
library(tidytuesdayR)
library(dplyr)
library(ggplot2)


RAWDAT <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-02-18/food_consumption.csv')

##explore a bit
RAWDAT %>% count(country) 


selected <- RAWDAT %>% filter(country %in% c('Argentina',"Canada","Pakistan"))


ggplot(data=RAWDAT, aes(y=consumption,x=factor(food_category),fill=co2_emmission))+
  theme_minimal() + 
  geom_boxplot() 



ui <- fluidPage(
  theme = shinythemes::shinytheme("sandstone"),
  titlePanel("Food Consumption and C02 Emmissions"),
  sidebarLayout(
    sidebarPanel(width=2),
    mainPanel(
      tabsetPanel( 
        tabPanel("All Countries","Displaying data for all countries",
                 br(),
                 br(),
                 dataTableOutput("alldat"),
                 plotOutput("plot",width="100%")
        ),
        tabPanel("Each Country","Please select a country from the drop down menu below:",
                 br(),
                 br(),
                 selectInput("country","Country",setNames(RAWDAT$country,RAWDAT$country)),
                 tableOutput("dat"),
                 plotOutput("countryplot")
        )
        
      )))
)

server <- function(input, output, session){
  selected <- reactive(RAWDAT %>% filter(country==input$country))
  
  output$alldat <- renderDataTable(RAWDAT,options=list(pageLength=5))
  
  output$dat <- renderTable(selected(),width="100%")
  
  output$plot <- renderPlot({
    ggplot(data=RAWDAT, aes(y=consumption,x=factor(food_category),fill=co2_emmission))+
      theme_minimal() + 
      geom_boxplot() +
      labs(x="food_category")},res=96)
  
  output$countryplot <- renderPlot({
    plot(selected()$consumption,selected()$co2_emmission)})
  #selected() %>%
  #ggplot(aes(y=consumption,x=factor(food_category),fill=co2_emmission))+
  #theme_minimal() + 
  #geom_point()
  #labs(x="food_category")},res=96)
}


shinyApp(ui = ui, server = server)



