#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
source('VS-lite_functions.R')
T <- read.csv("T.csv")
P <- read.csv("P.csv")

# Convert to a 12 x nyear data_frame
T <- pivot_wider(T, names_from = year, values_from = tmean)
P <- pivot_wider(P, names_from = year, values_from = ppt)

# Strip out the first row, which is a the month index
T <- as.data.frame(T[,-1])
P <- as.data.frame(P[ ,-1])
syear <- 1895
eyear <- 2019

phi <- 35.6

tout <- VSLite(syear, eyear, phi, T = T, P = P, T1 = 10, T2 = 22, M1 = .01, M2 = 1,
              I_0=1,I_f=12)

trw <- as.data.frame(read.csv("trw.csv"))



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("VS-Lite"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("T2",
                        "Tmax:",
                        min = 20,
                        max = 35,
                        value = 20, 
                        post=' C'),
            sliderInput("M2", 
                        "M2:",
                        min=.1, max=1,
                        value=.8),
            sliderInput("Mmax", 
                        "Mmax (holding capacity of soil):",
                        min=.1, max=1,
                        value=.76),
            sliderInput("rootd", 
                        "rootdepth",
                        min=500, max=1000,
                        value=1000)
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           textOutput("cor")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    out <- reactive({
        VSLite(syear, eyear, phi, T = T, P = P, T1 = 10, T2 = input$T2, M1 = .01, M2 = input$M2,
               I_0=1,I_f=12, Mmax=input$Mmax, rootd=input$rootd)
    })
    
    observe({print(out()$potEV)})

    output$distPlot <- renderPlot({
        ggplot(data=NULL, aes(x=1:12, y=rowMeans(out()$gT))) +
            geom_line() +
            geom_line(aes(x=1:12, y=rowMeans(out()$gM)), color='blue') + 
            geom_line(aes(x=1:12, y=rowMeans(out()$Gr)), color='red')
    })
    output$cor <- renderText({paste("correlation with trw: ", cor(as.vector(out()$trw), trw$trw_norm))})
}

# Run the application 
shinyApp(ui = ui, server = server)
