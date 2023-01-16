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

## Functions

# Define UI for application that draws a histogram
# Define UI for slider demo app ----
ui <- fluidPage(
    
    # App title ----
    titlePanel("Optimal and ESS LMA - Weng et al. 2017 model"),
    
    #plotOutput("essPlot"),
    
    # Output: Table summarizing the values entered ----
    
    plotOutput("essPlot"),
    
    hr(),
    
    fluidRow(
        column(4,
               h4("Leaf Trait Linkages "),
               # Input: Leaf Lifespan ----
               sliderInput("c", "Leaf lifespan (c):",
                           min = 0, max = 35,
                           value =28.57, step=1,
                           animate = TRUE),
               
               # Input: Functional Leaf N  ----
               sliderInput("A", "Functional Leaf N (A):",
                           min = 0, max = 3,
                           value = 1.5, step = 0.25,
                           animate = TRUE),
               
               
               # Input: Structural Leaf N ----
               sliderInput("B", "Structural Leaf N (B):",
                           min = 0, max = 16,
                           value = 8, step = 1,
                           animate = TRUE),
               
               # Input: Leaf respiration parameter ----
               sliderInput("r", "Leaf respiration (r):",
                           min = 0, max = 0.25,
                           value = 0.015, step = 0.005,
                           animate = TRUE)
               
        ),
        column(4,
               h4("Canopy Carbon Cycle"),
               
               # Input: Photosynthesis parameter  ----
               sliderInput("V", "Photosynthesis parameter (V):",
                           min = 0, max = 5,
                           value = 1.2, step = 0.1,
                           animate = TRUE),
               
               # Input: Saturate rate of GPP with LAI ----
               sliderInput("k", "Saturation rate of GPP with LAI (k):",
                           min = 0, max = 2,
                           value = 0.5, step = 0.1,
                           animate = TRUE),
               
               # Input: Growth Carbon Cost ----
               sliderInput("G", "Growth Carbon Cost (G):",
                           min = 0, max = 5,
                           value = 1.3333, step = 0.05,
                           animate = TRUE)
               

        ),
        column(4,
               h4("Nitrogen Cycle"),
               
               # Input: Soil N residence time parameter ----
               sliderInput("s", "Soil N residence time (s):",
                           min = 0, max = 1500,
                           value = 500, step = 0.05,
                           animate = TRUE),   
               
               # Input: Minimum LMA ----
               sliderInput("sigma_min", "Minimum LMA (sigma_min):",
                           min = 0, max = 1,
                           value = 0.02, step = 0.01,
                           animate = TRUE),     
               
                # Input: lma_min  ----
                sliderInput("lma_min", "lma - plot limit",
                            min = -2, max = 2,
                            value = c(0,1)),   
        
                # Input: lma_min  ----
                sliderInput("n_ref_min", "N_ref_min - plot limit",
                            min = -50, max = 50,
                            value = c(0,40))
               )
    )
)



# Define server logic required to draw a histogram
server <- function(input, output) {

    output$essPlot <- renderPlot({
        
        c = input$c
        A = input$A
        B = input$B
        r = input$r
        V = input$V
        k = input$k
        G = input$G
        s = input$s
        sigma_min = input$sigma_min
        
        # generate curves based on input$parms parametersfrom ui.R
        
        # ESS Thresholds
        sigma_max <- (sqrt((V-G/c)/(r*A)) - 1)*A/B
        sigma_ess_min <- (sqrt(V*A/((exp(1)^2*r)))- A)/B
        
        n_1 <- (A + B*sigma_min)/(k*c*sigma_min)*log(V/((A+B*sigma_min)^2*r/A+ G/c))
        n_2 <- A/(exp(1)*k*c*sigma_min)*sqrt(V/(r*A))*log(V/(V/exp(1)^2 + G/c))
        
        # LMA curve
        y <- seq(input$lma_min[1],input$lma_min[2],length = 20000) # sequence of LMA values
        
        fig_3_transpose <- tibble(lma = y,
                                  n_ref_ess = (A + B*y)/(k*c*sigma_min)*log(V/((A+B*y)^2*r/A + G/c)),
                                  n_ref_opt = (A + B*y)/(k*c*y)*log(V/((A+B*y)^2*r/A + G/c)))
        
        ## VIZ - draw the ESS lma and N relationship based on the parameters
        ## ESS LMA
        ggplot(fig_3_transpose) + 
            #thresholds
            geom_hline(yintercept = n_1, color = "blue2", linetype = 2) + 
            geom_hline(yintercept = n_2, color = "red", linetype = 2) + 
            geom_vline(xintercept = sigma_min, color = "green4") + 
            geom_vline(xintercept = sigma_max, color = "green4") + 
            geom_line(aes(x=lma, y= n_ref_ess), size = 1) +
            #geom_segment(aes(x = sigma_min, xend = sigma_min, y = n_1, yend = 50), size = 1) + 
            ylim(input$n_ref_min[1], input$n_ref_min[2]) +
            xlim(input$lma_min[1], input$lma_min[2]) + 
            coord_flip() + 
            # Formatting
            theme_bw() +
            ggtitle("Figure 3") +
            labs(y = "Reference N mineralization rate", x = "LMA")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
