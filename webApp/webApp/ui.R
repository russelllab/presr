library(shiny)
library(gapminder)
library(plotly)
library(ggplot2)
library(readr)
library(NGLVieweR)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(msaR)
library(DT)

data <- read_tsv("../scripts/ML/all_predictions.tsv")
#data <- read_tsv("all_predictions.tsv")
new_data <- data
new_data$Name = paste0(new_data$WT,new_data$Pos)

# Define UI for application that draws a histogram
shinyUI(navbarPage("PRESR: PREdictor of STXBP1 Related disorder",
        tabPanel("VISUALIZATION", icon = icon("chart-area"),
            fluidPage(
              tags$style(
                HTML('
                     #structure_panel {
                     display: flex;
                     align-items: center;
                     justify-content: center;
                     top: 50%;
                     font-size: 20px;
                     }
                     ')
                  ),
              fluidRow(
                #column(1),
                column(2,
                      align="center",
                      wellPanel(
                                useShinyjs(),
                                
                                style = "height: 300px;",
                                
                                div(style="display: inline-block;vertical-align:center; padding-bottom: 10%;",
                                  h1(id="big-heading", "Input"),
                                  tags$style(HTML("#big-heading{color: black; font-size: 25px; text-align: center; text-decoration: underline;}"))
                                  ),
                                
                                div(style="display: inline-block; padding-left: 3%;", actionButton("inputHelp", "", icon = icon("question-circle"))),
                                
                                div(style="display: inline-block; vertical-align:top; font-size: 18px;",
                                  selectInput("position", "MUNC18-1 position",
                                              #choices = levels(factor(new_data$`Name`)),
                                              choices = new_data$`Name`,
                                              selected = "V84")
                                  ),
   
                                div(style="vertical-align:top; width: 75px; font-size: 18px;",
                                  uiOutput("mutation2")
                                )
                            ),
                      wellPanel(
                        style = "height: 375px;",
                        div(style="display: inline-block;vertical-align:center; padding-bottom: 10%;",
                            h1(id="big-heading", "Output")
                        ),
                        div(style="display: inline-block; padding-left: 3%;", actionButton("outputHelp", "", icon = icon("question-circle"))),
                        div(style="font-size:18px;", tableOutput(outputId = 'otable')),
                        div(style="font-size:18px;", htmlOutput(outputId = 'otext')),
                      )
                      
                      ),
                column(4,
                       #style='border: 1px solid grey;',
                       align = "center",
                       wellPanel(
                         style = "height: 700px;",
                         div(style="display: inline-block;vertical-align:center; padding-bottom: 10px;",
                             h1(id="big-heading", "Structure")
                         ),
                        div(style="display: inline-block; padding-left: 3%;", actionButton("structureHelp", "", icon = icon("question-circle"))),
                        NGLVieweROutput(outputId = 'structure', height = "550px"),
                        div(style="display: inline-block; font-size: 15px; padding-top: 3%;",
                            switchInput(inputId = "spin", 'Spin', value = TRUE)
                        ),
                        div(style="display: inline-block; font-size: 15px; padding-top: 3%;",
                            actionButton(inputId = "snapshot", 'Snapshot', icon("fas fa-image"))
                        ),
                        div(style="display: inline-block; font-size: 15px; padding-top: 3%;",
                            actionButton(inputId = "fullscreen", 'Fullscreen', icon("fas fa-expand"))
                        )
                        #plotlyOutput(outputId = "plot3")
                       )
                ),
                column(3,
                       align="center",
                       #NGLVieweROutput(outputId = 'structure'),
                       wellPanel(
                         style = "height: 700px;",
                         div(style="display: inline-block;vertical-align:center; padding-bottom: 10px;",
                             h1(id="big-heading", "Sequence-based features")
                         ),
                         div(style="display: inline-block; padding-left: 3%;", actionButton("seqFeaturesHelp", "", icon = icon("question-circle"))),
                         #splitLayout(cellWidths = c("50%", "50%"),
                           plotlyOutput(outputId = "plot1", height="300px", width="400px"),
                           plotlyOutput(outputId = "plot2", height="300px", width="400px")
                         #)
                       )
                ),
                column(3,
                       align="center",
                       #NGLVieweROutput(outputId = 'structure'),
                       wellPanel(
                         style = "height: 700px;",
                         div(style="display: inline-block;vertical-align:center; padding-bottom: 10px;",
                             h1(id="big-heading", "Structure-based features")
                         ),
                         div(style="display: inline-block; padding-left: 3%;", actionButton("strucFeaturesHelp", "", icon = icon("question-circle"))),
                         #splitLayout(cellWidths = c("50%", "50%"),
                                     plotlyOutput(outputId = "plot3", height="300px", width="400px"),
                                     plotlyOutput(outputId = "plot4", height="300px", width="400px"),
                         #)
                       )
                )
                #column(1)
                #column(1,
                       #align="center",
                       #wellPanel(
                      #   style = "height: 700px;",
                       #  div(style="display: inline-block;vertical-align:center; padding-bottom: 10px;",
                        #     h1(id="big-heading", "Other panel")
                         #),
                         #div(style="display: inline-block; padding-left: 3%;", actionButton("inputhelp", "", icon = icon("question-circle"))),
                         #plotlyOutput(outputId = "plot3", height="300px"),
                         #plotlyOutput(outputId = "plot4", height="300px")

                       #)
                #)
              ),
              
            )
        ),
        tabPanel("DATA", icon = icon("table"),
                 fluidPage(theme = shinytheme("cerulean"),
                   tags$head(tags$style(HTML(".shiny-output-error-validation{color: mediumorchid;}"))),
                   fluidRow(
                     column(10, offset=1,
                      wellPanel(
                        align = "center",
                        div(style="padding-bottom: 20px;",
                             h1(id="big-heading", "Predicted probabilities")
                        ),
                        #mainPanel(DT::dataTableOutput("table"))
                        DT::dataTableOutput("table")
                      )
                    )
                   )
                 )
        ),
        tabPanel("ALIGNMENTS", icon = icon("fas fa-align-justify"),
                 fluidPage(theme = shinytheme("cerulean"),
                   tags$head(tags$style(HTML(".shiny-output-error-validation{color: mediumorchid;}"))),
                   fluidRow(
                     column(12,
                            
                      wellPanel(
                        align = "center",
                        div(style="padding-bottom: 20px;",
                            h1(id="big-heading", "Orthologs")
                        ),
                        msaROutput(outputId = 'msa1')
                      )
                      )
                     ),
                   fluidRow(
                     column(12,
                      wellPanel(
                        align = "center",
                        div(style="padding-bottom: 20px;",
                            h1(id="big-heading", "Paralogs (human-specific)")
                        ),
                        msaROutput(outputId = 'msa2')
                      )
                     )
                   )
                 )
        ),
        tabPanel("ABOUT", icon = icon("github", lib = "font-awesome"),
                 fluidPage(theme = shinytheme("cerulean"),
                           fluidRow(
                             column(6,
                                    wellPanel(
                                      align = 'center',
                                      useShinyjs(),
                                      htmlOutput("about"),
                                      actionButton(inputId = 'github',
                                                   label = '',
                                                   icon = icon("github"),
                                                   onclick ="window.open('https://github.com/russelllab/presr', '_blank')",
                                                   style='padding:4px; font-size:250%'
                                      )
                                    )
                             ),
                             column(6,
                                    imageOutput(outputId = "workflow")
                             )
                           )
                 ),
                 tags$head(tags$style(HTML(".shiny-output-error-validation{color: mediumorchid;}")))
        )
  )
)


# Run the application 
#shinyApp(ui = ui, server = server)