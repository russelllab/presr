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

#data <- read_tsv("../scripts/ML/all_predictions.tsv")
data <- read_tsv("all_predictions.tsv")
new_data <- data
new_data$Name = paste0(new_data$WT,new_data$Pos)

# Define UI for application that draws a histogram
ui <- navbarPage("PRESS",
        tabPanel("VISUALIZATION", icon = icon("chart-area"),
            fluidPage(
              tags$style(
                HTML('
                     #structure_panel {
                     display: flex;
                     align-items: center;
                     justify-content: center;
                     top: 50%;
                     }
                     ')
                  ),
              fluidRow(
                column(2,
                      align="center",
                      wellPanel(
                                useShinyjs(),
                                
                                style = "height: 300px;",
                                
                                div(style="display: inline-block;vertical-align:center; padding-bottom: 10%;",
                                  h1(id="big-heading", "Input panel"),
                                  tags$style(HTML("#big-heading{color: black; font-size: 25px; text-align: center; text-decoration: underline;}"))
                                  ),
                                
                                div(style="display: inline-block; padding-left: 3%;", actionButton("inputhelp", "", icon = icon("question-circle"))),
                                
                                div(style="display: inline-block; vertical-align:top; font-size: 18px;",
                                  selectInput("position", "Munc18-1 position",
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
                            h1(id="big-heading", "Output panel")
                        ),
                        div(style="display: inline-block; padding-left: 3%;", actionButton("outhelp", "", icon = icon("question-circle"))),
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
                             h1(id="big-heading", "Structure panel")
                         ),
                        div(style="display: inline-block; padding-left: 3%;", actionButton("inputhelp", "", icon = icon("question-circle"))),
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
                             h1(id="big-heading", "Homology panel")
                         ),
                         div(style="display: inline-block; padding-left: 3%;", actionButton("inputhelp", "", icon = icon("question-circle"))),
                         plotlyOutput(outputId = "plot1", height="300px"),
                         plotlyOutput(outputId = "plot2", height="300px")
                       )
                ),
                column(3,
                       #style='border: 1px solid grey',
                       align="center",
                       wellPanel(
                         style = "height: 700px;",
                         div(style="display: inline-block;vertical-align:center; padding-bottom: 10px;",
                             h1(id="big-heading", "Other panel")
                         ),
                         div(style="display: inline-block; padding-left: 3%;", actionButton("inputhelp", "", icon = icon("question-circle"))),
                         plotlyOutput(outputId = "plot3", height="300px")
                         #msaROutput(outputId = 'msa1'),
                         #msaROutput(outputId = 'msa2')
                       )
                )
              ),
              
            )
        ),
        tabPanel("DATA", icon = icon("table"),
                 fluidPage(theme = shinytheme("cerulean"),
                   tags$head(tags$style(HTML(".shiny-output-error-validation{color: red;}"))),
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
                   tags$head(tags$style(HTML(".shiny-output-error-validation{color: red;}"))),
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
                        height="350px;",
                        align = "center",
                        div(style="padding-bottom: 20px;",
                            h1(id="big-heading", "Paralogs")
                        ),
                        msaROutput(outputId = 'msa2')
                      )
                     )
                   )
                 )
        ),
        tabPanel("ABOUT", icon = icon("info-circle"),
                 fluidPage(theme = shinytheme("cerulean")),
                 tags$head(tags$style(HTML(".shiny-output-error-validation{color: red;}")))
        )
)

default_structure <- function(position){
    struc <- NGLVieweR("4JEU") %>%
      stageParameters(backgroundColor = "white") %>%
      addRepresentation("cartoon", param = list(name = "cartoon",
                                                sele = ":A",
                                                #sele = paste(substr(position, start = 2, stop = nchar(position)), ":A and .CA", sep=' '),
                                                #labelType = "format",
                                                #labelFormat = position, # or enter custom text
                                                #labelGrouping = "residue", # or "atom" (eg. sele = "20:A.CB")
                                                colorValue = "grey",
                                                colorScheme = "element"
                                                #colorScheme="residueindex"
                                                )) %>%
      setQuality("high") %>%
      setFocus(0) %>%
      setSpin(TRUE)
      
    return(struc)
}

make_color <- function(position, mut){
  pos <- substr(position, start = 2, stop = nchar(position))
  row_index <- match(c(pos), new_data$`Pos`)
  mut_val <- as.numeric(new_data[row_index, mut])
  print(mut_val)
  #return("green") 
  if (mut_val > 0.5) {
    return("red") 
  }
  else {
    return("green") 
  }
}

make_text <- function(position, mut) {
  pos <- substr(position, start = 2, stop = nchar(position))
  row_index <- match(c(pos), new_data$`Pos`)
  mut_val <- as.numeric(new_data[row_index, mut])
  print(mut_val)
  if (mut_val > 0.5) {
    return("<div style='font-size:20px; color: #c32607; display:inline;'><b><br>Diseased<br></b></div>") 
  }
  else {
    return("<div style='font-size:20px; color: #0e6c1a; display:inline;'><b><br>Neutral<br></b></div>") 
  }
}

make_table <- function(position, mut) {
  pos <- substr(position, start = 2, stop = nchar(position))
  row_index <- match(c(pos), new_data$`Pos`)
  mut_val <- as.numeric(new_data[row_index, mut])
  wt <- new_data$`WT`[row_index]
  wt_val <- as.numeric(new_data[row_index, wt])
  TYPE <- c('Wild Type', 'Variant')
  VAL <- c(wt_val, mut_val)
  #VAL <- c(wt_val, mut_val, mut_val-wt_val)
  AA <- c(wt, mut)
  #AA <- c(wt, mut, '\u0394')
  df <- data.frame(TYPE, AA, VAL)
  names(df)[3] <- "PROB"
  return(df)
}

make_msa <- function(x) {
  if (x==1) {
    #proteinseqfile <- "/home/gurdeep/projects/munc18-1/webApp/para.fasta"
    proteinseqfile <- "para.fasta"
  }
  else {
    #proteinseqfile <- "/home/gurdeep/projects/munc18-1/webApp/ortho.fasta"
    proteinseqfile <- "ortho.fasta"
  }
  proteins <- ape::read.FASTA(proteinseqfile, type="AA")
  return(msaR(proteins, overviewbox = F,  colorscheme = "clustal2"))
}

subsetting <- function (position, df) {
  value <- df %>% filter(`Mut` == 'D49H')
  print (value)
  return (as.double(value$`Score`))
}

plotting <- function(df, pos, mut, titlename) {
  mut_score <- (df %>% filter(`Mut` == paste0(pos, mut)))$Score
  wt_aa <-substr(pos, start = 1, stop = 1)
  wt_score <- (df %>% filter(`Mut` == paste0(pos, wt_aa)))$Score
  print (wt_score)
  
  p <- df %>% filter(`Label` != 'Unassigned')
  gg_diag <- ggplot(p, aes(x=`Score`, color=`Label`, fill=`Label`)) + geom_histogram(bins=100) + geom_density(alpha=0.6) +
    geom_vline(aes(xintercept=mut_score), color="red", linetype="dashed", size=0.5) +
    geom_vline(aes(xintercept=wt_score), color="black", linetype="dashed", size=0.5) +
    ggtitle(titlename) +
    ylab("Count") +
    theme(axis.text=element_text(size=10)) +
    theme(axis.title=element_text(size=10)) +
    theme(legend.text = element_text(size = 10)) +
    theme(plot.title = element_text(size = 10, face = "bold")) +
    theme(legend.title=element_blank()) +
    theme(legend.position="bottom")
  
  if (titlename != 'DynaMine'){
    gg_diag <- gg_diag + annotate(geom="text", x=wt_score+1.0, y=10, label="WT", color="black", size = 4) +
                         annotate(geom="text", x=mut_score-1.75, y=10, label=paste0(pos,mut), color="red",size = 4) +
                         xlim(-10, 2)
  }
  else {
    gg_diag <- gg_diag + annotate(geom="text", x=wt_score+0.075, y=20, label=pos, color="black", size=4)
    #label=substr(pos, start = 2, stop = nchar(pos))
  }
  
  return(ggplotly(gg_diag))
}

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$mutation2 <- renderUI({
      selectInput("mutation", "Mutation",
                #choices = levels(factor(new_data$`WT`)),
                #choices = substr(input$position, 1, 1)
                choices = levels(factor(new_data$`WT`))[levels(factor(new_data$`WT`))!=substr(input$position, 1, 1)]
                #selected = "P"
                )
  })
    observeEvent(input$inputhelp, {
      showModal(modalDialog(
        title = "Help",
        HTML("<kbd>Munc18-1 position</kbd>: Munch18-1 positions<br>
        <kbd>Mutation</kbd>: Substituion"
        )
      ))
    })
  observeEvent(input$outputhelp, {
    showModal(modalDialog(
      title = "Help",
      HTML("<kbd>TYPE</kbd>: Amino acid type<br>
        <kbd>AA</kbd>: Amino acid<br>
        <kbd>PROB</kbd>: Predicted probability<br>"
      )
    ))
  })
    #output$hist <- renderPlot({hist(rnorm(input$num))})
    #data <- eventReactive(input$submit, {hist(rnorm(input$num), col=input$color)})
    #data <- eventReactive(input$submit, {ggplotly(ggplot(gapminder, aes(gdpPercap, size = pop, color=continent)) + geom_histogram(bins=input$num))})

    #df1 <- read_tsv(file = "para.tsv")
    #output$plot1 <- renderPlotly({
    #  plotting(df1, input$position, input$mutation, 'Paralogs')
    #  })
    observeEvent(input$mutation, {df1 <- read_tsv(file = "para.tsv")
    output$plot1 <- renderPlotly({
      plotting(df1, input$position, input$mutation, 'Paralogs')
    })})
    
    #df2 <- read_tsv(file = "ortho.tsv")
    #output$plot2 <- renderPlotly({
    #  plotting(df2, input$position, input$mutation, 'Orthologs')
    #})
    observeEvent(input$mutation, {df2 <- read_tsv(file = "ortho.tsv")
    output$plot2 <- renderPlotly({
      plotting(df2, input$position, input$mutation, 'Orthologs')
    })})
    
    #df3 <- read_tsv(file = "dyna.tsv")
    #output$plot3 <- renderPlotly({
    #  plotting(df3, input$position, input$mutation, 'DynaMine')
    #})
    observeEvent(input$mutation, {df3 <- read_tsv(file = "dyna.tsv")
    output$plot3 <- renderPlotly({
      plotting(df3, input$position, input$mutation, 'DynaMine')
    })})
    
    datatable(new_data)%>%
      formatStyle(
        'A',
        color = styleInterval(c(3.4, 3.8), c('white', 'blue', 'red')),
        backgroundColor = styleInterval(3.4, c('gray', 'yellow'))
      )
    output$table <- DT::renderDataTable({
        DT::datatable(new_data,
                      extensions = list("ColReorder" = NULL,
                                                  "Buttons" = NULL,
                                                  "FixedColumns" = list(leftColumns=1)),
                      options = list(
                        dom = 'lBRrftpi',
                        autoWidth=TRUE,
                        pageLength = 15,
                        lengthMenu = list(c(5, 10, 15, 50, -1), c('5', '10', '15','50', 'All')),
                        ColReorder = TRUE,
                        buttons =
                          list(
                            'copy',
                            'print',
                            list(
                              extend = 'collection',
                              buttons = c('csv', 'excel', 'pdf'),
                              text = 'Download'
                            )
                          )
                      )
        )%>%
        formatStyle(
          c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'),
          #color = styleInterval(c(3.4, 3.8), c('white', 'blue', 'red')),
          backgroundColor = styleInterval(0.5, c('lightgreen', 'indianred'))
        )
    })
    #out_table <- make_table(input$position, input$mutation)
    #out_table <- eventReactive(input$submit, {make_table(input$position, input$mutation)})
    #output$otext <- renderText({paste('<div style="font-size: 17px">Predicted as ', make_text(input$position, input$mutation), 'variant</div>')})
    observeEvent(input$mutation,{output$otext <- renderText({paste('<div style="font-size: 18px">Predicted as ', make_text(input$position, input$mutation), 'variant</div>')})})
    output$msa1 <- renderMsaR({make_msa(1)})
    output$msa2 <- renderMsaR({make_msa(2)})
    #output$otable <- renderTable(align = 'c',{make_table(input$position, input$mutation)})
    observeEvent(input$mutation, output$otable <- renderTable(align = 'c',{make_table(input$position, input$mutation)}))
    #default_struc <- default_structure()
    output$structure <- renderNGLVieweR({default_structure(input$position)})
    observeEvent(input$spin, {NGLVieweR_proxy("structure") %>% updateSpin(input$spin)})
    observeEvent(input$fullscreen, {NGLVieweR_proxy("structure") %>% updateFullscreen()})
    observeEvent(input$snapshot, {NGLVieweR_proxy("structure") %>% snapShot("Snapshot",
                                                                            param = list(
                                                                              antialias = TRUE,
                                                                              trim = TRUE,
                                                                              transparent = TRUE,
                                                                              scale = 1
                                                                            )
    )})
    toListen <- reactive({
      list(input$mutation)
    })
    observeEvent(input$mutation,
                 {
                   NGLVieweR_proxy("structure") %>%
                     removeSelection("label1") %>%
                     addSelection("label",
                                          param = list(
                                            name = "label1",
                                            #sele = paste(input$position, ":A and .CA", sep=' '),
                                            sele = paste(substr(input$position, start = 2, stop = nchar(input$position)), ":A and .CA", sep=' '),
                                            labelType = "format",
                                            labelFormat = paste(input$position,input$mutation, sep=""), # or enter custom text
                                            labelGrouping = "residue", # or "atom" (eg. sele = "20:A.CB")
                                            #color = "red",
                                            color = make_color(input$position,input$mutation),
                                            fontFamiliy = "sans-serif",
                                            xOffset = 2,
                                            yOffset = 0,
                                            zOffset = 0,
                                            fixedSize = TRUE,
                                            radiusType = 1,
                                            radiusSize = 2.5, # Label size
                                            showBackground = TRUE,
                                            backgroundColor="white"
                                          )
                     ) %>%
                     updateZoomMove(
                       #center = paste(input$position, ":A and .C", sep=''),
                       center = paste(substr(input$position, start = 2, stop = nchar(input$position)), ":A and .CA", sep=' '),
                       #zoom = paste(input$position, ":A and .C", sep=''),
                       20,
                       z_offSet = 80,
                       duration = 2000
                     ) %>%
                     removeSelection("sel1") %>%
                     addSelection("spacefill",
                                  param = list(name="sel1",
                                        sele = paste(substr(input$position, start = 2, stop = nchar(input$position)), ":A and .CA", sep=' '),
                                        #colorValue="grey",
                                        color = make_color(input$position,input$mutation),
                                        colorScheme="element")
                                  )
                   
                 }
      )
}

# Run the application 
shinyApp(ui = ui, server = server)

#removeSelection("sel1") %>%
 # addSelection("spacefill", param = list(name="sel1",sele = paste(input$position, ":A and .CA", sep=' '),colorValue="yellow",colorScheme="element"))