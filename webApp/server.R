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
    return("mediumorchid") 
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
  if (mut_val >= 0.5) {
    return("<div style='font-size:20px; color: mediumorchid; display:inline;'><b><br>Disease<br></b></div>") 
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
  if (x==2) {
    #proteinseqfile <- "/home/gurdeep/projects/munc18-1/webApp/para.fasta"
    proteinseqfile <- "spec_para2.fasta"
  }
  else {
    #proteinseqfile <- "/home/gurdeep/projects/munc18-1/webApp/ortho.fasta"
    proteinseqfile <- "ortho2.fasta"
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
  
  p <- df %>% filter(`Variant` != 'Unassigned')
  gg_diag <- ggplot(p, aes(x=`Score`, color=`Variant`, fill=`Variant`)) +
    #geom_histogram(bins=100) +
    geom_density(alpha=0.6) +
    geom_vline(aes(xintercept=mut_score), color="mediumorchid", linetype="dashed", size=0.5) +
    geom_vline(aes(xintercept=wt_score), color="black", linetype="dashed", size=0.5) +
    ggtitle(titlename) +
    theme(axis.text=element_text(size=12)) +
    theme(axis.title=element_text(size=12)) +
    theme(legend.text = element_text(size = 10)) +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    theme(legend.position="bottom") +
    scale_fill_manual(values=c("mediumorchid", "darkgreen"))+
    scale_color_manual(values=c("mediumorchid", "darkgreen"))
  
  
  #if (titlename != 'DynaMine'){
  if (titlename != 'Protein backbone dynamics' && titlename != 'Side chain residues'){
    gg_diag <- gg_diag + annotate(geom="text", x=wt_score+1.0, y=0.4, label="WT", color="black", size = 4) +
                         annotate(geom="text", x=mut_score-1.0, y=0.4, label=paste0(pos,mut), color="mediumorchid",size = 4) +
                         xlim(-7.5, 2) +
                        xlab("Score") +
                        ylab("Fraction of positions")
  }
  else if (titlename == 'Side chain residues'){
    gg_diag <- gg_diag +
      annotate(geom="text", x=wt_score+1.0, y=0.4, label=pos, color="black", size = 4) +
      xlim(-2, 10) +
      xlab("Number of contacts") +
      ylab("Fraction of positions")
  }
  else {
    gg_diag <- gg_diag + annotate(geom="text", x=wt_score+0.075, y=10, label=pos, color="black", size=4) +
                        xlim(0.6, 1.15) +
                        xlab("Decreasing backbone dynamics") +
                        ylab("Fraction of positions")
    #label=substr(pos, start = 2, stop = nchar(pos))
  }
  
  return(ggplotly(gg_diag))
}

plotting2 <- function(df, pos, mut, titlename) {
  mut_score <- (df %>% filter(`Name` == paste0(pos, mut)))$MUT
  wt_score <- (df %>% filter(`Name` == paste0(pos, mut)))$WT
  #wt_aa <-substr(pos, start = 1, stop = 1)
  #wt_score <- (df %>% filter(`Name` == paste0(pos, wt_aa)))$Score
  print (wt_score)
  
  df_pos <- df %>% filter(`Pos` == pos)
  
  gg_diag <- ggplot(df_pos, aes(x=`MUT`)) + geom_density(alpha=0.6) +
    geom_vline(aes(xintercept=mut_score), color="mediumorchid", linetype="dashed", size=0.5) +
    geom_vline(aes(xintercept=wt_score), color="black", linetype="dashed", size=0.5) +
    ggtitle(titlename) +
    ylab("Count") +
    theme(axis.text=element_text(size=10)) +
    theme(axis.title=element_text(size=10)) +
    theme(legend.text = element_text(size = 10)) +
    theme(plot.title = element_text(size = 10, face = "bold")) +
    theme(legend.title=element_blank()) +
    theme(legend.position="bottom") +
    annotate(geom="text", x=wt_score+0.35, y=0.1, label="WT", color="black", size = 4) +
    annotate(geom="text", x=mut_score-0.5, y=0.1, label=paste0(pos,mut), color="mediumorchid",size = 4)
  
  return(ggplotly(gg_diag))
}

# Define server logic required to draw a histogram
shinyServer(
  function(input, output) {
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
    #observeEvent(input$mutation, {df1 <- read_tsv(file = "para.tsv")
    #output$plot1 <- renderPlotly({
    #  plotting(df1, input$position, input$mutation, 'Paralogs')
    #})})
    
    #df2 <- read_tsv(file = "ortho.tsv")
    #output$plot2 <- renderPlotly({
    #  plotting(df2, input$position, input$mutation, 'Orthologs')
    #})
    observeEvent(input$mutation, {df1 <- read_tsv(file = "ortho.tsv")
    output$plot1 <- renderPlotly({
      plotting(df1, input$position, input$mutation, 'Conservation across orthologs')
    })})
    
    output$workflow <- renderImage({
      list(src = "webAppWorkflow.png", contentType = 'image/png',
           width = 800,
           height = 650,
           alt = "Workflow")
    }, deleteFile = FALSE)
    #df3 <- read_tsv(file = "dyna.tsv")
    #output$plot3 <- renderPlotly({
    #  plotting(df3, input$position, input$mutation, 'DynaMine')
    #})
    observeEvent(input$inputHelp, {
      showModal(modalDialog(
        title = "Input Help",
        HTML("<div style='font-size: 15px'>Choose a mutation in Munc18-1 protein.Please note that this predictor was trained on the <a href='https://www.uniprot.org/uniprot/P61764'>canonical/isoform 1</a> of Munc18-1 reported in UniProt.</div>"
        )
      ))
    })
    
    observeEvent(input$outputHelp, {
      showModal(modalDialog(
        title = "Output Help",
        HTML("<div style='font-size: 15px'>Predicted probabilities of wild type and variant to be disease.<br><kbd>Type</kbd>: Wild type or the variant<br><kbd>AA</kbd>: Amino acid<br><kbd>PROB</kbd>: PRESS predicted probability to be a disease variant (neutral is <0.5, otherwise disease)</div>"
        )
      ))
    })
    
    observeEvent(input$structureHelp, {
      showModal(modalDialog(
        title = "Structure Help",
        HTML("<div style='font-size: 15px'>
             The selected MUNC18-1 position is mapped onto the 3D structure (PDB-ID: <a href='https://www.rcsb.org/structure/4JEU'>4JEU:A</a>).
             The given position is displayed in mediumorchid if the variant was predicted to be be disease associated, and in green otherwise (neutral associted).
             </div>"
        )
      ))
    })
    
    observeEvent(input$seqFeaturesHelp, {
      showModal(modalDialog(
        title = "Structure-based features",
        HTML("<div style='font-size: 15px'>
              <p>
              <kbd>Orthologs</kbd>:<br>
              The plot displays disitribution of difference in conservation scores (orthologs) between the variant and the wild type (WT) across known MUNC18-1 variants.
              The difference in conservation score is always zero for the WT, and is shown with a black dotted line.
              The difference in conservation score of any other variant at a given position is shown with a mediumorchid dotted line.
              Lower the difference in conservation score (more negative), higher is the conservation of the WT amino acid at the given position, and higher is the probability of the variant to be disease associated.
              </p>
              <kbd>Paralogs</kbd>:<br>
              The plot displays disitribution of difference in conservation scores (human-specific paralogs) between the variant and the wild type (WT) across known MUNC18-1 variants.
              The difference in conservation score is always zero for the WT, and is shown with a black dotted line.
              The difference in conservation score of any other variant at a given position is shown with a mediumorchid dotted line.
              Lower the difference in conservation score (more negative), higher is the conservation of the WT amino acid at the given position, and higher is the probability of the variant to be disease associated.
              </p>
              </div>
             ")
      ))
    })
    
    observeEvent(input$strucFeaturesHelp, {
      showModal(modalDialog(
        title = "Structure-based features",
        HTML("<div style='font-size: 15px'>
              <p>
              <kbd>Protein backbone dynamics</kbd>:<br>
              The plot displays disitribution of predicted backbone dynamics (using <a href='https://bio2byte.be/dynamine/#'>DynaMine</a>) of WT MUNC18-1 positions.
              A higher score (>0.8) refers to high rigidity, and a  lower score (<0.7) refers to low rigidity (high flexibility).
              The position selected in the input panel is shown with a black dotted line.
              Higher the rigidity (low flexibility), higher is the probability of the variant to be disease associated.
              </p>
              <p>
              <kbd>Side-chain residues</kbd>:<br>
              The plot displays disitribution of number of residues in contact with the side-chain of the WT MUNC18-1 (using <a href='https://www.uniprot.org/uniprot/P61764#structure'>AlphaFold2</a>).
              The position selected in the input panel is shown with a black dotted line.
              Higher the number of contacts the side-chain of the WT residue makes, higher is the probability of the variant to be disease associated.
              </p>
              </div>
             ")
      ))
    })
    
    getPage<-function() {
      return(includeHTML("about.html"))
    }
    output$about<-renderUI({getPage()})
    
    observeEvent(input$mutation, {df3 <- read_tsv(file = "dyna.tsv")
    #observeEvent(input$mutation, {df2 <- read_tsv(file = "spec_para.tsv")
    output$plot3 <- renderPlotly({
      plotting(df3, input$position, input$mutation, 'Protein backbone dynamics')
    })})
    
    observeEvent(input$mutation, {df2 <- read_tsv(file = "spec_para.tsv")
    output$plot2 <- renderPlotly({
      plotting(df2, input$position, input$mutation, 'Conservation across paralogs')
    })})
    
    observeEvent(input$mutation, {df4 <- read_tsv(file = "side_chain_residues.tsv")
    output$plot4 <- renderPlotly({
      plotting(df4, input$position, input$mutation, 'Side chain residues')
    })})
    
    datatable(new_data)%>%
      formatStyle(
        'A',
        color = styleInterval(c(3.4, 3.8), c('white', 'blue', 'lightmediumorchid')),
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
          #color = styleInterval(c(3.4, 3.8), c('white', 'blue', 'mediumorchid')),
          backgroundColor = styleInterval(0.5, c('lightgreen', 'mediumorchid'))
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
)
# Run the application 
#shinyApp(ui = ui, server = server)