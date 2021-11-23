library(shiny)
library(tidyverse)
library(MultiAssayExperiment)
library(plotly)
library(DT)
library(org.Hs.eg.db)
library(shinyWidgets)
library(officer)
library(rvg)
library(colourpicker)
library(randomcoloR)
library(shinycssloaders)
library(shinythemes)
options(stringsAsFactors = F)

# load CCLE data
ccledat <- readRDS("/poolio/public_data/CCLE/CCLE_MAE_processed.2019-10-27.rds")

# just make a summary table for the assay info
assay_info <- data.frame("Assay Type" = names(ccledat$mae),
                         "Sample #" = sapply(names(ccledat$mae), function(x) {
                           ncol(assay(ccledat$mae[,,x])) }),
                         "Feature #" = sapply(names(ccledat$mae), function(x) {
                           nrow(assay(ccledat$mae[,,x])) }),
                         check.names = F)

# organise tsne dim data
tsne_dims <- lapply(names(ccledat$dimred), function(x)  {
  y <- cbind(ccledat$dimred[[x]], colData(ccledat$mae)[ccledat$cc,]) %>%
    as.data.frame()
  y$assay_type <- x
  return(y)
}) %>% bind_rows()

# function to generate pptx
gen_pptx <- function(plot, file, height = 5, width = 5, left = 1, top = 1) {
  read_pptx() %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = dml(ggobj = plot), 
            location = ph_location(height = height, width = width,
                                   left = left, top = top),
            bg = "transparent") %>%
    print(target = file)
}

############### ui ############### 
ui <- fluidPage(
  
  title = "CCLE", 
  theme = shinytheme("flatly"),
  titlePanel(tags$h3(tags$a(
    imageOutput("icon", inline = TRUE),
                href="http://10.197.211.94:3838"), 
    "Cancer Cell Line Encyclopedia")),
  
  tabsetPanel(
    tabPanel("Sample Info",
             plotOutput("pheno_tsne_plot") %>% withSpinner(color="#0dc5c1"),
             br(),
             DT::DTOutput('sample_info_dt')),
    tabPanel("Assay Info",
             checkboxGroupButtons("plot_type", "Plot Element:", 
                                  choices = list("Violin" = "violin",
                                                 "Boxplot" = "boxplot",
                                                 "Dotplot" = "dotplot"),
                                  justified = T, 
                                  size = 'xs', width = '25%',
                                  selected = "violin"),
             plotOutput("feature_levels_plot") %>% withSpinner(color="#0dc5c1"),
             div(style="display:inline-block",
                 downloadButton("dl_feature_levels_plot_ppt", label = "PPT"),
                 downloadButton("dl_feature_levels_plot_png", label = "PNG")),
             br(), br(),
             br(), br(),
             DT::DTOutput('assay_info_dt'),
             br(), br(),
             DT::DTOutput('feature_info_dt'),
             br(), br(),
             DT::DTOutput('feature_levels_dt')),
    tabPanel("tSNE",
               mainPanel(width = "100%",
                 plotOutput("expr_tsne_plot") %>% withSpinner(color="#0dc5c1"),
                 br(),
                 div(style="display:inline-block",
                     downloadButton("dl_expr_tsne_png", label = "Download PNG"))))
  )
) 

############### server ############### 
server <- function(input, output) {
  
  # icon
  output$icon <- renderImage(list(src = "../hexagons/ccle.png",
                                  height = "86px", width = "72px"), 
                             deleteFile = F)
  
  # samples table
  output$sample_info_dt <- DT::renderDT(
    DT::datatable(as.data.frame(colData(ccledat$mae)),
                  escape = F, rownames = F,
                  editable = "cell",
                  extensions = 'Buttons',
                  options = list(dom = 'Bfrtip',
                                 buttons = list('pageLength',
                                                list(extend = 'copy',
                                                     title = NULL),
                                                list(extend = 'print',
                                                     title = "CCLE Sample Info"),
                                                list(extend = 'csv',
                                                     filename = paste0("CCLE_sample_info.", Sys.Date()),
                                                     title = NULL),
                                                list(extend = 'excel',
                                                     filename = paste0("CCLE_sample_info.", Sys.Date()),
                                                     title = NULL),
                                                list(extend = 'pdf',
                                                     filename = paste0("CCLE_sample_info.", Sys.Date()),
                                                     title = "CCLE Sample Info")),
                                 pagelength = 10,
                                 lengthMenu = list(c(10, 25, 100, -1),
                                                   c('10', '25', '100','All'))))
  )
  
  # assay summary table
  output$assay_info_dt <- DT::renderDT(
    DT::datatable(assay_info, 
                  escape = F, 
                  rownames = F,
                  selection = list(mode = "single", selected = 2),
                  options = list(dom = 't'))
  )
  
  # get assay type
  get_assay_type <- reactive({
    return(names(ccledat$mae)[input$assay_info_dt_rows_selected])
  })
  
  # features table
  output$feature_info_dt <- DT::renderDT({
    assay_type <- get_assay_type()
    if (is_empty(assay_type)) {
      return()
    } else {
      DT::datatable(ccledat$fdat[[assay_type]], rownames = F,
                    selection = list(mode = "single", selected = 1),
                    extensions = 'Buttons',
                    options = list(dom = 'Bfrtip',
                                   buttons = list('pageLength'),
                                   pagelength = 5,
                                   lengthMenu = list(c(5, 10, 20),
                                                     c('5', '10', '20'))))
    }
  })
  
  # get feature
  get_feature <- reactive({
    assay_type <- get_assay_type()
    if (is_empty(assay_type)) {
      return()
    } else {
      feat <- ccledat$fdat[[assay_type]][input$feature_info_dt_rows_selected,1]
      return(as.character(feat))
    }
  })
  
  # format feature name
  format_feature <- reactive({
    assay_type <- get_assay_type()
    feature <- get_feature()
    if (is_empty(feature)) {
      return()
    }
    if (assay_type == "RNAseq") {
      feature <- ccledat$fdat[[assay_type]][feature, "Gene_Symbol"]
    }
    feature <- gsub(" |/|\\(|\\)", "_", feature)
    gsub("_{2,}", "_", feature)
  })
  
  # get levels for selected feature
  get_feature_levels <- reactive({
    assay_type <- get_assay_type()
    feature <- get_feature()
    if (is_empty(assay_type) | is_empty(feature)) {
      return()
    }
    levelsdat <- tsne_dims[tsne_dims$assay_type == assay_type, c("V1", "V2", "Name", "Group")]
    colnames(levelsdat) <- c("Dim1","Dim2","Sample Name","Group")
    levelsdat[,"Feature Level"] <- as.numeric(assay(ccledat$mae[feature, ccledat$cc, assay_type]))
    levelsdat[,"Z-Score"] <- scale(levelsdat$'Feature Level')[,1]
    return(levelsdat)
  })
  
  # render table of feature levels
  output$feature_levels_dt <- DT::renderDT({
    assay_type <- get_assay_type()
    fname <- format_feature()
    levelsdat <- get_feature_levels()
    if (is_empty(levelsdat)) {
      return()
    } else {
      levelsdat[,c(1:2,5:6)] <- apply(levelsdat[,c(1:2,5:6)], 2, function(x) round(x, 3))
      tmptitle <- paste(c("CCLE", assay_type, fname, "levels"), collapse = " ")
      tmpfilepre <- paste(c("CCLE", assay_type, fname, "levels"), collapse = "_")
      tmpfilepre <- paste0(tmpfilepre, ".", Sys.Date())
      DT::datatable(levelsdat, 
                    rownames = F,
                    selection = "none",
                    extensions = 'Buttons',
                    options = list(dom = 'Bfrtip',
                                   buttons = list('pageLength',
                                                  list(extend = 'csv',
                                                       filename = tmpfilepre,
                                                       title = NULL),
                                                  list(extend = 'excel',
                                                       filename = tmpfilepre,
                                                       title = NULL)),
                                   pagelength = 10,
                                   lengthMenu = list(c(10, 25, -1),
                                                     c('10', '25', 'All'))))
    }
  })
  
  # plot feature levels
  plot_feature_levels <- reactive({
    assay_type <- get_assay_type()
    feature <- get_feature()
    fname <- format_feature()
    plot_types <- input$plot_type
    plotdat <- get_feature_levels()
    if (is_empty(assay_type) | is_empty(fname) | is_empty(plotdat)) {
      return()
    }
    plotdat_summary <- plotdat %>% 
      dplyr::group_by(Group) %>%
      dplyr::summarise(mean = mean(`Z-Score`)) %>%
      dplyr::arrange(-mean)
    plotdat$Group <- factor(plotdat$Group, levels = plotdat_summary$Group)
    a <- ggplot(plotdat, aes(x = Group, y = `Z-Score`)) +
      theme_bw(base_size = 18) + 
      theme(panel.grid = element_blank(),
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
    if (assay_type == "RNAseq") {
      a <- a + ylab(paste0(fname," Expression [TPM]"))
    } else {
      a <- a + ylab(paste0(feature, " Level [AU]"))
    }
    if ("dotplot" %in% plot_types) {
      a <- a + geom_point(aes(color = Group),
                          alpha = .5,
                          position = position_dodge2(width = .5))
    }
    if ("violin" %in% plot_types) {
      a <- a + geom_violin(aes(fill = Group))
    }
    if ("boxplot" %in% plot_types) {
      if ("dotplot" %in% plot_types & "violin" %in% plot_types) {
        a <- a + geom_boxplot(aes(fill = Group), width = 0.2, outlier.shape = NA)
      } else if ("dotplot" %in% plot_types) {
        a <- a + geom_boxplot(aes(fill = Group), width = 0.8, outlier.shape = NA)
      } else if ("violin" %in% plot_types) {
        a <- a + geom_boxplot(aes(fill = Group), width = 0.2)
      } else {
        a <- a + geom_boxplot(aes(fill = Group), width = 0.8)
      }
    }
    return(a)
  })
  
  # render feature levels plot
  output$feature_levels_plot <- renderPlot({
    plot_feature_levels()
  })
  
  # download feature levels plot as png
  output$dl_feature_levels_plot_png <- downloadHandler(
    filename = function() {
      paste0("CCLE_", get_assay_type(), "_",
             format_feature(), "_levels_distribution.", Sys.Date(), ".png")
    },
    content = function(file) {
      plot <- plot_feature_levels()
      ggsave(plot = plot, filename = file, height = 5, width = 10)
    }
  )
  
  # download feature levels plot as ppt
  output$dl_feature_levels_plot_ppt <- downloadHandler(
    filename = function() {
      paste0("CCLE_", get_assay_type(), "_",
             format_feature(), "_levels_distribution.", Sys.Date(), ".pptx")
    },
    content = function(file) {
      file_pptx <- tempfile(fileext = ".pptx")
      gen_pptx(plot_feature_levels(), file_pptx)
      file.rename(from = file_pptx, to = file)
    }
  )
  
  # plot tSNE of pheno dat
  output$pheno_tsne_plot <- renderPlot({
    palette <- distinctColorPalette(length(unique(tsne_dims$Group)))
    ggplot(tsne_dims, aes(x = V1, y = V2)) +
      geom_point(aes(fill = Group), col = "black", shape = 21, stroke = .2) +
      facet_grid(cols = vars(assay_type)) +
      scale_fill_manual(values = palette, name = "Sample Type") +
      theme_bw(base_size = 20) +
      theme(panel.grid = element_blank(),
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.title = element_blank())
  })
  
  # plot tSNE of expr dat
  get_expr_tsne_plot <- reactive({
    assay_type <- get_assay_type()
    feature <- get_feature()
    if (is_empty(assay_type) | is_empty(feature)) {
      return()
    }
    plotdat <- tsne_dims
    plotdat$expr <- scale(as.numeric(assay(ccledat$mae[feature,ccledat$cc,assay_type])))[,1]
    if (any(plotdat$expr > 3)) {
      plotdat[plotdat$expr > 3,]$expr <- 3
    }
    if (any(plotdat$expr < -3)) {
      plotdat[plotdat$expr < -3,]$expr <- -3
    }
    colnames(plotdat)[c(1,2)] <- c("dim1","dim2")
    ggplot(plotdat, aes(x = dim1, y = dim2, label = Name)) +
      geom_point(aes(fill = expr), col = "black", shape = 21, stroke = .2) +
      facet_grid(cols = vars(assay_type)) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           name = paste0(feature, " expr\n[z-score]"), limits = c(-3,3)) +
      theme_bw(base_size = 20) +
      theme(panel.grid = element_blank(),
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.title = element_blank())
  })
  
  # render tSNE of expr dat
  output$expr_tsne_plot <- renderPlot({
    get_expr_tsne_plot()
  })

  # download expression tSNE as png
  output$dl_expr_tsne_png <- downloadHandler(
    filename = function() {
      paste0("CCLE_tSNE_", get_feature(), "_levels.", Sys.Date(), ".png")
    },
    content = function(file) {
      plot <- get_expr_tsne_plot()
      ggsave(plot = plot, filename = file, height = 4, width = 18)
    }
  )
}

############### run ############### 
shinyApp(ui = ui, server = server) 
