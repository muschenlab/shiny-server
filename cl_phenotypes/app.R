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
library(reactlog)

#reactlog
options(shiny.reactlog=TRUE)

# initialise reactvals
reactvals <- reactiveValues(si = NULL)

# load CCLE data
ccledat <- readRDS("ccle_data/CCLE_MAE_processed.2019-10-27.rds")

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
    href="http://137.184.200.69:3838"), 
    "Cancer Cell Line Phenotypes")),
  
  tabsetPanel(
    tabPanel("Plot", align = "center",
             radioGroupButtons("plot_type", "Plot Type:",
                               choices = list("Distribution" = "distro",
                                              "tSNE" = "tsne"),
                               justified = T,
                               size = 'xs',
                               width = '25%',
                               status = "default",
                               selected = "distro"),
             uiOutput("distro_ele_ui"),
             uiOutput("fill_by_ui"),
             uiOutput("main_plot_ui") %>% withSpinner(color="#0dc5c1"),
             uiOutput("zscore_slider"),
             div(style="display:inline-block",
                 downloadButton("dl_main_plot_ppt", label = "PPT",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 downloadButton("dl_main_plot_png", label = "PNG",
                                style = "font-size:12px;height:30px;padding:5px;")),
             br(), br(),
             tags$h4("Summary Info:"),
             DT::DTOutput('summary_info_dt'),
             br(), br(),
             tags$h4("Assay Info:"),
             DT::DTOutput('assay_info_dt'),
             br(), br(),
             tags$h4("Feature Info:"),
             DT::DTOutput('feature_info_dt'),
             br(), br(),
             tags$h4("Sample Info:"),
             DT::DTOutput('sample_info_dt')
    )
  )
) 

############### server ############### 
server <- function(input, output, session) {
  
  # icon
  output$icon <- renderImage(list(src = "../assets/images/hex/cl_phenotypes.png",
                                  height = "90px", width = "85px"), 
                             deleteFile = F)
  
  # plot choice elements
  output$distro_ele_ui <- renderUI({
    if (input$plot_type == "distro") {
      checkboxGroupButtons("distro_elements", "Plot Elements:", 
                           choices = list("Violin" = "violin",
                                          "Boxplot" = "boxplot",
                                          "Dotplot" = "dotplot"),
                           justified = T, 
                           size = 'xs', width = '25%',
                           selected = "violin")
    } else if (input$plot_type == "tsne") {
      checkboxGroupButtons("cluster_ele", "Cluster by assay:", 
                           choices = list("RPPA" = "RPPA",
                                          "RNAseq" = "RNAseq",
                                          "Metabolites" = "Metabolites",
                                          "Chromatin" = "Chromatin",
                                          "Combined" = "Combined"),
                           justified = T, 
                           size = 'xs', width = '30%',
                           selected = c("RPPA","RNAseq","Metabolites","Chromatin","Combined"))
    } else return(NULL)
  })
  output$fill_by_ui <- renderUI({
    if (input$plot_type == "tsne") {
      radioGroupButtons("fill_by", "Fill by:",
                    choices = list("Cell type" = "ct",
                                   "Expression" = "expr"),
                    justified = T,
                    size = 'xs',
                    width = '30%',
                    status = "default",
                    selected = "ct")
    } else return(NULL)
  })
  
  # min/max
  output$zscore_slider <- renderUI({
    if (input$plot_type == "tsne") {
      sliderInput(inputId = "zscore_lims", label = "Upper limit", 
                  min = -10, max = 10, value = c(-3,3), step = 0.25)
    } else return(NULL)
  })
  
  # samples table
  output$sample_info_dt <- DT::renderDT(
    DT::datatable(as.data.frame(colData(ccledat$mae)),
                  escape = F, rownames = F,
                  selection = list(target = "column", selected = 3, mode = "single"),
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
                                                     title = NULL)),
                                 pagelength = 10,
                                 lengthMenu = list(c(10, 25, 100, -1),
                                                   c('10', '25', '100','All'))))
  )
  
  # add grouping variable
  observeEvent(input$sample_info_dt_columns_selected, {
    si <- colData(ccledat$mae)[ccledat$cc, ]
    si$group_name <- si[, (input$sample_info_dt_columns_selected+1)]
    reactvals$si <- si
  })
  
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
    print(assay_type)
    feature <- get_feature()
    if (is_empty(feature)) {
      return()
    }
    if (assay_type == "RNAseq") {
      print(ccledat$fdat[[assay_type]][feature, ])
      feature <- ccledat$fdat[[assay_type]][feature, "Gene_Symbol"]
    }
    feature <- gsub(" |/|\\(|\\)", "_", feature)
    gsub("_{2,}", "_", feature)
  })
  
  # get levels for selected feature
  get_feature_levels <- reactive({
    assay_type <- get_assay_type()
    feature <- get_feature()
    if (is_empty(assay_type) | is_empty(feature)) return(NULL)
    levelsdat <- tsne_dims[tsne_dims$assay_type == assay_type, c("V1", "V2", "Name")]
    colnames(levelsdat) <- c("dim1","dim2","sample_name")
    levelsdat$group_name <- reactvals$si$group_name
    levelsdat$feature_level <- as.numeric(assay(ccledat$mae[feature, ccledat$cc, assay_type]))
    levelsdat$zscore <- scale(levelsdat$feature_level)[,1]
    return(levelsdat)
  })
  
  # summary info
  get_summarydat <- reactive({
    assay_type <- get_assay_type()
    fname <- format_feature()
    levelsdat <- get_feature_levels()
    if (is_empty(levelsdat)) return(NULL)
    summarydat <- levelsdat %>% group_by(group_name) %>%
      summarise(mean_feature_level = round(mean(feature_level), 1),
                mean_zscore = round(mean(zscore), 3))
  })
  
  # render table of feature levels
  output$summary_info_dt <- DT::renderDT({
    summarydat <- get_summarydat()
    if (is_empty(summarydat)) return(NULL)
    fname <- format_feature()
    DT::datatable(summarydat, 
                  rownames = F,
                  selection = list(target = "row", mode = "multiple"),
                  extensions = 'Buttons',
                  colnames = Hmisc::capitalize(gsub("_", " ", colnames(summarydat))),
                  options = list(dom = 'Bfrtip',
                                 buttons = list('pageLength',
                                                list(extend = 'csv',
                                                     filename = paste0("CCLE_", fname, "_expr_summary"),
                                                     title = NULL),
                                                list(extend = 'excel',
                                                     filename = paste0("CCLE_", fname, "_expr_summary"),
                                                     title = paste0("CCLE ", fname, " expression summary"))),
                                 pagelength = 10,
                                 lengthMenu = list(c(10, 25, -1),
                                                   c('10', '25', 'All'))))
  })
  
  # observe group highlighting
  observeEvent(input$summary_info_dt_rows_selected, {
    if (length(input$summary_info_dt_rows_selected) == 0) {
      reactvals$sel_groups <- NULL
    } else {
      summarydat <- get_summarydat()
      reactvals$sel_groups <- summarydat[input$summary_info_dt_rows_selected, ]$group_name
    }
  })
  
  # plot feature levels
  plot_distro <- reactive({
    assay_type <- get_assay_type()
    feature <- get_feature()
    fname <- format_feature()
    plot_ele <- input$distro_elements
    plotdat <- get_feature_levels()
    if (is_empty(assay_type) | is_empty(fname) | is_empty(plotdat)) {
      return()
    }
    plotdat_summary <- plotdat %>% 
      dplyr::group_by(group_name) %>%
      dplyr::summarise(mean = mean(zscore)) %>%
      dplyr::arrange(-mean)
    plotdat$group_name <- factor(plotdat$group_name, levels = plotdat_summary$group_name)
    if (is_empty(reactvals$sel_groups)) {
      plotdat$color <- plotdat$group_name
    } else {
      plotdat$color <- ifelse(plotdat$group_name %in% reactvals$sel_groups, "B", "A")
    }
    print(head(plotdat))
    a <- ggplot(plotdat, aes(x = group_name, y = zscore)) +
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
    if ("dotplot" %in% plot_ele) {
      a <- a + geom_point(aes(color = color),
                          alpha = .5,
                          position = position_dodge2(width = .5))
    }
    if ("violin" %in% plot_ele) {
      a <- a + geom_violin(aes(fill = color))
    }
    if ("boxplot" %in% plot_ele) {
      if ("dotplot" %in% plot_ele & "violin" %in% plot_ele) {
        a <- a + geom_boxplot(aes(fill = color), width = 0.2, outlier.shape = NA)
      } else if ("dotplot" %in% plot_ele) {
        a <- a + geom_boxplot(aes(fill = color), width = 0.8, outlier.shape = NA)
      } else if ("violin" %in% plot_ele) {
        a <- a + geom_boxplot(aes(fill = color), width = 0.2)
      } else {
        a <- a + geom_boxplot(aes(fill = color), width = 0.8)
      }
    }
    if (!is_empty(reactvals$sel_groups)) {
      a <- a + scale_color_manual(values = c("grey","firebrick")) +
        scale_fill_manual(values = c("grey","firebrick"))
    } 
    return(a)
  })
  
  # plot tSNE of cell types
  plot_tsne_ct <- reactive({
    assay_type <- get_assay_type()
    plotdat <- tsne_dims
    plotdat$group_name <- reactvals$si$group_name
    plotdat <- plotdat[which(plotdat$assay_type %in% input$cluster_ele), ]
    sel_groups <- reactvals$sel_groups
    if (is_empty(sel_groups)) {
      plotdat$color <- plotdat$group_name
    } else {
      plotdat$color <- ifelse(plotdat$group_name %in% sel_groups, "B", "A")
    }
    colnames(plotdat)[c(1,2)] <- c("dim1","dim2")
    print(head(plotdat))
    p <- ggplot(plotdat, aes(x = dim1, y = dim2)) +
      geom_point(aes(fill = color), col = "black", shape = 21, stroke = .1, size = 1.8) +
      facet_grid(cols = vars(assay_type)) +
      theme_bw(base_size = 20, base_family = "Arial") +
      theme(panel.grid = element_blank(),
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.title = element_blank(),
            legend.title = element_blank(),
            legend.position = "none")
    if (!is_empty(sel_groups)) {
      p <- p + scale_fill_manual(values = c("grey", "firebrick")) 
    } 
    return(p)
  })
  
  # plot tSNE of expr dat
  plot_tsne <- reactive({
    assay_type <- get_assay_type()
    feature <- get_feature()
    if (is_empty(assay_type) | is_empty(feature)) {
      return()
    }
    plotdat <- tsne_dims
    plotdat$expr <- scale(as.numeric(assay(ccledat$mae[feature, ccledat$cc, assay_type])))[,1]
    plotdat <- plotdat[which(plotdat$assay_type %in% input$cluster_ele), ]
    if (any(plotdat$expr > input$zscore_lims[2])) {
      plotdat[plotdat$expr > input$zscore_lims[2],]$expr <- input$zscore_lims[2]
    }
    if (any(plotdat$expr < input$zscore_lims[1])) {
      plotdat[plotdat$expr < input$zscore_lims[1],]$expr <- input$zscore_lims[1]
    }
    colnames(plotdat)[c(1,2)] <- c("dim1","dim2")
    ggplot(plotdat, aes(x = dim1, y = dim2, label = Name)) +
      geom_point(aes(fill = expr), col = "black", shape = 21, stroke = .1, size = 1.8) +
      facet_grid(cols = vars(assay_type)) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           name = paste0(feature, " Expression [z-score]")) +
      guides(fill = guide_colourbar(barwidth = 8, barheight = 1.5)) +
      theme_bw(base_size = 20, base_family = "Arial") +
      theme(panel.grid = element_blank(),
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 16, vjust = 0.9),
            legend.position = "bottom")
  })
  
  # render main plot
  get_plotdat <- reactive({
    if (input$plot_type == "distro") {
      ngroups <- length(unique(reactvals$si$group_name))
      p <- plot_distro()
      return(list(p = p, h = 350, w = max(c(150, ngroups*30))))
    } else if (input$plot_type == "tsne") {
      nclust <- length(input$cluster_ele)
      if (input$fill_by == "expr") {
        p <- plot_tsne()
        return(list(p = p, h = 400, w = nclust*240))
      } else {
        p <- plot_tsne_ct()
        return(list(p = p, h = 320, w = nclust*240))
      }
    }
  })
  output$main_plot <- renderPlot({
    get_plotdat()$p
  })
  output$main_plot_ui <- renderUI({
    plotdat <- get_plotdat()
    plotOutput("main_plot", height = plotdat$h, width = plotdat$w)
  })
  
  # download feature levels plot as png
  output$dl_main_plot_png <- downloadHandler(
    filename = function() {
      paste0("CCLE_", get_assay_type(), "_",
             format_feature(), "_", 
             input$plot_type, ".png")
    },
    content = function(file) {
      plotdat <- get_plotdat()
      ggsave(plot = plotdat$p, filename = file,
             unit = "mm",
             height = plotdat$h / 3,
             width = plotdat$w / 3)
    }
  )
  
  # download feature levels plot as ppt
  output$dl_main_plot_ppt <- downloadHandler(
    filename = function() {
      paste0("CCLE_", get_assay_type(), "_",
             format_feature(), "_", 
             input$plot_type, ".pptx")
    },
    content = function(file) {
      file_pptx <- tempfile(fileext = ".pptx")
      plotdat <- get_plotdat()
      gen_pptx(plotdat$p, file_pptx,
               height = (plotdat$h/3)*0.039,
               width = (plotdat$w/3)*0.039)
      file.rename(from = file_pptx, to = file)
    }
  )
}

############### run ############### 
shinyApp(ui = ui, server = server) 
