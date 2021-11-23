library(shiny)
library(tidyverse)
library(dplyr)
library(shinyjs)
library(DT)
library(shinythemes)
library(shinyWidgets)
library(shinycssloaders)
library(org.Hs.eg.db)
library(sortable)
library(fgsea)
library(GEOquery)
library(colourpicker)
library(rvg)
library(officer)
options(stringsAsFactors = F)

################################### Setup #####################################

# set reactive vals
reactvals <- reactiveValues(gse = NULL,
                            es = NULL)

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

##################################### UI ######################################
ui <- fluidPage(
  
  # title
  title = "GEO", 
  theme = shinytheme("cosmo"),
  titlePanel(tags$h2(tags$a(
    imageOutput("icon", inline = TRUE),
    href="http://10.197.211.94:80"), "GEO")),
  
  # change color for column selection
  tags$style(HTML('table.dataTable tr.selected td{background-color: pink !important;}')),
  tags$style(HTML('table.dataTable td.selected {background-color: #3388ff88 !important;}')),
  
  
  sidebarLayout(
    
    #### Sidebar ####
    sidebarPanel(width = 2,
                 textInput("gse", "GSE Acession Number"),
                 div(actionButton("go_button", "GO", style = "font-size:14px;height:40px;", width = "100%")),
                 div(tags$small("Enter a GEO series accession ID (e.g. GSE38463) and press go.")),
                 hr(),
                 uiOutput("geo_link"),
                 hr(),
                 div(style = "font-size:13px;", uiOutput("gpl_choice")),
                 tags$small("Some data series may contain multiple different platform types,
                                 e.g. different microarrays used, select the platform 
                                 you would like to analyse (see GEO page for details)"),
                 hr(),
                 uiOutput("upload_text")
    ),
    
    #### Main ####
    mainPanel(width = 10, 
              tabsetPanel(
                tabPanel("Experiment Info",
                         htmlOutput("abstract_title"),
                         htmlOutput("abstract_text")
                ),
                tabPanel("Expression", align = "center",
                         uiOutput("distro_plot_ele"),
                         uiOutput("expr_plot_ui") %>%
                           withSpinner(color="pink"),
                         uiOutput("color_choices"),
                         div(style="display:inline-block;",
                             actionButton("log_expr", "Log",
                                          style = "font-size:12px;height:30px;padding:5px;"),
                             actionButton("linearize", "Linearize",
                                          style = "font-size:12px;height:30px;padding:5px;"),
                             actionButton("normalize", "Quantile normalize",
                                          style = "font-size:12px;height:30px;padding:5px;")),
                         br(),
                         div(style="display:inline-block;",
                             downloadButton("dl_expr_plot_ppt", label = "PPT",
                                            style = "font-size:12px;height:30px;padding:5px;"),
                             downloadButton("dl_expr_plot_png", label = "PNG",
                                            style = "font-size:12px;height:30px;padding:5px;")),
                         br(), br(),
                         tags$h3("Expression summary:"),
                         div(DT::dataTableOutput("expr_summary_dt"),
                             style = "font-size:90%"),
                         div(style="display:inline-block;",
                             actionButton("rank_toggle", "Re-order groups",
                                          style = "font-size:12px;height:30px;padding:5px;")),
                         div(style = "font-size:12px", uiOutput("group_ranks")),
                         br(), br(),
                         tags$h3("Feature info:"),
                         div(DT::dataTableOutput("fdat_dt"),
                             style = "font-size:90%"),
                         div(style="display:inline-block;",
                             checkboxInput("fdat_dt_select", "Add/remove all current genes", value = T),
                             h5(textOutput("n_features"))),
                         tags$h3("Sample info:"),
                         tags$h5("Select rows to change samples included and click columns names in footer to choose columns to group by."),
                         div(DT::dataTableOutput("pdat_dt"),
                             style = "font-size:90%"),
                         div(style="display:inline-block;",
                             checkboxInput("pdat_dt_select", "Add/remove all current samples"),
                             h5(textOutput("n_samples")))
                )#,
                # tabPanel("Clustering",
                #          radioGroupButtons("clust_plot_type", "Plot Type:",
                #                            choices = list("PCA" = "PCA",
                #                                           "Heatmap" = "heatmap"),
                #                            justified = T,
                #                            size = 'xs',
                #                            width = '25%',
                #                            status = "default",
                #                            selected = "pca"),
                #          uiOutput("clust_plot_ui") %>% 
                #            withSpinner(color="pink"),
                #          div(style="display:inline-block;",
                #              downloadButton("dl_clust_ppt", label = "PPT",
                #                             style = "font-size:12px;height:30px;padding:5px;"),
                #              downloadButton("dl_clust_png", label = "PNG",
                #                             style = "font-size:12px;height:30px;padding:5px;")),
                #          DT::dataTableOutput("groups_dt"),
                #          rank_list(text = "Drag groups to order",
                #            labels = list("one","two","three"),
                #            input_id = "rank_groups")
                # ),
                # tabPanel("DE",
                #          uiOutput("choose_pc1"),
                #          uiOutput("choose_pc2"),
                #          plotlyOutput("pca"),
                #          verbatimTextOutput("pca_hover")
                # ),
                # tabPanel("GSEA",
                #          uiOutput("choose_numer"),
                #          uiOutput("choose_denom"),
                #          selectInput("gs_type", "Gene set input type:",
                #                      c("MSigDB: H - hallmark" = "H",
                #                        "MSigDB: C1 - positional" = "c1",
                #                        "MSigDB: C2 - curated" = "c2",
                #                        "MSigDB: C3 - motif" = "c3",
                #                        "MSigDB: C4 - computational" = "c4",
                #                        "MsigDB: C5 - GO" = "c5",
                #                        "MSigDB: C6 - oncogenic" = "c6",
                #                        "MSigDB: C7 - immunological" = "c7",
                #                        "CUSTOM: GMT file" = "GMT")),
                #          uiOutput("gmt_upload"),
                #          uiOutput("choose_id_type"),
                #          selectInput("org", "Organism:",
                #                      c("Homo sapiens" = "Hs",
                #                        "Mus musculus" = "Mm")),
                #          plotOutput("gsea_plot"))
              )
    )
  )
)


################################### Server ####################################
server <- function(input, output) {
  
  # icon
  output$icon <- renderImage(list(src = "../hexagons/geo.png",
                                  height = "90px", width = "85px"), 
                             deleteFile = F)
  
  # get GSE data
  get_gse <- observeEvent(input$go_button, {
    # reactvals$gse <- NULL
    gse <- input$gse
    if (is.null(gse) | gse == "") return(NULL)
    withProgress(message = 'Fetching data from GEO', value = .3, { 
      gse <- GEOquery::getGEO(gse)
      reactvals$gse <- gse
    })
    reactvals$es <- gse[[1]]
    output$geo_link <- renderUI(tags$a(href=paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", input$gse),
                                       "Link to GEO record"))
  })
  
  # update GPL input choices
  output$gpl_choice <- renderUI({
    gse <- reactvals$gse
    if (is_empty(gse)) return(NULL)
    gpl_list <- unname(sapply(gse, annotation))
    names(gpl_list) <- unname(sapply(gpl_list, function(x) {
      gpl <- getGEO(x)
      paste0(x, ": ", gpl@header$title)
    }))
    selectInput("gpl_select", "Platform", gpl_list)
  })
  
  # select ES from list of platforms
  observeEvent(input$gpl_select, {
    gse <- reactvals$gse
    plat <- lapply(gse, function(x) x@annotation)
    if (is_empty(input$gpl_select)) return(NULL)
    reactvals$es <- gse[[grep(input$gpl_select, plat)]]
  })
  
  # experiment info output
  output$abstract_title <- renderUI({
    es <- reactvals$es
    if (is_empty(es)) {
      return(tags$h5("Enter GEO accession number in left panel"))
    } else {
      return(tags$h4("Abstract"))
    }
  })
  output$abstract_text <- renderUI({
    es <- reactvals$es
    print(es)
    if (is_empty(es)) return(NULL)
    tags$h5(experimentData(es)@abstract)
  })
  
  # pheno table
  output$pdat_dt <- DT::renderDataTable({
    es <- reactvals$es
    if (is_empty(es)) return(NULL)
    pdat <- pData(es)
    useless_cols <- apply(pdat, 2, function(x) length(unique(x)) == 1)
    pdat <- pdat[,!useless_cols]
    useless_cols <- colnames(pdat) %in% c("supplementary_file","data_row_count")
    pdat <- pdat[,!useless_cols]
    pData(reactvals$es) <- pdat
    selected <- grep("condition", colnames(pdat), ignore.case = T)-1
    selected <- ifelse(is_empty(selected), 1, selected[[1]])
    DT::datatable(
      data = pdat, 
      rownames = F,
      colnames = Hmisc::capitalize(gsub("[_\\.]", " ", colnames(pdat))),
      selection = list(mode = 'multiple', target = "row+column", 
                       selected = list(rows = 2:nrow(pdat), cols = selected)),
      # options = list(columnDefs = list(list(
      #   processing = F,
      #   targets = 0:(ncol(pdat)-1),
      #   render = JS(
      #     "function(data, type, row, meta) {",
      #     "return type === 'display' && data.length > 20 ?",
      #     "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
      #     "}")
      # ))),
      # callback = JS('table.page(3).draw(false);')
    )
  })
  
  # handle row selection/sample filtering
  pdat_dt_proxy <- DT::dataTableProxy("pdat_dt")
  observeEvent(input$pdat_dt_select, {
    if (isTRUE(input$pdat_dt_select)) {
      selected <- input$pdat_dt_rows_selected
      current <- input$pdat_dt_rows_all
      combined <- unique(c(selected, current))
      nsamples <- length(combined)
      DT::selectRows(pdat_dt_proxy, combined)
    } else {
      selected <- input$pdat_dt_rows_selected
      current <- input$pdat_dt_rows_all
      filtered <- selected[!selected %in% current]
      nsamples <- length(filtered)
      DT::selectRows(pdat_dt_proxy, NULL)
      DT::selectRows(pdat_dt_proxy, filtered)
    }
    output$n_samples <- renderText(paste0("Total selected samples: ", nsamples))
  })
  observeEvent(input$pdat_dt_rows_selected, {
    nsamples <- ncol(reactvals$es[, input$pdat_dt_rows_selected])
    output$n_samples <- renderText(paste0("Total selected samples: ", nsamples))
  })
  
  # feature info
  output$fdat_dt <- DT::renderDataTable({
    es <- reactvals$es
    if (is.null(es)) return(NULL)
    fdat <- fData(es)
    useless_cols <- apply(fdat, 2, function(x) length(unique(x)) == 1)
    fdat <- fdat[,!useless_cols]
    selected <- grep("gene|symbol", colnames(fdat), ignore.case = T)-1
    selected <- ifelse(is_empty(selected), 1, selected)
    fData(reactvals$es) <- fdat
    DT::datatable(
      data = fdat, 
      rownames = F,
      colnames = Hmisc::capitalize(gsub("[_\\.]", " ", colnames(fdat))),
      selection = list(mode = 'multiple', target = "row+column", 
                       selected = list(rows = c(1:nrow(fdat)), cols = selected)),
      options = list(columnDefs = list(list(
        targets = 0:(ncol(fdat)-1),
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 20 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
          "}")
      ))),
      callback = JS('table.page(3).draw(false);')
      )
  })
  
  # handle row selection from features table
  fdat_dt_proxy <- DT::dataTableProxy("fdat_dt")
  observeEvent(input$fdat_dt_select, {
    if (isTRUE(input$fdat_dt_select)) {
      selected <- input$fdat_dt_rows_selected
      current <- input$fdat_dt_rows_all
      combined <- unique(c(selected, current))
      nsamples <- length(combined)
      DT::selectRows(fdat_dt_proxy, combined)
    } else {
      selected <- input$fdat_dt_rows_selected
      current <- input$fdat_dt_rows_all
      filtered <- selected[!selected %in% current]
      nsamples <- length(filtered)
      DT::selectRows(fdat_dt_proxy, NULL)
      DT::selectRows(fdat_dt_proxy, filtered)
    }
    output$n_features <- renderText(paste0("Total selected samples: ", nfeatures))
  })
  observeEvent(input$fdat_dt_rows_selected, {
    nfeatures <- nrow(reactvals$es[input$fdat_dt_rows_selected, ])
    output$n_features <- renderText(paste0("Total selected samples: ", nfeatures))
  })
  
  # generate expr summary
  get_exprdat <- reactive({
    es <- reactvals$es
    if (is_empty(es)) return(NULL)
    fdat_col_idx <- input$fdat_dt_columns_selected+1
    if (length(fdat_col_idx) > 1) {
      fData(es)$gene_name <- Reduce(paste, fData(es)[, fdat_col_idx])
    } else {
      fData(es)$gene_name <- fData(es)[, fdat_col_idx]
    }
    pdat_col_idx <- input$pdat_dt_columns_selected+1
    if (length(pdat_col_idx) > 1) {
      pData(es)$group_name <- Reduce(paste, pData(es)[, pdat_col_idx])
    } else {
      pData(es)$group_name <- pData(es)[, pdat_col_idx]
    }
    samples_idx <- input$pdat_dt_rows_selected
    genes_idx <- input$fdat_dt_rows_selected
    es <- es[genes_idx, samples_idx]
    if (nrow(es) > 5000) {
      warning("Greater than 5000 features selected, just using subsample of 5000")
      es <- es[sample(1:nrow(es), 5000, replace = F), ]
    }
    if (length(unique(pData(es)$group_name)) > 20) {
      warning("Greater than 20 sample groups selected, just using subsample of 20")
      es <- es[, sample(1:ncol(es), 20, replace = F)]
    }
    exprdat <- exprs(es) %>%
      as.data.frame() %>%
      mutate(gene_name = fData(es)$gene_name) %>%
      gather("sample_id", "expr", -ncol(.)) 
    exprdat$group_name <- pData(es)[match(exprdat$sample_id, colnames(es)),]$group_name
    reactvals$exprdat <- exprdat
    reactvals$group_levels <- as.character(unique(exprdat$group_name))
    return(exprdat)
  })
  
  # summarise expression by group
  get_expr_summary <- reactive({  
    exprdat <- get_exprdat()
    exprdat <- reactvals$exprdat
    if (is_empty(exprdat)) return(NULL)
    exprdat %>%
      group_by(gene_name, group_name) %>%
      summarise(n = n(),
                mean = round(mean(expr, na.rm = T), 2),
                sd = round(sd(expr, na.rm = T), 2))
  })
  
  # output expr summary dt
  output$expr_summary_dt <- DT::renderDataTable({
    exprsum <- get_expr_summary()
    if (is_empty(exprsum)) return(NULL)
    DT::datatable(exprsum, 
                  rownames = F,
                  colnames = Hmisc::capitalize(gsub("_", " ", colnames(exprsum))),
                  editable = 'cell',
                  selection = list(target = 'row', mode = 'multiple', selected = NULL),
                  caption = "Select columns or edit cells to change ordering and labels")
  })
  
  # handle label editting/grouo highlighting
  summary_dt_proxy <- DT::dataTableProxy("expr_summary_dt")
  observeEvent(input$expr_summary_dt_cell_edit, {
    exprsum <- get_expr_summary()
    edited <- editData(exprsum, input$expr_summary_dt_cell_edit, 
                       proxy = 'summary_dt_proxy', rownames = F)
    exprdat <- reactvals$exprdat
    idx <- match(exprdat$group_name, exprsum$group_name)
    exprdat$group_name <- edited$group_name[idx]
    reactvals$exprdat <- exprdat
    levels <- as.character(reactvals$group_levels)
    idx <- match(levels, unique(exprsum$group_name))
    levels <- as.character(unique(edited$group_name)[idx])
    reactvals$group_levels <- levels
  })
  observeEvent(input$expr_summary_dt_rows_selected, {
    exprsum <- get_expr_summary()
    if (is_empty(input$expr_summary_dt_rows_selected)) {
      reactvals$selected_groups <- NULL
    } else {
      reactvals$selected_groups <- exprsum[input$expr_summary_dt_rows_selected,]$group_name
    }
  })
  observeEvent(input$group_rankings, {
    reactvals$group_levels <- as.character(unique(input$group_rankings))
  })
  
  # log/normalise data
  observeEvent(input$log_expr, {
    es <- reactvals$es
    if (is_empty(es) | input$log_expr == 0) return(NULL)
    if ((input$log_expr %% 2) != 0) {
      exprs(es) <- log2(exprs(es) + 0.01)
      reactvals$es <- es
    } else {
      exprs(es) <- 2^exprs(es)
      reactvals$es <- es
    }
  })
  observeEvent(input$linearize, {
    es <- reactvals$es
    if (is_empty(es) | input$linearize == 0) return(NULL)
    if ((input$linearize %% 2) != 0) {
      exprs(es) <- 2^exprs(es)
      reactvals$es <- es
    } else {
      exprs(es) <- log2(exprs(es) + 0.01)
      reactvals$es <- es
    }
  })
  observeEvent(input$normalize, {
    es <- reactvals$es
    if (is_empty(es) | input$log_expr == 0) return(NULL)
    if ((input$log_expr %% 2) != 0) {
      exprs(es) <- preprocessCore::normalize.quantiles(exprs(es))
      reactvals$es <- es
    } else return(NULL)
  })
  
  # group ranking
  output$group_ranks <- renderUI({
    if ((input$rank_toggle %% 2) != 0) {
      rank_list(text = "Drag groups to order",
                input_id = "group_rankings",
                labels = unique(reactvals$group_levels))
    } else return(NULL)
  })
  
  # group ranking
  output$distro_plot_ele <- renderUI({
    exprdat <- reactvals$exprdat
    if (length(unique(exprdat$gene_name)) == 1) {
      checkboxGroupButtons("plot_type", "Plot Element:",
                           choices = list("Violin" = "violin",
                                          "Boxplot" = "boxplot",
                                          "Dotplot" = "dotplot",
                                          "Range" = "range"),
                           justified = T,
                           size = 'xs', width = '25%',
                           selected = c("violin","boxplot","dotplot"))
    } else {
      return(NULL)
    }
  })
    
  # color elements
  output$color_choices <- renderUI({
    exprdat <- reactvals$exprdat
    selected_groups <- reactvals$selected_groups
    if (length(unique(exprdat$gene_name)) == 1 & !is_empty(selected_groups)) {
        div(style="display:inline-block;width:10%;", 
            colourInput("col_sel","Color#1", value = "#999999",
                         allowTransparent = TRUE),
            colourInput("col_unsel","Color#2", value = "#3498dbff",
                        allowTransparent = TRUE))
    } else {
      return(NULL)
    }
  })
  
  # generate distro plot
  plot_distribution <- reactive({
    exprdat <- reactvals$exprdat
    if (is_empty(exprdat)) return(NULL)
    print(reactvals$selected_groups)
    exprdat$group_name <- factor(exprdat$group_name, levels = unique(reactvals$group_levels))
    plot_types <- input$plot_type
    if (is_empty(reactvals$selected_groups)) {
      exprdat$group_col <- exprdat$group_name
    } else {
      exprdat$group_col <- ifelse(exprdat$group_name %in% reactvals$selected_groups, "Unselected", "Selected")
      exprdat$group_col <- factor(exprdat$group_col, levels = c("Selected", "Unselected"))
    }
    p <- ggplot(exprdat, aes(x = group_name, y = expr)) +
      theme_bw(base_size = 18) + 
      ylab(paste0(unique(exprdat$gene_symbol)," Expression [AU]")) +
      theme(panel.grid = element_blank(),
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
    if (!is_empty(reactvals$selected_groups)) {
      p <- p + scale_color_manual(values = c(input$col_sel, input$col_unsel)) +
        scale_fill_manual(values = c(input$col_sel, input$col_unsel))
    }
    if ("dotplot" %in% plot_types) {
      p <- p + geom_point(aes(color = group_col),
                          alpha = .5,
                          position = position_dodge2(width = .5))
    }
    if ("violin" %in% plot_types) {
      p <- p + geom_violin(aes(fill = group_col), alpha = .75)
    }
    if ("range" %in% plot_types) {
      p <- p + stat_summary(aes(color = group_col), 
                            geom = "pointrange", alpha = 1, 
                            fun.data = "mean_cl_boot")
    }
    if ("boxplot" %in% plot_types) {
      if ("dotplot" %in% plot_types & "violin" %in% plot_types) {
        p <- p + geom_boxplot(aes(fill = group_col), width = 0.2, outlier.shape = NA)
      } else if ("dotplot" %in% plot_types) {
        p <- p + geom_boxplot(aes(fill = group_col), width = 0.8, outlier.shape = NA)
      } else if ("violin" %in% plot_types) {
        p <- p + geom_boxplot(aes(fill = group_col), width = 0.2)
      } else {
        p <- p + geom_boxplot(aes(fill = group_col), width = 0.8)
      }
    }
    return(p)
  })
  
  # generate bivariate plot
  plot_bivariate <- reactive({
    exprdat <- reactvals$exprdat
    if (is_empty(exprdat)) return(NULL)
    if (is_empty(reactvals$selected_groups)) {
      exprdat$group_col <- exprdat$group_name
    } else {
      exprdat$group_col <- ifelse(exprdat$group_name %in% reactvals$selected_groups, "B", "A")
    }
    exprdat$group_name <- factor(exprdat$group_name, levels = unique(reactvals$group_levels))
    genes <- unique(exprdat$gene_name)
    exprdat <- spread(exprdat, "gene_name", "expr")
    genenames <-  colnames(exprdat)[c(4:5)]
    colnames(exprdat)[c(4:5)] <- c("geneA", "geneB")
    print(head(exprdat))
    p <- ggplot(exprdat, aes(x = geneA, y = geneB, color = group_col)) +
      theme_bw(base_size = 18) + 
      xlab(paste0(genenames[1]," Expression [AU]")) +
      ylab(paste0(genenames[2]," Expression [AU]")) +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(color = "black"),
            axis.text.x = element_text(color = "black")) +
      geom_point()
    return(p)
  })
  
  # generate bivariate plot
  plot_qc_box <- reactive({
    exprdat <- reactvals$exprdat
    if (is_empty(exprdat)) return(NULL)
    p <- ggplot(exprdat, aes(x = group_name, y = expr)) +
      geom_boxplot() +
      theme_bw(base_size = 18) + 
      ylab("Expression [AU]") +
      theme(panel.grid = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
    return(p)
  })
  
  # output distro plot
  sel_expr_plot <- reactive({
    exprsum <- get_expr_summary()
    if (is_empty(exprsum)) return(NULL)
    if (length(unique(exprsum$gene_name)) == 1) {
      plot_distribution()
    } else if (length(unique(exprsum$gene_name)) == 2) {
      plot_bivariate()
    } else {
      plot_qc_box()
    }
  })
  
  # output expr plot
  output$expr_plot <- renderPlot({
    sel_expr_plot()
  })
  get_dims <- reactive({
    exprsum <- get_expr_summary()
    ncond <- length(unique(exprsum$group_name))
    width <- min(70*ncond, 1000)
    list(w = width, h = 350)
  })
  output$expr_plot_ui <- renderUI({
    dims <- get_dims()
    plotOutput("expr_plot", height = dims$h, width = dims$w)
  })
  
  # download expr plot as png
  output$dl_expr_plot_png <- downloadHandler(
    filename = function() {
      exprsum <- get_expr_summary()
      if (length(unique(exprsum$gene_name)) > 2) {
        genes <- "overall"
      } else {
        genes <- paste0(unique(exprsum$gene_name), collapse = "_")
      }
      gseid <- input$gse
      paste0(gseid, "_", genes, "_expr.png")
    },
    content = function(file) {
      plot <- sel_expr_plot()
      dims <- get_dims()
      ggsave(plot = plot, filename = file, units = "mm",
             height = dims$h/3, width = dims$w/3)
    }
  )
  
  # download expr plot as ppt
  output$dl_expr_plot_ppt <- downloadHandler(
    filename = function() {
      exprsum <- get_expr_summary()
      if (length(unique(exprsum$gene_name)) > 2) {
        genes <- "overall"
      } else {
        genes <- paste0(unique(exprsum$gene_name), collapse = "_")
      }
      gseid <- input$gse
      paste0(gseid, "_", genes, "_expr.pptx")
    },
    content = function(file) {
      plot <- sel_expr_plot()
      dims <- get_dims()
      file_pptx <- tempfile(fileext = ".pptx")
      gen_pptx(plot, file_pptx,
               height = (dims$h/3)*0.039, 
               width = (dims$w/3)*0.039)
      file.rename(from = file_pptx, to = file)
    }
  )
  
  # # perform pca
  # get_pca_data <- reactive({
  #   gpl_data <- get_gpl_data()
  #   if (is.null(gpl_data)) {
  #     return()
  #   } else {
  #     pca_data <- t(na.omit(exprs(gpl_data)))
  #     pca_data <- as.data.frame(prcomp(pca_data, scale = T)$x)
  #     cond_levels <- get_cond_levels()[[1]] 
  #     pca_data$group <- as.character(cond_levels)
  #     return(pca_data)
  #   }
  # })
  
  # # input ui for choosing PCs
  # output$choose_pc1 <- renderUI({
  #   pca_data <- get_pca_data()
  #   if (is.null(pca_data)) {
  #     tags$h5("Enter GSE ID, GPL ID and groupings info before attempting PCA.")
  #   } else {
  #     selectInput('pc1', 'First PC:', colnames(pca_data), "PC1")
  #   }
  # })
  # output$choose_pc2 <- renderUI({
  #   pca_data <- get_pca_data()
  #   if (is.null(pca_data)) {
  #     return()
  #   } else {
  #     selectInput('pc2', 'Second PC:', colnames(pca_data), "PC2")
  #   }
  # })
  
  # # plot PCA
  # output$pca <- renderPlotly({
  #   pca_data <- get_pca_data()
  #   if (is.null(pca_data)) {
  #     return()
  #   } else {
  #     col_name <- "group"
  #     ggplot(data = pca_data, aes_string(x = input$pc1, y = input$pc2, 
  #                                        col = col_name)) +
  #       geom_point() +
  #       theme_bw(base_size = 18) +
  #       theme(legend.position = "bottom", 
  #             legend.title=element_blank(),
  #             panel.grid = element_blank()) +
  #       xlab(input$pc1) + ylab(input$pc2)
  #   }
  # })
  
  # # input ui for choosing contrast levels for GSEA
  # output$choose_numer <- renderUI({
  #   gpl_data <- get_gpl_data()
  #   if (is.null(gpl_data)) {
  #     tags$h5("Enter GSE ID, GPL ID and groupings info before attempting GSEA.")
  #   } else {
  #     cond_levels <- get_cond_levels()
  #     selectInput('numer', 'Numerator:', levels(cond_levels[[1]]))
  #   }
  # })
  # output$choose_denom <- renderUI({
  #   gpl_data <- get_gpl_data()
  #   if (is.null(gpl_data)) {
  #     return()
  #   } else {
  #     cond_levels <- get_cond_levels()[[1]]
  #     lvl2 <- ifelse(length(levels(cond_levels)) > 1,
  #                    levels(cond_levels)[2],
  #                    levels(cond_levels)[1])
  #     selectInput('denom', 'Denominator:', levels(cond_levels), lvl2)
  #   }
  # })
  # 
  # # file upload for custom
  # output$gmt_upload <- renderUI({
  #   gs_type <- input$gs_type
  #   if (gs_type == "GMT") {
  #     fileInput("gmt_file", "Choose GMT file:",
  #               accept=c('text/plain'))
  #   } else {
  #     return()
  #   }
  # })
  # 
  # # gmt id types
  # output$choose_id_type <- renderUI({
  #   gs_type <- input$gs_type
  #   if (gs_type == "GMT") {
  #     selectInput("gmt_id_type", "Choose ID type:",
  #                 keytypes(org.Hs.eg.db))
  #   } else {
  #     return()
  #   }
  # })
  # 
  # # match anno columns
  # output$match_anno1 <- renderUI({
  #   gs_type <- input$gs_type
  #   if (gs_type == "GMT") {
  #     selectInput("gmt_id_type", "Choose ID type:",
  #                 ketypes(org.Hs.eg.db))
  #   } else {
  #     return()
  #   }
  # })
  # 
  # # annotation data
  # anno <- reactive({
  #   gpl_data <- get_gpl_data()
  #   if (input$org == "Hs") {
  #     org_db <- "org.Hs.eg.db"
  #   } else if (input$org == "Mm") {
  #     org_db <- "org.Mm.eg.db"
  #   }
  #   bitr(fData(gpl_data)$'GENE SYMBOL', 
  #        fromType = "SYMBOL",
  #        toType = c("ENTREZID","ENSEMBL","REFSEQ",
  #                   "ACCNUM","UNIPROT"),
  #        OrgDb = org_db)
  # })
  # 
  # # getting gene set
  # get_gene_sets <- reactive({
  #   gs_type <- input$gs_type
  #   if (gs_type == "GMT") {
  #     gmt_file <- input$gmt_file
  #     if (is.null(gmt_file)) return()
  #     gene_set <- read_gmt(gmt_file$datapath)
  #     return(gene_set)
  #   } else {
  #     gene_set <- eval(parse(text = paste0(input$org, ".", gs_type)))
  #     return(gene_set)
  #   }
  # })
  # 
  # # temp
  # output$temp <- renderTable({
  #   gene_set <- get_gene_sets()
  #   if (is.null(gene_set)) return()
  #   summaryTable <- data.frame("Set name" = names(gene_set),
  #                              "N genes" = lengths(gene_set))
  #   summaryTable
  # })
  # 
  # # perform gsea test
  # gsea_test <- reactive({
  #   gpl_data <- get_gpl_data()
  #   if (is.null(gpl_data)) {
  #     return()
  #   } else {
  #     cond_levels <- get_cond_levels()
  #     gene_sets <- get_gene_sets()
  #     gene_sets <- term_to_gene(gene_sets)
  #     numer_idx <- which(cond_levels$cond_levels == input$numer)
  #     denom_idx <- which(cond_levels$cond_levels == input$denom)
  #     exprs_mat <- exprs(gpl_data)
  #     entrezids <- mapIds(org.Hs.eg.db,
  #                         keys = fData(gpl_data)[,'GENE SYMBOL'],
  #                         keytype = "SYMBOL",
  #                         column = "ENTREZID")
  #     rownames(exprs_mat) <- entrezids
  #     gm_values <- gm_ratios(exprs_mat, numer_idx, denom_idx)
  #     gm_values <- gm_values[order(-gm_values)]
  #     gsea_res <- GSEA(gm_values, TERM2GENE = gene_sets,
  #                      minGSSize = 10, pvalueCutoff = 0.1)
  #     return(gsea_res)
  #   }
  # })
  # 
  # # plot PCA
  # output$gsea_plot <- renderPlot({
  #   gsea_res <- gsea_test()
  #   if (is.null(gsea_res)) {
  #     return()
  #   } else {
  #     if (is.null(chosen_gene_set)) {
  #       chosen_gene_set <- gsea_res$ID[1]
  #     }
  #     gseaplot(gsea_res, geneSetID = chosen_gene_set)
  #   }
  # })
}

# Run the application 
shinyApp(ui = ui, server = server)
