library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(shinycssloaders)
library(colourpicker)
library(tidyverse)
library(DT)
library(rvg)
library(officer)
library(plotly)

################################### Setup #####################################

# load data
gene_summary_list <- list("B-cell vs solid" = readxl::read_excel("dep_data/bcell_vs_st_crispr_rnai.xlsx"),
                          "T-cell vs solid" = readxl::read_excel("dep_data/tcell_vs_st_crispr_rnai.xlsx"))

# set reactive vals
reactvals <- reactiveValues(gene_summary_df = NULL)

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
    title = "DEP", 
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
                     div(tags$small("Choose cell type comparison")),
                     div(style = "font-size:13px;", uiOutput("comparison_choice"))
        ),
        
        #### Main ####
        mainPanel(width = 10, 
                  tabsetPanel(
                      tabPanel("Genes",
                               tags$h5(paste0("")),
                               tags$h3("Gene-level dependency summary:"),
                               tags$h5("Select rows to plot gene dependency"),
                               div(DT::dataTableOutput("gene_summary_dt"),
                                   style = "font-size:90%")
                      ),
                      tabPanel("Drugs",
                               uiOutput("choose_pc1"),
                               uiOutput("choose_pc2"),
                               plotlyOutput("pca"),
                               verbatimTextOutput("pca_hover")
                      )
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
    
    # cell type comparison choice UI
    output$comparison_choice <- renderUI({
        selectInput('comparison', 'Comparison:', names(gene_summary_list), "B-cell vs solid")
    })
    
    # observe cell type comparison choice
    observeEvent(input$comparison, {
        if (is_empty(input$comparison)) return(NULL)
        reactvals$gene_summary_df <- gene_summary_list[[input$comparison]]
    })
    
    # feature info table for feature/row selection
    output$gene_summary_dt <- DT::renderDataTable({
        gene_summary_df <- reactvals$gene_summary_df
        if (is.null(gene_summary_df)) return(NULL)
        DT::datatable(
            data = gene_summary_df, 
            rownames = F,
            colnames = Hmisc::capitalize(gsub("[_\\.]", " ", colnames(gene_summary_df))),
            selection = list(mode = 'single', target = "row", 
                             selected = 1),
            options = list(columnDefs = list(list(
                targets = 0:(ncol(gene_summary_df)-1),
                render = JS(
                    "function(data, type, row, meta) {",
                    "return type === 'display' && data.length > 20 ?",
                    "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                    "}")
            ))),
            callback = JS('table.page(3).draw(false);')
        )
    })
    
    # # handle row selection from features table
    # fdat_rows_dt_proxy <- DT::dataTableProxy("fdat_rows_dt")
    # observeEvent(input$fdat_rows_dt_select, {
    #     if (isTRUE(input$fdat_rows_dt_select)) {
    #         selected <- input$fdat_rows_dt_rows_selected
    #         current <- input$fdat_rows_dt_rows_all
    #         combined <- unique(c(selected, current))
    #         nfeatures <- length(combined)
    #         DT::selectRows(fdat_rows_dt_proxy, combined)
    #     } else {
    #         selected <- input$fdat_rows_dt_rows_selected
    #         current <- input$fdat_rows_dt_rows_all
    #         filtered <- selected[!selected %in% current]
    #         nfeatures <- length(filtered)
    #         DT::selectRows(fdat_rows_dt_proxy, NULL)
    #         DT::selectRows(fdat_rows_dt_proxy, filtered)
    #     }
    #     output$n_features <- renderText(paste0("Total selected features: ", nfeatures))
    # })
    # observeEvent(input$fdat_rows_dt_rows_selected, {
    #     nfeatures <- nrow(reactvals$es[input$fdat_rows_dt_rows_selected, ])
    #     if (is_empty(nfeatures)) nfeatures <- 0
    #     output$n_features <- renderText(paste0("Total selected features: ", nfeatures))
    # })
    # 
    # # CRISPR data
    # get_crispr_dat <- reactive({
    #     genesdf <- reactvals$genesdf
    #     if (is_empty(genesdf)) return(NULL)
    #     
    #     return(crispr_dat)
    # })
    # 
    # # format es
    # get_es <- reactive({
    #     es <- reactvals$es
    #     useless_cols <- apply(pData(es), 2, function(x) length(unique(x)) == 1)
    #     useless_cols <- (useless_cols | colnames(pData(es)) %in% c("supplementary_file","data_row_count"))
    #     pData(es) <- pData(es)[,!useless_cols]
    #     useless_cols <- apply(fData(es), 2, function(x) length(unique(x)) == 1)
    #     fData(es) <- fData(es)[,!useless_cols]
    #     reactvals$es <- es
    #     return(es)
    # })
    # 
    # # generate expr summary
    # get_exprdat <- reactive({
    #     es <- reactvals$es
    #     if (is_empty(es)) return(NULL)
    #     fdat_col_idx <- input$fdat_cols_dt_columns_selected+1
    #     if (length(fdat_col_idx) < 1) {
    #         return(NULL)
    #     } else if (length(fdat_col_idx) > 1) {
    #         fData(es)$gene_name <- Reduce(paste, fData(es)[, fdat_col_idx])
    #     } else {
    #         fData(es)$gene_name <- fData(es)[, fdat_col_idx]
    #     }
    #     pdat_col_idx <- input$pdat_cols_dt_columns_selected+1
    #     if (length(pdat_col_idx) < 1) {
    #         return(NULL)
    #     } else if (length(pdat_col_idx) > 1) {
    #         pData(es)$group_name <- Reduce(paste, pData(es)[, pdat_col_idx])
    #     } else {
    #         pData(es)$group_name <- pData(es)[, pdat_col_idx]
    #     }
    #     samples_idx <- input$pdat_rows_dt_rows_selected
    #     genes_idx <- input$fdat_rows_dt_rows_selected
    #     if (length(samples_idx) < 1 | length(genes_idx) < 1) { return(NULL) }
    #     es <- es[genes_idx, samples_idx]
    #     if (nrow(es) > 5000) {
    #         warning("Greater than 5000 features selected, just using subsample of 5000")
    #         es <- es[sample(1:nrow(es), 5000, replace = F), ]
    #     }
    #     if (length(unique(pData(es)$group_name)) > 20) {
    #         warning("Greater than 20 sample groups selected, just using subsample of 20")
    #         es <- es[, sample(1:ncol(es), 20, replace = F)]
    #     }
    #     exprdat <- exprs(es) %>%
    #         as.data.frame() %>%
    #         mutate(gene_name = fData(es)$gene_name) %>%
    #         gather("sample_id", "expr", -ncol(.)) 
    #     exprdat$group_name <- pData(es)[match(exprdat$sample_id, colnames(es)),]$group_name
    #     reactvals$exprdat <- exprdat
    #     reactvals$group_levels <- as.character(unique(exprdat$group_name))
    #     return(exprdat)
    # })
    # 
    # # summarise expression by group
    # get_expr_summary <- reactive({  
    #     exprdat <- get_exprdat()
    #     if (is_null(exprdat)) return(NULL)
    #     exprdat %>%
    #         group_by(gene_name, group_name) %>%
    #         summarise(n = n(),
    #                   mean = round(mean(expr, na.rm = T), 2),
    #                   sd = round(sd(expr, na.rm = T), 2))
    # })
    # 
    # # output expr summary dt
    # output$expr_summary_dt <- DT::renderDataTable({
    #     exprsum <- get_expr_summary()
    #     if (is_empty(exprsum)) return(NULL)
    #     DT::datatable(exprsum, 
    #                   rownames = F,
    #                   colnames = Hmisc::capitalize(gsub("_", " ", colnames(exprsum))),
    #                   editable = 'cell',
    #                   selection = list(target = 'row', mode = 'multiple', selected = NULL),
    #                   caption = "Select columns or edit cells to change ordering and labels")
    # })
    # 
    # # handle label editting/group highlighting
    # summary_dt_proxy <- DT::dataTableProxy("expr_summary_dt")
    # observeEvent(input$expr_summary_dt_cell_edit, {
    #     exprsum <- get_expr_summary()
    #     edited <- editData(exprsum, input$expr_summary_dt_cell_edit, 
    #                        proxy = 'summary_dt_proxy', rownames = F)
    #     exprdat <- reactvals$exprdat
    #     idx <- match(exprdat$group_name, exprsum$group_name)
    #     exprdat$group_name <- edited$group_name[idx]
    #     reactvals$exprdat <- exprdat
    #     levels <- as.character(reactvals$group_levels)
    #     idx <- match(levels, unique(exprsum$group_name))
    #     levels <- as.character(unique(edited$group_name)[idx])
    #     reactvals$group_levels <- levels
    # })
    # observeEvent(input$expr_summary_dt_rows_selected, {
    #     exprsum <- get_expr_summary()
    #     if (is_empty(input$expr_summary_dt_rows_selected)) {
    #         reactvals$selected_groups <- NULL
    #     } else {
    #         reactvals$selected_groups <- exprsum[input$expr_summary_dt_rows_selected,]$group_name
    #     }
    # })
    # observeEvent(input$group_rankings, {
    #     reactvals$group_levels <- as.character(unique(input$group_rankings))
    # })
    # 
    # # log/normalise data
    # observeEvent(input$log_expr, {
    #     es <- reactvals$es
    #     if (is_empty(es) | input$log_expr == 0) return(NULL)
    #     if ((input$log_expr %% 2) != 0) {
    #         exprs(es) <- log2(exprs(es) + 0.01)
    #         reactvals$es <- es
    #     } else {
    #         exprs(es) <- 2^exprs(es)
    #         reactvals$es <- es
    #     }
    # })
    # observeEvent(input$linearize, {
    #     es <- reactvals$es
    #     if (is_empty(es) | input$linearize == 0) return(NULL)
    #     if ((input$linearize %% 2) != 0) {
    #         exprs(es) <- 2^exprs(es)
    #         reactvals$es <- es
    #     } else {
    #         exprs(es) <- log2(exprs(es) + 0.01)
    #         reactvals$es <- es
    #     }
    # })
    # observeEvent(input$normalize, {
    #     es <- reactvals$es
    #     if (is_empty(es) | input$log_expr == 0) return(NULL)
    #     if ((input$log_expr %% 2) != 0) {
    #         exprs(es) <- preprocessCore::normalize.quantiles(exprs(es))
    #         reactvals$es <- es
    #     } else return(NULL)
    # })
    # 
    # # group ranking
    # output$group_ranks <- renderUI({
    #     if ((input$rank_toggle %% 2) != 0) {
    #         rank_list(text = "Drag groups to order",
    #                   input_id = "group_rankings",
    #                   labels = unique(reactvals$group_levels))
    #     } else return(NULL)
    # })
    # 
    # # group ranking
    # output$distro_plot_ele <- renderUI({
    #     exprdat <- reactvals$exprdat
    #     if (is_empty(exprdat)) return(NULL) 
    #     if (length(unique(exprdat$gene_name)) == 1) {
    #         checkboxGroupButtons("plot_type", "Plot Element:",
    #                              choices = list("Violin" = "violin",
    #                                             "Boxplot" = "boxplot",
    #                                             "Dotplot" = "dotplot",
    #                                             "Range" = "range"),
    #                              justified = T,
    #                              size = 'xs', width = '25%',
    #                              selected = c("violin","boxplot","dotplot"))
    #     } else {
    #         return(NULL)
    #     }
    # })
    # 
    # # color elements
    # output$color_choice_1 <- renderUI({
    #     exprdat <- reactvals$exprdat
    #     selected_groups <- reactvals$selected_groups
    #     if (is_empty(selected_groups)) return(NULL)
    #     colourpicker::colourInput("col_sel", "Choose colors", value = "#88888888",
    #                               allowTransparent = TRUE)
    # })
    # output$color_choice_2 <- renderUI({
    #     exprdat <- reactvals$exprdat
    #     selected_groups <- reactvals$selected_groups
    #     if (is_empty(selected_groups)) return(NULL)
    #     colourpicker::colourInput("col_unsel", NULL, value = "#3498dbff",
    #                               allowTransparent = TRUE)
    # })
    # 
    # # generate distro plot
    # plot_distribution <- reactive({
    #     exprdat <- reactvals$exprdat
    #     if (is_empty(exprdat)) return(NULL)
    #     exprdat$group_name <- factor(exprdat$group_name, levels = unique(reactvals$group_levels))
    #     plot_types <- input$plot_type
    #     if (is_empty(reactvals$selected_groups)) {
    #         exprdat$group_col <- exprdat$group_name
    #     } else {
    #         exprdat$group_col <- ifelse(exprdat$group_name %in% reactvals$selected_groups, "Unselected", "Selected")
    #         exprdat$group_col <- factor(exprdat$group_col, levels = c("Selected", "Unselected"))
    #     }
    #     p <- ggplot(exprdat, aes(x = group_name, y = expr)) +
    #         theme_bw(base_size = 18) + 
    #         ylab(paste0(unique(exprdat$gene_symbol)," Expression [AU]")) +
    #         theme(panel.grid = element_blank(),
    #               legend.position = "none",
    #               axis.title.x = element_blank(),
    #               axis.ticks.x = element_blank(),
    #               axis.text.y = element_text(color = "black"),
    #               axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
    #     if (!is_empty(reactvals$selected_groups)) {
    #         p <- p + scale_color_manual(values = c(input$col_sel, input$col_unsel)) +
    #             scale_fill_manual(values = c(input$col_sel, input$col_unsel))
    #     }
    #     if ("dotplot" %in% plot_types) {
    #         p <- p + geom_point(aes(color = group_col),
    #                             position = position_dodge2(width = .5))
    #     }
    #     if ("violin" %in% plot_types) {
    #         p <- p + geom_violin(aes(fill = group_col))
    #     }
    #     if ("range" %in% plot_types) {
    #         p <- p + stat_summary(aes(color = group_col), 
    #                               geom = "pointrange",
    #                               fun.data = "mean_cl_boot")
    #     }
    #     if ("boxplot" %in% plot_types) {
    #         if ("dotplot" %in% plot_types & "violin" %in% plot_types) {
    #             p <- p + geom_boxplot(aes(fill = group_col), width = 0.15, outlier.shape = NA)
    #         } else if ("dotplot" %in% plot_types) {
    #             p <- p + geom_boxplot(aes(fill = group_col), width = 0.6, outlier.shape = NA)
    #         } else if ("violin" %in% plot_types) {
    #             p <- p + geom_boxplot(aes(fill = group_col), width = 0.15)
    #         } else {
    #             p <- p + geom_boxplot(aes(fill = group_col), width = 0.15)
    #         }
    #     }
    #     return(p)
    # })
    # 
    # # generate bivariate plot
    # plot_bivariate <- reactive({
    #     exprdat <- reactvals$exprdat
    #     if (is_empty(exprdat)) return(NULL)
    #     if (is_empty(reactvals$selected_groups)) {
    #         exprdat$group_col <- exprdat$group_name
    #     } else {
    #         exprdat$group_col <- ifelse(exprdat$group_name %in% reactvals$selected_groups, "B", "A")
    #     }
    #     exprdat$group_name <- factor(exprdat$group_name, levels = unique(reactvals$group_levels))
    #     genes <- unique(exprdat$gene_name)
    #     exprdat <- spread(exprdat, "gene_name", "expr")
    #     genenames <-  colnames(exprdat)[c(4:5)]
    #     colnames(exprdat)[c(4:5)] <- c("geneA", "geneB")
    #     p <- ggplot(exprdat, aes(x = geneA, y = geneB, color = group_col)) +
    #         theme_bw(base_size = 18) + 
    #         xlab(paste0(genenames[1]," Expression [AU]")) +
    #         ylab(paste0(genenames[2]," Expression [AU]")) +
    #         theme(panel.grid = element_blank(),
    #               axis.text.y = element_text(color = "black"),
    #               axis.text.x = element_text(color = "black")) +
    #         geom_point()
    #     return(p)
    # })
    # 
    # # generate boxplot plot
    # plot_qc_box <- reactive({
    #     exprdat <- reactvals$exprdat
    #     if (is_empty(exprdat)) return(NULL) 
    #     p <- ggplot(exprdat, aes(x = group_name, y = expr)) +
    #         geom_boxplot() +
    #         theme_bw(base_size = 18) + 
    #         ylab("Expression [AU]") +
    #         theme(panel.grid = element_blank(),
    #               axis.title.x = element_blank(),
    #               axis.ticks.x = element_blank(),
    #               axis.text.y = element_text(color = "black"),
    #               axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
    #     return(p)
    # })
    # 
    # # output distro plot
    # sel_expr_plot <- reactive({
    #     exprsum <- get_expr_summary()
    #     if (is_empty(exprsum)) {
    #         validate("Please select features & samples")
    #     }
    #     if (length(unique(exprsum$gene_name)) == 1) {
    #         plot_distribution()
    #     } else if (length(unique(exprsum$gene_name)) == 2) {
    #         plot_bivariate()
    #     } else {
    #         plot_qc_box()
    #     }
    # })
    # 
    # # output expr plot
    # output$expr_plot <- renderPlot({
    #     sel_expr_plot()
    # })
    # get_dims <- reactive({
    #     exprsum <- get_expr_summary()
    #     ncond <- length(unique(exprsum$group_name))
    #     width <- max(c(120, min(c(70*ncond, 1000))))
    #     list(w = width, h = 350)
    # })
    # output$expr_plot_ui <- renderUI({
    #     sel_expr_plot()
    #     dims <- get_dims()
    #     plotOutput("expr_plot", height = dims$h, width = dims$w)
    # })
    # 
    # # download expr plot as png
    # output$dl_expr_plot_png <- downloadHandler(
    #     filename = function() {
    #         exprsum <- get_expr_summary()
    #         if (length(unique(exprsum$gene_name)) > 2) {
    #             genes <- "overall"
    #         } else {
    #             genes <- paste0(unique(exprsum$gene_name), collapse = "_")
    #         }
    #         gseid <- input$gse
    #         paste0(gseid, "_", genes, "_expr.png")
    #     },
    #     content = function(file) {
    #         plot <- sel_expr_plot()
    #         dims <- get_dims()
    #         ggsave(plot = plot, filename = file, units = "mm",
    #                height = dims$h/3, width = dims$w/3)
    #     }
    # )
    # 
    # # download expr plot as ppt
    # output$dl_expr_plot_ppt <- downloadHandler(
    #     filename = function() {
    #         exprsum <- get_expr_summary()
    #         if (length(unique(exprsum$gene_name)) > 2) {
    #             genes <- "overall"
    #         } else {
    #             genes <- paste0(unique(exprsum$gene_name), collapse = "_")
    #         }
    #         gseid <- input$gse
    #         paste0(gseid, "_", genes, "_expr.pptx")
    #     },
    #     content = function(file) {
    #         plot <- sel_expr_plot()
    #         dims <- get_dims()
    #         file_pptx <- tempfile(fileext = ".pptx")
    #         gen_pptx(plot, file_pptx,
    #                  height = (dims$h/3)*0.039, 
    #                  width = (dims$w/3)*0.039)
    #         file.rename(from = file_pptx, to = file)
    #     }
    # )
}

# Run the application 
shinyApp(ui = ui, server = server)
