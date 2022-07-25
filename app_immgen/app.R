library(shiny)
library(plyr)
library(tidyverse)
library(shinyWidgets)
library(shinycssloaders)
library(shinythemes)
library(officer)
library(rvg)
library(SummarizedExperiment)
options(stringsAsFactors = F)

################################# Setup #######################################

# initialise reactive vals
reactvals <- reactiveValues(sel_samples = NULL)

# load data
counts_se <- readRDS("/srv/shiny-server/app_immgen/immgen_data/Immgen_RNAseq_norm_counts_SE.2020-01-10.rds")

# load plot template
plot_templ <- readRDS("/srv/shiny-server/app_immgen/immgen_data/Immgen_simple_lineage_ggplot_obj.rds")

# load celltype info
cell_info <- read.csv("/srv/shiny-server/app_immgen/immgen_data/Immgen_celltypes.csv")

# pre-process counts to average by group
groups <- plot_templ$data$names
avnorm <- lapply(groups, function(x) {
    idx <- which(colData(counts_se)$group_simple == x)
    if(length(idx) == 1) {
        return(assay(counts_se, "norm_counts")[,idx])
    } else {
        return(rowMeans(assay(counts_se, "norm_counts")[,idx]))
    }
}) %>% bind_cols() %>% as.data.frame()
colnames(avnorm) <- groups
rownames(avnorm) <- rownames(counts_se)

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

#################################### UI #######################################

ui <- fluidPage(

    # title
    title = "Immgen", 
    theme = shinytheme("cosmo"),
    
    titlePanel(tags$h3(tags$a(
        imageOutput("icon", inline = TRUE),
        href="http://10.197.211.94:80"), 
        "ImmGen: Immunological Genome Project", style = "color:black")),

    tabsetPanel(
        tabPanel("Expression Plot", width = "100%", align = "center",
                 br(),
                 radioGroupButtons("plot_type", "Plot Type:",
                                   choices = list("Lineage" = "lineage",
                                                  "Barplot" = "barplot",
                                                  "Boxplot" = "boxplot"),
                                   justified = T,
                                   size = 'xs',
                                   width = '20%',
                                   status = "default",
                                   selected = "lineage"),
                 uiOutput("main_plot_ui") %>% withSpinner(color="#ffffff"),
                 div(style="display:inline-block;",
                     downloadButton("dl_main_plot_ppt", label = "PPT",
                                    style = "font-size:12px;height:30px;padding:5px;"),
                     downloadButton("dl_main_plot_png", label = "PNG",
                                    style = "font-size:12px;height:30px;padding:5px;")),
                 br(), br(),
                 tags$h4("Gene Info"),
                 DT::DTOutput('gene_info_dt'),
                 br(), br(),
                 tags$h4("Cell Type Info"),
                 DT::DTOutput('cell_info_dt')),
        tabPanel("Expression Data", width = "100%",
                 DT::DTOutput('expr_level_dt'))
    )
)

################################# Server ######################################

server <- function(input, output, session) {

    # icon
    output$icon <- renderImage(list(src = "../hexagons/immgen.png",
                                    height = "90px", width = "85px"), 
                               deleteFile = F)
    
    # gene info table output
    output$gene_info_dt <- DT::renderDT(
        DT::datatable(as.data.frame(rowData(counts_se)),
                      escape = F, 
                      rownames = F,
                      colnames = Hmisc::capitalize(gsub("_", " ", colnames(rowData(counts_se)))),
                      selection = list(mode = "single"),
                      options = list(pagelength = 5,
                                     lengthMenu = c(5, 10, 15, 20)))
    )
    
    # cell info table output
    output$cell_info_dt <- DT::renderDT({
      cell_info <-   distinct(as.data.frame(colData(counts_se)))
      cell_info <- cell_info[order(cell_info$lineage, cell_info$cell_family),]
      reactvals$cell_info <- cell_info
      DT::datatable(cell_info, 
                    rownames = F,
                    selection = list(mode = 'multiple', target = "row", selected = c(3:4)),
                    colnames = Hmisc::capitalize(gsub("_", " ", colnames(cell_info))))
    })
    
    # get selected samples
    observeEvent(input$cell_info_dt_rows_selected, {
      cell_info <- reactvals$cell_info
      reactvals$sel_samples <- unique(cell_info[input$cell_info_dt_rows_selected,]$group_simple)
    })
    
    # get selected gene from dt
    get_gene_name <- reactive({
        return(rownames(counts_se)[input$gene_info_dt_rows_selected])
    })
    
    # get expr level info
    get_expr_level <- reactive({
        gene_name <- get_gene_name()
        if (is_empty(gene_name)) {
            return()
        } else {
            expr <- data.frame(sample_id = row.names(colData(counts_se)),
                               cell_type_id = colData(counts_se)$celltype,
                               lineage = colData(counts_se)$lineage,
                               cell_group = colData(counts_se)$group_simple,
                               counts = as.numeric(assay(counts_se, "counts")[gene_name,]),
                               normalised_counts = as.numeric(assay(counts_se, "norm_counts")[gene_name,]),
                               TPM = as.numeric(assay(counts_se, "iTPM")[gene_name,]),
                               sorting = colData(counts_se)$sorting_markers)
            return(expr)
        }
    })
    
    # get z-scores
    get_zscores <- reactive({
        gene_name <- get_gene_name()
        if (is_empty(gene_name)) {
            return()
        } else {
            feature_idx <- which(rownames(avnorm) == gene_name)
            zscores <- scale(as.numeric(2^avnorm[feature_idx, ]))[,1]
            return(zscores)
        }
    })
    
    # generate lineage plot
    gen_lineage_plot <- reactive({
        plot <- plot_templ
        gene_name <- get_gene_name()
        zscores <- get_zscores()
        if (is_empty(zscores)) {
          showNotification("Select a gene to plot from the table.")  
          return(plot)
        } else {
          plot$data$col <- zscores
          plot <- plot + 
              scale_fill_gradient2(low = "white", mid = "white", high = "firebrick", 
                                   name = paste0(gene_name, " expression\n[z-score]")) +
              theme(text = element_text(family = "Arial"),
                    legend.title = element_text(size = 12),
                    legend.key = element_rect(color = "black"),
                    legend.key.size = unit(.8, "cm"),
                    legend.text = element_text(size = 12)) +
              guides(fill = guide_colourbar(title.position="top", 
                                            frame.colour = "black",
                                            frame.linewidth = .8))
          return(plot)
        }
    })
    
    # generate boxplot
    gen_boxplot <- reactive({
      gene_name <- get_gene_name()
      plotdat <- na.omit(get_expr_level())
      validate(need(!is_empty(plotdat), "Select a gene from the table below"))
      plotdat$cell_group <- factor(plotdat$cell_group,
                                   levels = unique(na.omit(cell_info$group_simple)))
      n <- ddply(plotdat, .(lineage), function(x) length(unique(x$cell_group)))
      n$lineage <- factor(n$lineage, levels = unique(cell_info$lineage))
      n <- n[order(n$lineage),]
      ymax <- ceiling(max(plotdat$TPM))
      cat_lines <- data.frame(x = c(0.5, (cumsum(n$V1)[-length(cumsum(n$V1))])+.75),
                              xend = cumsum(n$V1)+.25,
                              y = ymax + (ymax/20),
                              yend = ymax + (ymax/20),
                              lineage = n$lineage)
      cat_labs <- data.frame(x = cat_lines$x+((cat_lines$xend-cat_lines$x)/2),
                             y = ymax + (ymax/10),
                             label = n$lineage)
      ggplot(plotdat, aes(x = cell_group, y = TPM)) +
        geom_boxplot(aes(fill = cell_group %in% c("FrBC","FrE"))) +
        geom_segment(data = cat_lines, lwd = 3, alpha = .5,
                     aes(x = x, y = y, xend = xend, yend = yend, color = lineage)) +
        coord_cartesian(clip = "off", expand = F, 
                        ylim = c(0, ceiling(max(plotdat$TPM))),
                        xlim = c(-.15, (length(unique(plotdat$cell_group))+0.25))) +
        geom_text(data = cat_labs, aes(x=x, y=y, label=label), size = 5) +
        theme_bw(base_size = 16) +
        scale_fill_manual(values = c("grey60","firebrick4")) +
        ylab(paste0(gene_name, " expression [TPM]")) +
        theme(panel.grid = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_line(size = .4, color = "black"),
              axis.line.y = element_line(size = .4, color = "black"),
              axis.line.x = element_line(size = .4, color = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(),
              axis.text = element_text(color = "black", family = "Arial"),
              text = element_text(color = "black", family = "Arial"),
              legend.position = "none",
              plot.margin = unit(c(4,1,2,1), "lines"))
    })
    
    # generate barplot
    gen_barplot <- reactive({
      gene_name <- get_gene_name()
      plotdat <- get_expr_level()
      validate(need(!is_empty(plotdat), "Select a gene from the table below"))
      plotdat$cell_group <- factor(plotdat$cell_group,
                                   levels = unique(na.omit(cell_info$group_simple)))
      plotdat <- na.omit(plotdat)
      n <- ddply(plotdat, .(lineage), function(x) length(unique(x$cell_group)))
      n$lineage <- factor(n$lineage, levels = unique(cell_info$lineage))
      n <- n[order(n$lineage),]
      ymax <- ceiling(max(plotdat$TPM))
      cat_lines <- data.frame(x = c(0.5, (cumsum(n$V1)[-length(cumsum(n$V1))])+.75),
                              xend = cumsum(n$V1)+.25,
                              y = ymax + (ymax/20),
                              yend = ymax + (ymax/20),
                              lineage = n$lineage)
      cat_labs <- data.frame(x = cat_lines$x+((cat_lines$xend-cat_lines$x)/2),
                             y = ymax + (ymax/10),
                             label = n$lineage)
      ggplot(plotdat, aes(x = cell_group, y = TPM)) +
        geom_bar(aes(fill = cell_group %in% reactvals$sel_samples),
                 stat = "identity", position = "dodge") +
        geom_segment(data = cat_lines, lwd = 3, alpha = .5,
                     aes(x = x, y = y, xend = xend, yend = yend, color = lineage)) +
        coord_cartesian(clip = "off", expand = F, 
                        ylim = c(0, ceiling(max(plotdat$TPM))),
                        xlim = c(-.15, (length(unique(plotdat$cell_group))+0.25))) +
        geom_text(data = cat_labs, aes(x=x, y=y, label=label), size = 5) +
        theme_bw(base_size = 16) +
        scale_fill_manual(values = c("grey60","firebrick4")) +
        ylab(paste0(gene_name, " expression [TPM]")) +
        theme(panel.grid = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_line(size = .4, color = "black"),
              axis.line.y = element_line(size = .4, color = "black"),
              axis.line.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(),
              axis.text = element_text(color = "black", family = "Arial"),
              text = element_text(color = "black", family = "Arial"),
              legend.position = "none",
              plot.margin = unit(c(4,1,2,1), "lines"))
    })
    
    # choose correct plot type
    gen_main_plot <- reactive({
      validate(need(!is_empty(input$plot_type), "Select a plot type to continue"))
      if (input$plot_type == "lineage") {
        gen_lineage_plot()
      } else if (input$plot_type == "boxplot") {
        gen_boxplot()
      } else if (input$plot_type == "barplot") {
        gen_barplot()
      } 
    })
    output$main_plot <- renderPlot({
      gen_main_plot()
    })
    output$main_plot_ui <- renderUI({
      if (input$plot_type == "lineage") {
        plotOutput("main_plot", height = 500, width = 700)
      } else {
        plotOutput("main_plot", height = 350, width = 800)
      }
    })
    
    # download main plot as png
    output$dl_main_plot_png <- downloadHandler(
        filename = function() {
            gene_name <- get_gene_name()
            if (is_empty(gene_name)) gene_name <- "cell_types"
            paste0("Immgen_", input$plot_type, "_", gene_name, ".png")
        },
        content = function(file) {
            plot <- gen_main_plot()
            if (input$plot_type == "lineage") {
              ggsave(plot = plot, filename = file, units = "mm",
                     height = 500/3, width = 700/3)
            } else {
              ggsave(plot = plot, filename = file, units = "mm",
                     height = 350/3, width = 800/3)
            }
        }
    )
    
    # download main plot as ppt
    output$dl_main_plot_ppt <- downloadHandler(
        filename = function() {
            gene_name <- get_gene_name()
            if (is_empty(gene_name)) gene_name <- "cell_types"
            paste0("Immgen_", input$plot_type, "_", gene_name, ".pptx")
        },
        content = function(file) {
            plot <- gen_main_plot()
            file_pptx <- tempfile(fileext = ".pptx")
            if (input$plot_type == "lineage") {
              gen_pptx(plot = plot, file = file,
                     height = (500/3)*0.039,
                     width = (700/3)*0.039)
            } else {
              gen_pptx(plot = plot, file = file, 
                     height = (350/3)*0.039,
                     width = (800/3)*0.039)
            }
            file.rename(from = file_pptx, to = file)
        }
    )
    
    # expression data table output
    output$expr_level_dt <- DT::renderDT(
        DT::datatable(as.data.frame(get_expr_level()), 
                      escape = F, 
                      rownames = F,
                      caption = paste0("Expression levels: ", get_gene_name()),
                      colnames = Hmisc::capitalize(gsub("_", " ", colnames(get_expr_level()))),
                      editable = "cell",
                      extensions = 'Buttons',
                      options = list(dom = 'frtipB',
                                     buttons = list('pageLength',
                                                    list(extend = 'copy',
                                                         title = NULL),
                                                    list(extend = 'print',
                                                         title = paste0("Immgen RNAseq expression levels: ", 
                                                                        get_gene_name())),
                                                    list(extend = 'csv',
                                                         filename = paste0("Immgen_RNAseq_expression_levels_", 
                                                                           get_gene_name()),
                                                         title = NULL),
                                                    list(extend = 'excel',
                                                         filename = paste0("Immgen_RNAseq_expression_levels_", 
                                                                           get_gene_name()),
                                                         title = NULL),
                                                    list(extend = 'pdf',
                                                         filename = paste0("Immgen_RNAseq_expression_levels_", 
                                                                           get_gene_name()),
                                                         title = paste0("Immgen RNAseq expression levels: ", 
                                                                        get_gene_name()))),
                                     pagelength = 10,
                                     lengthMenu = list(c(10, 25, 100, -1),
                                                       c('10', '25', '100','All')))) %>%
            DT::formatRound(c("normalised_counts", "TPM"), 3)
    )
    
    # gene-gene correlation
    reactive({
      geneA <- get_gene_name()
      if (is_empty(geneA)) return(NULL)
      ggmat <- lapply(rownames(counts_se), function(geneB) {
        cor.test(assay(counts_se, "TPM")[geneA,],
                 assay(counts_se, "TPM")[geneB,])
      })
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
