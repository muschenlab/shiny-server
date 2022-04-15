library(shiny)
library(tidyverse)
library(shinyWidgets)
library(shinycssloaders)
library(shinythemes)
library(officer)
library(rvg)
library(DT)

# ##  adding lineage & filter
# es <- readRDS("/poolio/expr_sets/es/MILE_leukemias_ES.rds")
# var <- apply(exprs(es), 1, var)
# es <- es[which(var > quantile(var, .3)), ]
# es$lineage <- NA
# pData(es)[es$disease %in% c("B-ALL","CLL","pre-B-ALL","Burkitt's"),]$lineage <- "B-cell"
# pData(es)[es$disease %in% c("T-ALL"),]$lineage <- "T-cell"
# pData(es)[es$disease %in% c("AML","CML","MDS"),]$lineage <- "Myeloid"
# pData(es)[es$disease %in% c("normal"),]$lineage <- "Normal"
# es <- es[,c(which(es$lineage == "Normal"),
#             which(es$disease == "pre-B-ALL"),
#             which(es$disease == "B-ALL"),
#             which(es$disease == "CLL"),
#             which(es$disease == "Burkitt's"),
#             which(es$lineage == "T-cell"),
#             which(es$lineage == "Myeloid"))]
# saveRDS(es, "MILE_leukemias_ES.rds")

################################# Setup #######################################

# load data
es <- readRDS("mile_data/MILE_leukemias_ES.rds")

# load plot templates
dendro_templ <- readRDS("mile_data/MILE_dendrogram_plot_outline.rds")
lineage_templ <- readRDS("mile_data/MILE_lineage_plot_outline.rds")

# set default factor levels
reactvals <- reactiveValues(grouplevels = unique(es$group),
                            lineagelevels = unique(es$lineage))
grouplevels = unique(es$group)
lineagelevels = unique(es$lineage)

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

# JS for substr
substr_js <- paste0("function(data, type, row, meta) {",
                    "return type === 'display' && data.length > 30 ?",
                    "'<span title=\"' + data + '\">' + data.substr(0, 30)",
                    "+ '...</span>' : data;}")

#################################### UI #######################################

ui <- fluidPage(
  
  # title
  title = "MILE", 
  theme = shinytheme("cosmo"),
  
  titlePanel(tags$h3(tags$a(
    imageOutput("icon", inline = TRUE),
    href="http://10.197.211.94:3838"), 
    "MILE", style = "color:black")),
  
  tabsetPanel(
    tabPanel("Expression", width = "100%",
             radioGroupButtons("plot_type", "Plot Element:",
                               choices = list("Lineage" = "lineage",
                                              "Barplot" = "barplot",
                                              "Boxplot" = "boxplot"),
                               justified = T,
                               size = 'xs',
                               width = '20%',
                               status = "default",
                               selected = "lineage"),
             fluidRow(column(width = 12, align="center",
                             uiOutput("main_plot_ui") %>%
                               withSpinner(color="#916f6f"))),
             div(style="display:inline-block;",
                 downloadButton("dl_main_plot_ppt", label = "PPT",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 downloadButton("dl_main_plot_png", label = "PNG",
                                style = "font-size:12px;height:30px;padding:5px;")),
             br(), br(),
             tags$h3("Gene Info:"),
             DT::DTOutput('gene_info_dt'),
             br(), br(),
             tags$h3("Expression Summary:"),
             DT::DTOutput('expr_summary_dt')),
    tabPanel("Expression Data", width = "100%",
             DT::DTOutput('expr_level_dt'),
             downloadButton("dl_expr_xls", label = "Excel",
                            style = "font-size:12px;height:30px;padding:5px;"))
  )
)

################################# Server ######################################

server <- function(input, output) {
  
  # icon
  output$icon <- renderImage(list(src = "../hexagons/mile.png",
                                  height = "90px", width = "85px"), 
                             deleteFile = F)
  
  # gene info table output
  output$gene_info_dt <- DT::renderDT({
    fdat <- as.data.frame(fData(es))
    colnames(fdat) <- Hmisc::capitalize(gsub("_", " ", colnames(fdat)))
    DT::datatable(data = fdat,
                  rownames = F,
                  escape = F,
                  options = list(
                    pagelength = 5,
                    lengthMenu = list(c(5, 10, 20),
                                      c('5','10','20')),
                    columnDefs = list(list(
                      targets = 2:4,
                      render = JS(
                        "function(data, type, row, meta) {",
                        "return type === 'display' && data.length > 15 ?",
                        "'<span title=\"' + data + '\">' + data.substr(0, 15) + '...</span>' : data;",
                        "}")))
                  ),
                  callback = JS('table.page(3).draw(false);'),
                  caption = "Select probes to plot.")
  })
  
  # get selected probes from dt
  get_features <- reactive({
    return(fData(es)[input$gene_info_dt_rows_selected,])
  })
  
  # get all expr dat 
  get_expr_dat <- reactive({
    finfo <- get_features()
    if (is_empty(finfo)) return(NULL)
    expr <- cbind(finfo[,c("probe_id", "gene_symbol")],
                  exprs(es[finfo$probe_id,])) %>% 
      gather("sample_id", "expr", -c("probe_id", "gene_symbol"))
    expr$group <- pData(es)[match(expr$sample_id, pData(es)$sample_id),]$group
    expr$lineage <- pData(es)[match(expr$sample_id, pData(es)$sample_id),]$lineage
    expr$group <- factor(expr$group, levels = reactvals$grouplevels)
    expr$lineage <- factor(expr$lineage, levels = reactvals$lineagelevels)
    return(expr)
  })
  
  # summarise expr data
  get_expr_summary <- reactive({
    expr <- get_expr_dat()
    if (is_empty(expr)) return(NULL)
    if (length(unique(expr$gene_symbol)) > 1) {
      showNotification("WARNING: multiple genes selected, average will be taken.",
                       type = "warning")
    }
    exprsum <- expr %>% group_by(group) %>%
      summarise(gene_symbol = paste0(unique(gene_symbol), collapse = ";"), 
                n = n(),
                average = mean(expr, na.rm = T),
                SD = sd(expr, na.rm = T),
                lineage = unique(lineage)) %>%
      group_by(gene_symbol) %>%
      mutate(zscore = scale(average)[,1])
    return(exprsum)
  })
  
  # expr summary table output
  output$expr_summary_dt <- DT::renderDT({
    exprdat <- as.data.frame(get_expr_summary())
    if (nrow(exprdat) < 1) return(NULL)
    exprdat[,c(4:5,7)] <- apply(exprdat[,c(4:5,7)], 2, function(x) round(x, 3))
    gene_names <- paste0(unique(exprdat$gene_symbol), collapse="_")
    colnames(exprdat) <- Hmisc::capitalize(gsub("_"," ",colnames(exprdat)))
    DT::datatable(exprdat,  
                  rownames = F,
                  extensions = 'Buttons',
                  options = list(dom = 'rtipB',
                                 buttons = list('pageLength',
                                                list(extend = 'copy',
                                                     title = NULL),
                                                list(extend = 'excel',
                                                     filename = paste0("MILE_expr_summary_", gene_names),
                                                     title = NULL)))) 
  })
  
  # generate lineage plot
  gen_lineage_plot <- reactive({
    plot <- lineage_templ +
      theme(plot.margin = margin(l=100))
    finfo <- get_features()
    if (nrow(finfo) < 1) {
      showNotification("Select a gene from the table to plot.")  
      return(list(p=plot, h=5, w=6.5))
    }
    gene_name <- unique(finfo$gene_symbol)
    exprdat <- get_expr_summary()
    if (is_empty(exprdat)) {
      return(NULL)
    } else {
      plot$data$col <- exprdat[match(plot$data$names, exprdat$group),]$zscore
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
      return(list(p=plot, h=5, w=6.5))
    }
  })
  
  # generate boxplot
  gen_boxplot <- reactive({
    plotdat <- na.omit(get_expr_dat())
    validate(need(nrow(plotdat) > 1, "Select a gene from the table below"))
    print(head(plotdat))
    gene_name <- unique(plotdat$gene_symbol)
    print(gene_name)
    nlin <- plotdat %>%
      group_by(lineage) %>%
      summarise(n = length(unique(group))) %>%
      mutate(lineage = factor(lineage, levels=levels(plotdat$lineage))) %>%
      arrange(lineage)
    print(head(nlin))
    ymax <- ceiling(max(plotdat$expr))
    cat_lines <- data.frame(x = c(0.5, (cumsum(nlin$n)[-length(cumsum(nlin$n))])+.75),
                            xend = cumsum(nlin$n)+.25,
                            y = ymax + (ymax/20),
                            yend = ymax + (ymax/20),
                            lineage = nlin$lineage)
    cat_labs <- data.frame(x = cat_lines$x+((cat_lines$xend-cat_lines$x)/2),
                           y = ymax + (ymax/10),
                           label = nlin$lineage)
    plot <- ggplot(plotdat, aes(x = group, y = expr)) +
      geom_boxplot(aes(fill = lineage), width = 0.75, alpha = 0.75) +
      geom_segment(data = cat_lines, lwd = 3, alpha = .5,
                   aes(x = x, y = y, xend = xend, yend = yend, color = lineage)) +
      coord_cartesian(clip = "off", expand = F, 
                      ylim = c(0, ymax),
                      xlim = c(-.15, (length(unique(plotdat$group))+0.25))) +
      geom_text(data = cat_labs, aes(x=x, y=y, label=label), size = 5) +
      theme_bw(base_size = 16) +
      scale_fill_manual(values = c("grey40","steelblue", "goldenrod", "firebrick4")) +
      scale_color_manual(values = c("grey40","steelblue", "goldenrod", "firebrick4")) +
      ylab(paste0(gene_name, " expression [AU]")) +
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
    return(list(p=plot, h=4, w=7))
  })
  
  # generate barplot
  gen_barplot <- reactive({
    plotdat <- get_expr_summary()
    validate(need(!is_empty(plotdat), "Select a gene from the table below"))
    gene_name <- unique(plotdat$gene_symbol)
    plotdat <- na.omit(plotdat)
    nlin <- plotdat %>%
      group_by(lineage) %>%
      summarise(n = length(unique(group))) %>%
      mutate(lineage = factor(lineage, levels=levels(plotdat$lineage))) %>%
      arrange(lineage)
    ymax <- ceiling(max(plotdat$average + plotdat$SD))
    cat_lines <- data.frame(x = c(0.5, (cumsum(nlin$n)[-length(cumsum(nlin$n))])+.75),
                            xend = cumsum(nlin$n)+.25,
                            y = ymax + (ymax/20),
                            yend = ymax + (ymax/20),
                            lineage = nlin$lineage)
    cat_labs <- data.frame(x = cat_lines$x+((cat_lines$xend-cat_lines$x)/2),
                           y = ymax + (ymax/10),
                           label = nlin$lineage)
    plot <- ggplot(plotdat, aes(x = group, y = average)) +
      geom_bar(aes(fill = lineage), alpha = 0.75,
               stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = average - SD, ymax = average + SD), width = 0.3) +
      geom_segment(data = cat_lines, lwd = 3, alpha = .5,
                   aes(x = x, y = y, xend = xend, yend = yend, color = lineage)) +
      coord_cartesian(clip = "off", expand = F, 
                      ylim = c(0, ymax),
                      xlim = c(-.15, (length(unique(plotdat$group))+0.25))) +
      geom_text(data = cat_labs, aes(x=x, y=y, label=label), size = 5) +
      theme_bw(base_size = 16) +
      scale_fill_manual(values = c("grey40","steelblue", "goldenrod", "firebrick4")) +
      scale_color_manual(values = c("grey40","steelblue", "goldenrod", "firebrick4")) +
      ylab(paste0(gene_name, " expression [AU]")) +
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
    return(list(p=plot, h=4, w=7))
  })
  
  # choose correct plot type
  gen_main_plot <- reactive({
    validate(need(!is_empty(input$plot_type), "Select a plot type to continue"))
    print(input$plot_type)
    if (input$plot_type == "lineage") {
      print("A")
      gen_lineage_plot()
    } else if (input$plot_type == "boxplot") {
      print("B")
      gen_boxplot()
    } else if (input$plot_type == "barplot") {
      print("C")
      gen_barplot()
    } else return(NULL)
  })
  
  # output main plot
  output$main_plot_ui <- renderUI ({
    plotdat <- gen_main_plot()
    output$main_plot <- renderPlot(plotdat$p)
    plotOutput("main_plot", width = plotdat$w*100, height = plotdat$h*100)
  })
  
  # download main plot as png
  output$dl_main_plot_png <- downloadHandler(
    filename = function() {
      gene_name <- get_expr_summary()$gene_symbol[1]
      if (is_empty(gene_name)) gene_name <- "diagram"
      paste0("MILE_", input$plot_type, "_", gene_name, ".png")
    },
    content = function(file) {
      plotdat <- gen_main_plot()
      ggsave(plot = plotdat$p, filename = file, 
             height = plotdat$h, width = plotdat$w)
    }
  )
  
  # download main plot as ppt
  output$dl_main_plot_ppt <- downloadHandler(
    filename = function() {
      gene_name <- get_expr_summary()$gene_symbol[1]
      if (is_empty(gene_name)) gene_name <- "diagram"
      paste0("MILE_", input$plot_type, "_", gene_name, ".pptx")
    },
    content = function(file) {
      plotdat <- gen_main_plot()
      file_pptx <- tempfile(fileext = ".pptx")
      gen_pptx(plotdat$p, file_pptx, 
               height = plotdat$h,
               width = plotdat$w)
      file.rename(from = file_pptx, to = file)
    }
  )
  
  # expression data table output
  output$expr_level_dt <- DT::renderDT({
    expr <- get_expr_dat()
    genes <- paste(unique(expr$gene_symbol), collapse = "_")
    DT::datatable(as.data.frame(expr), 
                  escape = F, 
                  rownames = F,
                  caption = paste0("Expression levels: ", genes),
                  colnames = Hmisc::capitalize(gsub("_", " ", colnames(expr))),
                  editable = "cell",
                  extensions = 'Buttons',
                  options = list(dom = 'frtipB',
                                 buttons = list('pageLength',
                                                list(extend = 'csv',
                                                     filename = paste0("MILE_expression_levels_", 
                                                                       genes, ".csv"),
                                                     title = NULL),
                                                list(extend = 'excel',
                                                     filename = paste0("MILE_expression_levels_", 
                                                                       genes, ".xlsx"),
                                                     title = NULL)),
                                 pagelength = 10,
                                 lengthMenu = list(c(10, 25, 100, -1),
                                                   c('10','25','100','All')))) %>%
      DT::formatRound(c("expr"), 3)
  })
  
  # download full expression info
  output$dl_expr_xls <- downloadHandler(
    filename = function() {
      expr <- get_expr_dat()
      genes <- paste(unique(expr$gene_symbol), collapse = "_")
      paste0("MILE_", genes, "_levels.xlsx")
    },
    content = function(file) {
      expr <- get_expr_dat()
      writexl::write_xlsx(expr, path=file)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
