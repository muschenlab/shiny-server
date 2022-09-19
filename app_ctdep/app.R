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
library(Biobase)

################################### Setup #####################################

# load data
load("dep_data/depdata.rda")

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

# reactive values
reactvals <- reactiveValues()

# boxplot function
gen_boxplot <- function(plotdat, x,  y, xlab = "", ylab) {
    if (is_empty(plotdat)) return(NULL)
    ggplot(plotdat, aes_string(x = x, y = y)) +
        geom_boxplot(aes_string(fill = x)) +
        coord_flip() +
        theme_bw(base_size=14) +
        theme(panel.grid = element_blank(),
              legend.position = "none") +
        xlab(xlab) + ylab(ylab) +
        scale_fill_manual(values=pickcols(jellypal, length(levels(plotdat[,x]))))
}

# density plot function
gen_densplot <- function(plotdat, x, y, xlab, ylab = "Frequency") {
    if (is_empty(plotdat)) return(NULL)
    ggplot(plotdat, aes_string(x = x)) +
        geom_density(aes_string(fill = y)) +
        theme_bw(base_size=14) +
        theme(panel.grid = element_blank(),
              legend.position = "bottom",
              legend.box = "horizontal",
              legend.title = element_blank()) +
        xlab(xlab) + ylab(ylab) +
        scale_fill_manual(values=pickcols(jellypal, length(levels(plotdat[,y]))))
}

# function just to print plain labels on log scale axis
plain <- function(x,...) {
    format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

# pallette function
jellypal <- c("#714481cc","#712e7ccc","#9b6295cc",
              "#b36f8ecc","#b9869fcc","#d49aa2cc",
              "#a2b2a8cc","#a7c3c6cc","#7d8b99cc",
              "#528199cc","#426284cc",
              "#eccad3cc","#eed9e0ff","#d0e6ea99")
    
pickcols <- function(pal, n) {
    idx <- c(round(seq(1, length(pal), length(pal)/(n-1))), length(pal))
    pal[idx]
}

##################################### UI ######################################
ui <- fluidPage(
    
    # title
    title = "ctDEP", 
    theme = shinytheme("cosmo"),
    titlePanel(tags$h2(tags$a(
        imageOutput("icon", inline = TRUE),
        href="http://10.197.211.94:80"), "ctDEP")),
    
    # change colors for DT row/column selection
    tags$style(HTML('table.dataTable tr.selected td{background-color: #B9869F99 !important;}')),
    tags$style(HTML('table.dataTable td.selected {background-color: #D0E6EA99 !important;}')),
    tags$style(HTML(".tabbable > .nav > li > a { color:#A7C4C6FF}")),
    tags$style(HTML(".tabbable > .nav > li[class=active] > a { color:black}")),
    
    # main UI
    tabsetPanel(
        tabPanel("Compounds",
                 
                 # celltype comparison choice
                 div(style = "font-size:13px;", uiOutput("comparison_choice")),
                 tags$h5(paste0("Select celltype comparison of interest, ",
                                "then select a row from the table below to explore sensitivity metrics")),
                 
                 # scatter plot of average
                 tags$h3("Average sensitivity per compound:"),
                 tags$h5("Click or drag a box around compounds to select"),
                 fluidRow(align="center",
                          plotOutput("sens_scatter", height = 450, width = 450, brush = "scatter_brush")),
                 
                 # metrics table
                 br(),br(),
                 tags$h3("Selected compounds:"),
                 tags$h5("Select a row to plot compound/target sensitivity accross cell lines below"),
                 div(DT::dataTableOutput("celltype_dt"),
                     style = "font-size:90%"),
                 
                 # rankings plot
                 tags$h4("Differential rankings:"),
                 tags$h5(paste0("Differential sensitivity ranks indicate how selectively sensitive the celltype of interest ",
                                "is relative to the control population (higher rank = more selective). Three primary metrics are ",
                                "considered - the CRISPR effect score, RNAi effect score, and drug sensitivity score (DSS).")),
                 fluidRow(align="center",
                          plotOutput("metrics_overview", height = 175, width = 400) %>% withSpinner(color = "#D0E6EA99")
                 ),
                 
                 # DSS plots
                 tags$h4("Drug sensitivity scores (DSS):"),
                 tags$h5(paste0("Drug sensitivity scores integrate the dose-response curve into a single metric - ",
                                "(higher score = more sensitive, range 0-100)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "20%", "22%"),
                                      plotOutput("dss_dens", height = 225), 
                                      plotOutput("dss_box_ds", height = 225), 
                                      plotOutput("dss_box_st", height = 225),
                                      plotOutput("ec50_plot", height = 225))),
                 # CRISPR plots
                 tags$h4("CRISPR:"),
                 tags$h5(paste0("CRISPR effect scores indicate the effect of gene knock-out on viability - ",
                                "(lower score = more sensitive)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "20%"),
                                      plotOutput("crispr_dens", height = 225), 
                                      plotOutput("crispr_box_ds", height = 225), 
                                      plotOutput("crispr_box_st", height = 225))),
                 # RNAi plots
                 tags$h4("RNAi:"),
                 tags$h5(paste0("CRISPR effect scores indicate the effect of gene knock-down on viability - ",
                                "(lower scores = more sensitive)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "20%"),
                                      plotOutput("rnai_dens", height = 225), 
                                      plotOutput("rnai_box_ds", height = 225), 
                                      plotOutput("rnai_box_st", height = 225))),
                 
                 # data for individual cell lines
                 br(),br(),
                 tags$h3("Cell-line level data:"),
                 tags$h4("Drugs:"),
                 div(DT::dataTableOutput("cl_drugdat"),
                     style = "font-size:90%"),
                 tags$h4("CRISPR:"),
                 div(DT::dataTableOutput("cl_crisprdat"),
                     style = "font-size:90%"),
                 tags$h4("RNAi:"),
                 div(DT::dataTableOutput("cl_rnaidat"),
                     style = "font-size:90%")
                 
        ),
        tabPanel("CRISPR",
                 tags$h3("under construction"))
    )
)


################################### Server ####################################
server <- function(input, output, session) {
    
    # icon
    output$icon <- renderImage(list(src = "../hexagons/jelly.png",
                                    height = "95px", width = "85px"), 
                               deleteFile = F)
    
    # cell type comparison choice UI
    output$comparison_choice <- renderUI({
        selectInput('comparison', 'Celltype comparison:', 
                    names(data$ct_diff),
                    "B-cell vs solid tumor")
    })

    # # observe cell type comparison choice
    # observeEvent(input$comparison, {
    #     if (is_empty(input$comparison)) return(NULL)
    #     reactvals$ct_summary_df <- data$ct_diff[[input$comparison]]
    #     reactvals$ct_summary_sel <- reactvals$ct_summary_df
    #     reactvals$ct_sel <- ifelse(grepl("B-cell", input$comparison),
    #                                "B-cell", "T-cell")
    #     print(input$comparison)
    #     print(reactvals$ct_sel)
    # })
    
    # lazy tmp way to do this before I fix custom ct
    get_ct <- reactive({
        ifelse(grepl("B-cell", input$comparison), "B-cell", "T-cell") 
    })
    
    # scatter plot of average sensitivty scores per compound
    output$sens_scatter <- renderPlot({
        selmetric <- "DSS3_avg" ## tmp
        ct <- get_ct()
        plotdat <- data$ct_diff[[input$comparison]] %>%
            select_if(grepl(selmetric, names(.))) %>%
            distinct()
        if (is_empty(plotdat)) return(NULL)
        colvar <- c(paste0("`", selmetric, "_", ct, "`"), 
                    paste0("`", selmetric, "_ST`"))
        plotdat %>%
            ggplot(aes_string(x=colvar[[1]], y=colvar[[2]])) +
            geom_point(color = "#71448199", size = 2) +
            geom_abline() +
            scale_x_continuous(name = paste0("Average ", selmetric, " score: ", ct)) +
            scale_y_continuous(name = paste0("Average ", selmetric, " score: ", "solid tumor")) +
            theme_bw(base_size = 18) +
            theme()
    })
    
    ### This for some reason causes an error whereby info is shared between user sessions
    ### doesn't appear to occur outside of observeEvent functions??
    # #  selecting summary table rows from scatter brushing
    # observeEvent(input$scatter_brush, {
    #     reactvals$ct_summary_sel <- reactvals$ct_summary_df %>%
    #         brushedPoints(input$scatter_brush)
    # })
    
    get_ct_dt_data <- reactive({
        if (is_empty(input$scatter_brush)) {
            data$ct_diff[[input$comparison]] %>%
                select_at(1:9) %>%
                mutate_if(is.numeric, round, digits = 3)
        } else {
            data$ct_diff[[input$comparison]] %>%
                brushedPoints(input$scatter_brush) %>%
                select_at(1:9) %>%
                mutate_if(is.numeric, round, digits = 3)
        }
    }) 
    
    # feature info table for feature/row selection
    output$celltype_dt <- DT::renderDataTable({
        ct_summary_sel <- get_ct_dt_data()
        if (is.null(ct_summary_sel)) return(NULL)
        DT::datatable(
            data = ct_summary_sel, 
            rownames = F,
            selection = list(mode = 'single', target = "row", selected = 1),
            options = list(columnDefs = list(list(
                targets = 0:(ncol(ct_summary_sel)-1),
                render = JS(
                    "function(data, type, row, meta) {",
                    "return type === 'display' && data.length > 30 ?",
                    "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                    "}")
            ))),
            callback = JS('table.page(3).draw(false);')
        )
    })
    
    ### This for some reason causes an error whereby info is shared between user sessions
    ### doesn't appear to occur outside of observeEvent functions??
    # # select gene
    # observeEvent(input$celltype_dt_rows_selected, {
    #     if (is_empty(input$celltype_dt_rows_selected)) return(NULL)
    #     df <- reactvals$ct_summary_sel[input$celltype_dt_rows_selected, ]
    #     reactvals$gene <- df$Gene
    #     print(reactvals$gene)
    # })
    
    get_ct_summary_df <- reactive({
        if (is_empty(input$scatter_brush)) {
            data$ct_diff[[input$comparison]] %>%
                slice(input$celltype_dt_rows_selected)
        } else {
            data$ct_diff[[input$comparison]] %>%
                brushedPoints(input$scatter_brush) %>%
                slice(input$celltype_dt_rows_selected)
        }
    })
    
    # metrics overview
    output$metrics_overview <- renderPlot({
        plotdat <- get_ct_summary_df() %>%
            select_at(1:9) %>%
            gather("metric", "rank", -c(1:3)) %>%
            mutate(rank = rank*100, 
                   metric = factor(metric, 
                                   levels = c("Avg_rank","dDSS1_rank","dDSS2_rank","dDSS3_rank","dRNAi_rank","dCRISPR_rank"),
                                   labels = c("overall","dDSS1","dDSS2","dDSS3","dRNAi","dCRISPR")))
        cpd <- unique(plotdat$Compound)
        gene <- unique(plotdat$Gene)
        ggplot(plotdat, aes(x=metric, y=rank)) +
            geom_point(size=4, color="#714481ff", shape=21) +
            coord_flip(clip="off") +
            scale_y_continuous(limits=c(0,100), expand=c(0,0), name="Differential rank") +
            xlab("Metric") +
            theme_bw(base_size = 14) +
            theme(panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.line.x = element_line(),
                  panel.grid.major.y = element_line(size = 2),
                  plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) +
            ggtitle(paste0("Compound: ", cpd, "; Gene: ", gene))
    })
    
    # get drug metrics for selected row
    get_drug_metrics <- reactive({
        ct_summary_df <- get_ct_summary_df()
        ct <- get_ct()
        if (nrow(ct_summary_df) < 1) return(NULL)
        metrics <- data$ctd_metrics %>%
            dplyr::filter(treatmentid == unique(ct_summary_df$Compound))
        if (is_empty(metrics)) return(NULL)
        if (ct == "B-cell") {
            fctlvls <- rev(c("Solid","B-cell","T-cell","Myeloma","Myeloid"))
            metrics <- metrics %>% 
                mutate(disease = factor(disease, levels = fctlvls))
        }
        if (ct == "T-cell") {
            fctlvls <- rev(c("Solid","T-cell","B-cell","Myeloma","Myeloid"))
            metrics <- metrics %>% 
                mutate(disease = factor(disease, levels = fctlvls))
        }
        return(metrics)
    })
    
    # subtype factor level formatting
    get_drug_subtype_metrics <- reactive({
        plotdat <- get_drug_metrics() 
        ct <- get_ct()
        if (is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            dplyr::filter(disease %in% c(ct, "Solid")) %>%
            mutate(st = ifelse(disease == ct,
                               as.character(subtype), as.character(disease))) %>%
            dplyr::filter(!is.na(st))
        if (ct == "B-cell") {
            plotdat <- plotdat %>%
                mutate(st = factor(st, levels = rev(c("Solid","B-ALL","CLL","Mantle cell",
                                                      "Burkitts","DLBCL","other NHL", "Hodgkin"))))
        }
        if (ct == "T-cell") {
            plotdat <- plotdat %>%
                mutate(st = factor(st, levels = rev(c("Solid","T-ALL","T-cell lymphoma"))))
        }
        return(plotdat)
    })
    
    # DSS1 density plot by disease
    output$dss_dens <- renderPlot({
        plotdat <- get_drug_metrics() 
        ct <- get_ct()
        if (is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            dplyr::filter(disease %in% c(ct, "Solid")) %>%
            droplevels()
        gen_densplot(plotdat, "DSS1", "disease", xlab = "DSS1 score")
    })
    
    # DSS1 boxplot by disease
    output$dss_box_ds <- renderPlot({
        plotdat <- get_drug_metrics() 
        if (is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            dplyr::filter(!is.na(disease)) %>%
            droplevels()
        gen_boxplot(plotdat, "disease", "DSS1", ylab = "DSS1 score") 
    })
    
    # DSS1 boxplot by subtype
    output$dss_box_st <- renderPlot({
        plotdat <- get_drug_subtype_metrics() %>%
            droplevels()
        if (is_empty(plotdat)) return(NULL)
        gen_boxplot(plotdat, "st", "DSS1", ylab = "DSS1 score") 
    })
    
    # EC50 boxplot by subtype
    output$ec50_plot <- renderPlot({
        plotdat <- get_drug_subtype_metrics() %>%
            droplevels()
        if (is_empty(plotdat)) return(NULL)
        gen_boxplot(plotdat, "st", "EC50", ylab = "EC50 (\U03BCM)") +
            scale_y_log10(labels=plain) #+
            # theme(axis.text.x = element_text(angle=45, hjust=1))
    })
    
    # get crispr data
    get_crispr_dat <- reactive({
        ct_summary_df <- get_ct_summary_df()
        gene <- unique(ct_summary_df$Gene)
        geneidx <- which(fData(data$crispr_es)$genesymbol == gene)
        df <- data.frame(depmap_id = pData(data$crispr_es)$DepMap_ID,
                         cell_line = pData(data$crispr_es)$stripped_cell_line_name,
                         disease = pData(data$crispr_es)$disease,
                         subtype = pData(data$crispr_es)$subtype,
                         crispr_effect = exprs(data$crispr_es)[geneidx, ]) %>%
            mutate(disease = factor(disease, levels = levels(data$crispr_es$disease)),
                   subtype = factor(subtype, levels = levels(data$crispr_es$subtype)))
        if (is_empty(df)) return(NULL)
        return(df)
    })
    
    # CRISPR density plot by disease
    output$crispr_dens <- renderPlot({
        ct <- get_ct()
        plotdat <- get_crispr_dat() %>%
            dplyr::filter(disease %in% c(ct, "Solid")) %>%
            droplevels()
        gen_densplot(plotdat, "crispr_effect", "disease", xlab = "CRISPR effect score")
    })
    
    # CRISPR boxplot by disease
    output$crispr_box_ds <- renderPlot({
        plotdat <- get_crispr_dat() %>%
            dplyr::filter(!is.na(disease)) %>%
            droplevels()
        gen_boxplot(plotdat, "disease", "crispr_effect", ylab = "CRISPR effect score")
    })
    
    # CRISPR boxplot by subtype
    output$crispr_box_st <- renderPlot({
        ct <- get_ct()
        plotdat <- get_crispr_dat() %>%
            dplyr::filter(disease %in% c(ct, "Solid")) %>%
            mutate(st = ifelse(disease == ct,
                               as.character(subtype), as.character(disease))) %>%
            mutate(st = factor(st, levels = c(levels(data$crispr_es$subtype), "Solid"))) %>%
            dplyr::filter(!is.na(st)) %>%
            droplevels()
        gen_boxplot(plotdat, "st", "crispr_effect", ylab = "CRISPR effect score")
    })
    
    # get RNAi data
    get_rnai_dat <- reactive({
        ct_summary_df <- get_ct_summary_df()
        gene <- unique(ct_summary_df$Gene)
        geneidx <- which(fData(data$rnai_es)$genesymbol == gene)
        data.frame(depmap_id = pData(data$rnai_es)$DepMap_ID,
                   cell_line = pData(data$rnai_es)$stripped_cell_line_name ,
                   disease = pData(data$rnai_es)$disease,
                   subtype = pData(data$rnai_es)$subtype,
                   rnai_effect = exprs(data$rnai_es)[geneidx, ]) %>%
            mutate(disease = factor(disease, levels = levels(data$rnai_es$disease)),
                   subtype = factor(subtype, levels = levels(data$rnai_es$subtype)))
    })
    
    # RNAi density plot by disease
    output$rnai_dens <- renderPlot({
        ct <- get_ct()
        plotdat <- get_rnai_dat() %>%
            dplyr::filter(disease %in% c(ct, "Solid")) %>%
            droplevels()
        gen_densplot(plotdat, "rnai_effect", "disease", xlab = "RNAi effect score")
    })
    
    # RNAi boxplot by disease
    output$rnai_box_ds <- renderPlot({
        plotdat <- get_rnai_dat() %>%
            dplyr::filter(!is.na(disease)) %>%
            droplevels()
        gen_boxplot(plotdat, "disease", "rnai_effect", ylab = "RNAi effect score")
    })
    
    # RNAi boxplot by subtype
    output$rnai_box_st <- renderPlot({
        ct <- get_ct()
        plotdat <- get_rnai_dat() %>%
            dplyr::filter(disease %in% c(ct, "Solid")) %>%
            mutate(st = ifelse(disease == ct,
                               as.character(subtype), as.character(disease))) %>%
            mutate(st = factor(st, levels = c(levels(data$rnai_es$subtype), "Solid"))) %>%
            dplyr::filter(!is.na(st)) %>%
            droplevels()
        gen_boxplot(plotdat, "st", "rnai_effect", ylab = "RNAi effect score")
    })
    
    # cell-line level drug data
    output$cl_drugdat <- DT::renderDataTable({
        ct_summary_df <- get_ct_summary_df()
        ct <- get_ct()
        dat <- data$ctd_metrics %>%
            dplyr::filter(treatmentid == unique(ct_summary_df$Compound)) %>%
            dplyr::filter(disease %in% c(ct)) %>%
            select(EC50, DSS1, DSS2, DSS3, minc, maxc, slope, sampleid, treatmentid, disease, subtype)
        if (is.null(dat)) return(NULL)
        DT::datatable(
            data = dat, 
            rownames = F,
            options = list(
                columnDefs = list(list(
                    targets = 0:(ncol(dat)-1),
                    render = JS(
                        "function(data, type, row, meta) {",
                        "return type === 'display' && data.length > 30 ?",
                        "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                        "}"))),
                lengthMenu = c(5, 10, 25), pageLength = 5),
            callback = JS('table.page(3).draw(false);')
        )
    })
    
    # cell-line level drug data
    output$cl_rnaidat <- DT::renderDataTable({
        ct <- get_ct()
        dat <- get_rnai_dat() %>%
            dplyr::filter(disease %in% c(ct)) 
        if (is.null(dat)) return(NULL)
        DT::datatable(
            data = dat, 
            rownames = F,
            options = list(
                columnDefs = list(list(
                    targets = 0:(ncol(dat)-1),
                    render = JS(
                        "function(data, type, row, meta) {",
                        "return type === 'display' && data.length > 30 ?",
                        "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                        "}"))),
                lengthMenu = c(5, 10, 25), pageLength = 5),
            callback = JS('table.page(3).draw(false);')
        )
    })
    
    # cell-line level drug data
    output$cl_crisprdat <- DT::renderDataTable({
        ct <- get_ct()
        dat <- get_crispr_dat() %>%
            dplyr::filter(disease %in% c(ct)) 
        if (is.null(dat)) return(NULL)
        DT::datatable(
            data = dat, 
            rownames = F,
            options = list(
                columnDefs = list(list(
                    targets = 0:(ncol(dat)-1),
                    render = JS(
                        "function(data, type, row, meta) {",
                        "return type === 'display' && data.length > 30 ?",
                        "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                        "}"))),
                lengthMenu = c(5, 10, 25), pageLength = 5),
            callback = JS('table.page(3).draw(false);')
        )
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
