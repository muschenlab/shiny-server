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
library(RMariaDB)

################################### Setup #####################################

### Database config
dbname <- "DDDB"
cnf <- list.files("data", pattern = paste0(dbname, ".cnf$"), full.names = T)

# list celltype comparison choices
ctcomparisons <- list.files("data/ctres", "diffsens.rds", full.names = T)
names(ctcomparisons) <- sub("^.+/(.+)_diffsens.rds$", "\\1", ctcomparisons)
names(ctcomparisons) <- gsub("_", " ", names(ctcomparisons))

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
reactvals <- reactiveValues(gene = NULL, 
                            ct1 = NULL,
                            ct2 = NULL,
                            cpdlvl = NULL,
                            si = NULL,
                            selcpd_metric = "DSS4")

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

# print plain labels on log scale axis
plain <- function(x,...) {
    format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

# palette function
jellypal <- c("#714481cc","#712e7ccc","#9b6295cc",
              "#b36f8ecc","#b9869fcc","#d49aa2cc",
              "#a2b2a8cc","#a7c3c6cc","#7d8b99cc",
              "#a2b2a8cc","#a7c3c6cc","#7d8b99cc", # temp
              "#528199cc","#426284cc",
              "#eccad3cc","#eed9e0ff","#d0e6ea99")
pickcols <- function(pal, n) {
    idx <- c(round(seq(1, length(pal), length(pal)/(n-1))), length(pal))
    pal[idx]
}

##################################### UI ######################################
ui <- fluidPage(
    
    # title
    title = "Drug Discovery", 
    theme = shinytheme("cosmo"),
    titlePanel(tags$h2(tags$a(
        imageOutput("icon", inline = TRUE),
        href="http://137.184.200.69:3838/"), "Drug Discovery")),
    
    # change colors for DT row/column selection
    tags$style(HTML('table.dataTable tr.selected td{background-color: #B9869F99 !important;}')),
    tags$style(HTML('table.dataTable td.selected {background-color: #D0E6EA99 !important;}')),
    tags$style(HTML(".tabbable > .nav > li > a { color:#A7C4C6FF}")),
    tags$style(HTML(".tabbable > .nav > li[class=active] > a { color:black}")),
    
    # celltype comparison choice
    tags$h3("Celltype comparison:"),
    div(style = "font-size:13px;", 
        fluidRow(column(width=3,uiOutput("comparison_choice_ui")))),
    
    # main UI
    tabsetPanel(
        tabPanel("Compounds",
                 
                 # scatter plot of average
                 tags$h3("Differential sensitivity per compound:"),
                 tags$h5("Drag a box around points to select compounds of interest"),
                 fluidRow(align="center", 
                          plotOutput("cpd_scatter", height = 415, width = 550, 
                                     brush = "cpd_scatter_brush", 
                                     hover = "cpd_scatter_hover") %>% 
                            withSpinner(color = "#D0E6EA99")),
                 fluidRow(align="center", 
                          verbatimTextOutput("cpd_scatter_hover_text")),
                 
                 # metrics table
                 br(),br(),
                 tags$h3("Selected compounds:"),
                 tags$h5("Select a row to plot compound/target sensitivity accross cell lines below"),
                 div(DT::dataTableOutput("cpd_summary_dt"),
                     style = "font-size:90%"),
                 downloadButton("dl_cpd_summary_xls", label = "Download compounds summary",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 
                 # selected compound ranking info
                 br(),br(),
                 tags$h3("Selected compound:"),
                 fluidRow(align="center", uiOutput("sel_cpd_text")),

                 # rankings plot
                 tags$h4("Differential ranking:"),
                 tags$h5(paste0("Differential sensitivity ranks indicate how selectively sensitive the celltype of interest ",
                                "is relative to the control population (higher rank = more selective). Three primary metrics are ",
                                "considered - the CRISPR effect score, RNAi effect score, and drug sensitivity score (DSS).")),
                 fluidRow(align="center",
                          plotOutput("metrics_overview", height = 175, width = 400)
                 ),

                 # Compound metrics plots
                 br(),
                 tags$h3("Drug sensitivity metrics:"),
                 tags$h5(paste0("Drug sensitivity scores integrate the dose-response curve into a single metric - ",
                                "(higher score = more sensitive, range 0-100)")),
                 # CTD
                 fluidRow(splitLayout(cellWidths = c("20%","18%","27%","27%"),
                                      uiOutput("ctd_dss_dens_ui"),
                                      uiOutput("ctd_dss_group_ui"),
                                      uiOutput("ctd_dss_st_ui"),
                                      uiOutput("ctd_ec50_st_ui"))),
                 # GDSC
                 fluidRow(splitLayout(cellWidths = c("20%","19%","26%","26%"),
                                      uiOutput("gdsc_dss_dens_ui"),
                                      uiOutput("gdsc_dss_group_ui"),
                                      uiOutput("gdsc_dss_st_ui"),
                                      uiOutput("gdsc_ec50_st_ui"))),
                 # data download
                 div(DT::dataTableOutput("cpd_metrics_dt"),
                     style = "font-size:90%"),
                 downloadButton("dl_cpd_metrics_xls", label = "Download compound metrics",
                                style = "font-size:12px;height:30px;padding:5px;"),

                 # CRISPR plots
                 br(), br(),
                 tags$h3("Target Dependency:"),
                 tags$h5(paste0("CRISPR effect scores indicate the effect of gene knock-out on viability - ",
                                "(lower score = more sensitive)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "28%"),
                                      plotOutput("cpd_crispr_dens", height = 225,
                                                 brush = "crispr_dens_brush"),
                                      plotOutput("cpd_crispr_box_ds", height = 225),
                                      plotOutput("cpd_crispr_box_st", height = 225))),

                 # RNAi plots
                 tags$h4("RNAi:"),
                 tags$h5(paste0("RNAi effect scores indicate the effect of gene knock-down on viability - ",
                                "(lower scores = more sensitive)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "28%"),
                                      plotOutput("cpd_rnai_dens", height = 225),
                                      plotOutput("cpd_rnai_box_ds", height = 225),
                                      plotOutput("cpd_rnai_box_st", height = 225))),
                 
                 # dependency metrics for compound target
                 div(DT::dataTableOutput("cpd_dep_metrics_dt"),
                     style = "font-size:90%"),
                 downloadButton("dl_cpd_dep_metrics_xls", label = "Download dependency metrics",
                                style = "font-size:12px;height:30px;padding:5px;")
        ),
        tabPanel("Dependency",

                 tags$h3("Differential gene dependency scores:"),
                 tags$h5("Drag a box around points to select genes of interest"),
                 fluidRow(align="center",
                          splitLayout(cellWidths = c("30%", "45%"),
                                      tags$h4("CRISPR"),
                                      tags$h4("RNAi"))),
                 fluidRow(align="center",
                          splitLayout(cellWidths = c("35%", "45%"),
                                      plotOutput("crispr_scatter", height = 425,
                                                 brush = "crispr_scatter_brush",
                                                 hover = "crispr_scatter_hover") %>%
                                        withSpinner(color = "#D0E6EA99"),
                                      plotOutput("rnai_scatter", height = 425,
                                                 brush = "rnai_scatter_brush",
                                                 hover = "rnai_scatter_hover") %>%
                                        withSpinner(color = "#D0E6EA99"))),
                 fluidRow(align="center",
                          splitLayout(cellWidths = c("50%", "50%"),
                                      verbatimTextOutput("crispr_scatter_hover_text"),
                                      verbatimTextOutput("rnai_scatter_hover_text"))),

                 # gene dependency summary table
                 br(),br(),
                 tags$h3("Selected genes:"),
                 tags$h5("Select a row to plot gene dependency accross cell lines below"),
                 div(DT::dataTableOutput("genedep_summary_dt"),
                     style = "font-size:90%"),
                 downloadButton("dl_genedep_summary_xls", label = "Download dependency summary",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 
                 # selected gene ranking info
                 br(),br(),
                 tags$h3("Selected gene:"),
                 fluidRow(align="center", uiOutput("sel_gene_text")),
                 
                 # rankings plot
                 tags$h4("Differential ranking:"),
                 tags$h5(paste0("Differential sensitivity ranks indicate how selectively dependent the celltype of interest ",
                                "is relative to the control population (higher rank = more selective). Three primary metrics are ",
                                "considered - the CRISPR effect score, RNAi effect score, and drug sensitivity score (DSS).")),
                 fluidRow(align="center",
                          plotOutput("metrics_overview_dep", height = 175, width = 400)
                 ),
                 
                 # CRISPR plots
                 br(), br(),
                 tags$h3("CRISPR metrics:"),
                 tags$h5(paste0("CRISPR effect scores indicate the effect of gene knock-out on viability - ",
                                "(more negative = more dependent on gene)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "28%"),
                                      plotOutput("crispr_dens_dep", height = 225,
                                                 brush = "crispr_dens_brush"),
                                      plotOutput("crispr_box_ds_dep", height = 225),
                                      plotOutput("crispr_box_st_dep", height = 225))),
                 # RNAi plots
                 br(), br(),
                 tags$h3("RNAi metrics:"),
                 tags$h5(paste0("RNAi scores indicate the effect of gene knock-out on viability - ",
                                "(more negative = more dependent on gene)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "28%"),
                                      plotOutput("rnai_dens_dep", height = 225,
                                                 brush = "crispr_dens_brush"),
                                      plotOutput("rnai_box_ds_dep", height = 225),
                                      plotOutput("rnai_box_st_dep", height = 225)))
                 
        )#,
        # tabPanel("Expression",
        #          
        #          # scatter plot of average
        #          tags$h3("Average protein expression:"),
        #          tags$h5("Drag a box around points to select genes of interest"),
        #          fluidRow(align="center",
        #                   plotOutput("prot_scatter", height = 450, width = 450,
        #                              brush = "prot_scatter_brush",
        #                              hover = "prot_scatter_hover") %>%
        #                       withSpinner(color = "#D0E6EA99")),
        #          verbatimTextOutput("prot_hover_text"),
        #          
        #          # expr table
        #          br(),br(),
        #          tags$h3("Selected genes/proteins:"),
        #          tags$h5("Select a row to plot expression accross cell lines below"),
        #          div(DT::dataTableOutput("expr_summary_dt"),
        #              style = "font-size:90%"),
        #          
        #          # expr plots
        #          uiOutput("prot_expr_dens_ui")
        # ),
        # tabPanel("Prognosis",
        #          
        #          # scatter plot 
        #          tags$h3("Gene hazard ratios:"),
        #          tags$h5("Drag a box around points to select genes of interest"),
        #          fluidRow(align="center",
        #                   plotOutput("prog_scatter", height = 450, width = 450,
        #                              brush = "prog_scatter_brush",
        #                              hover = "prog_scatter_hover") %>%
        #                       withSpinner(color = "#D0E6EA99")),
        #          verbatimTextOutput("prog_scatter_hover_text"),
        #          
        #          # expr table
        #          br(),br(),
        #          tags$h3("Selected genes:"),
        #          tags$h5("Select a row to plot prognosis effects per subtype/study"),
        #          div(DT::dataTableOutput("prog_summary_dt"),
        #              style = "font-size:90%"),
        # 
        #          # expr plots
        #          br(),br(),
        #          tags$h3("HR effects:"),
        #          uiOutput("prog_forest_plot_ui")
        #          
        # )
    )
)


################################### Server ####################################
server <- function(input, output, session) {
    
    # icon
    output$icon <- renderImage(list(src = "../assets/images/hex/jelly.png",
                                    height = "95px", width = "85px"), 
                               deleteFile = F)
    
    # cell type comparison choice UI
    output$comparison_choice_ui <- renderUI({
        selectInput('comparison_choice', "Selected cell type comparison", 
                    names(ctcomparisons), names(ctcomparisons)[[1]])
    })
    
    # load ctcomp res
    load_comparison <- reactive({
      print("load_comparison")
      if (is_empty(input$comparison_choice)) return(NULL)
      ctcomp <- readRDS(ctcomparisons[[input$comparison_choice]])
      # ct/sample info
      reactvals$ct1 <- ct1 <- ctcomp$ct1
      reactvals$ct2 <- ct2 <- ctcomp$ct2
      reactvals$si <- ctcomp$si
      reactvals$ct1_si <- ctcomp$si %>%
        filter(ct==ct1)
      reactvals$ct2_si <- ctcomp$si %>%
        filter(ct==ct2)
      # cpd-level
      selcpd_metric <- reactvals$selcpd_metric 
      cpdlvl <- ctcomp$cpdlvl
      reactvals$cpdlvl <- cpdlvl %>%
        mutate(signif = -log10(cpdlvl[,grep(paste0("^p", selcpd_metric, "$"), colnames(cpdlvl))]),
               score = cpdlvl[,grep(paste0("^r", selcpd_metric, "$"), colnames(cpdlvl))])
      # gene-level
      reactvals$genelvl <- ctcomp$genelvl %>%
        mutate(nl10pCRISPR = -log10(pCRISPR),
               nl10pRNAi = -log10(pRNAi))
      return(NULL)
    })
    
    ################################ Compounds ################################
    
    # scatter plot differential compound scores
    output$cpd_scatter <- renderPlot({
      print("cpd_scatter")
      tmp <- load_comparison()
      if(is_empty(reactvals$cpdlvl)) return(NULL)
      plotdat <- reactvals$cpdlvl %>% dplyr::filter(!is.na(cpd_name))
      ct1 <- reactvals$ct1
      ct2 <- reactvals$ct2
      plotdat <- plotdat %>%
        mutate(color = ifelse(signif < 1.3, "ns",
                              ifelse(score > 0, "sensitive", "resistant"))) %>%
        mutate(color = factor(color, levels= c("sensitive", "ns", "resistant"))) %>%
        filter(!is.na(color))
      plotdat %>%
        ggplot(aes(x=score, y=signif)) +
        geom_vline(xintercept = 0, lty=2, alpha=0.4) +
        geom_point(aes(color=color), size = 2) +
        scale_color_manual(values=c("#71448199","grey60","#528199cc")) +
        scale_x_continuous(name = paste0("\u0394", reactvals$selcpd_metric,
                                         " ", ct1, "/", ct2)) +
        scale_y_continuous(name = "Significance (-log10 p-value)") +
        theme_bw(base_size = 17) +
        theme(legend.title = element_blank(),
              legend.background = element_blank(),
              panel.grid = element_blank()) 
    })
    
    # brushed cpd
    get_sel_cpd <- reactive({
        print("get_sel_cpd")
        cpdsub <- reactvals$cpdlvl
        if (is_empty(cpdsub)) return(NULL)
        cpdsub <- cpdsub %>%
          mutate_if(is.numeric, round, digits = 3)
        if (is_empty(input$cpd_scatter_brush)) return(cpdsub)
        cpdsub %>% 
          brushedPoints(input$cpd_scatter_brush)
    }) 
    # hovered cpd
    output$cpd_scatter_hover_text <- renderText({
      if (is_empty(input$cpd_scatter_hover)) return(NULL)
      cpdsub <- reactvals$cpdlvl %>%
        nearPoints(input$cpd_scatter_hover)
      paste0("Compounds near cursor: ", paste0(unique(cpdsub$cpd_name), collapse=";"))
    })
   
    # cpd summary table
    output$cpd_summary_dt <- DT::renderDataTable({
      print("cpd_summary_dt")
      cpdsub <- get_sel_cpd()
      if (is.null(cpdsub)) return(NULL)
      cpdsub <- cpdsub %>%
        dplyr::select(cpd_name, pubchem_cid, cpd_score, target_genes,
                      DSS4_score, CRISPR_score, RNAi_score,
                      dDSS4, rDSS4, pDSS4, starts_with("avDSS")) %>%
        dplyr::rename(Compound = cpd_name,
                      `Target genes` = target_genes,
                      `Compound score` = cpd_score,
                      `PubChem CID` = pubchem_cid)
      DT::datatable(
        data = cpdsub,
        rownames = F,
        colnames = gsub("_", " ", colnames(cpdsub)),
        selection = list(mode = 'single', target = "row", selected = 1)
      )
    })
    output$dl_cpd_summary_xls <- downloadHandler(
      filename = function() {
        paste0("Compounds_summary_table_", reactvals$ct1, "_vs_", 
               reactvals$ct2, ".xlsx")
      },
      content = function(file) {
        drug_summary <- get_sel_cpd()
        if(is_empty(drug_summary)) return(NULL)
        writexl::write_xlsx(drug_summary, path=file)
      }
    )
    
    ######### score summary for selected compound  ########### 
    
    # text just highlighting the selected compound/gene
    output$sel_cpd_text <- renderUI({
        print("sel_cpd_text")
        plotdat <- get_sel_cpd()
        if(is_empty(plotdat)) return(NULL)
        if(is_empty(input$cpd_summary_dt_rows_selected)) return(NULL)
        plotdat <- plotdat %>%
            dplyr::slice(input$cpd_summary_dt_rows_selected)
        HTML(paste0("<h4>Selected compound: <b>", plotdat$cpd_name,
                    "</b>, Targets: <b>", plotdat$target_genes, "</b></h4>"))
    })
    # metrics overview
    output$metrics_overview <- renderPlot({
        print("metrics_overview")
        plotdat <- get_sel_cpd()
        if(is_empty(plotdat)) return(NULL)
        if(is_empty(input$cpd_summary_dt_rows_selected)) return(NULL)
        plotdat <- plotdat %>%
            dplyr::slice(input$cpd_summary_dt_rows_selected) %>%
            dplyr::select(cpd_name, cpd_score, DSS4_score, CRISPR_score, RNAi_score) %>%
            gather("metric", "percentile", -c(1)) %>%
            mutate(percentile = percentile*100,
                   metric = factor(metric,
                                   levels = c("RNAi_score","CRISPR_score","DSS4_score",
                                              "cpd_score"),
                                   labels = c("RNAi score","CRISPR score","DSS4 score",
                                              "Overall compound score")))
        ggplot(plotdat, aes(x=metric, y=percentile)) +
            geom_point(size=4, color="#714481ff", shape=21) +
            coord_flip(clip="off") +
            scale_y_continuous(limits=c(0,100), expand=c(0,0), name="Percentile") +
            xlab("Metric") +
            theme_bw(base_size = 14) +
            theme(panel.border = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.line.x = element_line(),
                  panel.grid.major.y = element_line(size = 2),
                  plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
    })

    # get selected compound
    get_sel_cpd_info <- reactive({
      print("get_sel_cpd_info")
      cpdsummary <- get_sel_cpd()
      if(is.null(input$cpd_summary_dt_rows_selected) |
         is.null(cpdsummary)) return(NULL)
      cpdsummary <- cpdsummary %>%
        dplyr::slice(input$cpd_summary_dt_rows_selected)
      if (is_empty(cpdsummary)) return(NULL)
      cpd_name <- unique(cpdsummary$cpd_name)
      print(cpd_name)
      reactvals$cpd_name <- cpd_name
      gene_symbols <- unique(unlist(strsplit(cpdsummary$target_genes, ";")))
      orgdb <- org.Hs.eg.db::org.Hs.eg.db
      geneinfo <- AnnotationDbi::select(orgdb, keys = gene_symbols,
                                        keytype = "SYMBOL", column = "ENTREZID")
      print(geneinfo)
      reactvals$target_gene_symbols <- gene_symbols
      reactvals$target_gene_ids <- unique(geneinfo$ENTREZID)
      return(cpd_name)
    })

    ######### sensitivity metrics for selected compound  ########### 
    
    # get drug metrics for selected row
    get_sel_cpd_metrics <- reactive({
      print("get_sel_cpd_metrics")
      si <- reactvals$si
      ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
      cpdnames <- get_sel_cpd_info() # selected compound names
      if (is_empty(cpdnames)) return(NULL)
      # query db for selected compounds
      db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
      query <- paste0('SELECT * FROM drug_metrics ',
                      'WHERE cpd_name IN ("',
                      paste0(cpdnames, collapse = '", "'),
                      '") AND sample_id IN ("',
                      paste0(si$sample_id, collapse = '", "'),
                      '");')
      queryres <- dbSendQuery(db, query)
      metrics <- dbFetch(queryres)
      dbClearResult(queryres)
      dbDisconnect(db)
      # limit the EC50
      metrics <- metrics %>%
        mutate(EC50 = ifelse(logEC50 > logmaxc, logmaxc*2, logEC50)) 
      # format sample levels
      matchidx <- match(metrics$sample_id, si$sample_id)
      si$st <- ifelse(si$ct==ct1, si$ds_subtype, ct2)
      stlvls <- as.character(unique(si$st))
      stlvls <- c(stlvls[!grepl(ct2, stlvls)], ct2)
      si$dg <- ifelse(si$ct==ct1, ct1, si$ds_group)
      dglvls <- unique(si$dg)
      dglvls <- c(dglvls[!grepl(ct1, dglvls)], ct1)
      metrics <- metrics %>%
        mutate(ct = si[matchidx,]$ct,
               st = si[matchidx,]$st,
               dg = si[matchidx,]$dg) %>%
        mutate(ct = factor(ct, levels=c(ct1, ct2)),
               st = factor(st, levels=c(stlvls)),
               dg = factor(dg, levels=c(dglvls))) %>%
        dplyr::filter(!is.na(ct) & !is.na(dg) & !dg=="Other")
      if (is_empty(metrics) | nrow(metrics) < 1) return(NULL)
      return(metrics)
    })

    # drug sensitivity density plot
    output$ctd_dss_dens_ui <- renderUI({
      print("ctd_dss_dens_ui")
      plotdat <- get_sel_cpd_metrics()
      if (is_empty(plotdat)) return(NULL)
      plotdat <- plotdat %>%
        filter(dataset == "CTD")
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_densplot(plotdat, reactvals$selcpd_metric,
                        "ct", xlab = paste0(reactvals$selcpd_metric, " score")) +
        ggtitle("CTD")
      output$ctd_dss_dens_plot <- renderPlot(p)
      plotOutput("ctd_dss_dens_plot", height = 225)
    })
    output$gdsc_dss_dens_ui <- renderUI({
      print("gdsc_dss_dens_ui")
      plotdat <- get_sel_cpd_metrics()
      if (is_empty(plotdat)) return(NULL)
      plotdat <- plotdat %>%
        filter(dataset %in% c("GDSC1","GDSC2"))
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_densplot(plotdat, reactvals$selcpd_metric,
                        "ct", xlab = paste0(reactvals$selcpd_metric, " score")) +
        ggtitle("GDSC")
      output$gdsc_dss_dens_plot <- renderPlot(p)
      plotOutput("gdsc_dss_dens_plot", height = 225)
    })
    
    # drug sensitivity boxplots - disease groups
    output$ctd_dss_group_ui <- renderUI({
      print("ctd_dss_group_ui")  
      plotdat <- get_sel_cpd_metrics()
      if (is_empty(plotdat)) return(NULL)
      plotdat <- plotdat %>%
        filter(dataset == "CTD")
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "dg", reactvals$selcpd_metric,
                       ylab = paste0(reactvals$selcpd_metric, " score"))
      output$ctd_dss_group_plot <- renderPlot(p)
      plotOutput("ctd_dss_group_plot", height = 225)
    })
    output$gdsc_dss_group_ui <- renderUI({
      print("gdsc_dss_group_ui")    
      plotdat <- get_sel_cpd_metrics()
      if (is_empty(plotdat)) return(NULL)
      plotdat <- plotdat %>%
        filter(dataset %in% c("GDSC1","GDSC2"))
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "dg", reactvals$selcpd_metric,
                       ylab = paste0(reactvals$selcpd_metric, " score"))
      output$gdsc_dss_group_plot <- renderPlot(p)
      plotOutput("gdsc_dss_group_plot", height = 225)
    })

    # drug sensitivity boxplots - subtypes
    output$ctd_dss_st_ui <- renderUI({
      print("ctd_dss_st_ui")    
      plotdat <- get_sel_cpd_metrics()
      if (is_empty(plotdat)) return(NULL)
      plotdat <- plotdat %>%
        filter(dataset == "CTD")
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "st", reactvals$selcpd_metric,
                       ylab = paste0(reactvals$selcpd_metric, " score")) +
        ggtitle(paste0(reactvals$ct1, " subtypes:"))
      output$ctd_dss_st_plot <- renderPlot(p)
      plotOutput("ctd_dss_st_plot", height = 225)
    })
    output$gdsc_dss_st_ui <- renderUI({
      print("gdsc_dss_st_ui") 
      plotdat <- get_sel_cpd_metrics()
      if (is_empty(plotdat)) return(NULL)
      plotdat <- plotdat %>%
        filter(dataset %in% c("GDSC1","GDSC2"))
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "st", reactvals$selcpd_metric,
                       ylab = paste0(reactvals$selcpd_metric, " score")) +
        ggtitle(paste0(reactvals$ct1, " subtypes:"))
      output$gdsc_dss_st_plot <- renderPlot(p)
      plotOutput("gdsc_dss_st_plot", height = 225)
    })

    # EC50 boxplots - subtypes
    output$ctd_ec50_st_ui <- renderUI({
      print("ctd_ec50_st_ui")
      plotdat <- get_sel_cpd_metrics()
      if (is_empty(plotdat)) return(NULL)
      plotdat <- plotdat %>%
        filter(dataset == "CTD")
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "st", "EC50", ylab = "EC50 (\U03BCM)") +
        ggtitle("")
      output$ctd_ec50_st_plot <- renderPlot(p)
      plotOutput("ctd_ec50_st_plot", height = 225)
    })
    output$gdsc_ec50_st_ui <- renderUI({
      print("gdsc_ec50_st_ui")
      plotdat <- get_sel_cpd_metrics()
      if (is_empty(plotdat)) return(NULL)
      plotdat <- plotdat %>%
        filter(dataset %in% c("GDSC1","GDSC2"))
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "st", "EC50", ylab = "EC50 (\U03BCM)")
      output$gdsc_ec50_st_plot <- renderPlot(p)
      plotOutput("gdsc_ec50_st_plot", height = 225)
    })

    # data table for selected cpd
    output$cpd_metrics_dt <- DT::renderDataTable({
        metrics <- get_sel_cpd_metrics()
        if(is_empty(metrics) | is_empty(reactvals$ct1)) return(NULL)
        metrics <- metrics %>%
          dplyr::select(sample_id, treatment_id, cpd_name, ## add cell line annotation here also
                        logEC50, DSS1, DSS2, DSS3, DSS4) %>%
        DT::datatable(data = metrics,
                      rownames = F,
                      options = list(lengthMenu = c(5, 10, 25),
                                     pageLength = 5))
    })
    # DT download buttons
    output$dl_cpd_metrics_xls <- downloadHandler(
        filename = function() {
            paste0("Compound_metrics_", reactvals$cpd_name,
                   "_", reactvals$ct1,  ".xlsx")
        },
        content = function(file) {
            metrics <- get_sel_cpd_metrics()
            if (is_empty(metrics) | is_empty(reactvals$ct1)) return(NULL)
            writexl::write_xlsx(metrics, path=file)
        })

    ######### dependency metrics for selected compound targets ########### 
    
    # get dependency metrics for selected compound targets
    get_dep_metrics_cpd <- reactive({
      print("get_dep_metrics_cpd")
      si <- reactvals$si
      ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
      entrezids <- reactvals$target_gene_ids
      if (is_empty(entrezids)) return(NULL)
      # query db for selected compounds
      db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
      query <- paste0('SELECT * FROM dep_data ',
                      'WHERE entrez_id IN ("',
                      paste0(entrezids, collapse = '", "'),
                      '") AND sample_id IN ("',
                      paste0(si$sample_id, collapse = '", "'),
                      '");')
      queryres <- dbSendQuery(db, query)
      depdata <- dbFetch(queryres)
      dbClearResult(queryres)
      dbDisconnect(db)
      # format sample levels
      depdata <- merge(depdata, si, by="sample_id")
      depdata$st <- ifelse(depdata$ct==ct1, depdata$ds_subtype, ct2)
      stlvls <- as.character(unique(depdata$st))
      stlvls <- c(stlvls[!grepl(ct2, stlvls)], ct2)
      depdata$dg <- ifelse(depdata$ct==ct1, ct1, depdata$ds_group)
      dglvls <- unique(depdata$dg)
      dglvls <- c(dglvls[!grepl(ct1, dglvls)], ct1)
      depdata <- depdata %>%
        mutate(ct = factor(ct, levels=c(ct1, ct2)),
               st = factor(st, levels=c(stlvls)),
               dg = factor(dg, levels=c(dglvls))) %>%
        dplyr::filter(!is.na(ct) & !is.na(dg) & !dg=="Other")
      if (is_empty(depdata) | nrow(depdata) < 1) return(NULL)
      return(depdata)
    })
    
    ### CRISPR
    # CRISPR density plot by disease
    output$cpd_crispr_dens <- renderPlot({
        print("crispr_dens")
        plotdat <- get_dep_metrics_cpd()
        validate(need(!is_empty(plotdat),
                      "      No CRISPR data associated with selected compound"))
        plotdat <- plotdat %>%
          dplyr::filter(assay_type == "CRISPR")
        validate(need(nrow(plotdat) > 1,
                      "      No CRISPR data associated with selected compound"))
        gen_densplot(plotdat, "score", "ct", xlab = "CRISPR effect score")
    })
    # CRISPR boxplot by disease
    output$cpd_crispr_box_ds <- renderPlot({
      print("crispr_box_ds")
      plotdat <- get_dep_metrics_cpd()
      validate(need(!is_empty(plotdat),
                    "      No CRISPR data associated with selected compound"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "CRISPR")
      validate(need(nrow(plotdat) > 1,
                    "      No CRISPR data associated with selected compound"))
      gen_boxplot(plotdat, "dg", "score", ylab = "CRISPR effect score")
    })
    # CRISPR boxplot by subtype
    output$cpd_crispr_box_st <- renderPlot({
      print("crispr_box_ds")
      plotdat <- get_dep_metrics_cpd()
      validate(need(!is_empty(plotdat),
                    "      No CRISPR data associated with selected compound"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "CRISPR")
      gen_boxplot(plotdat, "st", "score", ylab = "CRISPR effect score")
    })
    
    ### RNAi
    # RNAi density plot by disease
    output$cpd_rnai_dens <- renderPlot({
      print("crispr_dens")
      plotdat <- get_dep_metrics_cpd()
      validate(need(!is_empty(plotdat),
                    "      No CRISPR data associated with selected compound"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "RNAi")
      validate(need(nrow(plotdat) > 1,
                    "      No RNAi data associated with selected compound"))
      gen_densplot(plotdat, "score", "ct", xlab = "RNAi effect score")
    })
    # RNAi boxplot by disease
    output$cpd_rnai_box_ds <- renderPlot({
      print("rnai_box_ds")
      plotdat <- get_dep_metrics_cpd()
      validate(need(!is_empty(plotdat),
                    "      No RNAi data associated with selected compound"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "RNAi")
      validate(need(nrow(plotdat) > 1,
                    "      No RNAi data associated with selected compound"))
      gen_boxplot(plotdat, "dg", "score", ylab = "RNAi effect score")
    })
    # RNAi boxplot by subtype
    output$cpd_rnai_box_st <- renderPlot({
      print("rnai_box_st")
      plotdat <- get_dep_metrics_cpd()
      validate(need(!is_empty(plotdat),
                    "      No RNAi data associated with selected compound"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "RNAi")
      gen_boxplot(plotdat, "st", "score", ylab = "RNAi effect score")
    })
    
    # data table of dependency data
    output$cpd_dep_metrics_dt <- DT::renderDataTable({
        dat <- get_dep_metrics_cpd() %>%
          dplyr::select_at(c(1:5,8,10:17))
        if (is_empty(dat) | is_empty(reactvals$ct1)) return(NULL)
        # if (!is_empty(input$crispr_dens_brush)) {
        #     brushinfo <- input$crispr_dens_brush
        #     dat <- dat %>% dplyr::filter(score > brushinfo$xmin &
        #                                  score < brushinfo$xmax)
        # }
        DT::datatable(
            data = dat,
            rownames = F,
            selection = "none",
            options = list(lengthMenu = c(5, 10, 25), pageLength = 5))
    })
    output$dl_cpd_dep_metrics_xls <- downloadHandler(
      filename = function() {
        paste0("CRISPR_metrics_", 
               paste0(reactvals$target_gene_symbols, collapse="_"),
               "_", reactvals$ct1,  ".xlsx")
      },
      content = function(file) {
        metrics <- get_dep_metrics_cpd()
        if (is_empty(metrics)) return(NULL)
        writexl::write_xlsx(metrics, path=file)
      }
    )
    
    ############################## Dependency ####################################
     
    # scatter plot of differential CRISPR scores
    output$crispr_scatter <- renderPlot({
      print("crispr_scatter")
      tmp <- load_comparison()
      if(is_empty(reactvals$genelvl)) return(NULL)
      plotdat <- reactvals$genelvl
      ct1 <- reactvals$ct1
      ct2 <- reactvals$ct2
      plotdat <- plotdat %>%
        mutate(color = ifelse(nl10pCRISPR < 1.3, "ns",
                              ifelse(dCRISPR > 0, "dependency", "tolerance"))) %>%
        # filter(color != "ns") %>%
        mutate(color = factor(color, levels= c("dependency", "ns", "tolerance")))
      plotdat %>%
        ggplot(aes(x=dCRISPR, y=nl10pCRISPR)) +
        geom_vline(xintercept = 0, lty=2, alpha=0.4) +
        geom_point(aes(color=color), size = 2) +
        scale_color_manual(values=c("#71448199","grey80","#528199cc")) +
        scale_x_continuous(name = paste0("\u0394 CRISPR score ", ct1, "/", ct2)) +
        scale_y_continuous(name = "Significance (-log10 p-value)") +
        theme_bw(base_size = 17) +
        theme(legend.position = "none",
              panel.grid = element_blank())
    })
    # scatter plot of differential RNAi scores
    output$rnai_scatter <- renderPlot({
      print("rnai_scatter")
      tmp <- load_comparison()
      if(is_empty(reactvals$genelvl)) return(NULL)
      plotdat <- reactvals$genelvl
      ct1 <- reactvals$ct1
      ct2 <- reactvals$ct2
      plotdat <- plotdat %>%
        mutate(color = ifelse(nl10pRNAi < 1.3, "ns",
                              ifelse(dRNAi > 0, "dependency", "tolerance"))) %>%
        # filter(color != "ns") %>%
        mutate(color = factor(color, levels= c("dependency", "ns", "tolerance")))
      plotdat %>%
        ggplot(aes(x=dRNAi, y=nl10pRNAi)) +
        geom_vline(xintercept = 0, lty=2, alpha=0.4) +
        geom_point(aes(color=color), size = 2) +
        scale_color_manual(values=c("#71448199","grey80","#528199cc"), name=paste0(ct1, "-specific")) +
        scale_x_continuous(name = paste0("\u0394 RNAi score ", ct1, "/", ct2)) +
        scale_y_continuous(name = "Significance (-log10 p-value)") +
        theme_bw(base_size = 17) +
        theme(legend.background = element_blank(),
              panel.grid = element_blank())
    })

    # handle brushing
    brushed_crispr <- reactive({
      print("brushed_crispr")
      if (is_empty(input$crispr_scatter_brush)) return(NULL)
      reactvals$rnai_scatter_brush <- F
      reactvals$crispr_scatter_brush <- T
    })
    brushed_rnai <- reactive({
      print("brushed_rnai")
      if (is_empty(input$rnai_scatter_brush)) return(NULL)
      reactvals$rnai_scatter_brush <- T
      reactvals$crispr_scatter_brush <- F
    })

    # get all selected compounds
    get_sel_genedep <- reactive({
      print("get_sel_genedep")
      genelvl <- reactvals$genelvl
      if(is_empty(genelvl)) return(NULL)
      genesub <- genelvl
      crispr_scatter_brush <- brushed_crispr()
      rnai_scatter_brush <- brushed_rnai()
      if (!is_empty(reactvals$crispr_scatter_brush)) {
        if (reactvals$rnai_scatter_brush == T) {
          genesub <- genesub %>%
            brushedPoints(input$rnai_scatter_brush)
        } else if (reactvals$crispr_scatter_brush == T) {
          genesub <- genesub %>%
            brushedPoints(input$crispr_scatter_brush)
        }
      } # else {
      #   genesub <- genesub %>%
      #     filter((gene_rank > .9 | gene_rank < .9) |
      #             (dCRISPR_rank > .9 | dCRISPR_rank < .9) |
      #              (dRNAi_rank > .9 | dRNAi_rank < .9))
      # }
      genesub %>%
        select_at(c(1:2, 4,6,7,5,8, 13:14, 16:17, 19:20, 21:22)) %>%
        mutate_at(c(3:5,7:9,12:15), round, digits=3)
    })

    # output text for hovered point - crispr
    output$crispr_hover_text <- renderText({
        if (is_empty(input$crispr_scatter_hover)) return(NULL)
        subdat <- get_sel_genedep() %>%
            nearPoints(input$crispr_scatter_hover)
        paste0("Genes near cursor: ", paste0(unique(subdat$gene_symbol), collapse=";"))
    })
    output$rnai_hover_text <- renderText({
      if (is_empty(input$rnai_scatter_hover)) return(NULL)
      subdat <- get_sel_genedep() %>%
        nearPoints(input$rnai_scatter_hover)
      paste0("Genes near cursor: ", paste0(unique(subdat$gene_symbol), collapse=";"))
    })

    # dependency info table
    output$genedep_summary_dt <- DT::renderDataTable({
        print("genedep_summary_dt")
        sel_genedep <- get_sel_genedep()
        if (is.null(sel_genedep)) return(NULL)
        DT::datatable(
            data = sel_genedep,
            rownames = F,
            selection = list(mode = 'single', target = "row", selected = 1),
            options = list(columnDefs = list(list(
              targets = 0:1,
              render = JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 30 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                "}")
            ))), callback = JS('table.page(3).draw(false);')
        )
    })
    output$dl_genedep_summary_xls <- downloadHandler(
      filename = function() {
        paste0("Gene_dependency_summary_table_", reactvals$ct1, "_vs_", 
               reactvals$ct2, ".xlsx")
      },
      content = function(file) {
        gene_summary <- get_sel_genedep()
        if(is_empty(gene_summary)) return(NULL)
        writexl::write_xlsx(gene_summary, path=file)
      }
    )
    
    ######### score summary for selected compound  ########### 
    
    # text just highlighting the selected gene
    output$sel_gene_text <- renderUI({
      print("sel_gene_text")
      plotdat <- get_sel_genedep()
      if(is_empty(plotdat)) return(NULL)
      if(is_empty(input$genedep_summary_dt_rows_selected)) return(NULL)
      plotdat <- plotdat %>%
        dplyr::slice(input$genedep_summary_dt_rows_selected)
      HTML(paste0("<h4>Selected gene: <b>", plotdat$gene_symbol,
                  "</b>, Compounds targetting: <b>", plotdat$compounds, "</b></h4>"))
    })
    # metrics overview
    output$metrics_overview_dep <- renderPlot({
      print("metrics_overview")
      plotdat <- get_sel_genedep()
      if(is_empty(plotdat)) return(NULL)
      if(is_empty(input$genedep_summary_dt_rows_selected)) return(NULL)
      plotdat <- plotdat %>%
        dplyr::slice(input$genedep_summary_dt_rows_selected) %>%
        dplyr::select(gene_symbol, gene_score, CRISPR_score, RNAi_score, DSS4_score) %>%
        gather("metric", "percentile", -c(1)) %>%
        mutate(percentile = percentile*100,
               metric = factor(metric,
                               levels = c("DSS4_score","RNAi_score","CRISPR_score",
                                          "gene_score"),
                               labels = c("DSS4 score","RNAi score","CRISPR score",
                                          "Overall gene score")))
      ggplot(plotdat, aes(x=metric, y=percentile)) +
        geom_point(size=4, color="#714481ff", shape=21) +
        coord_flip(clip="off") +
        scale_y_continuous(limits=c(0,100), expand=c(0,0), name="Percentile") +
        xlab("Metric") +
        theme_bw(base_size = 14) +
        theme(panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.x = element_line(),
              panel.grid.major.y = element_line(size = 2),
              plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
    })
    
    # get selected gene
    get_sel_gene_info <- reactive({
      print("get_sel_gene_info")
      depsummary <- get_sel_genedep()
      if(is.null(input$genedep_summary_dt_rows_selected) |
         is.null(depsummary)) { print("line 992"); return(NULL) }
      depsummary <- depsummary %>%
        dplyr::slice(input$genedep_summary_dt_rows_selected)
      if (is_empty(depsummary)) { print("line 996"); return(NULL) }
      print(head(depsummary))
      entrez_id <- unique(depsummary$entrez_id)
      gene_symbol <- unique(depsummary$gene_symbol)
      compounds <- unique(unlist(strsplit(depsummary$compounds, ";")))
      reactvals$dep_compounds <- compounds
      reactvals$dep_gene_symbol <- gene_symbol
      reactvals$dep_gene_id <- entrez_id
      return(entrez_id)
    })

    ######### dependency metrics for selected gene #############

    # get dependency metrics for selected compound targets
    get_dep_metrics_gene <- reactive({
      print("get_dep_metrics_gene")
      si <- reactvals$si
      ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
      entrezid <- get_sel_gene_info()
      if (is_empty(entrezid)) { print("line 1014"); return(NULL) }
      # query db for selected compounds
      db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
      query <- paste0('SELECT * FROM dep_data ',
                      'WHERE entrez_id IN ("',
                      paste0(entrezid, collapse = '", "'),
                      '") AND sample_id IN ("',
                      paste0(si$sample_id, collapse = '", "'),
                      '");')
      queryres <- dbSendQuery(db, query)
      depdata <- dbFetch(queryres)
      dbClearResult(queryres)
      dbDisconnect(db)
      # format sample levels
      depdata <- merge(depdata, si, by="sample_id")
      depdata$st <- ifelse(depdata$ct==ct1, depdata$ds_subtype, ct2)
      stlvls <- as.character(unique(depdata$st))
      stlvls <- c(stlvls[!grepl(ct2, stlvls)], ct2)
      depdata$dg <- ifelse(depdata$ct==ct1, ct1, depdata$ds_group)
      dglvls <- unique(depdata$dg)
      dglvls <- c(dglvls[!grepl(ct1, dglvls)], ct1)
      depdata <- depdata %>%
        mutate(ct = factor(ct, levels=c(ct1, ct2)),
               st = factor(st, levels=c(stlvls)),
               dg = factor(dg, levels=c(dglvls))) %>%
        dplyr::filter(!is.na(ct) & !is.na(dg) & !dg=="Other")
      if (is_empty(depdata) | nrow(depdata) < 1) return(NULL)
      return(depdata)
    })

    ### CRISPR
    # CRISPR density plot by disease
    output$crispr_dens_dep <- renderPlot({
      print("crispr_dens_dep")
      plotdat <- get_dep_metrics_gene()
      validate(need(!is_empty(plotdat),
                    "      No CRISPR data associated with selected compound"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "CRISPR")
      validate(need(nrow(plotdat) > 1,
                    "      No CRISPR data associated with selected compound"))
      gen_densplot(plotdat, "score", "ct", xlab = "CRISPR effect score")
    })
    # CRISPR boxplot by disease
    output$crispr_box_ds_dep <- renderPlot({
      print("crispr_box_ds_dep")
      plotdat <- get_dep_metrics_gene()
      validate(need(!is_empty(plotdat),
                    "      No CRISPR data associated with selected compound"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "CRISPR")
      validate(need(nrow(plotdat) > 1,
                    "      No CRISPR data associated with selected compound"))
      gen_boxplot(plotdat, "dg", "score", ylab = "CRISPR effect score")
    })
    # CRISPR boxplot by subtype
    output$crispr_box_st_dep <- renderPlot({
      print("crispr_box_st_dep")
      plotdat <- get_dep_metrics_gene()
      validate(need(!is_empty(plotdat),
                    "      No CRISPR data associated with selected gene"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "CRISPR")
      gen_boxplot(plotdat, "st", "score", ylab = "CRISPR effect score")
    })

    ### RNAi
    # RNAi density plot by disease
    output$rnai_dens_dep <- renderPlot({
      print("rnai_dens_dep")
      plotdat <- get_dep_metrics_gene()
      validate(need(!is_empty(plotdat),
                    "      No CRISPR data associated with selected gene"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "RNAi")
      validate(need(nrow(plotdat) > 1,
                    "      No RNAi data associated with selected gene"))
      gen_densplot(plotdat, "score", "ct", xlab = "RNAi effect score")
    })
    # RNAi boxplot by disease
    output$rnai_box_ds_dep <- renderPlot({
      print("rnai_box_ds_dep")
      plotdat <- get_dep_metrics_gene()
      validate(need(!is_empty(plotdat),
                    "      No RNAi data associated with selected gene"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "RNAi")
      validate(need(nrow(plotdat) > 1,
                    "      No RNAi data associated with selected gene"))
      gen_boxplot(plotdat, "dg", "score", ylab = "RNAi effect score")
    })
    # RNAi boxplot by subtype
    output$rnai_box_st_dep <- renderPlot({
      print("rnai_box_st_dep")
      plotdat <- get_dep_metrics_gene()
      validate(need(!is_empty(plotdat),
                    "      No RNAi data associated with selected gene"))
      plotdat <- plotdat %>%
        dplyr::filter(assay_type == "RNAi")
      gen_boxplot(plotdat, "st", "score", ylab = "RNAi effect score")
    })
    # 
    # # data table of dependency data
    # output$cpd_dep_metrics_dt <- DT::renderDataTable({
    #   dat <- get_dep_metrics_cpd() %>%
    #     dplyr::select_at(c(1:5,8,10:17))
    #   if (is_empty(dat) | is_empty(reactvals$ct1)) return(NULL)
    #   # if (!is_empty(input$crispr_dens_brush)) {
    #   #     brushinfo <- input$crispr_dens_brush
    #   #     dat <- dat %>% dplyr::filter(score > brushinfo$xmin &
    #   #                                  score < brushinfo$xmax)
    #   # }
    #   DT::datatable(
    #     data = dat,
    #     rownames = F,
    #     selection = "none",
    #     options = list(lengthMenu = c(5, 10, 25), pageLength = 5))
    # })
    # output$dl_cpd_dep_metrics_xls <- downloadHandler(
    #   filename = function() {
    #     paste0("CRISPR_metrics_", 
    #            paste0(reactvals$target_gene_symbols, collapse="_"),
    #            "_", reactvals$ct1,  ".xlsx")
    #   },
    #   content = function(file) {
    #     metrics <- get_dep_metrics_cpd()
    #     if (is_empty(metrics)) return(NULL)
    #     print(head(metrics))
    #     writexl::write_xlsx(metrics, path=file)
    #   }
    # )
    # 
    
    
    # ######### sensitivity metrics for selected compound  ########### 
    # 
    # # get drug metrics for selected row
    # get_sel_cpd_metrics <- reactive({
    #   print("get_sel_cpd_metrics")
    #   si <- reactvals$si
    #   ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
    #   cpdnames <- get_sel_cpd_info() # selected compound names
    #   if (is_empty(cpdnames)) return(NULL)
    #   # query db for selected compounds
    #   db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
    #   query <- paste0('SELECT * FROM drug_metrics ',
    #                   'WHERE cpd_name IN ("',
    #                   paste0(cpdnames, collapse = '", "'),
    #                   '") AND sample_id IN ("',
    #                   paste0(si$sample_id, collapse = '", "'),
    #                   '");')
    #   queryres <- dbSendQuery(db, query)
    #   metrics <- dbFetch(queryres)
    #   dbClearResult(queryres)
    #   dbDisconnect(db)
    #   # limit the EC50
    #   metrics <- metrics %>%
    #     mutate(EC50 = ifelse(logEC50 > logmaxc, logmaxc*2, logEC50)) 
    #   # format sample levels
    #   matchidx <- match(metrics$sample_id, si$sample_id)
    #   si$st <- ifelse(si$ct==ct1, si$ds_subtype, ct2)
    #   stlvls <- as.character(unique(si$st))
    #   stlvls <- c(stlvls[!grepl(ct2, stlvls)], ct2)
    #   si$dg <- ifelse(si$ct==ct1, ct1, si$ds_group)
    #   dglvls <- unique(si$dg)
    #   dglvls <- c(dglvls[!grepl(ct1, dglvls)], ct1)
    #   metrics <- metrics %>%
    #     mutate(ct = si[matchidx,]$ct,
    #            st = si[matchidx,]$st,
    #            dg = si[matchidx,]$dg) %>%
    #     mutate(ct = factor(ct, levels=c(ct1, ct2)),
    #            st = factor(st, levels=c(stlvls)),
    #            dg = factor(dg, levels=c(dglvls))) %>%
    #     dplyr::filter(!is.na(ct) & !is.na(dg) & !dg=="Other")
    #   if (is_empty(metrics) | nrow(metrics) < 1) return(NULL)
    #   return(metrics)
    # })
    # 
    # # drug sensitivity density plot
    # output$ctd_dss_dens_ui <- renderUI({
    #   print("ctd_dss_dens_ui")
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset == "CTD")
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_densplot(plotdat, reactvals$selcpd_metric,
    #                     "ct", xlab = paste0(reactvals$selcpd_metric, " score")) +
    #     ggtitle("CTD")
    #   output$ctd_dss_dens_plot <- renderPlot(p)
    #   plotOutput("ctd_dss_dens_plot", height = 225)
    # })
    # output$gdsc_dss_dens_ui <- renderUI({
    #   print("gdsc_dss_dens_ui")
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset %in% c("GDSC1","GDSC2"))
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_densplot(plotdat, reactvals$selcpd_metric,
    #                     "ct", xlab = paste0(reactvals$selcpd_metric, " score")) +
    #     ggtitle("GDSC")
    #   output$gdsc_dss_dens_plot <- renderPlot(p)
    #   plotOutput("gdsc_dss_dens_plot", height = 225)
    # })
    # 
    # # drug sensitivity boxplots - disease groups
    # output$ctd_dss_group_ui <- renderUI({
    #   print("ctd_dss_group_ui")  
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset == "CTD")
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "dg", reactvals$selcpd_metric,
    #                    ylab = paste0(reactvals$selcpd_metric, " score"))
    #   output$ctd_dss_group_plot <- renderPlot(p)
    #   plotOutput("ctd_dss_group_plot", height = 225)
    # })
    # output$gdsc_dss_group_ui <- renderUI({
    #   print("gdsc_dss_group_ui")    
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset %in% c("GDSC1","GDSC2"))
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "dg", reactvals$selcpd_metric,
    #                    ylab = paste0(reactvals$selcpd_metric, " score"))
    #   output$gdsc_dss_group_plot <- renderPlot(p)
    #   plotOutput("gdsc_dss_group_plot", height = 225)
    # })
    # 
    # # drug sensitivity boxplots - subtypes
    # output$ctd_dss_st_ui <- renderUI({
    #   print("ctd_dss_st_ui")    
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset == "CTD")
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "st", reactvals$selcpd_metric,
    #                    ylab = paste0(reactvals$selcpd_metric, " score")) +
    #     ggtitle(paste0(reactvals$ct1, " subtypes:"))
    #   output$ctd_dss_st_plot <- renderPlot(p)
    #   plotOutput("ctd_dss_st_plot", height = 225)
    # })
    # output$gdsc_dss_st_ui <- renderUI({
    #   print("gdsc_dss_st_ui") 
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset %in% c("GDSC1","GDSC2"))
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "st", reactvals$selcpd_metric,
    #                    ylab = paste0(reactvals$selcpd_metric, " score")) +
    #     ggtitle(paste0(reactvals$ct1, " subtypes:"))
    #   output$gdsc_dss_st_plot <- renderPlot(p)
    #   plotOutput("gdsc_dss_st_plot", height = 225)
    # })
    # 
    # # EC50 boxplots - subtypes
    # output$ctd_ec50_st_ui <- renderUI({
    #   print("ctd_ec50_st_ui")
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset == "CTD")
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "st", "EC50", ylab = "EC50 (\U03BCM)") +
    #     ggtitle("")
    #   output$ctd_ec50_st_plot <- renderPlot(p)
    #   plotOutput("ctd_ec50_st_plot", height = 225)
    # })
    # output$gdsc_ec50_st_ui <- renderUI({
    #   print("gdsc_ec50_st_ui")
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset %in% c("GDSC1","GDSC2"))
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "st", "EC50", ylab = "EC50 (\U03BCM)")
    #   output$gdsc_ec50_st_plot <- renderPlot(p)
    #   plotOutput("gdsc_ec50_st_plot", height = 225)
    # })
    # 
    # # data table for selected cpd
    # output$cpd_metrics_dt <- DT::renderDataTable({
    #   metrics <- get_sel_cpd_metrics()
    #   if(is_empty(metrics) | is_empty(reactvals$ct1)) return(NULL)
    #   metrics <- metrics %>%
    #     dplyr::select(sample_id, treatment_id, cpd_name, ## add cell line annotation here also
    #                   logEC50, DSS1, DSS2, DSS3, DSS4) %>%
    #     DT::datatable(data = metrics,
    #                   rownames = F,
    #                   options = list(lengthMenu = c(5, 10, 25),
    #                                  pageLength = 5))
    # })
    # # DT download buttons
    # output$dl_cpd_metrics_xls <- downloadHandler(
    #   filename = function() {
    #     paste0("Compound_metrics_", reactvals$cpd_name,
    #            "_", reactvals$ct1,  ".xlsx")
    #   },
    #   content = function(file) {
    #     metrics <- get_sel_cpd_metrics()
    #     if (is_empty(metrics) | is_empty(reactvals$ct1)) return(NULL)
    #     writexl::write_xlsx(metrics, path=file)
    #   })
    # 
    

    # ############################## Expression ####################################
    # 
    # # av proteomics data
    # get_prot_av <- reactive({
    #     print("get_prot_av")
    #     ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
    #     prot <- data$sanger_prot_se
    #     idx1 <- which(!is.na(prot$Sanger_Model_ID) & prot$Sanger_Model_ID %in% reactvals$ct1_si$Sanger_Model_ID)
    #     idx2 <- which(!is.na(prot$Sanger_Model_ID) & prot$Sanger_Model_ID %in% reactvals$ct2_si$Sanger_Model_ID)
    #     av <- data.frame(genesymbol = rowData(prot)$name,
    #                      protein_id = rowData(prot)$ID,
    #                      avExpr_1 = rowMeans(assay(prot, "imputed")[, idx1]),
    #                      avExpr_2 = rowMeans(assay(prot, "imputed")[, idx2])) %>%
    #         mutate(dExpr = avExpr_1 - avExpr_2)
    #     colnames(av) <- sub("_1$", paste0("_",ct1), colnames(av))
    #     colnames(av) <- sub("_2$", paste0("_",ct2), colnames(av))
    #     return(av)
    # })
    # get_prot_dat <- reactive({
    #     print("get_prot_dat")
    #     ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
    #     prot <- data$sanger_prot_se
    #     idx1 <- which(!is.na(prot$Sanger_Model_ID) & prot$Sanger_Model_ID %in% reactvals$ct1_si$Sanger_Model_ID)
    #     idx2 <- which(!is.na(prot$Sanger_Model_ID) & prot$Sanger_Model_ID %in% reactvals$ct2_si$Sanger_Model_ID)
    #     prot$ct <- NA; prot$ct[idx2] <- ct2; prot$ct[idx1] <- ct1
    #     av <- get_prot_av()
    #     prot <- prot[, c(idx1, idx2)]
    #     protsub <- prot[abs(av$dExpr) > quantile(abs(av$dExpr), .5, na.rm=T), ]
    #     protdat <- assay(protsub, "imputed") %>%
    #         as.data.frame() %>%
    #         mutate(genesymbol = rowData(protsub)$name, .before=1) %>%
    #         gather("cl_id", "value", -1) %>%
    #         mutate(ct = colData(protsub)[match(.$cl_id, protsub$Sanger_Model_ID),]$ct) %>%
    #         dplyr::filter(!is.na(value))
    #     return(protdat)
    # })
    # get_prot_summary <- reactive({
    #     print("get_prot_summary")
    #     ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
    #     protdat <- get_prot_dat()
    #     av <- get_prot_av()
    #     av <- av[, -which(colnames(av)=="dExpr")]
    #     stats <- protdat %>%
    #         select(genesymbol, ct, value) %>%
    #         mutate(ct = factor(ct, levels=c(ct1,ct2))) %>%
    #         nest(data = -genesymbol) %>%
    #         mutate(dExpr = map_dbl(data, ~filt_ttest(.x, metric="value", stat="statistic")),
    #                pExpr = map_dbl(data, ~filt_ttest(.x, metric="value", stat="p.value"))) %>%
    #         mutate(qExpr = p.adjust(pExpr)) %>%
    #         mutate(dExpr_rank = rank(dExpr)/length(dExpr),
    #                .before=3) %>%
    #         select(-data)
    #     summary <- merge(stats, av, by = c("genesymbol"))
    #     summary <- summary %>% arrange(-dExpr_rank)
    #     return(summary)
    # })
    # 
    # # scatter plot of average sensitivity scores per compound
    # output$prot_scatter <- renderPlot({
    #     plotdat <- get_prot_summary()
    #     if (is_empty(plotdat)) return(NULL)
    #     tmp <- c(paste0("`avExpr_", reactvals$ct1, "`"),
    #              paste0("`avExpr_", reactvals$ct2, "`"))
    #     plotdat %>%
    #         ggplot(aes_string(x=tmp[[1]], y=tmp[[2]])) +
    #         geom_point(aes(color=ifelse(qExpr > 0.001, "B", "A")), size = 2) +
    #         geom_abline() +
    #         scale_x_continuous(name = paste0("Average protein level: ", reactvals$ct1)) +
    #         scale_y_continuous(name = paste0("Average protein level: ", reactvals$ct2)) +
    #         scale_color_manual(values=c("#71448199","#D0E6EA66")) +
    #         theme_bw(base_size = 18) +
    #         theme(legend.position = "none")
    # })
    # 
    # # get selected points from scatter
    # get_sel_expr <- reactive({
    #     print("get_sel_expr")
    #     subdat <- get_prot_summary() %>%
    #         mutate_if(grepl("avExpr|dExpr", names(.)), round, digits = 3)
    #     if (!is_empty(input$prot_scatter_brush)) {
    #         subdat <- subdat %>%
    #             brushedPoints(input$prot_scatter_brush)
    #     }
    #     return(subdat)
    # })
    # 
    # # output text for hovered
    # output$prot_hover_text <- renderText({
    #     if (is_empty(input$prot_scatter_hover)) return(NULL)
    #     subdat <- get_prot_summary() %>%
    #         nearPoints(input$prot_scatter_hover)
    #     paste0("Genes near cursor: ", paste0(unique(subdat$genesymbol), collapse=";"))
    # })
    # 
    # # expr info table
    # output$expr_summary_dt <- DT::renderDataTable({
    #     expr_summary_sel <- get_sel_expr() %>%
    #         dplyr::select(-"pExpr")
    #     if (is.null(expr_summary_sel)) return(NULL)
    #     DT::datatable(
    #         data = expr_summary_sel,
    #         rownames = F,
    #         selection = list(mode = 'multiple', target = "row", selected = 1),
    #         options = list(rowCallback = JS(
    #             "function(row, data, displayNum, index){",
    #             "  var x = data[3];",
    #             "  $('td:eq(3)', row).html(x.toExponential(2));",
    #             "}"
    #         ))
    #     )
    # })
    # 
    # # get selected gene(s)
    # get_sel_genes_expr <- reactive({
    #     print("get_sel_genes_expr")
    #     if(is_empty(input$expr_summary_dt_rows_selected)) return(NULL)
    #     expr_summary <- get_sel_expr() %>%
    #         dplyr::slice(input$expr_summary_dt_rows_selected)
    #     if (is_empty(expr_summary)) return(NULL)
    #     sel_genes <- unique(expr_summary$genesymbol)
    #     reactvals$sel_genes <- sel_genes
    #     return(sel_genes)
    # })
    # 
    # # get expr values for selected genes
    # get_sel_prot_levels <- reactive({
    #     print("get_sel_prot_levels")
    #     sel_genes <- get_sel_genes_expr()
    #     if (is_empty(sel_genes)) return(NULL)
    #     prot <- data$sanger_prot_se
    #     prot <- prot[which(rowData(prot)$name %in% sel_genes),]
    #     exprdat  <- assay(prot, "imputed") %>%
    #         as.data.frame() %>%
    #         mutate(genesymbol = rowData(prot)$name, .before=1) %>%
    #         gather("sampleid", "expr", -1)
    #     if (is_empty(exprdat) | nrow(exprdat) < 1) return(NULL)
    #     return(exprdat)
    # })
    # 
    # # ct filtering and formatting
    # format_prot_levels <- reactive({
    #     print("format_prot_levels_ct")
    #     exprdat <- get_sel_prot_levels()
    #     if (is_empty(exprdat)) return(NULL)
    #     ct1 <- exprdat %>%
    #         dplyr::filter(!is.na(sampleid) & sampleid %in% reactvals$ct1_si$Sanger_Model_ID) %>%
    #         mutate(ct = reactvals$ct1)
    #     ct2 <- exprdat %>%
    #         dplyr::filter(!is.na(sampleid) & sampleid %in% reactvals$ct2_si$Sanger_Model_ID) %>%
    #         mutate(ct = reactvals$ct2)
    #     plotdat <- rbind(ct1, ct2) %>%
    #         mutate(ct = factor(ct, levels=c(reactvals$ct1, reactvals$ct2)))
    #     matchidx <- match(plotdat$sampleid, data$si$Sanger_Model_ID)
    #     plotdat <- plotdat %>%
    #         mutate(ds_group = data$si[matchidx,]$ds_group,
    #                ds_subtype = data$si[matchidx,]$ds_subtype) %>%
    #         mutate(st = ifelse(ct == reactvals$ct1,
    #                            as.character(ds_subtype),
    #                            as.character(ct)))
    #     lvls <- unique(plotdat$st)
    #     lvls <- c(lvls[!grepl(reactvals$ct2, lvls)], reactvals$ct2)
    #     plotdat <- plotdat %>%
    #         mutate(st = factor(st, levels=lvls))
    #     return(plotdat)
    # })
    # 
    # # drug sensitivity density plot
    # output$prot_expr_dens_ui <- renderUI({
    #     print("prot_expr_dens_ui")
    #     plotdat <- format_prot_levels()
    #     if (is_empty(plotdat) | nrow(plotdat) < 1) return(NULL)
    #     sel_genes <- get_sel_genes_expr()
    #     if (length(sel_genes) == 1) {
    #         p <- gen_densplot(plotdat, "expr", "ct", xlab = "Protein levels [AU]") +
    #             ggtitle(paste0(sel_genes, " expression"))
    #         output$prot_expr_dens_plot <- renderPlot(p)
    #         fluidRow(splitLayout(cellWidths = c("20%"),
    #                              plotOutput("prot_expr_dens_plot", height = 225)))
    #     } else if (length(sel_genes) < 101) {
    #         ctav <- plotdat %>% 
    #             group_by(genesymbol, ct) %>%
    #             summarise(av=mean(expr, na.rm=T)) %>%
    #             ungroup() %>%
    #             pivot_wider(names_from = ct,
    #                         values_from = av) %>%
    #             mutate(diff = .[,2] - .[,3]) %>%
    #             arrange(diff)
    #         p <- plotdat %>% 
    #             group_by(sampleid, genesymbol) %>% 
    #             summarise(ct = ct,
    #                       av = mean(expr, na.rm=T)) %>% 
    #             group_by(genesymbol) %>% 
    #             mutate(z = scale(av)[,1]) %>%
    #             ungroup %>%
    #             mutate(genesymbol = factor(genesymbol, levels=unique(ctav$genesymbol))) %>%
    #             ggplot(aes(x=sampleid, y=genesymbol, fill=z)) + 
    #             geom_tile() +
    #             facet_wrap(~ct, nrow=1, scales="free_x") +
    #             theme_bw() +
    #             theme(axis.text.x = element_blank(),
    #                   axis.ticks = element_blank(),
    #                   panel.grid = element_blank(),
    #                   panel.background = element_blank(),
    #                   panel.border = element_rect(fill=NA),
    #                   strip.background = element_blank(),
    #                   strip.text = element_text(size=12)) +
    #             scale_fill_gradient2(low="#426284cc", mid="#d0e6ea99", high="#714481",
    #                                  name="Relative expression [AU]") +
    #             scale_y_discrete(expand=c(0,0)) +
    #             xlab("Cell lines") + ylab("")
    #         output$prot_expr_hm_plot <- renderPlot(p)
    #         plotheight <- 200 + (length(sel_genes)*15)
    #         fluidRow(plotOutput("prot_expr_hm_plot", height = plotheight, width = 1200))
    #     } else {
    #         tags$h5("Too many genes selected to plot")
    #     }
    # })
    # 
    # ############################# Prognosis #####################################
    # 
    # # subset prognosis data
    # get_prog_dat <- reactive({
    #     print("get_prog_dat")
    #     prog <- data$prognosis_metrics
    #     ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
    #     ds <- unique(prog$ds_type); st <- unique(prog$ds_subtype)
    #     
    #     # if (ct2 == "Solid tumor") {
    #     #     ct2dat <- prog %>% filter(dataset!="meta" & !ds_group %in% c("Lymphoid", "Myeloid")) %>% mutate(group = "ct2")
    #     # } else {
    #     #     if (is.na(cat2)) return(NULL)
    #     #     ct2dat <- prog %>% filter(dataset!="meta" & (!!as.symbol(cat2)) == ct2) %>% mutate(group = "ct2")
    #     # }
    #     cat1 <- ifelse(ct1 %in% ds, "ds_type", ifelse(ct1 %in% st, "ds_subtype", NA))
    #     cat2 <- ifelse(ct2 %in% ds, "ds_type", ifelse(ct2 %in% st, "ds_subtype", NA))
    #     if (is.na(cat1) | is.na(cat2)) return(NULL)
    #     ct1dat <- prog %>% filter(dataset!="meta" & (!!as.symbol(cat1)) == ct1) %>% mutate(group = "ct1")
    #     ct1meta <- prog %>% filter(dataset=="meta" & (!!as.symbol(cat1)) == ct1) %>% mutate(group = "ct1")
    #     ct2meta <- prog %>% filter(dataset=="meta" & (!!as.symbol(cat2)) == ct2) %>% mutate(group = "ct2")
    #     progdat <- rbind(ct1dat, ct1meta, ct2meta) %>% distinct()
    #     return(progdat)
    # })
    # 
    # # get prognosis summary data
    # get_prog_summary <- reactive({
    #     print("get_prog_summary")
    #     progdat <- get_prog_dat()
    #     if (is_empty(progdat)) return(NULL)
    #     progsummary <- progdat %>% 
    #         filter(dataset == "meta" & is.na(ds_subtype)) %>% 
    #         select(gene, n, coef, se_coef, group) %>%
    #         pivot_wider(names_from = group, values_from = c(n, coef, se_coef)) %>%
    #         group_by(gene) %>%
    #         summarise(coef_ct1 = coef_ct1,
    #                   coef_ct2 = coef_ct2,
    #                   se_coef_ct1 = se_coef_ct1,
    #                   se_coef_ct2 = se_coef_ct2,
    #                   dPrognosis = (coef_ct1 - coef_ct2) / sum(se_coef_ct1+se_coef_ct2)) %>%
    #         filter(!is.na(dPrognosis)) %>%
    #         mutate(prog_rank = rank(-dPrognosis)) %>%
    #         arrange(-prog_rank)
    #     return(progsummary)
    # })
    # 
    # # differential prognosis scatter
    # output$prog_scatter <- renderPlot({
    #     print("prog_scatter")
    #     progsum <- get_prog_summary()
    #     if (is_empty(progsum)) return(NULL)
    #     progsum %>%
    #         ggplot(aes(x=coef_ct1, y=coef_ct2)) +
    #         geom_point(aes(color=ifelse(abs(dPrognosis) > 2, "A", "B")), size = 2) +
    #         geom_abline() +
    #         scale_x_continuous(name = paste0("Gene hazard ratio (log): ", reactvals$ct1)) +
    #         scale_y_continuous(name = paste0("Gene hazard ratio (log): ", reactvals$ct2)) +
    #         scale_color_manual(values=c("#71448199","#D0E6EA66")) +
    #         theme_bw(base_size = 18) +
    #         theme(legend.position = "none")
    # })
    # 
    # # get all selected genes - prognosis
    # get_sel_progsummary <- reactive({
    #     print("get_sel_progsummary")
    #     progsum <- get_prog_summary()
    #     if(is_empty(progsum)) return(NULL)
    #     if (!is_empty(reactvals$gdsc_scatter_brush)) {
    #         progsum <- progsum %>%
    #             brushedPoints(input$prog_scatter_brush)
    #     }
    #     return(progsum)
    # })
    # 
    # # output text for hovered - prognosis
    # output$prog_scatter_hover_text <- renderText({
    #     if (is_empty(input$prog_scatter_hover)) return(NULL)
    #     progsum <- get_prog_summary() %>%
    #         nearPoints(input$prog_scatter_hover)
    #     paste0("Genes near cursor: ", paste0(unique(progsum$gene), collapse=";"))
    # })
    # 
    # # DT table - prognosis summary
    # output$prog_summary_dt <- DT::renderDataTable({
    #     print("prog_summary_dt")
    #     progsum <- get_sel_progsummary()
    #     if (is.null(progsum)) return(NULL)
    #     progsum <- progsum %>%
    #         select_at(c(1,7,6,2,4,3,5)) %>%
    #         mutate_at(c(3:7), ~round(.,3)) %>%
    #         dplyr::rename(Gene_symbol = gene,
    #                       dPrognosis_rank = prog_rank,
    #                       dPrognosis_score = dPrognosis,
    #                       logHR_ct1 = coef_ct1,
    #                       SE_ct1 = se_coef_ct1,
    #                       logHR_ct2 = coef_ct2,
    #                       SE_ct2 = se_coef_ct2)
    #     colnames(progsum) <- sub("ct1", reactvals$ct1, colnames(progsum))
    #     colnames(progsum) <- sub("ct2", reactvals$ct2, colnames(progsum))
    #     DT::datatable(
    #         data = progsum,
    #         rownames = F,
    #         colnames = gsub("_", " ", colnames(progsum)),
    #         selection = list(mode = 'single', target = "row", selected = 1))
    # })
    # # output$dl_prog_summary_xls <- downloadHandler(
    # #     filename = function() {
    # #         paste0("Prognosis_summary_table_", reactvals$ct1, "_vs_",
    # #                reactvals$ct2, ".xlsx")
    # #     },
    # #     content = function(file) {
    # #         drug_summary <- get_sel_progsummary()
    # #         if(is_empty(drug_summary)) return(NULL)
    # #         writexl::write_xlsx(drug_summary, path=file)
    # #     })
    # 
    # # get selected gene(s)
    # get_sel_gene_prog <- reactive({
    #     print("get_sel_gene_prog")
    #     if(is_empty(input$prog_summary_dt_rows_selected)) return(NULL)
    #     progsum <- get_sel_progsummary() %>%
    #         dplyr::slice(input$prog_summary_dt_rows_selected)
    #     if (is_empty(progsum)) return(NULL)
    #     selegene <- unique(progsum$gene)
    #     return(selegene)
    # })
    # 
    # # Forest plot of selected gene
    # output$prog_forest_plot_ui <- renderUI({
    #     progdat <- get_prog_dat()
    #     selgene <- get_sel_gene_prog()
    #     ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
    #     if (is_empty(progdat) | is_empty(selgene)) return(NULL)
    #     group_levels <- rev(c(ct2, ct1, unique(progdat$ds_subtype)))
    #     plotdat <- progdat %>%
    #         filter(gene == selgene) %>%
    #         mutate(col = ifelse(dataset == "meta", "meta", "study")) %>%
    #         mutate(grouping = ifelse(is.na(ds_subtype), ds_type, ds_subtype)) %>%
    #         mutate(grouping = factor(grouping, group_levels)) %>%
    #         mutate(label = paste0("HR:", round(exp_coef, 2), "\npval:", round(pval, 3)))
    #     textdat <- plotdat %>%
    #         filter(dataset == "meta")
    #     print(plotdat)
    #     xpos <- max(log(plotdat$upper_95ci)) + (max(log(plotdat$upper_95ci))/5)
    #     print(max(log(plotdat$upper_95ci)))
    #     p <- ggplot(plotdat) +
    #         geom_vline(xintercept = 0, lty=2, alpha=0.25) +
    #         geom_pointrange(data = filter(plotdat, group == "ct1" | dataset == "meta"),
    #                         aes(y = grouping, group = dataset, x = coef, color = col,
    #                             xmin = log(lower_95ci), xmax = log(upper_95ci)),
    #                         position = position_dodge(0.8), size = 0.9, linewidth = 0.9) +
    #         geom_text(data = textdat, x = xpos,
    #                   aes(y = grouping, label = label), 
    #                   size = 4.5, hjust = 0) +
    #         coord_cartesian(clip="off") +
    #         scale_color_manual(values = c( "#71448199", "#528199cc")) +
    #         ylab("") + xlab(paste0("logHR: high ", selgene, " expression")) +
    #         theme_bw(base_size = 16) +
    #         theme(panel.grid = element_blank(),
    #               legend.position = "none",
    #               plot.margin=unit(c(10,100,10,10), units = "pt"),
    #               panel.border = element_rect(linewidth=0.25),
    #               axis.ticks.y = element_blank())
    #     output$prog_forest_plot <- renderPlot(p)
    #     pheight <- 200 + (length(group_levels) * 40)
    #     fluidRow(align="center",
    #              plotOutput("prog_forest_plot", height = pheight, width = 600))
    # })

}

# Run the application 
shinyApp(ui = ui, server = server)
