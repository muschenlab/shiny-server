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

# dbconfig
db_name <- "DDDB"
cnf_file <- list.files("data", pattern = paste0(db_name, ".cnf$"), full.names = T)

##### functions
toptail <- function(x, n=100) {
  a <- head(x, n)
  b <- tail(x, n)
  rbind(a,b)
}
gen_pptx <- function(plot, file, height = 5, width = 5, left = 1, top = 1) {
  read_pptx() %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = dml(ggobj = plot), 
            location = ph_location(height = height, width = width,
                                   left = left, top = top),
            bg = "transparent") %>%
    print(target = file)
}
gen_boxplot <- function(plotdat, x,  y, color = "ct", xlab = "", ylab) {
  if (is_empty(plotdat)) return(NULL)
  ggplot(plotdat, aes_string(x = x, y = y)) +
    geom_boxplot(aes_string(fill = color)) +
    coord_flip() +
    theme_bw(base_size=14) +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    xlab(xlab) + ylab(ylab) +
    scale_fill_manual(values=c("#994444", "grey80"))
    # scale_fill_manual(values=pickcols(jellypal, length(levels(plotdat[,x]))))
}
gen_densplot <- function(plotdat, x, y, color = "ct", xlab, ylab = "Frequency") {
  if (is_empty(plotdat)) return(NULL)
  ggplot(plotdat, aes_string(x = x)) +
    geom_density(aes_string(fill = color)) +
    theme_bw(base_size=14) +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.title = element_blank()) +
    xlab(xlab) + ylab(ylab) +
    scale_fill_manual(values=c("#994444", "grey80"))
    # scale_fill_manual(values=pickcols(jellypal, length(levels(plotdat[,y]))))
}
unlog_axis <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}
query_db <- function(query, cnffile = cnf_file, dbname = db_name) {
  db <- dbConnect(RMariaDB::MariaDB(), default.file = cnffile, group = dbname)
  queryres <- dbSendQuery(db, query)
  res <- dbFetch(queryres)
  dbClearResult(queryres)
  dbDisconnect(db)
  return(res)
}

################################### Setup #####################################

## load contrast options
ctcompdf <- query_db('SELECT DISTINCT comp_id, ct, primaryct FROM comparisons;') %>%
  group_by(comp_id) %>% 
  arrange(-primaryct) %>% 
  summarise(ct1 = ct[1], ct2 = ct[2]) %>%
  ungroup() %>% 
  filter(ct1 %in% c("B-cell", "T-cell", "Bone", "Breast",
                    "Colorectal", "Lung", "Melanoma",
                    "CNS", "Pancreas", "Myeloid", "Kidney")) %>%
  mutate(ct1_label = factor(ct1, 
                            levels = c("B-cell", "T-cell", "Myeloid", "Breast",
                                       "Colorectal", "Lung", "Melanoma",
                                       "CNS",  "Bone", "Pancreas", "Kidney"),
                            labels = c("B-cell", "T-cell", "Myeloid", "Breast",
                                       "Colorectal", "Lung", "Melanoma",
                                       "Neuronal",  "Osteoblastic", "Pancreatic", "Renal"))) %>%
  arrange(ct1_label) %>%
  filter(ct2 == "Solid tumor") %>%
  mutate(ct2_label = "Other lineages") %>%
  mutate(compname = paste0(ct1_label, " vs ", ct2_label))
# idxorder <- unique(c(grep("B-cell vs Solid tumor", ctcompdf$compname),
#                      grep("B-cell", ctcompdf$ct1), 1:nrow(ctcompdf)))
# ctcompdf <- ctcompdf[idxorder,]

## pre-load gene and compound info
gi <- query_db('SELECT * FROM gene_info;')
di <- query_db('SELECT * FROM drug_info;')
dilong <- di %>%
  tidyr::separate_rows(target_genes) %>%
  dplyr::rename(target_gene = target_genes) %>%
  dplyr::filter(!is.na(target_gene))

####### TEMP - put in db
pathinfo <- readRDS("data/wp_hs_curated_pathinfo.2025-08-26.rds")

# reactive values
reactvals <- reactiveValues(gene = NULL, 
                            ct1 = NULL,
                            ct2 = NULL,
                            si = NULL,
                            selcpd_metric = "CSS")

##################################### UI ######################################
ui <- fluidPage(
    
    # title
    title = "Delineate", 
    theme = shinytheme("cosmo"),
    titlePanel(tags$h2(tags$a(
        imageOutput("icon", inline = TRUE),
        href="http://137.184.200.69:3838/"), "Delineate: Dependency Lineage Target Explorer")),
    
    # change colors for DT row/column selection
    tags$style(HTML('table.dataTable tr.selected td{background-color: #B9869F99 !important;}')),
    tags$style(HTML('table.dataTable td.selected {background-color: #D0E6EA99 !important;}')),
    tags$style(HTML(".tabbable > .nav > li > a { color:#A7C4C6FF}")),
    tags$style(HTML(".tabbable > .nav > li[class=active] > a { color:black}")),
    
    # lineage comparison choice
    fluidRow(align="center",
             div(style = "display:inline-block;font-size:13px;", 
                 uiOutput("comparison_choice_ui")),
             div(style = "display:inline-block;", 
                 actionButton(style = "font-size:13px;height:30px;padding:5px;text-align:center;", #line-height:10px;border-radius:4px;
                              "loadcomp", label = "Load comparison", icon = shiny::icon("refresh")))),
        
    
    # main UI
    tabsetPanel(id = "maintabs",
      tabPanel("Pathway dependency",
               # dotplot
               tags$h3("Differential pathway sensitivity:"),
               fluidRow(align="center",
                        plotOutput("pathlvl_dotplot", height = 400, width = 700,
                                   hover = "pathlvl_dotplot_hover") %>%
                          withSpinner(color = "#D0E6EA99")),
               fluidRow(align="center",
                        verbatimTextOutput("pathlvl_dotplot_hover_text")),
               br(),br(),
               # plots for 3-different levels
               fluidRow(align="center",
                        splitLayout(cellWidths = c("25%", "25%", "25%"),
                                    tags$h4("Protein dependency"),
                                    tags$h4("Compound sensitivity"),
                                    tags$h4("Metabolite abundance"))),
               fluidRow(align="center",
                        splitLayout(cellWidths = c("25%", "25%", "25%"),
                                    plotOutput("pathcomp_gene_volc", height = 375) %>%
                                      withSpinner(color = "#D0E6EA99"),
                                    plotOutput("pathcomp_cpd_volc", height = 375) %>%
                                      withSpinner(color = "#eed9e0ff"),
                                    plotOutput("pathcomp_metab_volc", height = 375) %>%
                                      withSpinner(color = "#9b629588"))),
               # summary table
               tags$h3("Pathway info:"),
               div(DT::dataTableOutput("pathlvl_summary_dt"),
                   style = "font-size:90%"),
               downloadButton("dl_pathlvl_summary_xls", label = "Download pathways summary",
                              style = "font-size:12px;height:30px;padding:5px;")
               
      ),
        tabPanel("Compound sensitivity",

                 # scatter plot of average
                 tags$h3("Differential sensitivity per compound:"),
                 tags$h5("Drag a box around points to select compounds of interest"),
                 fluidRow(align="center",
                          plotOutput("cpd_volcano", height = 420, width = 660,
                                     brush = "cpd_volcano_brush",
                                     hover = "cpd_volcano_hover") %>%
                            withSpinner(color = "#D0E6EA99")),
                 fluidRow(align="center",
                          verbatimTextOutput("cpd_volcano_hover_text")),

                 # metrics table
                 br(),br(),
                 tags$h3("Selected compounds:"),
                 tags$h5("Select a row to plot compound/target sensitivity accross cell lines below"),
                 div(DT::dataTableOutput("cpd_summary_dt"),
                     style = "font-size:90%"),
                 downloadButton("dl_cpd_summary_xls", label = "Download compounds summary",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 actionButton("cpd_summary_dt_reset", "Reset selection",
                              style = "font-size:12px;height:30px;padding:5px;"),

                 # selected compound ranking info
                 br(),br(),
                 tags$h3("Selected compound:"),
                 fluidRow(align="center", uiOutput("sel_cpd_text")),

                 # rankings plot
                 tags$h4("Differential ranking:"),
                 tags$h5(paste0("Differential sensitivity ranks indicate how selectively sensitive the lineage of interest ",
                                "is relative to the control population (higher rank = more selective). Three primary metrics are ",
                                "considered - the CRISPR effect score, RNAi effect score, and compound sensitivity score (CSS).")),
                 fluidRow(align="center",
                          plotOutput("selcpd_metrics_overview", height = 175, width = 400)
                 ),

                 # Compound metrics plots
                 br(),
                 tags$h3("Drug sensitivity metrics:"),
                 tags$h5(paste0("Drug sensitivity scores integrate the dose-response curve into a single metric - ",
                                "(higher score = more sensitive, range 0-100)")),
                 # CTD
                 fluidRow(splitLayout(cellWidths = c("20%","18%","27%","27%"),
                                      uiOutput("ctd_css_dens_ui"),
                                      uiOutput("ctd_css_group_ui"),
                                      uiOutput("ctd_css_st_ui"),
                                      uiOutput("ctd_ec50_st_ui"))),
                 # GDSC
                 fluidRow(splitLayout(cellWidths = c("20%","19%","26%","26%"),
                                      uiOutput("gdsc_css_dens_ui"),
                                      uiOutput("gdsc_css_group_ui"),
                                      uiOutput("gdsc_css_st_ui"),
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
                 # fluidRow(uiOutput("cpd_crispr_corr_ui")),
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
        tabPanel("Protein dependency",

                 tags$h3("Differential gene dependency scores:"),
                 tags$h5("Drag a box around points to select genes of interest"),
                 fluidRow(align="center",
                          splitLayout(cellWidths = c("30%", "40%"),
                                      tags$h4("CRISPR"),
                                      tags$h4("RNAi"))),
                 fluidRow(align="center",
                          splitLayout(cellWidths = c("35%", "40%"),
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
                 
                 # downloadButton("dl_crispr_volc_png", label = "PNG",
                 #                style = "font-size:12px;height:30px;padding:5px;"),

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
                 tags$h5(paste0("Differential sensitivity ranks indicate how selectively dependent the lineage of interest ",
                                "is relative to the control population (higher rank = more selective). Three primary metrics are ",
                                "considered - the CRISPR effect score, RNAi effect score, and compound sensitivity score (CSS).")),
                 fluidRow(align="center",
                          plotOutput("seldep_metrics_overview", height = 175, width = 400)
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
                 
        ),
        tabPanel("Metabolite abundance",
                 
                 tags$h3("Differential metabolite levels:"),
                 tags$h5("Drag a box around points to select metabolites of interest"),
                 fluidRow(align="center",
                          plotOutput("metab_volcano", height = 420, width = 660,
                                     brush = "metab_volcano_brush",
                                     hover = "metab_volcano_hover") %>%
                            withSpinner(color = "#D0E6EA99")),
                 fluidRow(align="center",
                          verbatimTextOutput("metab_volcano_hover_text")),
                 
                 # gene dependency summary table
                 br(),br(),
                 tags$h3("Selected metabolites:"),
                 tags$h5("Select a row to plot metabolite levels accross cell lines below"),
                 div(DT::dataTableOutput("metab_summary_dt"),
                     style = "font-size:90%"),
                 downloadButton("dl_metab_summary_xls", label = "Download metabolites summary",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 
                 # metabolite dist plots
                 br(), br(),
                 tags$h3("Metabolite distribution:"),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "28%"),
                                      plotOutput("metab_dens", height = 225,
                                                 brush = "crispr_dens_brush"),
                                      plotOutput("metab_box_ds", height = 225),
                                      plotOutput("metab_box_st", height = 225)))

                 
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
    
    # check for parameters in url
    get_url_parameters <- reactive({
      query <- parseQueryString(session$clientData$url_search)
      reactvals$selcompid <- ifelse(is_empty(query$compid), 
                                    ctcompdf$comp_id[1], query$compid)
      if (!is_empty(query$compid)) load_comparison()
      if (!is_empty(query$tab)) updateTabsetPanel(inputId = "maintabs",
                                                  selected = query$tab)
    })
    
    # cell type comparison choice UI
    output$comparison_choice_ui <- renderUI({
      get_url_parameters()
      reactvals$selcompid <- ifelse(is_empty(reactvals$selcompid), 
                                    ctcompdf$comp_id[1], reactvals$selcompid)
      selcompname <- ctcompdf[which(ctcompdf$comp_id==reactvals$selcompid),]$compname
      selectInput('comparison_choice', "Selected lineage comparison", 
                  ctcompdf$compname, selcompname)
    })
    
    # react to cell type comparison choice/load button
    observeEvent(input$loadcomp, {
      print("loadbuttonclicked")
      reactvals$selcompid <- ctcompdf[which(ctcompdf$compname==input$comparison_choice),]$comp_id
      load_comparison()
    })
    # eventReactive(input$comparison_choice, {
    #   
    #   load_comparison()
    # })
    
    # load lineage comparison data
    load_comparison <- reactive({
      print("load_comparison")
      # if (is_empty(input$comparison_choice)) return(NULL)
      # set comparison info
      compidx <- which(ctcompdf$comp_id == reactvals$selcompid)
      reactvals$comp_id <- comp_id <- ctcompdf[compidx,]$comp_id
      reactvals$ct1 <- ct1 <- ctcompdf[compidx,]$ct1
      reactvals$ct2 <- ct2 <- ctcompdf[compidx,]$ct2
      reactvals$ct1_label <- ct1_label <- ctcompdf[compidx,]$ct1_label
      reactvals$ct2_label <- ct1_label <- ctcompdf[compidx,]$ct2_label
      # fetch sample info
      query <- paste0('SELECT sample_info.*, comparisons.ct ',
                      'FROM sample_info INNER JOIN comparisons ON ',
                      'sample_info.sample_id=comparisons.sample_id ',
                      'WHERE comparisons.comp_id = "', comp_id, '";')
      reactvals$si <- query_db(query) %>%
        mutate(ct = factor(ct, levels=c(ct1, ct2)))
      # chemgen res
      reactvals$pathlvl <- readRDS(paste0("data/ctres/WP_GLMres_",
                                          ct1, "_vs_", ct2, ".2025-08-26.rds")) %>%
        mutate(pathway_rank = rank(coef, ties.method = "random"))
      # fetch diff cpd sens
      query <- paste0('SELECT * FROM diff_data INNER JOIN drug_info ON ',
                      'diff_data.feature_id=drug_info.cpd_id ',
                      'WHERE diff_data.comp_id = "', comp_id, '"',
                      'AND diff_data.metric IN ("CSS","AAC","EC50")', ';')
      reactvals$cpdlvl <- query_db(query) %>% arrange(-d)
      # fetch diff dependency
      query <- paste0('SELECT * FROM diff_data INNER JOIN gene_info ON ',
                      'diff_data.feature_id=gene_info.entrez_id ',
                      'WHERE diff_data.comp_id = "', comp_id, '"',
                      'AND diff_data.metric IN ("CRISPR","RNAi")', ';')
      reactvals$genelvl <- query_db(query) %>% 
        mutate(d = -d,  r= -r) %>% # invert scores so high = more sens
        arrange(-d)
      # fetch diff metab
      reactvals$metablvl <- readRDS(paste0("data/ctres/",
                                    ct1, "_vs_", ct2, "_dmetab.rds")) %>%
        mutate(nl10p = -log10(p + 1e-25))
      return(NULL)
    })

    ################################# Pathways #################################
    # plot
    output$pathlvl_dotplot <- renderPlot({
      print("pathlvl_dotplot")
      if(is_empty(reactvals$pathlvl)) return(NULL)
      plotdat <- reactvals$pathlvl
      ct1 <- reactvals$ct1
      ct2 <- reactvals$ct2
      plotdat <- plotdat %>%
        mutate(color = ifelse(coef < -0.01, "resistant",
                              ifelse(coef > 0.01, "sensitive", "ns"))) %>%
        mutate(color = factor(color, levels= c("sensitive", "ns", "resistant"))) %>%
        filter(!is.na(color))
      plotdat %>%
        ggplot(aes(x=pathway_rank, y=coef, size=coef)) +
        geom_hline(yintercept = 0, lty=2, alpha=0.4) +
        geom_point(aes(color=color)) +
        scale_color_manual(values=c("#994444","grey80","#528199cc")) +
        scale_x_continuous(name = paste0("Pathways (n=", nrow(plotdat), ")")) +
        scale_y_continuous(name = "\u0394 Pathway score \n(Coefficient)") +
        theme_bw(base_size = 17) +
        theme(legend.title = element_blank(),
              legend.background = element_blank(),
              panel.grid = element_blank())
    })
    # brushed
    get_sel_path <- reactive({
      print("get_sel_path")
      selpath <- reactvals$pathlvl
      if (is_empty(selpath)) return(NULL)
      selpath <- selpath %>%
        mutate_if(is.numeric, round, digits = 3)
      if (is_empty(input$pathlvl_dotplot_click)) return(selpath)
      selpath %>%
        nearPoints(input$pathlvl_dotplot_click)
    })
    # hovered
    output$pathlvl_dotplot_hover_text <- renderText({
      if (is_empty(reactvals$pathlvl) | is_empty(input$pathlvl_dotplot_hover)) return(NULL)
      pathsub <- reactvals$pathlvl %>%
        nearPoints(input$pathlvl_dotplot_hover)
      paste0("Pathways near cursor: ", paste0(unique(pathsub$pathway_name), collapse=";"))
    })
    # summary table
    output$pathlvl_summary_dt <- DT::renderDataTable({
      print("pathlvl_summary_dt")
      selpath <- get_sel_path()
      if (is.null(selpath)) return(NULL)
      selpath <- selpath %>%
        mutate_if(is.numeric, round, digits = 3) %>%
        mutate(topfeat = ifelse(coef < 0, botf, topf)) %>%
        dplyr::select(pathway_id:coef, topfeat) %>%
        dplyr::rename(`Pathway ID` = pathway_id,
                      `Pathway Name` = pathway_name,
                      `Score (coef)` = coef,
                      `Top diff features` = topfeat)
      DT::datatable(
        data = selpath,
        rownames = F,
        options = list(pageLength = 25),
        selection = list(mode = 'single', target = "row", selected = 2))
    })
    output$dl_pathlvl_summary_xls <- downloadHandler(
      filename = function() {
        paste0("Chemgen_pathways_summary_table_", reactvals$ct1, "_vs_",
               reactvals$ct2, ".xlsx")
      },
      content = function(file) {
        dat <- get_sel_path()
        if(is_empty(dat)) return(NULL)
        writexl::write_xlsx(dat, path=file)
      }
    )
    
    # select corresponding path/gene/cpds for selected row
    get_sel_components <- reactive({
      print("get_sel_components")
      if(is_empty(input$pathlvl_summary_dt_rows_selected)) return(NULL)
      selrow <- reactvals$pathlvl %>%
        dplyr::slice(input$pathlvl_summary_dt_rows_selected)
      pisub <- pathinfo %>% filter(pathway_id %in% selrow$pathway_id)
      selgenes <- unique(pisub[which(!is.na(pisub$gene_symbol)),]$gene_symbol)
      selcpds <- unique(pisub[which(!is.na(pisub$chembl_id)),]$chembl_id)
      selmetab <- unique(pisub[which(!is.na(pisub$ccle_name)),]$ccle_name)
      reactvals$selgenes <- c("",selgenes)
      reactvals$selcpds <- c("",selcpds)
      reactvals$selmetab <- c("",selmetab)
      return(NULL)
    })
    
    ### plots
    output$pathcomp_gene_volc <- renderPlot({
      print("pathcomp_gene_volc")
      get_sel_components()
      if(is_empty(reactvals$genelvl)) return(NULL)
      plotdat <- reactvals$genelvl %>%
        filter(metric == "CRISPR") %>%
        mutate(selected = gene_symbol %in% reactvals$selgenes)
      ggplot(data=plotdat, aes(x=d, y=p, label=gene_symbol)) +
        geom_vline(xintercept = 0, lty=2, alpha=0.4) +
        geom_point(data=dplyr::filter(plotdat, selected==F), color="#88888811", size=1.5) +
        geom_point(data=dplyr::filter(plotdat, selected==T), color="#994444", size=4) +
        ggrepel::geom_text_repel(data=dplyr::filter(plotdat, selected),
                                 size=6, min.segment.length = 0) +
        scale_x_continuous(name = paste0("\u0394 Protein dependency ",
                                         reactvals$ct1, "/", reactvals$ct2)) +
        scale_y_continuous(name = "Significance (-log10 p-value)") +
        theme_bw(base_size = 17) +
        theme(legend.title = element_blank(),
              legend.background = element_blank(),
              panel.grid = element_blank(),
              legend.position = "none")
    })
    output$pathcomp_cpd_volc <- renderPlot({
      print("pathcomp_cpd_volc")
      if(is_empty(reactvals$cpdlvl)) return(NULL)
      plotdat <- reactvals$cpdlvl %>%
        dplyr::filter(metric == reactvals$selcpd_metric) %>%
        mutate(selected = chembl_id %in% reactvals$selcpds)
      ggplot(data=plotdat, aes(x=r, y=p, label=cpd_name)) +
        geom_vline(xintercept = 0, lty=2, alpha=0.4) +
        geom_point(data=dplyr::filter(plotdat, selected==F), color="#88888833", size=1.8) +
        geom_point(data=dplyr::filter(plotdat, selected==T), color="#994444", size=4) +
        ggrepel::geom_text_repel(data=dplyr::filter(plotdat, selected),
                                 size=6, min.segment.length = 0) +
        scale_x_continuous(name = paste0("\u0394 Compound sensitivity ",
                                         reactvals$ct1, "/", reactvals$ct2)) +
        scale_y_continuous(name = "Significance (-log10 p-value)") +
        theme_bw(base_size = 17) +
        theme(legend.title = element_blank(),
              legend.background = element_blank(),
              panel.grid = element_blank())
    })
    output$pathcomp_metab_volc <- renderPlot({
      print("pathcomp_metab_volc")
      if(is_empty(reactvals$metablvl)) return(NULL)
      plotdat <- reactvals$metablvl %>%
        mutate(selected = feature %in% reactvals$selmetab)
      ggplot(data=plotdat, aes(x=d, y=nl10p, label=feature)) +
        geom_vline(xintercept = 0, lty=2, alpha=0.4) +
        geom_point(data=dplyr::filter(plotdat, selected==F), color="#88888844", size=1.8) +
        geom_point(data=dplyr::filter(plotdat, selected==T), color="#994444", size=4) +
        ggrepel::geom_text_repel(data=dplyr::filter(plotdat, selected),
                                 size=6, min.segment.length = 0) +
        scale_x_continuous(name = paste0("\u0394 Metabolite level ",
                                         reactvals$ct1, "/", reactvals$ct2)) +
        scale_y_continuous(name = "Significance (-log10 p-value)") +
        theme_bw(base_size = 17) +
        theme(legend.title = element_blank(),
              legend.background = element_blank(),
              panel.grid = element_blank())
    })
    
    ################################# Compounds ################################
    
    # scatter plot differential compound scores
    output$cpd_volcano <- renderPlot({
      print("cpd_volcano")
      if (is_empty(reactvals$cpdlvl)) return(NULL)
      selcpdrow <- get_sel_cpdrow_info()
      if (is_empty(reactvals$selcpd_compounds)) {
        plotdat <- reactvals$cpdlvl %>%
          dplyr::filter(metric == reactvals$selcpd_metric) %>%
          mutate(color = ifelse(p < 2, "ns", ifelse(d > 0, "sensitive", "resistant"))) %>%
          mutate(color = factor(color, levels= c("ns","sensitive","resistant"))) %>%
          filter(!is.na(color))
        plotdat %>%
          ggplot(aes(x=r, y=p)) +
          geom_vline(xintercept = 0, lty=2, alpha=0.4) +
          geom_point(aes(color=color), size = 2) +
          scale_color_manual(values=c("grey80","#994444","#528199cc")) +
          scale_x_continuous(name = paste0("\u0394 CSS ", reactvals$ct1,
                                           "/", reactvals$ct2)) +
          scale_y_continuous(name = "Significance (-log10 p-value)") +
          theme_bw(base_size = 17) +
          theme(legend.title = element_blank(),
                legend.background = element_blank(),
                panel.grid = element_blank())
      } else {
        plotdat <- reactvals$cpdlvl %>%
          dplyr::filter(metric == reactvals$selcpd_metric) %>%
          mutate(color = ifelse(cpd_name %in% reactvals$selcpd_compounds, 
                                ifelse(d > 0, "sensitive", "resistant"), "unselected")) %>%
          mutate(color = factor(color, levels= c("unselected","sensitive","resistant"))) %>%
          filter(!is.na(color))
        plotdat %>%
          ggplot(aes(x=r, y=p, label=cpd_name)) +
          geom_vline(xintercept = 0, lty=2, alpha=0.4) +
          geom_point(data=filter(plotdat, color=="unselected"),
                     aes(color=color), size = 2, color="grey80") +
          geom_point(data=filter(plotdat, color=="sensitive"),
                     aes(color=color), size = 5, color="#994444") +
          geom_point(data=filter(plotdat, color=="resistant"),
                     aes(color=color), size = 5, color="#447799") +
          ggrepel::geom_text_repel(data=filter(plotdat, color!="unselected"),
                                   size=6, min.segment.length = 0) +
          scale_color_manual(values=c("grey80","#994444","#447799")) +
          scale_x_continuous(name = paste0("\u0394 CSS ", reactvals$ct1,
                                           "/", reactvals$ct2)) +
          scale_y_continuous(name = "Significance (-log10 p-value)",
                             limits = c(0,12)) +
          theme_bw(base_size = 17) +
          theme(legend.title = element_blank(),
                legend.background = element_blank(),
                panel.grid = element_blank())
      }
    })

    # brushed cpd
    get_sel_cpd <- reactive({
        print("get_sel_cpd")
        cpdsub <- reactvals$cpdlvl %>%
          dplyr::filter(metric == reactvals$selcpd_metric)
        if (is_empty(cpdsub)) return(NULL)
        if (is_empty(input$cpd_volcano_brush)) return(cpdsub)
        cpdsub %>%
          brushedPoints(input$cpd_volcano_brush)
    })
    # hovered cpd
    output$cpd_volcano_hover_text <- renderText({
      if (is_empty(input$cpd_volcano_hover)) return(NULL)
      cpdsub <- reactvals$cpdlvl %>%
        dplyr::filter(metric == reactvals$selcpd_metric)
        nearPoints(input$cpd_volcano_hover)
      paste0("Compounds near cursor: ", paste0(unique(cpdsub$cpd_name), collapse=";"))
    })

    # cpd summary table
    output$cpd_summary_dt <- DT::renderDataTable({
      print("cpd_summary_dt")
      cpdsub <- get_sel_cpd()
      if (is.null(cpdsub)) return(NULL)
      cpdsub <- cpdsub %>%
        dplyr::select(cpd_name, datasets, metric, r:n_ct2, pubchem_cid:aliases) %>%
        dplyr::rename(Compound = cpd_name,
                      `PubChem CID` = pubchem_cid,
                      `Target genes` = target_genes,
                      `Datasets` = datasets,
                      `Signif (-log10P)` = p)
      colnames(cpdsub) <- gsub("_", " ", colnames(cpdsub))
      colnames(cpdsub) <- gsub("ct1", reactvals$ct1_label, colnames(cpdsub))
      colnames(cpdsub) <- gsub("ct2", reactvals$ct2_label, colnames(cpdsub))
      DT::datatable(
        data = cpdsub,
        rownames = F,
        options = list(pageLength = 25),
        selection = list(mode = 'multiple', target = "row")
      )
    })
    output$dl_cpd_summary_xls <- downloadHandler(
      filename = function() {
        paste0("Compounds_summary_table_", reactvals$ct1_label, "_vs_",
               reactvals$ct2_label, ".xlsx")
      },
      content = function(file) {
        drug_summary <- get_sel_cpd()
        if(is_empty(drug_summary)) return(NULL)
        writexl::write_xlsx(drug_summary, path=file)
      }
    )

    # cpd row selection reset
    cpd_summary_dt_proxy <- dataTableProxy('cpd_summary_dt')
    observeEvent(input$cpd_summary_dt_reset, {
      print("cpd_summary_dt_reset")
      cpd_summary_dt_proxy %>% selectRows(NULL)
    })
    
    # gte info of selected row(s)
    get_sel_cpdrow_info <- reactive({
      print("get_sel_cpdrow")
      cpdsummary <- get_sel_cpd()
      if(is.null(input$cpd_summary_dt_rows_selected) |
         is.null(cpdsummary)) return(NULL)
      cpdsummary <- cpdsummary %>%
        dplyr::slice(input$cpd_summary_dt_rows_selected)
      print(head(cpdsummary))
      if (is_empty(cpdsummary)) return(NULL)
      compounds <- unique(cpdsummary$cpd_name)
      target_gene_symb <- c(unlist(unique(strsplit(cpdsummary$target_genes, split=", "))))
      target_gene_eid <- unique(gi[which(gi$gene_symbol %in% target_gene_symb),]$entrez_id)
      reactvals$selcpd_compounds <- compounds
      reactvals$selcpd_gene_symbols <- target_gene_symb
      reactvals$selcpd_gene_ids <- target_gene_eid
      return(cpdsummary)
    })
    
    ### score summary for selected compound  ####
    # text just highlighting the selected compound/gene
    output$sel_cpd_text <- renderUI({
        print("sel_cpd_text")
        plotdat <- get_sel_cpdrow_info()
        if(is_empty(plotdat)) return(NULL)
        if(is_empty(input$cpd_summary_dt_rows_selected)) return(NULL)
        HTML(paste0("<h4>Selected compound: <b>", plotdat$cpd_name,
                    "</b>, Targets: <b>", plotdat$target_genes, "</b></h4>"))
    })
    # metrics overview
    output$selcpd_metrics_overview <- renderPlot({
        print("selcpd_metrics_overview")
        if(is_empty(reactvals$selcpd_compounds)) return(NULL)
      print(reactvals$selcpd_gene_symbols)  
      crispr <- reactvals$genelvl %>% 
          filter(metric == "CRISPR") %>%
          mutate(prcntl = rank(d)/n()) %>% 
          filter(gene_symbol %in% reactvals$selcpd_gene_symbols) %>%
          dplyr::select(metric, prcntl)
        rnai <- reactvals$genelvl %>% 
          filter(metric == "RNAi") %>%
          mutate(prcntl = rank(d)/n()) %>% 
          filter(gene_symbol %in% reactvals$selcpd_gene_symbols) %>%
          dplyr::select(metric, prcntl)
        cpds <- reactvals$cpdlvl %>%
          dplyr::filter(metric==reactvals$selcpd_metric) %>%
          dplyr::filter(!is.na(d)) %>%
          mutate(prcntl = rank(d)/n()) %>%
          filter(cpd_name %in% reactvals$selcpd_compounds) %>%
          mutate(metric = "CSS") %>% dplyr::select(metric, prcntl)
        plotdat <- rbind(crispr, rnai, cpds) %>%
          mutate(percentile = prcntl*100,
                 metric = factor(metric,
                                 levels = c("CSS","RNAi","CRISPR"),
                                 labels = c("Cpd sensitivity","RNAi score","CRISPR score")))
        ggplot(plotdat, aes(x=metric, y=percentile)) +
          geom_point(size=5, color="#994444", shape=21) +
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


    ######### sensitivity metrics for selected compound  ###########

    # get drug metrics for selected row
    get_sel_cpd_metrics <- reactive({
      print("get_sel_cpd_metrics")
      si <- reactvals$si
      ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
      if (is_empty(reactvals$selcpd_compounds)) return(NULL)
      # query db for selected compounds
      query <- paste0('SELECT * FROM drug_metrics ',
                      'WHERE cpd_name IN ("',
                      paste0(reactvals$selcpd_compounds, collapse = '", "'),
                      '") AND sample_id IN ("',
                      paste0(si$sample_id, collapse = '", "'),
                      '");')
      drug_metrics <- query_db(query)
      drug_metrics <- merge(drug_metrics, si, by="sample_id")
      if (is_empty(drug_metrics) | nrow(drug_metrics) < 1) return(NULL)
      # format sample levels
      matchidx <- match(drug_metrics$sample_id, si$sample_id)
      si$st <- ifelse(si$ct==ct1, si$ds_subtype, ct2)
      stlvls <- as.character(unique(si$st))
      stlvls <- c(stlvls[!grepl(ct2, stlvls)], ct2)
      si$dg <- ifelse(si$ct==ct1, ct1, si$ds_group)
      dglvls <- unique(si$dg)
      dglvls <- c(dglvls[!grepl(ct1, dglvls)], ct1)
      drug_metrics <- drug_metrics %>%
        mutate(ct = si[matchidx,]$ct,
               st = si[matchidx,]$st,
               dg = si[matchidx,]$dg) %>%
        mutate(ct = factor(ct, levels=c(ct1, ct2)),
               st = factor(st, levels=c(stlvls)),
               dg = factor(dg, levels=c(dglvls))) %>%
        dplyr::filter(!is.na(ct) & !is.na(dg) & !dg=="Other")
      if (is_empty(drug_metrics) | nrow(drug_metrics) < 1) return(NULL)
      reactvals$selcpd_cpd_metrics <- drug_metrics
      return(drug_metrics)
    })

    # CSS density plot
    output$ctd_css_dens_ui <- renderUI({
      print("ctd_css_dens_ui")
      plotdat <- get_sel_cpd_metrics()
      if (is_empty(plotdat)) return(NULL)
      plotdat <- plotdat %>%
        filter(dataset == "CTD")
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_densplot(plotdat, "CSS",
                        "ct", xlab = paste0(reactvals$selcpd_metric, " score")) +
        ggtitle("CTD")
      output$ctd_css_dens_plot <- renderPlot(p)
      plotOutput("ctd_css_dens_plot", height = 225)
    })
    output$gdsc_css_dens_ui <- renderUI({
      print("gdsc_css_dens_ui")
      if (is_empty(reactvals$selcpd_cpd_metrics)) return(NULL)
      plotdat <- reactvals$selcpd_cpd_metrics %>%
        filter(dataset %in% c("GDSC1","GDSC2"))
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_densplot(plotdat, reactvals$selcpd_metric,
                        "ct", xlab = paste0(reactvals$selcpd_metric, " score")) +
        ggtitle("GDSC")
      output$gdsc_css_dens_plot <- renderPlot(p)
      plotOutput("gdsc_css_dens_plot", height = 225)
    })

    # drug sensitivity boxplots - disease groups
    output$ctd_css_group_ui <- renderUI({
      print("ctd_css_group_ui")
      if (is_empty(reactvals$selcpd_cpd_metrics)) return(NULL)
      plotdat <- reactvals$selcpd_cpd_metrics %>%
        filter(dataset == "CTD")
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "dg", reactvals$selcpd_metric,
                       ylab = paste0(reactvals$selcpd_metric, " score"))
      output$ctd_css_group_plot <- renderPlot(p)
      plotOutput("ctd_css_group_plot", height = 225)
    })
    output$gdsc_css_group_ui <- renderUI({
      print("gdsc_css_group_ui")
      if (is_empty(reactvals$selcpd_cpd_metrics)) return(NULL)
      plotdat <- reactvals$selcpd_cpd_metrics %>%
        filter(dataset %in% c("GDSC1","GDSC2"))
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "dg", reactvals$selcpd_metric,
                       ylab = paste0(reactvals$selcpd_metric, " score"))
      output$gdsc_css_group_plot <- renderPlot(p)
      plotOutput("gdsc_css_group_plot", height = 225)
    })

    # drug sensitivity boxplots - subtypes
    output$ctd_css_st_ui <- renderUI({
      print("ctd_css_st_ui")
      if (is_empty(reactvals$selcpd_cpd_metrics)) return(NULL)
      plotdat <- reactvals$selcpd_cpd_metrics %>%
        filter(dataset == "CTD")
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "st", reactvals$selcpd_metric,
                       ylab = paste0(reactvals$selcpd_metric, " score")) +
        ggtitle(paste0(reactvals$ct1, " subtypes:"))
      output$ctd_css_st_plot <- renderPlot(p)
      plotOutput("ctd_css_st_plot", height = 225)
    })
    output$gdsc_css_st_ui <- renderUI({
      print("gdsc_css_st_ui")
      if (is_empty(reactvals$selcpd_cpd_metrics)) return(NULL)
      plotdat <- reactvals$selcpd_cpd_metrics %>%
        filter(dataset %in% c("GDSC1","GDSC2"))
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "st", reactvals$selcpd_metric,
                       ylab = paste0(reactvals$selcpd_metric, " score")) +
        ggtitle(paste0(reactvals$ct1, " subtypes:"))
      output$gdsc_css_st_plot <- renderPlot(p)
      plotOutput("gdsc_css_st_plot", height = 225)
    })

    # EC50 boxplots - subtypes
    output$ctd_ec50_st_ui <- renderUI({
      print("ctd_ec50_st_ui")
      if (is_empty(reactvals$selcpd_cpd_metrics)) return(NULL)
      plotdat <- reactvals$selcpd_cpd_metrics %>%
        filter(dataset == "CTD")
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "st", "logEC50", ylab = "EC50 (log \U03BCM)") +
        ggtitle("")
      output$ctd_ec50_st_plot <- renderPlot(p)
      plotOutput("ctd_ec50_st_plot", height = 225)
    })
    output$gdsc_ec50_st_ui <- renderUI({
      print("gdsc_ec50_st_ui")
      if (is_empty(reactvals$selcpd_cpd_metrics)) return(NULL)
      plotdat <- reactvals$selcpd_cpd_metrics %>%
        filter(dataset %in% c("GDSC1","GDSC2"))
      if (nrow(plotdat) < 1) return(NULL)
      p <- gen_boxplot(plotdat, "st", "logEC50", ylab = "EC50 (log \U03BCM)")
      output$gdsc_ec50_st_plot <- renderPlot(p)
      plotOutput("gdsc_ec50_st_plot", height = 225)
    })

    # data table for selected cpd
    output$cpd_metrics_dt <- DT::renderDataTable({
        if(is_empty(reactvals$selcpd_cpd_metrics)) return(NULL)
        metrics <- reactvals$selcpd_cpd_metrics %>%
          dplyr::select(cpd_name, dataset, sample_name, rrid, logEC50, CSS, AAC_obs,
                 maxr_rel, ds_type, ds_subtype) %>%
          dplyr::rename("Compound name" = cpd_name,
                        "Dataset" = dataset,
                        "Sample name" = sample_name,
                        "RRID" = rrid,
                        "EC50 (log \U03BCM)" = logEC50,
                        "AAC" = AAC_obs,
                        "Max response" = maxr_rel,
                        "Cell type" = ds_type,
                        "Subtype" = ds_subtype)
        DT::datatable(data = metrics,
                      rownames = F,
                      options = list(pageLength = 10))
    })
    # DT download buttons
    output$dl_cpd_metrics_xls <- downloadHandler(
        filename = function() {
            paste0("Compound_metrics_", reactvals$sel_cpd_names,
                   "_", reactvals$ct1,  ".xlsx")
        },
        content = function(file) {
            metrics <- reactvals$selcpd_cpd_metrics
            if (is_empty(metrics) | is_empty(reactvals$ct1)) return(NULL)
            writexl::write_xlsx(metrics, path=file)
        })

    ######### dependency metrics for selected compound targets ###########

    # get dependency metrics for selected compound targets
    get_dep_metrics_cpd <- reactive({
      print("get_dep_metrics_cpd")
      si <- reactvals$si
      ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
      entrezids <- reactvals$selcpd_gene_ids
      if (is_empty(entrezids)) return(NULL)
      # query db for selected compounds
      query <- paste0('SELECT * FROM dep_data ',
                      'WHERE entrez_id IN ("',
                      paste0(entrezids, collapse = '", "'),
                      '") AND sample_id IN ("',
                      paste0(si$sample_id, collapse = '", "'),
                      '");')
      depdata <- query_db(query)
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
      reactvals$selcpd_depmetrics <- depdata
      return(depdata)
    })

    ## CRISPR
    # # correlation with each target
    # output$cpd_crispr_corr_ui <- renderUI({
    #   print("crispr_cpd_corr")
    #   plotdat <- get_dep_metrics_cpd()
    #   validate(need(!is_empty(plotdat),
    #                 "      No CRISPR data associated with selected compound"))
    #   cpdsub <- get_sel_cpd_metrics()
    #   plotdat <- plotdat %>%
    #     dplyr::filter(assay_type == "CRISPR")
    #   plotdat <- merge(plotdat, cpdsub)
    #   print(head(plotdat))
    #   validate(need(nrow(plotdat) > 1,
    #                 "      No CRISPR data associated with selected compound"))
    #   p <- ggplot(plotdat, aes(x=score, y=CSS)) +
    #     geom_point() +
    #     geom_smooth(method="lm") +
    #     facet_wrap(vars(entrez_id))
    #   output$cpd_crispr_corr <- renderPlot(p)
    #   plotOutput("cpd_crispr_corr",
    #              height = 225*(ceiling(length(unique(plotdat$entrez_id))/3)),
    #              width = min(600*length(unique(plotdat$entrez_id))/3, 1200))
    # })
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
      validate(need(!is_empty(reactvals$selcpd_depmetrics),""))
      plotdat <- reactvals$selcpd_depmetrics %>%
        dplyr::filter(assay_type == "CRISPR")
      validate(need(!is_empty(plotdat),""))
      gen_boxplot(plotdat, "dg", "score", ylab = "CRISPR effect score")
    })
    # CRISPR boxplot by subtype
    output$cpd_crispr_box_st <- renderPlot({
      print("crispr_box_ds")
      validate(need(!is_empty(reactvals$selcpd_depmetrics),""))
      plotdat <- reactvals$selcpd_depmetrics %>%
        dplyr::filter(assay_type == "CRISPR")
      validate(need(!is_empty(plotdat),""))
      gen_boxplot(plotdat, "st", "score", ylab = "CRISPR effect score")
    })

    ### RNAi
    # RNAi density plot by disease
    output$cpd_rnai_dens <- renderPlot({
      print("crispr_dens")
      validate(need(!is_empty(reactvals$selcpd_depmetrics),
                    "      No RNAi data associated with selected compound"))
      plotdat <- reactvals$selcpd_depmetrics %>%
        dplyr::filter(assay_type == "RNAi")
      validate(need(nrow(plotdat) > 1,
                    "      No RNAi data associated with selected compound"))
      gen_densplot(plotdat, "score", "ct", xlab = "RNAi effect score")
    })
    # RNAi boxplot by disease
    output$cpd_rnai_box_ds <- renderPlot({
      print("rnai_box_ds")
      validate(need(!is_empty(reactvals$selcpd_depmetrics),""))
      plotdat <- reactvals$selcpd_depmetrics %>%
        dplyr::filter(assay_type == "RNAi")
      validate(need(!is_empty(plotdat),""))
      gen_boxplot(plotdat, "dg", "score", ylab = "RNAi effect score")
    })
    # RNAi boxplot by subtype
    output$cpd_rnai_box_st <- renderPlot({
      print("rnai_box_st")
      validate(need(!is_empty(reactvals$selcpd_depmetrics),""))
      plotdat <- reactvals$selcpd_depmetrics %>%
        dplyr::filter(assay_type == "RNAi")
      validate(need(!is_empty(plotdat),""))
      gen_boxplot(plotdat, "st", "score", ylab = "RNAi effect score")
    })

    # data table of dependency data
    output$cpd_dep_metrics_dt <- DT::renderDataTable({
        dat <- get_dep_metrics_cpd()
        if (is_empty(dat) | is_empty(reactvals$ct1)) return(NULL)
        dat <- dat %>%
          dplyr::select_at(c(1:5,8,10:17))
        if (nrow(dat) < 1) return(NULL)
        # if (!is_empty(input$crispr_dens_brush)) {
        #     brushinfo <- input$crispr_dens_brush
        #     dat <- dat %>% dplyr::filter(score > brushinfo$xmin &
        #                                  score < brushinfo$xmax)
        # }
        DT::datatable(
            data = dat,
            rownames = F,
            selection = "none",
            options = list(pageLength = 10))
    })
    output$dl_cpd_dep_metrics_xls <- downloadHandler(
      filename = function() {
        paste0("CRISPR_metrics_",
               paste0(reactvals$sel_cpd_target_symbols, collapse="_"),
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
    gen_crispr_volc <- reactive({
      if(is_empty(reactvals$genelvl)) return(NULL)
      plotdat <- reactvals$genelvl %>%
        filter(metric == "CRISPR") %>%
        mutate(color = ifelse(p < 1.3, "ns", 
                              ifelse(d > 0, "dependency", "tolerance"))) %>%
        mutate(color = factor(color, levels= c("dependency", "ns", "tolerance")))
      plotdat %>%
        ggplot(aes(x=d, y=p)) +
        geom_vline(xintercept = 0, lty=2, alpha=0.4) +
        geom_point(aes(color=color), size = 2) +
        scale_color_manual(values=c("#994444","grey80","#447799")) +
        scale_x_continuous(name = paste0("\u0394 CRISPR score ", 
                                         reactvals$ct1, "/", reactvals$ct2)) +
        scale_y_continuous(name = "Significance (-log10 p-value)") +
        theme_bw(base_size = 17) +
        theme(legend.position = "none",
              panel.grid = element_blank())
    })
    output$crispr_scatter <- renderPlot({
      print("crispr_scatter")
      gen_crispr_volc()
    })
    # scatter plot of differential RNAi scores
    output$rnai_scatter <- renderPlot({
      print("rnai_scatter")
      if(is_empty(reactvals$genelvl)) return(NULL)
      plotdat <- reactvals$genelvl %>% 
        filter(metric == "RNAi") %>%
        mutate(color = ifelse(p < 1.3, "ns", 
                              ifelse(d > 0, "dependency", "tolerance"))) %>%
        mutate(color = factor(color, levels= c("dependency", "ns", "tolerance")))
      plotdat %>%
        ggplot(aes(x=d, y=p)) +
        geom_vline(xintercept = 0, lty=2, alpha=0.4) +
        geom_point(aes(color=color), size = 2) +
        scale_color_manual(values=c("#994444","grey80","#447799"),
                           name=paste0(reactvals$ct1, "-specific")) +
        scale_x_continuous(name = paste0("\u0394 RNAi score ", 
                                         reactvals$ct1, "/", reactvals$ct2)) +
        scale_y_continuous(name = "Significance (-log10 p-value)") +
        theme_bw(base_size = 17) +
        theme(legend.background = element_blank(),
              panel.grid = element_blank())
    })

    # download CRISPR volc as png
    output$dl_crispr_volc_png <- downloadHandler(
      filename = function() {
        paste0("CRISPR_diff_", reactvals$ct1, "_vs_", reactvals$ct2, ".png")
      },
      content = function(file) {
        p <- gen_crispr_volc()
        ggsave(plot = p, filename = file,
               height = 4, width = 6)
      }
    )
    
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

    # get all selected genes
    get_sel_genedep <- reactive({
      print("get_sel_genedep")
      genesub <- reactvals$genelvl
      if(is_empty(genesub)) return(NULL)
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
      } else {
        genesub <- genesub #%>%
          #toptail(n=250)
      }
      return(genesub)
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
        sel_genedep <- sel_genedep %>%
          dplyr::select(entrez_id:gene_name, metric:n_ct2) %>%
          mutate_if(is.numeric, round, digits=2) %>%
          arrange(-d) %>%
          dplyr::rename("Entrez ID" = entrez_id,
                        'Gene symbol' = gene_symbol,
                        'Gene name' = gene_name,
                        'Effect size' = r,
                        'Difference' = d,
                        'Significance (-log10P)' = p)
        colnames(sel_genedep) <- sub("_", " ", colnames(sel_genedep))
        colnames(sel_genedep) <- sub("ct1", reactvals$ct1_label, colnames(sel_genedep))
        colnames(sel_genedep) <- sub("ct2", reactvals$ct2_label, colnames(sel_genedep))
        DT::datatable(
            data = sel_genedep,
            rownames = F,
            options = list(pageLength = 25),
            selection = list(mode = 'single', target = "row", selected = 1)
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
    
    ## score summary for selected gene ##
    # get selected gene
    get_sel_deprow <- reactive({
      print("get_sel_deprow")
      depsummary <- get_sel_genedep()
      if(is.null(input$genedep_summary_dt_rows_selected) |
         is.null(depsummary)) return(NULL)
      depsummary <- depsummary %>%
        dplyr::slice(input$genedep_summary_dt_rows_selected)
      if (is_empty(depsummary)) return(NULL)
      entrez_id <- unique(depsummary$entrez_id)
      gene_symbol <- unique(depsummary$gene_symbol)
      compounds <- unique(dilong[which(dilong$target_gene %in% gene_symbol),]$cpd_name)
      reactvals$dep_compounds <- compounds
      reactvals$dep_gene_symbol <- gene_symbol
      reactvals$dep_gene_id <- entrez_id
      return(depsummary)
    })
    # text just highlighting the selected gene
    output$sel_gene_text <- renderUI({
      print("sel_gene_text")
      depselrow <- get_sel_deprow()
      if(is_empty(depselrow)) return(NULL)
      HTML(paste0("<h4>Selected protein: <b>", reactvals$dep_gene_symbol,
                  "</b>, Compounds targetting: <b>", 
                  paste0(reactvals$dep_compounds, collapse=", "), "</b></h4>"))
    })
    # metrics overview
    output$seldep_metrics_overview <- renderPlot({
      print("seldep_metrics_overview")
      seldeprow <- get_sel_deprow()
      if(is_empty(seldeprow)) return(NULL)
      crispr <- reactvals$genelvl %>% filter(metric == "CRISPR") %>%
        mutate(prcntl = rank(d)/n()) %>% 
        filter(gene_symbol %in% reactvals$dep_gene_symbol) %>%
        dplyr::select(metric, prcntl)
      rnai <- reactvals$genelvl %>% filter(metric == "RNAi") %>%
        mutate(prcntl = rank(d)/n()) %>% 
        filter(gene_symbol %in% reactvals$dep_gene_symbol) %>%
        dplyr::select(metric, prcntl)
      cpds <- reactvals$cpdlvl %>%
        dplyr::filter(metric==reactvals$selcpd_metric) %>%
        dplyr::filter(!is.na(d)) %>%
        mutate(prcntl = rank(d)/n()) %>%
        filter(cpd_name %in% reactvals$dep_compounds) %>%
        mutate(metric = "CSS") %>% dplyr::select(metric, prcntl)
      plotdat <- rbind(crispr, rnai, cpds) %>%
        mutate(percentile = prcntl*100,
               metric = factor(metric,
                               levels = c("CSS","RNAi","CRISPR"),
                               labels = c("Cpd sensitivity","RNAi score","CRISPR score")))
      ggplot(plotdat, aes(x=metric, y=percentile)) +
        geom_point(size=5, color="#994444", shape=21) +
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
    
    ######### dependency metrics for selected gene #############

    # get dependency metrics for selected compound targets
    get_dep_metrics_gene <- reactive({
      print("get_dep_metrics_gene")
      si <- reactvals$si
      ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
      deprow <- get_sel_deprow()
      entrezid <- reactvals$dep_gene_id
      if (is_empty(entrezid)) return(NULL)
      # query db for selected compounds
      query <- paste0('SELECT * FROM dep_data ',
                      'WHERE entrez_id IN ("',
                      paste0(entrezid, collapse = '", "'),
                      '") AND sample_id IN ("',
                      paste0(si$sample_id, collapse = '", "'),
                      '");')
      depdata <- query_db(query)
      # format sample levels
      depdata <- merge(depdata, si, by="sample_id")
      depdata$st <- ifelse(depdata$ct==ct1, depdata$ds_subtype, ct2)
      stlvls <- as.character(unique(depdata$st))
      stlvls <- c(stlvls[!grepl(ct2, stlvls)], ct2)
      depdata$dg <- ifelse(depdata$ct==ct1, ct1, depdata$ds_group)
      dglvls <- unique(depdata$dg)
      dglvls <- c(dglvls[!grepl(ct1, dglvls)], ct1)
      depdata <- depdata  %>%
        mutate(ct = factor(ct, levels=c(ct1, ct2)),
               st = factor(st, levels=c(stlvls)),
               dg = factor(dg, levels=c(dglvls))) %>%
        dplyr::filter(!is.na(ct) & !is.na(dg) & !dg=="Other")
      if (is_empty(depdata) | nrow(depdata) < 1) return(NULL)
      reactvals$seldeprow_depmetrics <- depdata
      return(NULL)
    })

    ### CRISPR
    # CRISPR density plot by disease
    output$crispr_dens_dep <- renderPlot({
      print("crispr_dens_dep")
      get_dep_metrics_gene()
      validate(need(!is_empty(reactvals$seldeprow_depmetrics),
                    "      No CRISPR data associated with selected compound"))
      plotdat <- reactvals$seldeprow_depmetrics %>%
        dplyr::filter(assay_type == "CRISPR")
      validate(need(nrow(plotdat) > 1,
                    "      No CRISPR data associated with selected compound"))
      gen_densplot(plotdat, "score", "ct", xlab = "CRISPR effect score")
    })
    # CRISPR boxplot by disease
    output$crispr_box_ds_dep <- renderPlot({
      print("crispr_box_ds_dep")
      validate(need(!is_empty(reactvals$seldeprow_depmetrics),
                    "      No CRISPR data associated with selected compound"))
      plotdat <- reactvals$seldeprow_depmetrics %>%
        dplyr::filter(assay_type == "CRISPR")
      validate(need(nrow(plotdat) > 1,
                    "      No CRISPR data associated with selected compound"))
      gen_boxplot(plotdat, "dg", "score", ylab = "CRISPR effect score")
    })
    # CRISPR boxplot by subtype
    output$crispr_box_st_dep <- renderPlot({
      print("crispr_box_st_dep")
      validate(need(!is_empty(reactvals$seldeprow_depmetrics),
                    "      No CRISPR data associated with selected protein"))
      plotdat <- reactvals$seldeprow_depmetrics %>%
        dplyr::filter(assay_type == "CRISPR")
      gen_boxplot(plotdat, "st", "score", ylab = "CRISPR effect score")
    })

    ### RNAi
    # RNAi density plot by disease
    output$rnai_dens_dep <- renderPlot({
      print("rnai_dens_dep")
      validate(need(!is_empty(reactvals$seldeprow_depmetrics),
                    "      No CRISPR data associated with selected protein"))
      plotdat <- reactvals$seldeprow_depmetrics %>%
        dplyr::filter(assay_type == "RNAi")
      validate(need(nrow(plotdat) > 1,
                    "      No RNAi data associated with selected protein"))
      gen_densplot(plotdat, "score", "ct", xlab = "RNAi effect score")
    })
    # RNAi boxplot by disease
    output$rnai_box_ds_dep <- renderPlot({
      print("rnai_box_ds_dep")
      validate(need(!is_empty(reactvals$seldeprow_depmetrics),
                    "      No RNAi data associated with selected protein"))
      plotdat <- reactvals$seldeprow_depmetrics %>%
        dplyr::filter(assay_type == "RNAi")
      validate(need(nrow(plotdat) > 1,
                    "      No RNAi data associated with selected protein"))
      gen_boxplot(plotdat, "dg", "score", ylab = "RNAi effect score")
    })
    # RNAi boxplot by subtype
    output$rnai_box_st_dep <- renderPlot({
      print("rnai_box_st_dep")
      validate(need(!is_empty(reactvals$seldeprow_depmetrics),
                    "      No RNAi data associated with selected protein"))
      plotdat <- reactvals$seldeprow_depmetrics %>%
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
    #            paste0(reactvals$sel_cpd_target_symbols, collapse="_"),
    #            "_", reactvals$ct1,  ".xlsx")
    #   },
    #   content = function(file) {
    #     metrics <- get_dep_metrics_cpd()
    #     if (is_empty(metrics)) return(NULL)
    #     print(head(metrics))
    #     writexl::write_xlsx(metrics, path=file)
    #   }
    # )
    
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
    
    ############################## Metabolites ####################################
    
    # volcano plot of differential metabolite levels
    output$metab_volcano <- renderPlot({
      print("metab_volcano")
      if(is_empty(reactvals$metablvl)) return(NULL)
      plotdat <- reactvals$metablvl
      ct1 <- reactvals$ct1
      ct2 <- reactvals$ct2
      plotdat <- plotdat %>%
        mutate(color = ifelse(nl10p < 1.3, "ns", ifelse(d > 0, "enriched", "depleted"))) %>%
        mutate(color = factor(color, levels= c("enriched", "ns", "depleted")))
      plotdat %>%
        ggplot(aes(x=d, y=nl10p)) +
        geom_vline(xintercept = 0, lty=2, alpha=0.4) +
        geom_point(aes(color=color), size = 2) +
        scale_color_manual(values=c("#994444","grey80","#447799")) +
        scale_x_continuous(name = paste0("\u0394 Metabolite level ", ct1, "/", ct2)) +
        scale_y_continuous(name = "Significance (-log10 p-value)") +
        theme_bw(base_size = 17) +
        theme(legend.title = element_blank(),
              legend.background = element_blank(),
              panel.grid = element_blank())
    })
    
    # get summary dat for selected - metab
    get_sel_metablvl <- reactive({
      print("get_sel_metablvl")
      metabsub <- reactvals$metablvl
      if(is_empty(metabsub)) return(NULL)
      if (!is_empty(reactvals$metab_volcano_brush)) {
        metabsub <- metabsub %>%
          brushedPoints(input$metab_volcano_brush)
      } 
      return(metabsub)
    })
    
    # output text for hovered point - metab
    output$metab_hover_text <- renderText({
      if (is_empty(input$metab_volcano_hover)) return(NULL)
      metabsub <- get_sel_metablvl() %>%
        nearPoints(input$metab_volcano_hover)
      paste0("Metab near cursor: ", paste0(unique(metabsub$metabolite), collapse=";"))
    })
    
    # summary table - metab
    output$metab_summary_dt <- DT::renderDataTable({
      print("metab_summary_dt")
      metabsub <- get_sel_metablvl()
      if (is.null(metabsub)) return(NULL)
      metabsub <- metabsub %>%
        arrange(d) %>%
        dplyr::select(-"p") %>%
        mutate(prcntl = prcntl*100) %>%
        mutate_if(is.numeric, round, digits=2) %>%
        dplyr::rename("Metabolite" = feature,
                      'FC' = d,
                      'Effect size' = r,
                      'Rank %' = prcntl,
                      'NL10 P-val' = nl10p)
      colnames(metabsub) <- sub("_", " ", colnames(metabsub))
      colnames(metabsub) <- sub("ct1", reactvals$ct1_label, colnames(metabsub))
      colnames(metabsub) <- sub("ct2", reactvals$ct2_label, colnames(metabsub))
      DT::datatable(
        data = metabsub,
        rownames = F,
        options = list(pageLength = 25),
        selection = list(mode = 'multiple', target = "row", selected = 1)
      )
    })
    output$dl_metab_summary_xls <- downloadHandler(
      filename = function() {
        paste0("Metabolites_summary_table_", reactvals$ct1, "_vs_", 
               reactvals$ct2, ".xlsx")
      },
      content = function(file) {
        metabsub <- get_sel_metablvl()
        if(is_empty(metabsub)) return(NULL)
        writexl::write_xlsx(metabsub, path=file)
      }
    )
    
    ######### metrics for selected metabolite #############
    
    # get selected metab
    get_selrow_metab <- reactive({
      print("get_selrow_metab")
      selmetablvl <- get_sel_metablvl()
      if(is.null(input$metab_summary_dt_rows_selected) |
         is.null(selmetablvl)) return(NULL)
      selmetablvl <- selmetablvl %>%
        dplyr::slice(input$metab_summary_dt_rows_selected)
      if (is_empty(selmetablvl)) return(NULL)
      reactvals$metab_feature <- unique(selmetablvl$feature)
      return(selmetablvl)
    })
    
    # get dependency metrics for selected compound targets
    get_metrics_selmetab <- reactive({
      print("get_metrics_selmetab")
      si <- reactvals$si
      ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
      metabrow <- get_selrow_metab()
      feature <- reactvals$metab_feature
      if (is_empty(feature)) return(NULL)
      # query db for selected compounds
      query <- paste0('SELECT * FROM metab_data ',
                      'WHERE feature IN ("',
                      paste0(feature, collapse = '", "'),
                      '") AND sample_id IN ("',
                      paste0(si$sample_id, collapse = '", "'),
                      '");')
      metabdat <- query_db(query)
      # format sample levels
      metabdat <- merge(metabdat, si, by="sample_id")
      metabdat$st <- ifelse(metabdat$ct==ct1, metabdat$ds_subtype, ct2)
      stlvls <- as.character(unique(metabdat$st))
      stlvls <- c(stlvls[!grepl(ct2, stlvls)], ct2)
      metabdat$dg <- ifelse(metabdat$ct==ct1, ct1, metabdat$ds_group)
      dglvls <- unique(metabdat$dg)
      dglvls <- c(dglvls[!grepl(ct1, dglvls)], ct1)
      metabdat <- metabdat %>%
        mutate(ct = factor(ct, levels=c(ct1, ct2)),
               st = factor(st, levels=c(stlvls)),
               dg = factor(dg, levels=c(dglvls))) %>%
        dplyr::filter(!is.na(ct) & !is.na(dg) & !dg=="Other")
      if (is_empty(metabdat) | nrow(metabdat) < 1) return(NULL)
      return(metabdat)
    })
    
    ### distplots - metab
    # density
    output$metab_dens <- renderPlot({
      print("crispr_dens_dep")
      plotdat <- get_metrics_selmetab()
      validate(need(!is_empty(plotdat),
                    "      No data associated with selected metabolite"))
      gen_densplot(plotdat, "level", "ct", xlab = "Metabolite level")
    })
    # boxplot by disease
    output$metab_box_ds <- renderPlot({
      print("crispr_box_ds_dep")
      plotdat <- get_metrics_selmetab()
      validate(need(!is_empty(plotdat),
                    "      No data associated with selected metabolite"))
      gen_boxplot(plotdat, "dg", "level", ylab = "Metabolite level")
    })
    # boxplot by subtype
    output$metab_box_st <- renderPlot({
      print("crispr_box_st_dep")
      plotdat <- get_metrics_selmetab()
      validate(need(!is_empty(plotdat),
                    "      No data associated with selected metabolite"))
      gen_boxplot(plotdat, "st", "level", ylab = "Metabolite level")
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
    #            paste0(reactvals$sel_cpd_target_symbols, collapse="_"),
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
    # output$ctd_css_dens_ui <- renderUI({
    #   print("ctd_css_dens_ui")
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset == "CTD")
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_densplot(plotdat, reactvals$selcpd_metric,
    #                     "ct", xlab = paste0(reactvals$selcpd_metric, " score")) +
    #     ggtitle("CTD")
    #   output$ctd_css_dens_plot <- renderPlot(p)
    #   plotOutput("ctd_css_dens_plot", height = 225)
    # })
    # output$gdsc_css_dens_ui <- renderUI({
    #   print("gdsc_css_dens_ui")
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset %in% c("GDSC1","GDSC2"))
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_densplot(plotdat, reactvals$selcpd_metric,
    #                     "ct", xlab = paste0(reactvals$selcpd_metric, " score")) +
    #     ggtitle("GDSC")
    #   output$gdsc_css_dens_plot <- renderPlot(p)
    #   plotOutput("gdsc_css_dens_plot", height = 225)
    # })
    # 
    # # drug sensitivity boxplots - disease groups
    # output$ctd_css_group_ui <- renderUI({
    #   print("ctd_css_group_ui")  
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset == "CTD")
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "dg", reactvals$selcpd_metric,
    #                    ylab = paste0(reactvals$selcpd_metric, " score"))
    #   output$ctd_css_group_plot <- renderPlot(p)
    #   plotOutput("ctd_css_group_plot", height = 225)
    # })
    # output$gdsc_css_group_ui <- renderUI({
    #   print("gdsc_css_group_ui")    
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset %in% c("GDSC1","GDSC2"))
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "dg", reactvals$selcpd_metric,
    #                    ylab = paste0(reactvals$selcpd_metric, " score"))
    #   output$gdsc_css_group_plot <- renderPlot(p)
    #   plotOutput("gdsc_css_group_plot", height = 225)
    # })
    # 
    # # drug sensitivity boxplots - subtypes
    # output$ctd_css_st_ui <- renderUI({
    #   print("ctd_css_st_ui")    
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset == "CTD")
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "st", reactvals$selcpd_metric,
    #                    ylab = paste0(reactvals$selcpd_metric, " score")) +
    #     ggtitle(paste0(reactvals$ct1, " subtypes:"))
    #   output$ctd_css_st_plot <- renderPlot(p)
    #   plotOutput("ctd_css_st_plot", height = 225)
    # })
    # output$gdsc_css_st_ui <- renderUI({
    #   print("gdsc_css_st_ui") 
    #   plotdat <- get_sel_cpd_metrics()
    #   if (is_empty(plotdat)) return(NULL)
    #   plotdat <- plotdat %>%
    #     filter(dataset %in% c("GDSC1","GDSC2"))
    #   if (nrow(plotdat) < 1) return(NULL)
    #   p <- gen_boxplot(plotdat, "st", reactvals$selcpd_metric,
    #                    ylab = paste0(reactvals$selcpd_metric, " score")) +
    #     ggtitle(paste0(reactvals$ct1, " subtypes:"))
    #   output$gdsc_css_st_plot <- renderPlot(p)
    #   plotOutput("gdsc_css_st_plot", height = 225)
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
    #                   logEC50, CSS, AAC_pred) %>%
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
    #         dplyr::select(genesymbol, ct, value) %>%
    #         mutate(ct = factor(ct, levels=c(ct1,ct2))) %>%
    #         nest(data = -genesymbol) %>%
    #         mutate(dExpr = map_dbl(data, ~filt_ttest(.x, metric="value", stat="statistic")),
    #                pExpr = map_dbl(data, ~filt_ttest(.x, metric="value", stat="p.value"))) %>%
    #         mutate(qExpr = p.adjust(pExpr)) %>%
    #         mutate(dExpr_rank = rank(dExpr)/length(dExpr),
    #                .before=3) %>%
    #         dplyr::select(-data)
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
    #         dplyr::select(gene, n, coef, se_coef, group) %>%
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
