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
library(Biobase)
library(SummarizedExperiment)

################################### Setup #####################################

# load data
load("dep_data/depdata.rda")

### tmp - quick fixes for data issues
data$si <- data$si %>% 
    dplyr::rename(COSMIC_ID = "COSMIC.identifier") %>%
    filter(!ds_type %in% c("Testicular", "Embryonal","NOS","Unknown",
                           "Engineered","Other","Eye","Gallbladder"))
data$gdsc_metrics <- data$gdsc_metrics %>%
    mutate(TARGET = ifelse(is.na(TARGET), "", TARGET))

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

# cell types
celltypes1 <- unique(c("B-cell",
                       "B-cell leukemia",
                       "B-cell lymphoma",
                       "Mantle cell lymphoma (MCL)",
                       "Plasma cell",
                       "T-cell",
                       "T-cell leukemia",
                       "T-cell lyphoma",
                       "Myeloid",
                       data$si$ds_type))
celltypes2 <- unique(c("Solid tumor",
                       "B-cell",
                       "B-cell leukemia",
                       "B-cell lymphoma",
                       "Mantle cell lymphoma (MCL)",
                       "Plasma cell",
                       "T-cell",
                       "T-cell leukemia",
                       "T-cell lyphoma",
                       "Myeloid",
                       data$si$ds_type))
solid_si <- data$si[which(!data$si$ds_type %in% c("B-cell",
                                                  "T-cell",
                                                  "Myeloid",
                                                  "Plasma cell")),]
haem_si <- data$si[which(data$si$ds_type %in% c("B-cell",
                                                "T-cell",
                                                "Myeloid",
                                                "Plasma cell")),]
bcell_leuk_si <- data$si[grep("B-ALL", data$si$ds_subtype),]
bcell_lymph_si <- data$si[which(data$si$ds_type == "B-cell" &
                                grepl("lymphoma", data$si$ds_subtype)),]
tcell_leuk_si <- data$si[grep("T-ALL", data$si$ds_subtype),]
tcell_lymph_si <- data$si[which(data$si$ds_type == "T-cell" &
                               grepl("lymphoma", data$si$ds_subtype)),]
mcl_si <- data$si[grep("MCL", data$si$ds_subtype),]

# run t-test only if sufficient observations
filt_ttest <- function(x, metric, stat = "statistic", grouping_var = "ct") {
    x <- x[, c(grouping_var, metric)] %>% na.omit()
    if(any(table(x[,grouping_var]) < 2) | length(unique(x[,metric,drop=T])) < 3) return(NA)
    tres <- t.test(reformulate(grouping_var, metric), x)
    return(tres[[stat]])
}

# # for testing
# reactvals <- list()
# reactvals$ct1_si <- data$si[grep("MCL", data$si$ds_subtype),]
# reactvals$ct2_si <- data$si[which(data$si$ds_type == "Colorectal"),]
# reactvals$ct1 <- "Mantle cell lymphoma (MCL)"
# reactvals$ct2 <- "Colorectal"

# reactive values
reactvals <- reactiveValues(drug_summary = data.frame(),
                            gene = "", 
                            ct1 = NULL, ct2 = NULL,
                            ct1_si = NULL, ct2_si = NULL)

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

# # JS formatting for DT
# js_substring <- c(
#     "function(data, type, row, meta) {",
#     "return type === 'display' && data.length > 30 ?",
#     "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
#     "}")
# js_scientific <- c(
#     "function(row, data, displayNum, index){",
#     "  var x = data[1];",
#     "  $('td:eq(1)', row).html(x.toExponential(3));",
#     "}")

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
    
    # celltype comparison choice
    tags$h3("Celltype comparison:"),
    div(style = "font-size:13px;", 
        fluidRow(column(width=3,uiOutput("ct1_choice_ui")),
                 column(width=3,uiOutput("ct2_choice_ui")))),
    
    # main UI
    tabsetPanel(
        tabPanel("Summary",
                 div(DT::dataTableOutput("gene_summary_dt"),
                     style = "font-size:90%")),
        tabPanel("Compounds",
                 
                 # scatter plot of average
                 tags$h3("Average sensitivity per protein/compound:"),
                 tags$h5("Drag a box around points to select compounds of interest"),
                 fluidRow(align="center", 
                          splitLayout(cellWidths = c("50%", "50%"),
                                      tags$h4("CTD"),
                                      tags$h4("GDSC"))),
                 fluidRow(align="center", 
                          splitLayout(cellWidths = c("50%", "50%"),
                                      plotOutput("ctd_scatter", height = 425, width = 425, 
                                                 brush = "ctd_scatter_brush", 
                                                 hover = "ctd_scatter_hover") %>% 
                                          withSpinner(color = "#D0E6EA99"),
                                      plotOutput("gdsc_scatter", height = 425, width = 425, 
                                                 brush = "gdsc_scatter_brush", 
                                                 hover = "gdsc_scatter_hover") %>% 
                                          withSpinner(color = "#D0E6EA99"))),
                 fluidRow(align="center", 
                          splitLayout(cellWidths = c("50%", "50%"),
                                      verbatimTextOutput("ctd_scatter_hover_text"),
                                      verbatimTextOutput("gdsc_scatter_hover_text"))),
                 
                 # metrics table
                 br(),br(),
                 tags$h3("Selected compounds:"),
                 tags$h5("Select a row to plot compound/target sensitivity accross cell lines below"),
                 div(DT::dataTableOutput("drug_summary_dt"),
                     style = "font-size:90%"),
                 downloadButton("dl_drug_summary_xls", label = "XLS",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 
                 # selected compound info
                 br(),br(),
                 tags$h3("Cell-line metrics for selected compound/target:"),
                 fluidRow(align="center", uiOutput("sel_cpd_text")),
                 
                 # rankings plot
                 tags$h4("Differential ranking:"),
                 tags$h5(paste0("Differential sensitivity ranks indicate how selectively sensitive the celltype of interest ",
                                "is relative to the control population (higher rank = more selective). Three primary metrics are ",
                                "considered - the CRISPR effect score, RNAi effect score, and drug sensitivity score (DSS).")),
                 fluidRow(align="center",
                          plotOutput("metrics_overview", height = 175, width = 400) 
                 ),
                 
                 # DSS plots
                 tags$h4("Drug sensitivity scores (DSS):"),
                 tags$h5(paste0("Drug sensitivity scores integrate the dose-response curve into a single metric - ",
                                "(higher score = more sensitive, range 0-100)")),
                 # CTD
                 fluidRow(splitLayout(cellWidths = c("20%","18%","27%","27%"),
                                      uiOutput("ctd_dss_dens_ui"),
                                      uiOutput("ctd_dss_group_ui"),
                                      uiOutput("ctd_dss_st_ui"),
                                      uiOutput("ctd_ec50_st_ui"))),
                 div(DT::dataTableOutput("cl_ctd_dt"),
                     style = "font-size:90%"),
                 downloadButton("dl_cl_ctd_xls", label = "XLS",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 # GDSC
                 fluidRow(splitLayout(cellWidths = c("20%","19%","26%","26%"),
                                      uiOutput("gdsc_dss_dens_ui"),
                                      uiOutput("gdsc_dss_group_ui"),
                                      uiOutput("gdsc_dss_st_ui"),
                                      uiOutput("gdsc_ec50_st_ui"))),
                 div(DT::dataTableOutput("cl_gdsc_dt"),
                     style = "font-size:90%"),
                 downloadButton("dl_cl_gdsc_xls", label = "XLS",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 
                 # CRISPR plots
                 tags$h4("CRISPR:"),
                 tags$h5(paste0("CRISPR effect scores indicate the effect of gene knock-out on viability - ",
                                "(lower score = more sensitive)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "28%"),
                                      plotOutput("crispr_dens", height = 225,
                                                 brush = "crispr_dens_brush"),
                                      plotOutput("crispr_box_ds", height = 225),
                                      plotOutput("crispr_box_st", height = 225))),
                 div(DT::dataTableOutput("cl_crisprdat"),
                     style = "font-size:90%"),
                 downloadButton("dl_cl_crispr_xls", label = "XLS",
                                style = "font-size:12px;height:30px;padding:5px;"),

                 # RNAi plots
                 tags$h4("RNAi:"),
                 tags$h5(paste0("RNAi effect scores indicate the effect of gene knock-down on viability - ",
                                "(lower scores = more sensitive)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "28%"),
                                      plotOutput("rnai_dens", height = 225),
                                      plotOutput("rnai_box_ds", height = 225),
                                      plotOutput("rnai_box_st", height = 225))),
                 div(DT::dataTableOutput("cl_rnaidat"),
                     style = "font-size:90%"),
                 downloadButton("dl_cl_rnai_xls", label = "XLS",
                                style = "font-size:12px;height:30px;padding:5px;")
        ),
        tabPanel("Dependency",

                 # scatter plot of average
                 tags$h3("Average CRISPR score per gene:"),
                 tags$h5("Drag a box around points to select genes of interest"),
                 fluidRow(align="center",
                          plotOutput("crispr_scatter", height = 450, width = 450,
                                     brush = "crispr_scatter_brush",
                                     hover = "crispr_scatter_hover") %>%
                              withSpinner(color = "#D0E6EA99")),
                 verbatimTextOutput("crispr_hover_text"),

                 # metrics table
                 br(),br(),
                 tags$h3("Selected genes:"),
                 tags$h5("Select a row to plot gene dependency accross cell lines below"),
                 div(DT::dataTableOutput("crispr_summary_dt"),
                     style = "font-size:90%")

        ),
        tabPanel("Expression",
                 
                 # scatter plot of average
                 tags$h3("Average protein expression:"),
                 tags$h5("Drag a box around points to select genes of interest"),
                 fluidRow(align="center",
                          plotOutput("prot_scatter", height = 450, width = 450,
                                     brush = "prot_scatter_brush",
                                     hover = "prot_scatter_hover") %>%
                              withSpinner(color = "#D0E6EA99")),
                 verbatimTextOutput("prot_hover_text"),
                 
                 # expr table
                 br(),br(),
                 tags$h3("Selected genes/proteins:"),
                 tags$h5("Select a row to plot expression accross cell lines below"),
                 div(DT::dataTableOutput("expr_summary_dt"),
                     style = "font-size:90%"),
                 
                 # expr plots
                 uiOutput("prot_expr_dens_ui")
        ),
        tabPanel("Prognosis",
                 
                 # scatter plot 
                 tags$h3("Gene hazard ratios:"),
                 tags$h5("Drag a box around points to select genes of interest"),
                 fluidRow(align="center",
                          plotOutput("prog_scatter", height = 450, width = 450,
                                     brush = "prog_scatter_brush",
                                     hover = "prog_scatter_hover") %>%
                              withSpinner(color = "#D0E6EA99")),
                 verbatimTextOutput("prog_scatter_hover_text"),
                 
                 # expr table
                 br(),br(),
                 tags$h3("Selected genes:"),
                 tags$h5("Select a row to plot prognosis effects per subtype/study"),
                 div(DT::dataTableOutput("prog_summary_dt"),
                     style = "font-size:90%"),

                 # expr plots
                 br(),br(),
                 tags$h3("HR effects:"),
                 uiOutput("prog_forest_plot_ui")
                 
        )
    )
)


################################### Server ####################################
server <- function(input, output, session) {
    
    # icon
    output$icon <- renderImage(list(src = "../hexagons/jelly.png",
                                    height = "95px", width = "85px"), 
                               deleteFile = F)
    
    # cell type comparison choice UI
    output$ct1_choice_ui <- renderUI({
        selectInput('ct1_choice', "Selected cell type", 
                    celltypes1, "B-cell")
    })
    output$ct2_choice_ui <- renderUI({
        selectInput('ct2_choice', "Control population", 
                    celltypes2, "Solid tumor")
    })
    
    # get sample info
    get_si <- reactive({
        print("get_si")
        ct1 <- input$ct1_choice
        ct2 <- input$ct2_choice
        if (is_empty(ct1) | is_empty(ct2)) return(NULL)
        if (ct1 == "Hematological") { ct1_si <- haem_si }
        else if (ct1 == "B-cell leukemia") { ct1_si <- bcell_leuk_si }
        else if (ct1 == "T-cell leukemia") { ct1_si <- tcell_leuk_si }
        else if (ct1 == "B-cell lymphoma") { ct1_si <- bcell_lymph_si }
        else if (ct1 == "T-cell lymphoma") { ct1_si <- tcell_lymph_si }
        else if (ct1 == "Mantle cell lymphoma (MCL)") { ct1_si <- mcl_si }
        else { ct1_si <- data$si[which(data$si$ds_type == ct1), ] }
        if (ct2 == "Solid tumor") { ct2_si <- solid_si }
        else if (ct2 == "Hematological") { ct2_si <- haem_si }
        else if (ct2 == "B-cell leukemia") { ct2_si <- bcell_leuk_si }
        else if (ct2 == "T-cell leukemia") { ct2_si <- tcell_leuk_si }
        else if (ct2 == "B-cell lymphoma") { ct2_si <- bcell_lymph_si }
        else if (ct2 == "T-cell lymphoma") { ct2_si <- tcell_lymph_si }
        else if (ct2 == "Mantle cell lymphoma (MCL)") { ct2_si <- mcl_si }
        else { ct2_si <- data$si[which(data$si$ds_type == ct2), ] }
        reactvals$ct1 <- ct1
        reactvals$ct2 <- ct2
        reactvals$ct1_si <- ct1_si
        reactvals$ct2_si <- ct2_si
        return(T)
    })
    
    ###################### Overall ########################
    
    concat_uniq <- function(x) { paste0(unique(x), collapse=";") }
    get_gene_summary <- reactive({
      cpd_summary <- get_drug_dep_summary() 
      if(is.null(cpd_summary)) return(NULL)
      gene_summary <- cbind(cpd_summary[,1:10],
                            data$gi[match(cpd_summary$genesymbol, data$gi$SYMBOL), ])
      gene_summary %>%
        dplyr::filter(Gene != "" & !is.na(Gene)) %>%
        group_by(genesymbol) %>%
        summarise(Gene_symbol = unique(genesymbol),
                  Gene_name = concat_uniq(GENENAME),
                  Entrez_ID = concat_uniq(ENTREZID),
                  Ensembl_ID = concat_uniq(ENSEMBL),
                  Drug_screens = concat_uniq(dataset),
                  Compounds = concat_uniq(treatmentid),
                  Overall_rank = mean(av_rank, na.rm=T),
                  dEC50_rank = mean(dEC50_rank, na.rm=T),
                  dDSS3_rank = mean(dDSS3_rank, na.rm=T),
                  dCRISPR_rank = mean(dCRISPR_rank, na.rm=T),
                  dRNAi_rank = mean(dRNAi_rank, na.rm=T),
                  CGC_tier = concat_uniq(CGC_tier),
                  CGC_cancers = concat_uniq(c(CGC_cancers_somatic, 
                                              CGC_cancers_germline)),
                  CGC_role = concat_uniq(CGC_role),
                  Transloc_partners = concat_uniq(CGC_translocation_partner),
                  Known_driver = concat_uniq(Mutability_cancers),
                  Mutfreq = concat_uniq(Mutability_mutfreq)) %>%
        mutate(across(8:12, round, 2))
    })
    output$gene_summary_dt <- DT::renderDataTable({
      print("gene_summary_dt")
      gene_summary <- get_gene_summary()
      if (is.null(gene_summary)) return(NULL)
      print(head(gene_summary))
      DT::datatable(
        data = gene_summary,
        rownames = F,
        colnames = gsub("_", " ", colnames(gene_summary)),
        selection = list(mode = 'single', target = "row", selected = 1),
        options = list(columnDefs = list(list(
          targets = c(1:6,12:16),
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 30 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
            "}")
        ))), callback = JS('table.page(3).draw(false);')
      )
    })
    
    ################################ Compounds ################################
    
    # subset CTD metrics
    get_ctd_metrics <- reactive({
        print("get_ctd_metrics")
        rtn <- get_si()
        if (is_empty(rtn)) return(NULL)
        ct1 <- data$ctd_metrics %>%
            dplyr::filter(!is.na(sampleid) &
                       sampleid %in% reactvals$ct1_si$sampleid) %>%
            mutate(ct = reactvals$ct1)
        ct2 <- data$ctd_metrics %>%
            dplyr::filter(!is.na(sampleid) & 
                       sampleid %in% reactvals$ct2_si$sampleid) %>%
            mutate(ct = reactvals$ct2)
        if (nrow(ct1) < 3 | nrow(ct2) < 3) return(NULL)
        ctd_metrics <- rbind(ct1, ct2)
        return(ctd_metrics)
    })
    
    # generate CTD summary data
    get_ctd_summary <- reactive({
        print("get_ctd_summary")
        ctd_metrics <- get_ctd_metrics()
        if (is_empty(ctd_metrics)) return(NULL)
        av <- ctd_metrics %>%
            group_by(ct, treatmentid) %>%
            summarise(avEC50 = mean(EC50, na.rm=T),
                      avDSS1 = mean(DSS1, na.rm=T),
                      avDSS2 = mean(DSS2, na.rm=T),
                      avDSS3 = mean(DSS3, na.rm=T),
                      genesymbol = unique(genesymbol)) %>%
            pivot_wider(names_from = ct,
                        values_from = c(avEC50, avDSS1, avDSS2, avDSS3))
        stats <- ctd_metrics %>%
            select(treatmentid, genesymbol, ct, DSS1, DSS2, DSS3, EC50) %>%
            group_by(ct) %>%
            mutate(zDSS1 = scale(DSS1)[,1],
                   zDSS2 = scale(DSS2)[,1],
                   zDSS3 = scale(DSS3)[,1]) %>%
            ungroup() %>%
            mutate(ct = factor(ct, levels=c(reactvals$ct1, reactvals$ct2))) %>%
            nest(data = -c(treatmentid, genesymbol)) %>%
            mutate(dEC50 = map_dbl(data, ~filt_ttest(.x, metric="EC50", stat="statistic")),
                   dDSS1 = map_dbl(data, ~filt_ttest(.x, metric="zDSS1", stat="statistic")),
                   dDSS2 = map_dbl(data, ~filt_ttest(.x, metric="zDSS2", stat="statistic")),
                   dDSS3 = map_dbl(data, ~filt_ttest(.x, metric="zDSS3", stat="statistic")),
                   pEC50 = map_dbl(data, ~filt_ttest(.x, metric="EC50", stat="p.value")),
                   pDSS1 = map_dbl(data, ~filt_ttest(.x, metric="zDSS1", stat="p.value")),
                   pDSS2 = map_dbl(data, ~filt_ttest(.x, metric="zDSS2", stat="p.value")),
                   pDSS3 = map_dbl(data, ~filt_ttest(.x, metric="zDSS3", stat="p.value"))) %>%
            mutate_at(vars(dEC50:dDSS3), ~ifelse(is.na(.), 0, .)) %>%
            mutate(pEC50 = p.adjust(pEC50),
                   pDSS1 = p.adjust(pDSS1),
                   pDSS2 = p.adjust(pDSS2),
                   pDSS3 = p.adjust(pDSS3)) %>%
            mutate(dataset = "CTD",
                   dEC50_rank = rank(-dEC50)/length(dEC50),
                   dDSS1_rank = rank(dDSS1)/length(dDSS1),
                   dDSS2_rank = rank(dDSS2)/length(dDSS2),
                   dDSS3_rank = rank(dDSS3)/length(dDSS3),
                   .before=3) %>%
            select(-data)
        ctd_summary <- merge(stats, av, by = c("treatmentid","genesymbol"))
        return(ctd_summary)
    })
    
    # subset GDSC metrics
    get_gdsc_metrics <- reactive({
        print("get_gdsc_metrics")
        rtn <- get_si()
        if (is.null(rtn)) return(NULL)
        ct1 <- data$gdsc_metrics %>%
            dplyr::filter(!is.na(COSMIC_ID) & COSMIC_ID %in% reactvals$ct1_si$COSMIC_ID) %>%
            mutate(ct = reactvals$ct1)
        ct2 <- data$gdsc_metrics %>%
            dplyr::filter(!is.na(COSMIC_ID) & COSMIC_ID %in% reactvals$ct2_si$COSMIC_ID) %>%
            mutate(ct = reactvals$ct2)
        gdsc_metrics <- rbind(ct1, ct2) %>% 
            dplyr::rename(treatmentid = DRUG_NAME) 
        return(gdsc_metrics)
    })
    
    # generate GDSC summary data
    get_gdsc_summary <- reactive({
        print("get_gdsc_summary")
        gdsc_metrics <- get_gdsc_metrics()
        if (is_empty(gdsc_metrics)) return(NULL)
        av <- gdsc_metrics %>%
            group_by(ct, treatmentid) %>%
            summarise(avEC50 = mean(EC50, na.rm=T),
                      avDSS1 = mean(DSS1, na.rm=T),
                      avDSS2 = mean(DSS2, na.rm=T),
                      avDSS3 = mean(DSS3, na.rm=T),
                      genesymbol = unique(TARGET)) %>%
            pivot_wider(names_from = ct,
                        values_from = c(avEC50, avDSS1, avDSS2, avDSS3))
        stats <- gdsc_metrics %>%
            dplyr::rename(genesymbol = TARGET) %>%
            select(treatmentid, genesymbol, ct, DSS1, DSS2, DSS3, EC50) %>%
            group_by(ct) %>%
            mutate(zDSS1 = scale(DSS1)[,1],
                   zDSS2 = scale(DSS2)[,1],
                   zDSS3 = scale(DSS3)[,1]) %>%
            ungroup() %>%
            mutate(ct = factor(ct, levels=c(reactvals$ct1, reactvals$ct2))) %>%
            nest(data = -c(treatmentid, genesymbol)) %>%
            mutate(dEC50 = map_dbl(data, ~filt_ttest(.x, metric="EC50", stat="statistic")),
                   dDSS1 = map_dbl(data, ~filt_ttest(.x, metric="zDSS1", stat="statistic")),
                   dDSS2 = map_dbl(data, ~filt_ttest(.x, metric="zDSS2", stat="statistic")),
                   dDSS3 = map_dbl(data, ~filt_ttest(.x, metric="zDSS3", stat="statistic")),
                   pEC50 = map_dbl(data, ~filt_ttest(.x, metric="EC50", stat="p.value")),
                   pDSS1 = map_dbl(data, ~filt_ttest(.x, metric="zDSS1", stat="p.value")),
                   pDSS2 = map_dbl(data, ~filt_ttest(.x, metric="zDSS2", stat="p.value")),
                   pDSS3 = map_dbl(data, ~filt_ttest(.x, metric="zDSS3", stat="p.value"))) %>%
            mutate_at(vars(dEC50:dDSS3), ~ifelse(is.na(.), 0, .)) %>%
            mutate(pEC50 = p.adjust(pEC50),
                   pDSS1 = p.adjust(pDSS1),
                   pDSS2 = p.adjust(pDSS2),
                   pDSS3 = p.adjust(pDSS3)) %>%
            mutate(dataset = "GDSC",
                   dEC50_rank = rank(-dEC50)/length(dEC50),
                   dDSS1_rank = rank(dDSS1)/length(dDSS1),
                   dDSS2_rank = rank(dDSS2)/length(dDSS2),
                   dDSS3_rank = rank(dDSS3)/length(dDSS3),
                   .before=3) %>%
            select(-data)
        gdsc_summary <- merge(stats, av, by = c("treatmentid","genesymbol"))
        return(gdsc_summary)
    })
    
    # add CRISPR/RNAi to drug metrics
    get_drug_dep_summary <- reactive({
        print("get_drug_dep_summary")
        ctd_summary <- get_ctd_summary()
        gdsc_summary <- get_gdsc_summary()
        if (is_empty(ctd_summary) & is_empty(gdsc_summary)) {
          return(NULL)
        } else if (is_empty(ctd_summary)) {
          drug_summary <- gdsc_summary %>%
            separate_rows(genesymbol, sep=";|, ") %>%
            distinct()
        } else if (is_empty(gdsc_summary)) {
          drug_summary <- ctd_summary %>%
            separate_rows(genesymbol, sep=";|, ") %>%
            distinct()
        } else {
          drug_summary <- rbind(ctd_summary, gdsc_summary) %>%
            separate_rows(genesymbol, sep=";|, ") %>%
            distinct()
        }
        genes <- unique(drug_summary$genesymbol)
        genes <- genes[!is.na(genes)]
        crispr <- data$crispr_es[which(fData(data$crispr_es)$genesymbol %in% genes), ]
        ct1idx <- which(crispr$DepMap_ID %in% reactvals$ct1_si$DepMap_ID)
        ct2idx <- which(crispr$DepMap_ID %in% reactvals$ct2_si$DepMap_ID)
        crispr <- data.frame(genesymbol = fData(crispr)$genesymbol,
                             avCRISPR_1 = rowMeans(exprs(crispr)[, ct1idx], na.rm=T),
                             avCRISPR_2 = rowMeans(exprs(crispr)[, ct2idx], na.rm=T)) %>%
            mutate(dCRISPR = avCRISPR_1 - avCRISPR_2) %>%
            mutate(dCRISPR_rank = rank(-dCRISPR)/length(dCRISPR))
        rnai <- data$rnai_es[which(fData(data$rnai_es)$genesymbol %in% genes), ]
        ct1idx <- which(rnai$CCLE_Name %in% reactvals$ct1_si$CCLE_Name)
        ct2idx <- which(rnai$CCLE_Name %in% reactvals$ct2_si$CCLE_Name)
        rnai <- data.frame(genesymbol = fData(rnai)$genesymbol,
                           avRNAi_1 = rowMeans(exprs(rnai)[, ct1idx], na.rm=T),
                           avRNAi_2 = rowMeans(exprs(rnai)[, ct2idx], na.rm=T)) %>%
            mutate(dRNAi = avRNAi_1 - avRNAi_2) %>%
            mutate(dRNAi_rank = rank(-dRNAi)/length(dRNAi))
        comb <- merge(crispr, rnai, by=c("genesymbol"), all=T)
        colnames(comb) <- sub("1", reactvals$ct1, colnames(comb))
        colnames(comb) <- sub("2", reactvals$ct2, colnames(comb))
        comb <- merge(drug_summary, comb, by=c("genesymbol"), all=T)
        comb <- comb[, c(grep("treatment|gene|dataset", colnames(comb)),
                         grep("_rank", colnames(comb)),
                         grep("^d[A-Za-z0-9]+$", colnames(comb)), 
                         grep("^p[A-Za-z0-9]+$", colnames(comb)), 
                         grep("^av", colnames(comb)))] 
        comb <- comb %>% 
            group_by(treatmentid, genesymbol) %>%
            mutate(av_rank = mean(c(dDSS3_rank, dCRISPR_rank, dRNAi_rank), na.rm=T), 
                        .before = 4) %>%
            ungroup() %>%
            arrange(-av_rank)
        return(comb)
    })
    
    # scatter plot of average CTD sensitivty scores per compound
    output$ctd_scatter <- renderPlot({
        print("ctd_scatter")
        selmetric <- "DSS3" ## tmp
        plotdat <- get_drug_dep_summary()
        if(is_empty(plotdat)) return(plotdat)
        sig <- plotdat[,grep(paste0("p",selmetric), colnames(plotdat))]
        plotdat <- plotdat %>%
            mutate(significance = ifelse(sig < 0.05, "signif","ns")) %>%
            filter(dataset=="CTD") %>%
            select_if(!grepl("genesymbol|RNAi|CRISPR", names(.))) %>%
            distinct()
        if (is_empty(plotdat)) return(NULL)
        tmp <- c(paste0("`av", selmetric, "_", reactvals$ct1, "`"),
                 paste0("`av", selmetric, "_", reactvals$ct2, "`"))
        plotdat %>%
            ggplot(aes_string(x=tmp[[1]], y=tmp[[2]])) +
            geom_point(aes(color=significance), size = 2) +
            geom_smooth(method="glm", se=F, color="grey40", lty=2) +
            scale_color_manual(values=c("#D0E6EA66","#71448199")) +
            scale_x_continuous(name = paste0(selmetric, " score: ", reactvals$ct1)) +
            scale_y_continuous(name = paste0(selmetric, " score: ", reactvals$ct2)) +
            theme_bw(base_size = 18) +
            theme(legend.position = "none") 
    })
    # scatter plot of average GDSC sensitivty scores per compound
    output$gdsc_scatter <- renderPlot({
        print("gdsc_scatter")
        selmetric <- "DSS3" ## tmp
        plotdat <- get_drug_dep_summary()
        if(is_empty(plotdat)) return(plotdat)
        sig <- plotdat[,grep(paste0("p",selmetric), colnames(plotdat))]
        plotdat <- plotdat %>%
            mutate(significance = ifelse(sig < 0.05, "signif","ns")) %>%
            filter(dataset=="GDSC") %>%
            select_if(!grepl("genesymbol|RNAi|CRISPR", names(.))) %>%
            distinct()
        if (is_empty(plotdat)) return(NULL)
        tmp <- c(paste0("`av", selmetric, "_", reactvals$ct1, "`"),
                 paste0("`av", selmetric, "_", reactvals$ct2, "`"))
        plotdat %>%
            ggplot(aes_string(x=tmp[[1]], y=tmp[[2]])) +
            geom_point(aes(color=significance), size = 2) +
            geom_smooth(method="glm", se=F, color="grey40", lty=2) +
            scale_color_manual(values=c("#D0E6EA66","#71448199")) +
            scale_x_continuous(name = paste0(selmetric, " score: ", reactvals$ct1)) +
            scale_y_continuous(name = paste0(selmetric, " score: ", reactvals$ct2)) +
            theme_bw(base_size = 18) +
            theme(legend.position = "none")
    })
    
    # handle brushing
    brushed_ctd <- reactive({
        print("brushed_ctd")
        if (is_empty(input$ctd_scatter_brush)) return(NULL)
        reactvals$gdsc_scatter_brush <- F
        reactvals$ctd_scatter_brush <- T
    })
    brushed_gdsc <- reactive({
        print("brushed_gdsc")
        if (is_empty(input$gdsc_scatter_brush)) return(NULL)
        reactvals$gdsc_scatter_brush <- T
        reactvals$ctd_scatter_brush <- F
    })
    
    # get all selected compounds 
    get_sel_drugsum <- reactive({
        print("get_sel_drugsum")
        subdat <- get_drug_dep_summary()
        if(is_empty(subdat)) return(NULL)
        subdat <- subdat %>%
            mutate_if(is.numeric, round, digits = 3)
        ctd_scatter_brush <- brushed_ctd()
        gdsc_scatter_brush <- brushed_gdsc()
        if (!is_empty(reactvals$gdsc_scatter_brush)) {
            if (reactvals$ctd_scatter_brush == T) {
                subdat <- subdat %>% 
                    filter(dataset == "CTD") %>%
                    brushedPoints(input$ctd_scatter_brush)
            } else if (reactvals$gdsc_scatter_brush == T) {
                subdat <- subdat %>% 
                    filter(dataset == "GDSC") %>%
                    brushedPoints(input$gdsc_scatter_brush)
            }
        } 
        return(subdat)
    }) 

    # output text for hovered
    output$ctd_scatter_hover_text <- renderText({
        if (is_empty(input$ctd_scatter_hover)) return(NULL)
        subdat <- get_drug_dep_summary() %>%
            filter(dataset=="CTD") %>%
            nearPoints(input$ctd_scatter_hover)
        paste0("Compounds near cursor: ", paste0(unique(subdat$treatmentid), collapse=";"))
    })
    output$gdsc_scatter_hover_text <- renderText({
        if (is_empty(input$gdsc_scatter_hover)) return(NULL)
        subdat <- get_drug_dep_summary() %>%
            filter(dataset=="GDSC") %>%
            nearPoints(input$gdsc_scatter_hover)
        paste0("Compounds near cursor: ", paste0(unique(subdat$treatmentid), collapse=";"))
    })
    
    # drug info table
    output$drug_summary_dt <- DT::renderDataTable({
        print("drug_summary_dt")
        drug_summary_sel <- get_sel_drugsum()
        if (is.null(drug_summary_sel)) return(NULL)
        drug_summary_sel <- drug_summary_sel %>%
            dplyr::rename(Compound = treatmentid,
                          `Target gene` = genesymbol,
                          `Data set` = dataset,
                          `Overall rank` = av_rank)
        DT::datatable(
            data = drug_summary_sel,
            rownames = F,
            colnames = gsub("_", " ", colnames(drug_summary_sel)),
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
    output$dl_drug_summary_xls <- downloadHandler(
        filename = function() {
            paste0("Drug_summary_table_", reactvals$ct1, "_vs_", 
                   reactvals$ct2, ".xlsx")
        },
        content = function(file) {
            drug_summary <- get_sel_drugsum()
            if(is_empty(drug_summary)) return(NULL)
            writexl::write_xlsx(drug_summary, path=file)
        })
    
    # text just highlighting the selected compound/gene
    output$sel_cpd_text <- renderUI({
        print("sel_cpd_text")
        plotdat <- get_sel_drugsum()
        if(is_empty(plotdat)) return(NULL)
        if(is_empty(input$drug_summary_dt_rows_selected)) return(NULL)
        plotdat <- plotdat %>%
            dplyr::slice(input$drug_summary_dt_rows_selected)
        cpd <- unique(plotdat$treatmentid)
        gene <- unique(plotdat$genesymbol)
        gene <- ifelse(gene=="", "none", gene)
        HTML(paste0("<h4>Selected compound: <b>", cpd,
                    "</b>, Targets: <b>", gene, "</b></h4>"))
    })

    # metrics overview
    output$metrics_overview <- renderPlot({
        print("metrics_overview")
        plotdat <- get_sel_drugsum()
        if(is_empty(plotdat)) return(NULL)
        if(is_empty(input$drug_summary_dt_rows_selected)) return(NULL)
        plotdat <- plotdat %>%
            dplyr::slice(input$drug_summary_dt_rows_selected) %>%
            select_at(1:10) %>%
            gather("metric", "rank", -c(1:3)) %>%
            mutate(rank = rank*100,
                   metric = factor(metric,
                                   levels = c("dEC50_rank","dDSS1_rank","dDSS2_rank","dDSS3_rank",
                                              "dCRISPR_rank","dRNAi_rank","av_rank"),
                                   labels = c("dEC50","dDSS1","dDSS2","dDSS3","dCRISPR",
                                              "dRNAi","average")))
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
                  plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
    })

    # get selected compound
    get_sel_cpd <- reactive({
      print("get_sel_cpd")
      drug_summary <- get_sel_drugsum()
      if(is.null(input$drug_summary_dt_rows_selected) |
         is.null(drug_summary)) return(NULL)
      drug_summary <- drug_summary %>%
        dplyr::slice(input$drug_summary_dt_rows_selected)
      if (is_empty(drug_summary)) return(NULL)
      cpd_id <- unique(drug_summary$treatmentid)
      formatted_cpd <- sub("(\\(.+\\))", "", cpd_id)
      formatted_cpd <- sub("\\:", "_", formatted_cpd)
      reactvals$cpd_id <- cpd_id
      reactvals$formatted_cpd <- formatted_cpd
      gene <- unique(drug_summary$genesymbol)
      reactvals$gene <- gene
      return(cpd_id)
    })
    
    # get drug metrics for selected row
    get_sel_ctd_metrics <- reactive({
        print("get_sel_ctd_metrics")
        cpd_id <- get_sel_cpd()
        if (is_empty(cpd_id)) return(NULL)
        ctd_metrics <- data$ctd_metrics %>%
            dplyr::filter(treatmentid == cpd_id)
        if (is_empty(ctd_metrics) | nrow(ctd_metrics) < 1) return(NULL)
        return(ctd_metrics)
    })
    get_sel_gdsc_metrics <- reactive({
        print("get_sel_gdsc_metrics")
        cpd_id <- get_sel_cpd()
        if (is_empty(cpd_id)) return(NULL)
        gdsc_metrics <- data$gdsc_metrics %>%
            dplyr::filter(DRUG_NAME == cpd_id)
        if (is_empty(gdsc_metrics) | nrow(gdsc_metrics) < 1) return(NULL)
        return(gdsc_metrics)
    })
    
    # ct filtering and formatting
    format_ct_ctd_metrics <- reactive({
        print("format_ct_ctd_metrics")
        metrics <- get_sel_ctd_metrics()
        if (is_empty(metrics)) return(NULL)
        ct1 <- metrics %>%
            dplyr::filter(!is.na(sampleid) & sampleid %in% reactvals$ct1_si$sampleid) %>%
            mutate(ct = reactvals$ct1)
        ct2 <- metrics %>%
            dplyr::filter(!is.na(sampleid) & sampleid %in% reactvals$ct2_si$sampleid) %>%
            mutate(ct = reactvals$ct2)
        plotdat <- rbind(ct1, ct2) %>%
            mutate(ct = factor(ct, levels=c(reactvals$ct1, reactvals$ct2)))
        return(plotdat)
    })
    format_ct_gdsc_metrics <- reactive({
        print("format_ct_gdsc_metrics")
        metrics <- get_sel_gdsc_metrics()
        if (is_empty(metrics)) return(NULL)
        ct1 <- metrics %>%
            dplyr::filter(!is.na(sampleid) & sampleid %in% reactvals$ct1_si$sampleid) %>%
            mutate(ct = reactvals$ct1)
        ct2 <- metrics %>%
            dplyr::filter(!is.na(sampleid) & sampleid %in% reactvals$ct2_si$sampleid) %>%
            mutate(ct = reactvals$ct2)
        plotdat <- rbind(ct1, ct2) %>%
            mutate(ct = factor(ct, levels=c(reactvals$ct1, reactvals$ct2)))
        return(plotdat)
    })
    
    # subtype formatting
    format_st_ctd_metrics <- reactive({
        print("format_st_ctd_metrics")
        metrics <- format_ct_ctd_metrics() 
        if (is_empty(metrics)) return(NULL)
        plotdat <- metrics %>%
            mutate(st = ifelse(ct == reactvals$ct1,
                               as.character(ds_subtype),
                               as.character(ct)))
        lvls <- unique(plotdat$st)
        lvls <- c(lvls[!grepl(reactvals$ct2, lvls)], reactvals$ct2)
        plotdat <- plotdat %>%
            mutate(st = factor(st, levels=lvls))
        return(plotdat)
    })
    format_st_gdsc_metrics <- reactive({
        print("format_st_gdsc_metrics")
        metrics <- format_ct_gdsc_metrics() 
        if (is_empty(metrics)) return(NULL)
        plotdat <- metrics %>%
            mutate(st = ifelse(ct == reactvals$ct1,
                               as.character(ds_subtype),
                               as.character(ct)))
        lvls <- unique(plotdat$st)
        lvls <- c(lvls[!grepl(reactvals$ct2, lvls)], reactvals$ct2)
        plotdat <- plotdat %>%
            mutate(st = factor(st, levels=lvls))
        return(plotdat)
    })
     
    # drug sensitivity density plot
    output$ctd_dss_dens_ui <- renderUI({
        print("ctd_dens_ui")
        plotdat <- format_ct_ctd_metrics()
        if (is_empty(plotdat)) return(NULL)
        p <- gen_densplot(plotdat, "DSS3", "ct", xlab = "DSS3 score") + 
            ggtitle("CTD")
        output$ctd_dss_dens_plot <- renderPlot(p)
        plotOutput("ctd_dss_dens_plot", height = 225)
    })
    output$gdsc_dss_dens_ui <- renderUI({
        print("gdsc_dens_ui")
        plotdat <- format_ct_gdsc_metrics()
        if (is_empty(plotdat)) return(NULL)
        p <- gen_densplot(plotdat, "DSS3", "ct", xlab = "DSS3 score") + 
            ggtitle("GDSC")
        output$gdsc_dss_dens_plot <- renderPlot(p)
        plotOutput("gdsc_dss_dens_plot", height = 225)
    })
    
    # drug sensitivity boxplots - disease groups
    output$ctd_dss_group_ui <- renderUI({
        metrics <- get_sel_ctd_metrics()
        if (is_empty(metrics)) return(NULL)
        plotdat <- metrics %>%
            dplyr::filter(!is.na(ds_group) & !ds_group %in% c("Unknown")) %>%
            mutate(ds_group = factor(ds_group))
        if (is_empty(plotdat) | nrow(plotdat) < 1) return(NULL)
        p <- gen_boxplot(plotdat, "ds_group", "DSS3", ylab = "DSS3 score")
        output$ctd_dss_group_plot <- renderPlot(p)
        plotOutput("ctd_dss_group_plot", height = 225)
    })
    output$gdsc_dss_group_ui <- renderUI({
        metrics <- get_sel_gdsc_metrics()
        if (is_empty(metrics)) return(NULL)
        plotdat <- metrics %>%
            dplyr::filter(!is.na(ds_group) & !ds_group %in% c("Unknown")) %>%
            mutate(ds_group = factor(ds_group))
        if (is_empty(plotdat) | nrow(plotdat) < 1) return(NULL)
        p <- gen_boxplot(plotdat, "ds_group", "DSS3", ylab = "DSS3 score")
        output$gdsc_dss_group_plot <- renderPlot(p)
        plotOutput("gdsc_dss_group_plot", height = 225)
    })

    # drug sensitivity boxplots - subtypes
    output$ctd_dss_st_ui <- renderUI({
        plotdat <- format_st_ctd_metrics()
        if (is_empty(plotdat)) return(NULL)
        p <- gen_boxplot(plotdat, "st", "DSS3", ylab = "DSS3 score") +
            ggtitle(paste0(reactvals$ct1, " subtypes:"))
        output$ctd_dss_st_plot <- renderPlot(p)
        plotOutput("ctd_dss_st_plot", height = 225)
    })
    output$gdsc_dss_st_ui <- renderUI({
        plotdat <- format_st_gdsc_metrics()
        if (is_empty(plotdat)) return(NULL)
        p <- gen_boxplot(plotdat, "st", "DSS3", ylab = "DSS3 score") +
            ggtitle(paste0(reactvals$ct1, " subtypes:"))
        output$gdsc_dss_st_plot <- renderPlot(p)
        plotOutput("gdsc_dss_st_plot", height = 225)
    })
    
    # EC50 boxplots - subtypes
    output$ctd_ec50_st_ui <- renderUI({
        plotdat <- format_st_ctd_metrics()
        if (is_empty(plotdat)) return(NULL)
        p <- gen_boxplot(plotdat, "st", "EC50", ylab = "EC50 (\U03BCM)") +
            ggtitle("")
        output$ctd_ec50_st_plot <- renderPlot(p)
        plotOutput("ctd_ec50_st_plot", height = 225)
    })
    output$gdsc_ec50_st_ui <- renderUI({
        plotdat <- format_st_gdsc_metrics()
        if (is_empty(plotdat)) return(NULL)
        p <- gen_boxplot(plotdat, "st", "EC50", ylab = "EC50 (\U03BCM)")
        output$gdsc_ec50_st_plot <- renderPlot(p)
        plotOutput("gdsc_ec50_st_plot", height = 225)
    })

    # data table for selected cpd
    output$cl_ctd_dt <- DT::renderDataTable({
        metrics <- format_ct_ctd_metrics() 
        if(is_empty(metrics) | is_empty(reactvals$ct1)) return(NULL)
        metrics <- metrics %>%
            select(sampleid, ds_group, ds_type, ds_subtype,
                   treatmentid, cpd_name, genesymbol,
                   EC50, DSS1, DSS2, DSS3, ec50_published, minc, maxc) %>%
            mutate_if(is.numeric, round, digits=3) %>%
            dplyr::rename("Cell line name" = sampleid,
                          "Disease group" = ds_group,
                          "Disease type" = ds_type,
                          "Disease subtype" = ds_subtype,
                          "Treatment ID" = treatmentid,
                          "Compound name" = cpd_name,
                          "Target gene" = genesymbol,
                          "Published EC50" = ec50_published,
                          "Min conc." = minc,
                          "Max conc." = maxc)
        DT::datatable(data = metrics,
                      rownames = F,
                      options = list(lengthMenu = c(5, 10, 25),
                                     pageLength = 5))
    })
    output$cl_gdsc_dt <- DT::renderDataTable({
        metrics <- format_ct_gdsc_metrics()
        if(is_empty(metrics) | is_empty(reactvals$ct1)) return(NULL)
        metrics <- metrics %>%
            select(CELL_LINE_NAME, ds_group, ds_type, ds_subtype,
                   DRUG_NAME, TARGET, TARGET_PATHWAY,
                   EC50, DSS1, DSS2, DSS3, minc, maxc) %>%
            mutate_if(is.numeric, round, digits=3) %>%
            dplyr::rename("Cell line name" = CELL_LINE_NAME,
                          "Disease group" = ds_group,
                          "Disease type" = ds_type,
                          "Disease subtype" = ds_subtype,
                          "Compound name" = DRUG_NAME,
                          "Target" = TARGET,
                          "Target pathway" = TARGET_PATHWAY,
                          "Min conc." = minc,
                          "Max conc." = maxc)
        DT::datatable(data = metrics,
                      rownames = F,
                      options = list(lengthMenu = c(5, 10, 25),
                                     pageLength = 5))
    })
    
    # DT download buttons
    output$dl_cl_ctd_xls <- downloadHandler(
        filename = function() {
            paste0("CTD_drug_metrics_", reactvals$formatted_cpd,
                   "_", reactvals$ct1,  ".xlsx")
        },
        content = function(file) {
            metrics <- format_ct_ctd_metrics()
            if (is_empty(metrics) | is_empty(reactvals$ct1)) return(NULL)
            writexl::write_xlsx(metrics, path=file)
        })
    output$dl_cl_gdsc_xls <- downloadHandler(
        filename = function() {
            paste0("GDSC_drug_metrics_", reactvals$formatted_cpd,
                   "_", reactvals$ct1,  ".xlsx")
        },
        content = function(file) {
            metrics <- format_ct_gdsc_metrics()
            if (is_empty(metrics) | is_empty(reactvals$ct1)) return(NULL)
            writexl::write_xlsx(metrics, path=file)
        })

    # get crispr data
    get_cpd_crispr_dat <- reactive({
        print("get_cpd_crispr_dat")
        gene <- reactvals$gene
        if(is_empty(gene) | gene=="") return(NULL)
        geneidx <- which(fData(data$crispr_es)$genesymbol == gene)
        df <- cbind(crispr_effect = exprs(data$crispr_es)[geneidx, ],
                    pData(data$crispr_es))
        if (is_empty(df)) return(NULL)
        return(df)
    })

    # filter CRISPR to just selected cell types
    get_filt_crispr_dat <- reactive({
        print("get_filt_crispr_dat")
        metrics <- get_cpd_crispr_dat()
        if (is_empty(metrics)) return(NULL)
        ct1 <- metrics %>%
            dplyr::filter(!is.na(DepMap_ID) & DepMap_ID %in% reactvals$ct1_si$DepMap_ID) %>%
            mutate(ct = reactvals$ct1)
        ct2 <- metrics %>%
            dplyr::filter(!is.na(DepMap_ID) & DepMap_ID %in% reactvals$ct2_si$DepMap_ID) %>%
            mutate(ct = reactvals$ct2)
        filt <- rbind(ct1, ct2) %>%
            mutate(ct = factor(ct, levels = c(reactvals$ct1, reactvals$ct2))) %>%
            mutate(st = ifelse(ct == reactvals$ct1,
                               as.character(ds_subtype),
                               as.character(ct)))
        lvls <- unique(filt$st)
        lvls <- c(lvls[!grepl(reactvals$ct2, lvls)], reactvals$ct2)
        filt %>%
            mutate(st = factor(st, levels=lvls))
    })

    # CRISPR density plot by disease
    output$crispr_dens <- renderPlot({
        print("crispr_dens")
        plotdat <- get_filt_crispr_dat()
        validate(need(!is_empty(plotdat),
                      "      No CRISPR data associated with selected compound"))
        gen_densplot(plotdat, "crispr_effect", "ct", xlab = "CRISPR effect score")
    })

    # CRISPR boxplot by disease
    output$crispr_box_ds <- renderPlot({
        print("crispr_box_ds")
        metrics <- get_cpd_crispr_dat()
        if(is_empty(metrics)) return(NULL)
        plotdat <- metrics %>%
            dplyr::filter(!is.na(DepMap_ID) & !ds_type %in% c("Unknown")) %>%
            mutate(ds_group = factor(ds_group))
        gen_boxplot(plotdat, "ds_group", "crispr_effect", ylab = "CRISPR effect score")
    })

    # CRISPR boxplot by subtype
    output$crispr_box_st <- renderPlot({
        print("crispr_box_st")
        plotdat <- get_filt_crispr_dat()
        if(is_empty(plotdat)) return(NULL)
        gen_boxplot(plotdat, "st", "crispr_effect", ylab = "CRISPR effect score")
    })

    # data table of CRISPR data
    output$cl_crisprdat <- DT::renderDataTable({
        dat <- get_cpd_crispr_dat()
        if (is_empty(dat) | is_empty(reactvals$ct1)) return(NULL)
        dat %>%
            dplyr::filter(!is.na(DepMap_ID) & DepMap_ID %in% reactvals$ct1_si$DepMap_ID)
        if (!is_empty(input$crispr_dens_brush)) {
            brushinfo <- input$crispr_dens_brush
            dat <- dat %>% dplyr::filter(crispr_effect > brushinfo$xmin &
                                  crispr_effect < brushinfo$xmax)
        }
        DT::datatable(
            data = dat,
            rownames = F,
            selection = "none",
            options = list(lengthMenu = c(5, 10, 25), pageLength = 5))
    })
    output$dl_cl_crispr_xls <- downloadHandler(
        filename = function() {
            ct <- reactvals$ct1
            dat <- get_cpd_crispr_dat()
            gene <- unique(dat$genesymbol)
            paste0("CRISPR_metrics_", gene, "_", ct,  ".xlsx")
        },
        content = function(file) {
            ct <- reactvals$ct1
            dat <- get_cpd_crispr_dat()
            if (is_empty(dat) | is(empty(ct))) return(NULL)
            dat %>%
                dplyr::filter(!is.na(DepMap_ID) & DepMap_ID %in% reactvals$ct1_si$DepMap_ID)
            writexl::write_xlsx(dat, path=file)
        })

    # get RNAi data
    get_rnai_dat <- reactive({
        print("get_rnai_dat")
        gene <- reactvals$gene
        if(is_empty(gene) | gene=="") return(NULL)
        geneidx <- which(fData(data$rnai_es)$genesymbol == gene)
        if(is_empty(geneidx)) return(NULL)
        cbind(rnai_effect = exprs(data$rnai_es)[geneidx, ],
              pData(data$rnai_es))
    })

    # filter RNAi to just selected cell types
    get_filt_rnai_dat <- reactive({
        print("get_filt_rnai_dat")
        metrics <- get_rnai_dat()
        if (is_empty(metrics)) return(NULL)
        ct1 <- metrics %>%
            dplyr::filter(!is.na(CCLE_Name) & CCLE_Name %in% reactvals$ct1_si$CCLE_Name) %>%
            mutate(ct = reactvals$ct1)
        ct2 <- metrics %>%
            dplyr::filter(!is.na(CCLE_Name) & CCLE_Name %in% reactvals$ct2_si$CCLE_Name) %>%
            mutate(ct = reactvals$ct2)
        filt <- rbind(ct1, ct2) %>%
            mutate(ct = factor(ct, levels = c(reactvals$ct1, reactvals$ct2))) %>%
            mutate(st = ifelse(ct == reactvals$ct1,
                               as.character(ds_subtype),
                               as.character(ct)))
        lvls <- unique(filt$st)
        lvls <- c(lvls[!grepl(reactvals$ct2, lvls)], reactvals$ct2)
        filt %>%
            mutate(st = factor(st, levels=lvls))
    })

    # RNAi density plot by disease
    output$rnai_dens <- renderPlot({
        print("rnai_dens")
        plotdat <- get_filt_rnai_dat()
        validate(need(!is_empty(plotdat),
                      "      No RNAi data associated with selected compound"))
        gen_densplot(plotdat, "rnai_effect", "ct", xlab = "RNAi effect score")
    })

    # RNAi boxplot by disease
    output$rnai_box_ds <- renderPlot({
        print("rnai_box_ds")
        metrics <- get_rnai_dat()
        if(is_empty(metrics)) return(NULL)
        plotdat <- metrics %>%
            dplyr::filter(!is.na(CCLE_Name) & !ds_type %in% c("Unknown")) %>%
            mutate(ds_group = factor(ds_group))
        gen_boxplot(plotdat, "ds_group", "rnai_effect", ylab = "RNAi effect score")
    })

    # RNAi boxplot by subtype
    output$rnai_box_st <- renderPlot({
        print("rnai_box_st")
        plotdat <- get_filt_rnai_dat()
        if(is_empty(plotdat)) return(NULL)
        gen_boxplot(plotdat, "st", "rnai_effect", ylab = "RNAi effect score")
    })

    # cell-line level RNAi data
    output$cl_rnaidat <- DT::renderDataTable({
        dat <- get_rnai_dat()
        if(is_empty(dat) | is_empty(reactvals$ct1)) return(NULL)
        dat %>%
            dplyr::filter(!is.na(CCLE_Name) & CCLE_Name %in% reactvals$ct1_si$CCLE_Name)
        DT::datatable(
            data = dat,
            rownames = F,
            selection = "none",
            options = list(lengthMenu = c(5, 10, 25), pageLength = 5))
    })
    output$dl_cl_rnai_xls <- downloadHandler(
        filename = function() {
            ct <- reactvals$ct1
            dat <- get_rnai_dat()
            gene <- unique(dat$genesymbol)
            paste0("RNAi_metrics_", gene, "_", ct,  ".xlsx")
        },
        content = function(file) {
            ct <- reactvals$ct1
            dat <- get_rnai_dat()
            if(is_empty(dat) | is_empty(ct)) return(NULL)
            dat %>%
                dplyr::filter(!is.na(CCLE_Name) & CCLE_Name %in% reactvals$ct1_si$CCLE_Name)
            writexl::write_xlsx(dat, path=file)
        })

    ############################## CRISPR ####################################

    # av crispr data
    get_crispr_av <- reactive({
        print("get_crispr_av")
        ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
        crispr <- data$crispr_es
        idx1 <- which(!is.na(crispr$DepMap_ID) & crispr$DepMap_ID %in% reactvals$ct1_si$DepMap_ID)
        idx2 <- which(!is.na(crispr$DepMap_ID) & crispr$DepMap_ID %in% reactvals$ct2_si$DepMap_ID)
        av <- data.frame(gene = fData(crispr)$genesymbol,
                         avCRISPR_1 = rowMeans(exprs(crispr)[, idx1]),
                         avCRISPR_2 = rowMeans(exprs(crispr)[, idx2])) %>%
            mutate(dCRISPR = avCRISPR_1 - avCRISPR_2)
        colnames(av) <- sub("_1$", paste0("_",ct1), colnames(av))
        colnames(av) <- sub("_2$", paste0("_",ct2), colnames(av))
        return(av)
    })
    get_crispr_dat <- reactive({
        print("get_crispr_dat")
        ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
        crispr <- data$crispr_es
        idx1 <- which(!is.na(crispr$DepMap_ID) & crispr$DepMap_ID %in% reactvals$ct1_si$DepMap_ID)
        idx2 <- which(!is.na(crispr$DepMap_ID) & crispr$DepMap_ID %in% reactvals$ct2_si$DepMap_ID)
        crispr$ct <- NA; crispr$ct[idx2] <- ct2; crispr$ct[idx1] <- ct1
        av <- get_crispr_av()
        crispr <- crispr[, c(idx1, idx2)]
        crisprsub <- crispr[abs(av$dCRISPR) > quantile(abs(av$dCRISPR), .85, na.rm=T), ]
        print("gather start")
        crisprdat <- exprs(crisprsub) %>%
            as.data.frame() %>%
            mutate(gene = fData(crisprsub)$genesymbol, .before=1) %>%
            gather("cl_id", "value", -1) %>%
            mutate(ct = pData(crisprsub)[match(.$cl_id, crisprsub$DepMap_ID),]$ct) %>%
            dplyr::filter(!is.na(value))
        print("gather end")
        return(crisprdat)
    })
    get_crispr_summary <- reactive({
        print("get_crispr_summary")
        ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
        crisprdat <- get_crispr_dat()
        av <- get_crispr_av()
        av <- av[, -which(colnames(av)=="dCRISPR")]
        print("stats start")
        stats <- crisprdat %>%
            select(gene, ct, value) %>%
            mutate(ct = factor(ct, levels=c(ct1,ct2))) %>%
            nest(data = -gene) %>%
            mutate(dCRISPR = map_dbl(data, ~filt_ttest(.x, metric="value", stat="statistic")),
                   pCRISPR = map_dbl(data, ~filt_ttest(.x, metric="value", stat="p.value"))) %>%
            mutate(qCRISPR = p.adjust(pCRISPR)) %>%
            mutate(dCRISPR_rank = rank(-dCRISPR)/length(dCRISPR),
                   .before=3) %>%
            select(-data)
        print("stats end")
        summary <- merge(stats, av, by = c("gene"))
        summary %>% arrange(-dCRISPR_rank)
    })

    # scatter plot of average sensitivty scores per compound
    output$crispr_scatter <- renderPlot({
        print("crispr_scatter")
        plotdat <- get_crispr_summary()
        if (is_empty(plotdat)) return(NULL)
        tmp <- c(paste0("`avCRISPR_", reactvals$ct1, "`"),
                 paste0("`avCRISPR_", reactvals$ct2, "`"))
        plotdat %>%
            ggplot(aes_string(x=tmp[[1]], y=tmp[[2]])) +
            geom_point(aes(color=ifelse(qCRISPR > 0.05, "B", "A")), size = 2) +
            geom_abline() +
            scale_x_continuous(name = paste0("avCRISPR score: ", reactvals$ct1)) +
            scale_y_continuous(name = paste0("avCRISPR score: ", reactvals$ct2)) +
            scale_color_manual(values=c("#71448199","#D0E6EA66")) +
            theme_bw(base_size = 18) +
            theme(legend.position = "none")
    })

    # get selected points from scatter
    get_sel_crispr <- reactive({
        print("get_sel_crispr")
        subdat <- get_crispr_summary() %>%
            mutate_if(grepl("avCRISPR|dCRISPR", names(.)), round, digits = 3)
        if (!is_empty(input$crispr_scatter_brush)) {
            subdat <- subdat %>%
                brushedPoints(input$crispr_scatter_brush)
        }
        return(subdat)
    })

    # output text for hovered
    output$crispr_hover_text <- renderText({
        if (is_empty(input$crispr_scatter_hover)) return(NULL)
        subdat <- get_crispr_summary() %>%
            nearPoints(input$crispr_scatter_hover)
        paste0("Genes near cursor: ", paste0(unique(subdat$gene), collapse=";"))
    })

    # dependency info table
    output$crispr_summary_dt <- DT::renderDataTable({
        print("crispr_summary_dt")
        crispr_summary_sel <- get_sel_crispr() %>%
            select(-pCRISPR)
        if (is.null(crispr_summary_sel)) return(NULL)
        DT::datatable(
            data = crispr_summary_sel,
            rownames = F,
            selection = list(mode = 'single', target = "row", selected = 1),
            options = list(rowCallback = JS(
                "function(row, data, displayNum, index){",
                "  var x = data[1];",
                "  $('td:eq(3)', row).html(x.toExponential(2));",
                "}"
            ))
        )
    })
    
    
    ############################## Expression ####################################
    
    # av proteomics data
    get_prot_av <- reactive({
        print("get_prot_av")
        ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
        prot <- data$sanger_prot_se
        idx1 <- which(!is.na(prot$Sanger_Model_ID) & prot$Sanger_Model_ID %in% reactvals$ct1_si$Sanger_Model_ID)
        idx2 <- which(!is.na(prot$Sanger_Model_ID) & prot$Sanger_Model_ID %in% reactvals$ct2_si$Sanger_Model_ID)
        av <- data.frame(genesymbol = rowData(prot)$name,
                         protein_id = rowData(prot)$ID,
                         avExpr_1 = rowMeans(assay(prot, "imputed")[, idx1]),
                         avExpr_2 = rowMeans(assay(prot, "imputed")[, idx2])) %>%
            mutate(dExpr = avExpr_1 - avExpr_2)
        colnames(av) <- sub("_1$", paste0("_",ct1), colnames(av))
        colnames(av) <- sub("_2$", paste0("_",ct2), colnames(av))
        return(av)
    })
    get_prot_dat <- reactive({
        print("get_prot_dat")
        ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
        prot <- data$sanger_prot_se
        idx1 <- which(!is.na(prot$Sanger_Model_ID) & prot$Sanger_Model_ID %in% reactvals$ct1_si$Sanger_Model_ID)
        idx2 <- which(!is.na(prot$Sanger_Model_ID) & prot$Sanger_Model_ID %in% reactvals$ct2_si$Sanger_Model_ID)
        prot$ct <- NA; prot$ct[idx2] <- ct2; prot$ct[idx1] <- ct1
        av <- get_prot_av()
        prot <- prot[, c(idx1, idx2)]
        protsub <- prot[abs(av$dExpr) > quantile(abs(av$dExpr), .5, na.rm=T), ]
        protdat <- assay(protsub, "imputed") %>%
            as.data.frame() %>%
            mutate(genesymbol = rowData(protsub)$name, .before=1) %>%
            gather("cl_id", "value", -1) %>%
            mutate(ct = colData(protsub)[match(.$cl_id, protsub$Sanger_Model_ID),]$ct) %>%
            dplyr::filter(!is.na(value))
        return(protdat)
    })
    get_prot_summary <- reactive({
        print("get_prot_summary")
        ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
        protdat <- get_prot_dat()
        av <- get_prot_av()
        av <- av[, -which(colnames(av)=="dExpr")]
        stats <- protdat %>%
            select(genesymbol, ct, value) %>%
            mutate(ct = factor(ct, levels=c(ct1,ct2))) %>%
            nest(data = -genesymbol) %>%
            mutate(dExpr = map_dbl(data, ~filt_ttest(.x, metric="value", stat="statistic")),
                   pExpr = map_dbl(data, ~filt_ttest(.x, metric="value", stat="p.value"))) %>%
            mutate(qExpr = p.adjust(pExpr)) %>%
            mutate(dExpr_rank = rank(dExpr)/length(dExpr),
                   .before=3) %>%
            select(-data)
        summary <- merge(stats, av, by = c("genesymbol"))
        summary <- summary %>% arrange(-dExpr_rank)
        return(summary)
    })
    
    # scatter plot of average sensitivity scores per compound
    output$prot_scatter <- renderPlot({
        plotdat <- get_prot_summary()
        if (is_empty(plotdat)) return(NULL)
        tmp <- c(paste0("`avExpr_", reactvals$ct1, "`"),
                 paste0("`avExpr_", reactvals$ct2, "`"))
        plotdat %>%
            ggplot(aes_string(x=tmp[[1]], y=tmp[[2]])) +
            geom_point(aes(color=ifelse(qExpr > 0.001, "B", "A")), size = 2) +
            geom_abline() +
            scale_x_continuous(name = paste0("Average protein level: ", reactvals$ct1)) +
            scale_y_continuous(name = paste0("Average protein level: ", reactvals$ct2)) +
            scale_color_manual(values=c("#71448199","#D0E6EA66")) +
            theme_bw(base_size = 18) +
            theme(legend.position = "none")
    })
    
    # get selected points from scatter
    get_sel_expr <- reactive({
        print("get_sel_expr")
        subdat <- get_prot_summary() %>%
            mutate_if(grepl("avExpr|dExpr", names(.)), round, digits = 3)
        if (!is_empty(input$prot_scatter_brush)) {
            subdat <- subdat %>%
                brushedPoints(input$prot_scatter_brush)
        }
        return(subdat)
    })
    
    # output text for hovered
    output$prot_hover_text <- renderText({
        if (is_empty(input$prot_scatter_hover)) return(NULL)
        subdat <- get_prot_summary() %>%
            nearPoints(input$prot_scatter_hover)
        paste0("Genes near cursor: ", paste0(unique(subdat$genesymbol), collapse=";"))
    })
    
    # expr info table
    output$expr_summary_dt <- DT::renderDataTable({
        expr_summary_sel <- get_sel_expr() %>%
            dplyr::select(-"pExpr")
        if (is.null(expr_summary_sel)) return(NULL)
        DT::datatable(
            data = expr_summary_sel,
            rownames = F,
            selection = list(mode = 'multiple', target = "row", selected = 1),
            options = list(rowCallback = JS(
                "function(row, data, displayNum, index){",
                "  var x = data[3];",
                "  $('td:eq(3)', row).html(x.toExponential(2));",
                "}"
            ))
        )
    })
    
    # get selected gene(s)
    get_sel_genes_expr <- reactive({
        print("get_sel_genes_expr")
        if(is_empty(input$expr_summary_dt_rows_selected)) return(NULL)
        expr_summary <- get_sel_expr() %>%
            dplyr::slice(input$expr_summary_dt_rows_selected)
        if (is_empty(expr_summary)) return(NULL)
        sel_genes <- unique(expr_summary$genesymbol)
        reactvals$sel_genes <- sel_genes
        return(sel_genes)
    })
    
    # get expr values for selected genes
    get_sel_prot_levels <- reactive({
        print("get_sel_prot_levels")
        sel_genes <- get_sel_genes_expr()
        if (is_empty(sel_genes)) return(NULL)
        prot <- data$sanger_prot_se
        prot <- prot[which(rowData(prot)$name %in% sel_genes),]
        exprdat  <- assay(prot, "imputed") %>%
            as.data.frame() %>%
            mutate(genesymbol = rowData(prot)$name, .before=1) %>%
            gather("sampleid", "expr", -1)
        if (is_empty(exprdat) | nrow(exprdat) < 1) return(NULL)
        return(exprdat)
    })
    
    # ct filtering and formatting
    format_prot_levels <- reactive({
        print("format_prot_levels_ct")
        exprdat <- get_sel_prot_levels()
        if (is_empty(exprdat)) return(NULL)
        ct1 <- exprdat %>%
            dplyr::filter(!is.na(sampleid) & sampleid %in% reactvals$ct1_si$Sanger_Model_ID) %>%
            mutate(ct = reactvals$ct1)
        ct2 <- exprdat %>%
            dplyr::filter(!is.na(sampleid) & sampleid %in% reactvals$ct2_si$Sanger_Model_ID) %>%
            mutate(ct = reactvals$ct2)
        plotdat <- rbind(ct1, ct2) %>%
            mutate(ct = factor(ct, levels=c(reactvals$ct1, reactvals$ct2)))
        matchidx <- match(plotdat$sampleid, data$si$Sanger_Model_ID)
        plotdat <- plotdat %>%
            mutate(ds_group = data$si[matchidx,]$ds_group,
                   ds_subtype = data$si[matchidx,]$ds_subtype) %>%
            mutate(st = ifelse(ct == reactvals$ct1,
                               as.character(ds_subtype),
                               as.character(ct)))
        lvls <- unique(plotdat$st)
        lvls <- c(lvls[!grepl(reactvals$ct2, lvls)], reactvals$ct2)
        plotdat <- plotdat %>%
            mutate(st = factor(st, levels=lvls))
        return(plotdat)
    })
    
    # drug sensitivity density plot
    output$prot_expr_dens_ui <- renderUI({
        print("prot_expr_dens_ui")
        plotdat <- format_prot_levels()
        if (is_empty(plotdat) | nrow(plotdat) < 1) return(NULL)
        sel_genes <- get_sel_genes_expr()
        if (length(sel_genes) == 1) {
            p <- gen_densplot(plotdat, "expr", "ct", xlab = "Protein levels [AU]") +
                ggtitle(paste0(sel_genes, " expression"))
            output$prot_expr_dens_plot <- renderPlot(p)
            fluidRow(splitLayout(cellWidths = c("20%"),
                                 plotOutput("prot_expr_dens_plot", height = 225)))
        } else if (length(sel_genes) < 101) {
            ctav <- plotdat %>% 
                group_by(genesymbol, ct) %>%
                summarise(av=mean(expr, na.rm=T)) %>%
                ungroup() %>%
                pivot_wider(names_from = ct,
                            values_from = av) %>%
                mutate(diff = .[,2] - .[,3]) %>%
                arrange(diff)
            p <- plotdat %>% 
                group_by(sampleid, genesymbol) %>% 
                summarise(ct = ct,
                          av = mean(expr, na.rm=T)) %>% 
                group_by(genesymbol) %>% 
                mutate(z = scale(av)[,1]) %>%
                ungroup %>%
                mutate(genesymbol = factor(genesymbol, levels=unique(ctav$genesymbol))) %>%
                ggplot(aes(x=sampleid, y=genesymbol, fill=z)) + 
                geom_tile() +
                facet_wrap(~ct, nrow=1, scales="free_x") +
                theme_bw() +
                theme(axis.text.x = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid = element_blank(),
                      panel.background = element_blank(),
                      panel.border = element_rect(fill=NA),
                      strip.background = element_blank(),
                      strip.text = element_text(size=12)) +
                scale_fill_gradient2(low="#426284cc", mid="#d0e6ea99", high="#714481",
                                     name="Relative expression [AU]") +
                scale_y_discrete(expand=c(0,0)) +
                xlab("Cell lines") + ylab("")
            output$prot_expr_hm_plot <- renderPlot(p)
            plotheight <- 200 + (length(sel_genes)*15)
            fluidRow(plotOutput("prot_expr_hm_plot", height = plotheight, width = 1200))
        } else {
            tags$h5("Too many genes selected to plot")
        }
    })
    
    ############################# Prognosis #####################################
    
    # subset prognosis data
    get_prog_dat <- reactive({
        print("get_prog_dat")
        prog <- data$prognosis_metrics
        ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
        ds <- unique(prog$ds_type); st <- unique(prog$ds_subtype)
        
        # if (ct2 == "Solid tumor") {
        #     ct2dat <- prog %>% filter(dataset!="meta" & !ds_group %in% c("Lymphoid", "Myeloid")) %>% mutate(group = "ct2")
        # } else {
        #     if (is.na(cat2)) return(NULL)
        #     ct2dat <- prog %>% filter(dataset!="meta" & (!!as.symbol(cat2)) == ct2) %>% mutate(group = "ct2")
        # }
        cat1 <- ifelse(ct1 %in% ds, "ds_type", ifelse(ct1 %in% st, "ds_subtype", NA))
        cat2 <- ifelse(ct2 %in% ds, "ds_type", ifelse(ct2 %in% st, "ds_subtype", NA))
        if (is.na(cat1) | is.na(cat2)) return(NULL)
        ct1dat <- prog %>% filter(dataset!="meta" & (!!as.symbol(cat1)) == ct1) %>% mutate(group = "ct1")
        ct1meta <- prog %>% filter(dataset=="meta" & (!!as.symbol(cat1)) == ct1) %>% mutate(group = "ct1")
        ct2meta <- prog %>% filter(dataset=="meta" & (!!as.symbol(cat2)) == ct2) %>% mutate(group = "ct2")
        progdat <- rbind(ct1dat, ct1meta, ct2meta) %>% distinct()
        return(progdat)
    })
    
    # get prognosis summary data
    get_prog_summary <- reactive({
        print("get_prog_summary")
        progdat <- get_prog_dat()
        if (is_empty(progdat)) return(NULL)
        progsummary <- progdat %>% 
            filter(dataset == "meta" & is.na(ds_subtype)) %>% 
            select(gene, n, coef, se_coef, group) %>%
            pivot_wider(names_from = group, values_from = c(n, coef, se_coef)) %>%
            group_by(gene) %>%
            summarise(coef_ct1 = coef_ct1,
                      coef_ct2 = coef_ct2,
                      se_coef_ct1 = se_coef_ct1,
                      se_coef_ct2 = se_coef_ct2,
                      dPrognosis = (coef_ct1 - coef_ct2) / sum(se_coef_ct1+se_coef_ct2)) %>%
            filter(!is.na(dPrognosis)) %>%
            mutate(prog_rank = rank(-dPrognosis)) %>%
            arrange(-prog_rank)
        return(progsummary)
    })
    
    # differential prognosis scatter
    output$prog_scatter <- renderPlot({
        print("prog_scatter")
        progsum <- get_prog_summary()
        if (is_empty(progsum)) return(NULL)
        progsum %>%
            ggplot(aes(x=coef_ct1, y=coef_ct2)) +
            geom_point(aes(color=ifelse(abs(dPrognosis) > 2, "A", "B")), size = 2) +
            geom_abline() +
            scale_x_continuous(name = paste0("Gene hazard ratio (log): ", reactvals$ct1)) +
            scale_y_continuous(name = paste0("Gene hazard ratio (log): ", reactvals$ct2)) +
            scale_color_manual(values=c("#71448199","#D0E6EA66")) +
            theme_bw(base_size = 18) +
            theme(legend.position = "none")
    })
    
    # get all selected genes - prognosis
    get_sel_progsummary <- reactive({
        print("get_sel_progsummary")
        progsum <- get_prog_summary()
        if(is_empty(progsum)) return(NULL)
        if (!is_empty(reactvals$gdsc_scatter_brush)) {
            progsum <- progsum %>%
                brushedPoints(input$prog_scatter_brush)
        }
        return(progsum)
    })

    # output text for hovered - prognosis
    output$prog_scatter_hover_text <- renderText({
        if (is_empty(input$prog_scatter_hover)) return(NULL)
        progsum <- get_prog_summary() %>%
            nearPoints(input$prog_scatter_hover)
        paste0("Genes near cursor: ", paste0(unique(progsum$gene), collapse=";"))
    })

    # DT table - prognosis summary
    output$prog_summary_dt <- DT::renderDataTable({
        print("prog_summary_dt")
        progsum <- get_sel_progsummary()
        if (is.null(progsum)) return(NULL)
        progsum <- progsum %>%
            select_at(c(1,7,6,2,4,3,5)) %>%
            mutate_at(c(3:7), ~round(.,3)) %>%
            dplyr::rename(Gene_symbol = gene,
                          dPrognosis_rank = prog_rank,
                          dPrognosis_score = dPrognosis,
                          logHR_ct1 = coef_ct1,
                          SE_ct1 = se_coef_ct1,
                          logHR_ct2 = coef_ct2,
                          SE_ct2 = se_coef_ct2)
        colnames(progsum) <- sub("ct1", reactvals$ct1, colnames(progsum))
        colnames(progsum) <- sub("ct2", reactvals$ct2, colnames(progsum))
        DT::datatable(
            data = progsum,
            rownames = F,
            colnames = gsub("_", " ", colnames(progsum)),
            selection = list(mode = 'single', target = "row", selected = 1))
    })
    # output$dl_prog_summary_xls <- downloadHandler(
    #     filename = function() {
    #         paste0("Prognosis_summary_table_", reactvals$ct1, "_vs_",
    #                reactvals$ct2, ".xlsx")
    #     },
    #     content = function(file) {
    #         drug_summary <- get_sel_progsummary()
    #         if(is_empty(drug_summary)) return(NULL)
    #         writexl::write_xlsx(drug_summary, path=file)
    #     })

    # get selected gene(s)
    get_sel_gene_prog <- reactive({
        print("get_sel_gene_prog")
        if(is_empty(input$prog_summary_dt_rows_selected)) return(NULL)
        progsum <- get_sel_progsummary() %>%
            dplyr::slice(input$prog_summary_dt_rows_selected)
        if (is_empty(progsum)) return(NULL)
        selegene <- unique(progsum$gene)
        return(selegene)
    })

    # Forest plot of selected gene
    output$prog_forest_plot_ui <- renderUI({
        progdat <- get_prog_dat()
        selgene <- get_sel_gene_prog()
        ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
        if (is_empty(progdat) | is_empty(selgene)) return(NULL)
        group_levels <- rev(c(ct2, ct1, unique(progdat$ds_subtype)))
        plotdat <- progdat %>%
            filter(gene == selgene) %>%
            mutate(col = ifelse(dataset == "meta", "meta", "study")) %>%
            mutate(grouping = ifelse(is.na(ds_subtype), ds_type, ds_subtype)) %>%
            mutate(grouping = factor(grouping, group_levels)) %>%
            mutate(label = paste0("HR:", round(exp_coef, 2), "\npval:", round(pval, 3)))
        textdat <- plotdat %>%
            filter(dataset == "meta")
        print(plotdat)
        xpos <- max(log(plotdat$upper_95ci)) + (max(log(plotdat$upper_95ci))/5)
        print(max(log(plotdat$upper_95ci)))
        p <- ggplot(plotdat) +
            geom_vline(xintercept = 0, lty=2, alpha=0.25) +
            geom_pointrange(data = filter(plotdat, group == "ct1" | dataset == "meta"),
                            aes(y = grouping, group = dataset, x = coef, color = col,
                                xmin = log(lower_95ci), xmax = log(upper_95ci)),
                            position = position_dodge(0.8), size = 0.9, linewidth = 0.9) +
            geom_text(data = textdat, x = xpos,
                      aes(y = grouping, label = label), 
                      size = 4.5, hjust = 0) +
            coord_cartesian(clip="off") +
            scale_color_manual(values = c( "#71448199", "#528199cc")) +
            ylab("") + xlab(paste0("logHR: high ", selgene, " expression")) +
            theme_bw(base_size = 16) +
            theme(panel.grid = element_blank(),
                  legend.position = "none",
                  plot.margin=unit(c(10,100,10,10), units = "pt"),
                  panel.border = element_rect(linewidth=0.25),
                  axis.ticks.y = element_blank())
        output$prog_forest_plot <- renderPlot(p)
        pheight <- 200 + (length(group_levels) * 40)
        fluidRow(align="center",
                 plotOutput("prog_forest_plot", height = pheight, width = 600))
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
