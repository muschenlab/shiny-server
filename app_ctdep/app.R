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
celltypes <- unique(c("Solid tumor",
                      "Hematological",
                      data$ctd_metrics$disease))
celltypes <- celltypes[celltypes != "Solid" & !is.na(celltypes)]
reactvals <- reactiveValues(drug_summary = data.frame(),
                            gene = "",
                            ct1 = "", ct2 = "",
                            celltypes = celltypes)

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
              "#528199cc","#426284cc",
              "#eccad3cc","#eed9e0ff","#d0e6ea99")
pickcols <- function(pal, n) {
    idx <- c(round(seq(1, length(pal), length(pal)/(n-1))), length(pal))
    pal[idx]
}

# run t-test only if sufficient observations
filt_ttest <- function(x, metric, stat = "statistic", grouping_var = "disease") {
    if(any(table(x$disease) < 2)) return(NA)
    tres <- t.test(reformulate(grouping_var, metric), x)
    return(tres[[stat]])
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
    
    # celltype comparison choice
    tags$h3("Celltype comparison:"),
    div(style = "font-size:13px;", 
        fluidRow(column(width=3,uiOutput("ct1_choice_ui")),
                 column(width=3,uiOutput("ct2_choice_ui")))),
    
    # main UI
    tabsetPanel(
        tabPanel("Compounds",
                 
                 # scatter plot of average
                 tags$h3("Average sensitivity per protein/compound:"),
                 tags$h5("Drag a box around points to select compounds of interest"),
                 fluidRow(align="center",
                          plotOutput("sens_scatter", height = 450, width = 450, 
                                     brush = "scatter_brush") %>% 
                              withSpinner(color = "#D0E6EA99")) ,
                 
                 # metrics table
                 br(),br(),
                 tags$h3("Selected compounds:"),
                 tags$h5("Select a row to plot compound/target sensitivity accross cell lines below"),
                 div(DT::dataTableOutput("drug_summary_dt"),
                     style = "font-size:90%"),
                 
                 # headers
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
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "20%", "22%"),
                                      plotOutput("dss_dens", height = 225), 
                                      plotOutput("dss_box_ds", height = 225), 
                                      plotOutput("dss_box_st", height = 225),
                                      plotOutput("ec50_plot", height = 225))),
                 div(DT::dataTableOutput("cl_drugdat"),
                     style = "font-size:90%"),
                 downloadButton("dl_cl_drug_xls", label = "XLS",
                                style = "font-size:12px;height:30px;padding:5px;"),
                 
                 # CRISPR plots
                 tags$h4("CRISPR:"),
                 tags$h5(paste0("CRISPR effect scores indicate the effect of gene knock-out on viability - ",
                                "(lower score = more sensitive)")),
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "20%"),
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
                 fluidRow(splitLayout(cellWidths = c("25%", "20%", "20%"),
                                      plotOutput("rnai_dens", height = 225), 
                                      plotOutput("rnai_box_ds", height = 225), 
                                      plotOutput("rnai_box_st", height = 225))),
                 div(DT::dataTableOutput("cl_rnaidat"),
                     style = "font-size:90%"),
                 downloadButton("dl_cl_rnai_xls", label = "XLS",
                                style = "font-size:12px;height:30px;padding:5px;")
        ),
        tabPanel("CRISPR",
                 
                 # scatter plot of average
                 tags$h3("Average CRISPR score per gene:"),
                 tags$h5("Drag a box around points to select genes of interest"),
                 fluidRow(align="center",
                          plotOutput("crispr_scatter", height = 450, width = 450, 
                                     brush = "crispr_scatter_brush") %>% 
                              withSpinner(color = "#D0E6EA99")),
                 
                 # metrics table
                 br(),br(),
                 tags$h3("Selected genes:"),
                 tags$h5("Select a row to plot gene dependency accross cell lines below"),
                 div(DT::dataTableOutput("crispr_summary_dt"),
                     style = "font-size:90%")
                 
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
                    reactvals$celltypes,
                    "B-cell")
    })
    output$ct2_choice_ui <- renderUI({
        selectInput('ct2_choice', "Control population", 
                    reactvals$celltypes,
                    "Solid tumor")
    })
    
    ################################ Compounds ################################
    
    # subset drug metrics
    get_drug_metrics <- reactive({
        print("get_drug_metrics")
        ct1 <- input$ct1_choice
        ct2 <- input$ct2_choice
        ct1 <- ifelse(ct1 == "Solid tumor", "Solid", ct1)
        ct2 <- ifelse(ct2 == "Solid tumor", "Solid", ct2)
        drug_metrics <- data$ctd_metrics %>%
            filter(disease %in% c(ct1,ct2))
        reactvals$ct1 <- ct1
        reactvals$ct2 <- ct2
        return(drug_metrics)
    })
    
    # generate summary data
    get_drug_summary <- reactive({
        print("get_drug_summary")
        drug_metrics <- get_drug_metrics() 
        if (is_empty(drug_metrics)) return(NULL)
        av <- drug_metrics %>%
            group_by(disease, treatmentid) %>%
            summarise(avEC50 = mean(EC50, na.rm=T),
                      avDSS1 = mean(DSS1, na.rm=T),
                      avDSS2 = mean(DSS2, na.rm=T),
                      avDSS3 = mean(DSS3, na.rm=T),
                      genesymbol = unique(gene_symbol_of_protein_target)) %>%
            pivot_wider(names_from = disease,
                        values_from = c(avEC50, avDSS1, avDSS2, avDSS3))
        stats <- drug_metrics %>%
            select(treatmentid, gene_symbol_of_protein_target, disease, DSS1, DSS2, DSS3, EC50) %>%
            rename(genesymbol = gene_symbol_of_protein_target) %>%
            mutate(disease = factor(disease, levels=c(reactvals$ct1,
                                                      reactvals$ct2))) %>%
            nest(-c(treatmentid, genesymbol)) %>%
            mutate(dEC50 = map_dbl(data, ~filt_ttest(.x, metric="EC50", stat="statistic")),
                   dDSS1 = map_dbl(data, ~filt_ttest(.x, metric="DSS1", stat="statistic")),
                   dDSS2 = map_dbl(data, ~filt_ttest(.x, metric="DSS2", stat="statistic")),
                   dDSS3 = map_dbl(data, ~filt_ttest(.x, metric="DSS3", stat="statistic"))) %>%
            mutate(dEC50_rank = rank(-dEC50)/length(dEC50),
                   dDSS1_rank = rank(dDSS1)/length(dDSS1),
                   dDSS2_rank = rank(dDSS2)/length(dDSS2),
                   dDSS3_rank = rank(dDSS3)/length(dDSS3),
                   .before=3) %>%
            select(-data)
        drug_summary <- merge(stats, av, by = c("treatmentid","genesymbol"))
        print(head(drug_summary))
        reactvals$drug_summary <- drug_summary
        return(drug_summary)
    })
    
    # add CRISPR/RNAi to drug metrics
    get_drug_dep_summary <- reactive({
        print("get_drug_dep_summary")
        drug_summary <- get_drug_summary()
        if (is_empty(drug_summary)) return(NULL)
        drug_summary <- drug_summary %>%
            separate_rows(genesymbol, sep=";") %>%
            distinct()
        genes <- unique(drug_summary$genesymbol)
        genes <- genes[!is.na(genes)]
        crispr <- data$crispr_es[which(fData(data$crispr_es)$genesymbol %in% genes), ]
        ct1idx <- which(crispr$disease == reactvals$ct1)
        ct2idx <- which(crispr$disease == reactvals$ct2)
        crispr <- data.frame(genesymbol = fData(crispr)$genesymbol,
                             avCRISPR_1 = rowMeans(exprs(crispr)[, ct1idx], na.rm=T),
                             avCRISPR_2 = rowMeans(exprs(crispr)[, ct2idx], na.rm=T)) %>%
            mutate(dCRISPR = avCRISPR_1 - avCRISPR_2) %>%
            mutate(dCRISPR_rank = rank(-dCRISPR)/length(dCRISPR))
        rnai <- data$rnai_es[which(fData(data$rnai_es)$genesymbol %in% genes), ]
        ct1idx <- which(rnai$disease == reactvals$ct1)
        ct2idx <- which(rnai$disease == reactvals$ct2)
        rnai <- data.frame(genesymbol = fData(rnai)$genesymbol,
                           avRNAi_1 = rowMeans(exprs(rnai)[, ct1idx], na.rm=T),
                           avRNAi_2 = rowMeans(exprs(rnai)[, ct2idx], na.rm=T)) %>%
            mutate(dRNAi = avRNAi_1 - avRNAi_2) %>%
            mutate(dRNAi_rank = rank(-dRNAi)/length(dRNAi))
        comb <- merge(crispr, rnai, by=c("genesymbol"), all=T)
        colnames(comb) <- sub("1", reactvals$ct1, colnames(comb))
        colnames(comb) <- sub("2", reactvals$ct2, colnames(comb))
        comb <- merge(drug_summary, comb, by=c("genesymbol"), all=T)
        comb <- comb[, c(grep("treat|gene", colnames(comb)), grep("_rank", colnames(comb)),
                         grep("^d[A-Za-z0-9]+$", colnames(comb)), grep("^av", colnames(comb)))]
        comb %>% 
            group_by(treatmentid, genesymbol) %>%
            mutate(av_rank = mean(c(dDSS3_rank, dCRISPR_rank, dRNAi_rank), na.rm=T), 
                        .before = 3) %>%
            ungroup() %>%
            arrange(-av_rank)
    })
    
    # scatter plot of average sensitivty scores per compound
    output$sens_scatter <- renderPlot({
        print("sens_scatter")
        selmetric <- "avDSS3" ## tmp
        plotdat <- get_drug_dep_summary()
        plotdat <- plotdat %>%
            select_if(grepl(selmetric, names(.))) %>%
            distinct()
        if (is_empty(plotdat)) return(NULL)
        tmp <- c(paste0("`", selmetric, "_", reactvals$ct1, "`"),
                 paste0("`", selmetric, "_", reactvals$ct2, "`"))
        plotdat %>%
            ggplot(aes_string(x=tmp[[1]], y=tmp[[2]])) +
            geom_point(color = "#71448199", size = 2) +
            geom_abline() +
            scale_x_continuous(name = paste0(selmetric, " score: ", reactvals$ct1)) +
            scale_y_continuous(name = paste0(selmetric, " score: ", reactvals$ct2)) +
            theme_bw(base_size = 18) +
            theme()
    })
    
    # get selected points from scatter
    get_sel_drugsum <- reactive({
        print("get_sel_drugsum")
        if (is_empty(input$scatter_brush)) {
            get_drug_dep_summary() %>%
                mutate_if(is.numeric, round, digits = 3)
        } else {
            get_drug_dep_summary() %>%
                brushedPoints(input$scatter_brush) %>%
                mutate_if(is.numeric, round, digits = 3)
        }
    }) 

    # drug info table
    output$drug_summary_dt <- DT::renderDataTable({
        print("drug_summary_dt")
        drug_summary_sel <- get_sel_drugsum() %>%
            select_at(1:10) %>%
            rename(Compound = treatmentid,
                   `Target gene` = genesymbol)
        if (is.null(drug_summary_sel)) return(NULL)
        DT::datatable(
            data = drug_summary_sel,
            rownames = F,
            selection = list(mode = 'single', target = "row", selected = 1),
            options = list(columnDefs = list(list(
                targets = 0:3,
                render = JS(
                    "function(data, type, row, meta) {",
                    "return type === 'display' && data.length > 30 ?",
                    "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                    "}")
            ))), callback = JS('table.page(3).draw(false);')
        )
    })
    
    # text just highlighting the selected compound/gene
    output$sel_cpd_text <- renderUI({
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
            select_at(1:9) %>%
            gather("metric", "rank", -c(1:2)) %>%
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
    
    # get drug metrics for selected row
    get_selcpd_metrics <- reactive({
        print("get_selcpd_metrics")
        if(is_empty(input$drug_summary_dt_rows_selected)) return(NULL)
        drug_summary <- get_sel_drugsum() %>%
            slice(input$drug_summary_dt_rows_selected)
        if (is_empty(drug_summary)) return(NULL)
        metrics <- data$ctd_metrics %>%
            filter(treatmentid == drug_summary$treatmentid)
        if (is_empty(metrics)) return(NULL)
        if (reactvals$ct1 == "B-cell") {
            fctlvls <- rev(c("Solid","B-cell","T-cell","Myeloma","Myeloid"))
            metrics <- metrics %>%
                mutate(disease = factor(disease, levels = fctlvls))
        }
        if (reactvals$ct1 == "T-cell") {
            fctlvls <- rev(c("Solid","T-cell","B-cell","Myeloma","Myeloid"))
            metrics <- metrics %>%
                mutate(disease = factor(disease, levels = fctlvls))
        }
        return(metrics)
    })
     
    # subtype factor level formatting
    get_selcpd_st_metrics <- reactive({
        print("get_selcpd_st_metrics")
        plotdat <- get_selcpd_metrics()
        if (is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            filter(disease %in% c(reactvals$ct1, reactvals$ct2)) %>%
            mutate(st = ifelse(disease == reactvals$ct1,
                               as.character(subtype), as.character(disease))) %>%
            mutate(st = ifelse(is.na(st), "NOS", st)) 
        # lvls <- unique(plotdat$st)
        # lvls <- c(lvls[!grepl(reactvals$ct2, lvls)], reactvals$ct2) ## add factors to MAE and remove this
        # plotdat <- plotdat %>%
        #     mutate(st = factor(st, levels=lvls))
        # B/T cell specific stuff can be removed once fixed factor levels in MAE
        if (reactvals$ct1 == "B-cell") {
            plotdat <- plotdat %>%
                mutate(st = factor(st, levels = rev(c("Solid","B-ALL","CLL","Mantle cell",
                                                      "Burkitts","DLBCL","other NHL", "Hodgkin","NOS"))))
        }
        if (reactvals$ct1 == "T-cell") {
            plotdat <- plotdat %>%
                mutate(st = factor(st, levels = rev(c("Solid","T-ALL","T-cell lymphoma","NOS"))))
        }
        return(plotdat)
    })

    # DSS1 density plot by disease
    output$dss_dens <- renderPlot({
        plotdat <- get_selcpd_metrics()
        if (is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            filter(disease %in% c(reactvals$ct1, "Solid")) %>%
            droplevels()
        gen_densplot(plotdat, "DSS3", "disease", xlab = "DSS3 score")
    })

    # DSS1 boxplot by disease
    output$dss_box_ds <- renderPlot({
        plotdat <- get_selcpd_metrics()
        if (is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            filter(!is.na(disease)) %>%
            droplevels()
        gen_boxplot(plotdat, "disease", "DSS3", ylab = "DSS3 score")
    })

    # DSS1 boxplot by subtype
    output$dss_box_st <- renderPlot({
        plotdat <- get_selcpd_st_metrics() %>%
            droplevels()
        if (is_empty(plotdat)) return(NULL)
        gen_boxplot(plotdat, "st", "DSS3", ylab = "DSS3 score")
    })

    # EC50 boxplot by subtype
    output$ec50_plot <- renderPlot({
        plotdat <- get_selcpd_st_metrics() %>%
            droplevels()
        if (is_empty(plotdat)) return(NULL)
        gen_boxplot(plotdat, "st", "EC50", ylab = "EC50 (\U03BCM)") +
            scale_y_log10(labels=plain) #+
            # theme(axis.text.x = element_text(angle=45, hjust=1))
    })
    
    # data table for selected cpd
    output$cl_drugdat <- DT::renderDataTable({
        selcpd_metrics <- get_selcpd_metrics()
        dat <- data$ctd_metrics %>%
            dplyr::filter(treatmentid == unique(selcpd_metrics$treatmentid)) %>%
            dplyr::filter(disease %in% c(reactvals$ct1)) %>%
            select(sampleid, treatmentid, disease, subtype, EC50, DSS1, DSS2, DSS3, minc, maxc, slope)
        if (is.null(dat)) return(NULL)
        DT::datatable(
            data = dat, 
            rownames = F,
            options = list(lengthMenu = c(5, 10, 25), pageLength = 5))
    })
    output$dl_cl_drug_xls <- downloadHandler(
        filename = function() {
            selcpd_metrics <- get_selcpd_metrics()
            ct <- reactvals$ct1
            cpd <- unique(selcpd_metrics$treatmentid)
            cpd <- sub("(\\(.+\\))", "", cpd)
            cpd <- sub("\\:", "_", cpd)
            paste0("Drug_metrics_", cpd, "_", ct,  ".xlsx")
        },
        content = function(file) {
            selcpd_metrics <- get_selcpd_metrics() %>%
                select(EC50, DSS1, DSS2, DSS3, minc, maxc, slope, sampleid, treatmentid, disease, subtype)  %>%
                dplyr::filter(disease %in% c(reactvals$ct1))
            writexl::write_xlsx(dat, path=file)
        })

    # get crispr data
    get_cpd_crispr_dat <- reactive({
        cpdsum <- get_sel_drugsum()
        if(is_empty(cpdsum)) return(NULL)
        cpdsum <- cpdsum %>%
            dplyr::slice(input$drug_summary_dt_rows_selected)
        gene <- unique(cpdsum$genesymbol)
        if(is_empty(gene) | gene=="") return(NULL)
        geneidx <- which(fData(data$crispr_es)$genesymbol == gene)
        print(gene); print(geneidx)
        df <- data.frame(crispr_effect = exprs(data$crispr_es)[geneidx, ],
                         disease = pData(data$crispr_es)$disease,
                         subtype = pData(data$crispr_es)$subtype) %>%
            mutate(disease = factor(disease, levels = levels(data$crispr_es$disease)),
                   subtype = factor(subtype, levels = levels(data$crispr_es$subtype)))
        if (is_empty(df)) return(NULL)
        return(df)
    })

    # CRISPR density plot by disease
    output$crispr_dens <- renderPlot({
        plotdat <- get_cpd_crispr_dat()
        validate(need(!is_empty(plotdat), 
                      "      No CRISPR data associated with selected compound"))
        plotdat <- plotdat %>%
            filter(disease %in% c(reactvals$ct1, reactvals$ct2)) %>%
            droplevels()
        print(head(plotdat))
        gen_densplot(plotdat, "crispr_effect", "disease", xlab = "CRISPR effect score")
    })

    # CRISPR boxplot by disease
    output$crispr_box_ds <- renderPlot({
        plotdat <- get_cpd_crispr_dat() 
        if(is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            filter(!is.na(disease)) %>%
            droplevels()
        gen_boxplot(plotdat, "disease", "crispr_effect", ylab = "CRISPR effect score")
    })

    # CRISPR boxplot by subtype
    output$crispr_box_st <- renderPlot({
        plotdat <- get_cpd_crispr_dat() 
        if(is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            filter(disease %in% c(reactvals$ct1, reactvals$ct2)) %>%
            mutate(st = ifelse(disease == reactvals$ct1,
                               as.character(subtype), as.character(disease))) %>%
            mutate(st = factor(st, levels = c(levels(data$crispr_es$subtype), "Solid"))) %>%
            filter(!is.na(st)) %>%
            droplevels()
        gen_boxplot(plotdat, "st", "crispr_effect", ylab = "CRISPR effect score")
    })
    
    # data table of CRISPR data
    output$cl_crisprdat <- DT::renderDataTable({
        dat <- get_cpd_crispr_dat() %>%
            dplyr::filter(disease %in% reactvals$ct1) 
        if (is_empty(dat)) return(NULL)
        if (!is_empty(input$crispr_dens_brush)) {
            brushinfo <- input$crispr_dens_brush
            dat <- dat %>% filter(crispr_effect > brushinfo$xmin &
                                  crispr_effect < brushinfo$xmax)
        }
        DT::datatable(
            data = dat, 
            rownames = F,
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
            dat <- get_cpd_crispr_dat() %>%
                dplyr::filter(disease %in% c(ct, "Solid")) 
            writexl::write_xlsx(dat, path=file)
        })

    # get RNAi data
    get_rnai_dat <- reactive({
        cpdsum <- get_sel_drugsum()
        if(is_empty(cpdsum)) return(NULL)
        cpdsum <- cpdsum %>%
            slice(input$drug_summary_dt_rows_selected)
        gene <- unique(cpdsum$genesymbol)
        if(is_empty(gene) | gene=="") return(NULL)
        geneidx <- which(fData(data$rnai_es)$genesymbol == gene)
        if(is_empty(geneidx)) return(NULL)
        # validate(need(!is_empty(geneidx), "Selected compound has no annotated target"))
        data.frame(rnai_effect = exprs(data$rnai_es)[geneidx, ],
                   disease = pData(data$rnai_es)$disease,
                   subtype = pData(data$rnai_es)$subtype) %>%
            mutate(disease = factor(disease, levels = levels(data$rnai_es$disease)),
                   subtype = factor(subtype, levels = levels(data$rnai_es$subtype)))
    })

    # RNAi density plot by disease
    output$rnai_dens <- renderPlot({
        plotdat <- get_rnai_dat()
        validate(need(!is_empty(plotdat), 
                      "      No RNAi data associated with selected compound"))
        plotdat <- plotdat %>%
            filter(disease %in% c(reactvals$ct1, reactvals$ct2)) %>%
            droplevels()
        print(head(plotdat))
        gen_densplot(plotdat, "rnai_effect", "disease", xlab = "RNAi effect score")
    })

    # RNAi boxplot by disease
    output$rnai_box_ds <- renderPlot({
        plotdat <- get_rnai_dat() 
        if(is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            filter(!is.na(disease)) %>%
            droplevels()
        gen_boxplot(plotdat, "disease", "rnai_effect", ylab = "RNAi effect score")
    })

    # RNAi boxplot by subtype
    output$rnai_box_st <- renderPlot({
        plotdat <- get_rnai_dat() 
        if(is_empty(plotdat)) return(NULL)
        plotdat <- plotdat %>%
            filter(disease %in% c(reactvals$ct1, reactvals$ct2)) %>%
            mutate(st = ifelse(disease == reactvals$ct1,
                               as.character(subtype), as.character(disease))) %>%
            mutate(st = factor(st, levels = c(levels(data$rnai_es$subtype), "Solid"))) %>%
            filter(!is.na(st)) %>%
            droplevels()
        gen_boxplot(plotdat, "st", "rnai_effect", ylab = "RNAi effect score")
    })
    
    # cell-line level RNAi data
    output$cl_rnaidat <- DT::renderDataTable({
        dat <- get_rnai_dat() %>%
            dplyr::filter(disease %in% c(reactvals$ct1)) 
        if (is.null(dat)) return(NULL)
        DT::datatable(
            data = dat, 
            rownames = F,
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
            dat <- get_rnai_dat() %>%
                dplyr::filter(disease %in% c(ct, "Solid")) 
            writexl::write_xlsx(dat, path=file)
        })
    
    ############################## CRISPR ####################################
    
    # summarise crispr data
    get_crispr_dat <- reactive({
        ct1 <- reactvals$ct1; ct2 <- reactvals$ct2
        crispr <- data$crispr_es[, which(data$crispr_es$disease %in% c(ct1,ct2))] 
        var <- apply(exprs(crispr), 1, var)
        crispr <- crispr[var > quantile(var, .5, na.rm=T), ] 
        crisprdat <- exprs(crispr) %>%
            as.data.frame() %>%
            mutate(gene = fData(crispr)$genesymbol, .before=1) %>%
            gather("cl_id", "value", -1) %>%
            mutate(disease = pData(crispr)[match(.$cl_id, crispr$DepMap_ID),]$disease)
        return(crisprdat)
    })
    get_crispr_summary <- reactive({
        crisprdat <- get_crispr_dat()
        av <- crisprdat %>%
            group_by(disease, gene) %>%
            summarise(avCRISPR = mean(value, na.rm=T)) %>%
            pivot_wider(names_from = disease,
                        values_from = avCRISPR)
        stats <- crisprdat %>%
            filter(!is.na(value)) %>%
            select(gene, disease, value) %>%
            mutate(disease = factor(disease, levels=c(ct1,ct2))) %>%
            nest(-c(gene)) %>%
            mutate(dCRISPR = map_dbl(data, ~filt_ttest(.x, metric="value", stat="statistic")),
                   pCRISPR = map_dbl(data, ~filt_ttest(.x, metric="value", stat="p.value"))) %>%
            mutate(dCRISPR_rank = rank(-dCRISPR)/length(dCRISPR),
                   .before=3) %>%
            select(-data)
        summary <- merge(stats, av, by = c("gene"))
        summary %>% arrange(-dCRISPR_rank)
    })
    
    # scatter plot of average sensitivty scores per compound
    output$crispr_scatter <- renderPlot({
        plotdat <- get_crispr_summary()
        print(head(plotdat))
        print(tail(plotdat))
        if (is_empty(plotdat)) return(NULL)
        tmp <- c(paste0("`", reactvals$ct1, "`"),
                 paste0("`", reactvals$ct2, "`"))
        plotdat %>%
            ggplot(aes_string(x=tmp[[1]], y=tmp[[2]])) +
            geom_point(color = "#71448199", size = 2) +
            geom_abline() +
            scale_x_continuous(name = paste0("avCRISPR score: ", reactvals$ct1)) +
            scale_y_continuous(name = paste0("avCRISPR score: ", reactvals$ct2)) +
            theme_bw(base_size = 18) +
            theme()
    })
    
    # get selected points from scatter
    get_sel_crispr <- reactive({
        if (is_empty(input$crispr_scatter_brush)) {
            get_crispr_summary() %>%
                mutate_if(is.numeric, round, digits = 3)
        } else {
            get_crispr_summary() %>%
                brushedPoints(input$crispr_scatter_brush) %>%
                mutate_if(is.numeric, round, digits = 3)
        }
    }) 
    
    # drug info table
    output$crispr_summary_dt <- DT::renderDataTable({
        crispr_summary_sel <- get_sel_crispr() 
        if (is.null(crispr_summary_sel)) return(NULL)
        DT::datatable(
            data = crispr_summary_sel,
            rownames = F,
            selection = list(mode = 'single', target = "row", selected = 1)
        )
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
