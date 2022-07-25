library(shiny)
library(tidyverse)
library(fgsea)
library(shinythemes)
library(shinycssloaders)
library(shinyWidgets)
library(DT)
library(DBI)
library(RMariaDB)

################################### Setup #####################################

# function to generate pptx
library(rvg)
library(officer)
gen_pptx <- function(plot, file, height = 5, width = 5, left = 1, top = 1) {
   read_pptx() %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = dml(ggobj = plot), 
              location = ph_location(height = height, width = width,
                                     left = left, top = top),
              bg = "transparent") %>%
      print(target = file)
}

# db info
dbname <- "expr"
cnf <- list.files("/srv/shiny-server/app_expr/expr_data", pattern = paste0(dbname, ".cnf$"), full.names = T)
print(cnf)

# get sample info
db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
query <- "SELECT * FROM sample_info;"
queryres <- dbSendQuery(db, query)
sample_info <- dbFetch(queryres)
dbClearResult(queryres)
dbDisconnect(db)

# summarise sample info
si_summary <- sample_info[,c(1,3:8,10:11)] %>%
   group_by(tissue_type, disease_state, cell_type, cell_subtype, sample_type, group_name) %>%
   summarize(n = n(),
             studies = paste0(unique(na.omit(study_accession)), collapse = ";"),
             references = paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed/",
                                 unique(na.omit(pmid)),
                                 "'>",unique(na.omit(pmid)),"</a>",
                                 collapse = ";"),
             platforms = paste0(unique(na.omit(platform)), collapse = ";")) %>%
   as.data.frame()

# get feature info 
db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
query <- paste0("SELECT * FROM feature_info;")
queryres <- dbSendQuery(db, query)
feature_info <- dbFetch(queryres)
dbClearResult(queryres)
dbDisconnect(db)

# set reactive values
reactvals <- reactiveValues(selgenes = c("PAX5"),
                            si_summary = si_summary)

##################################### UI ######################################

ui <- fluidPage(
   
   # title
   title = "Expression Explorer", 
   theme = shinytheme("cosmo"),
   titlePanel(tags$h4(tags$a(
      imageOutput("icon", inline = TRUE),
      href="http://10.197.211.94:80"),
      "Expression Explorer")),
   
   # change color for column selection
   tags$style(HTML('table.dataTable tr.selected td{background-color: #88ddee77 !important;}')),
   tags$style(HTML('table.dataTable td.selected {background-color: #bbbbbb !important;}')),
   
   # main
   tabsetPanel(
      tabPanel("Expression Level", align = "center",
               checkboxGroupButtons("plot_type", "Plot Element:", 
                                    choices = list("Violin" = "violin",
                                                   "Boxplot" = "boxplot",
                                                   "Dotplot" = "dotplot"),
                                    justified = T, 
                                    size = 'xs', width = '25%',
                                    selected = c("violin","boxplot","dotplot")),
               uiOutput("expr_plot_ui") %>% 
                  withSpinner(color="#88ddee88"),
               div(style="display:inline-block;",
                   downloadButton("dl_expr_plot_ppt", label = "PPT",
                                  style = "font-size:12px;height:30px;padding:5px;"),
                   downloadButton("dl_expr_plot_png", label = "PNG",
                                  style = "font-size:12px;height:30px;padding:5px;"),
                   downloadButton("dl_expr_dat_xls", label = "Data XLSX",
                                  style = "font-size:12px;height:30px;padding:5px;")),
               br(), br(),
               tags$h3("Expression summary:"),
               DT::dataTableOutput("expr_summary_dt"),
               #div(style="display:inline-block;",
               #    actionButton("expand_levels", "Expand",
               #             style = "font-size:12px;height:30px;padding:5px;")),
	       #div(style = "font-size:12px", uiOutput("group_ranks")),
	       div(style="display:inline-block;",
                   actionButton("reorder_levels", "Reorder",
                                style = "font-size:12px;height:30px;padding:5px;")),
               br(), br(),
               tags$h3("Sample info:"),
               DT::dataTableOutput("samples_dt"),
               div(style="display:inline-block;",
                  checkboxInput("samples_dt_sel", "Select/De-select all")),
               div(style="display:inline-block;",
                  h5(textOutput("n_samples_sel"))),
               br(), br(),
               tags$h3("Feature info:"),
               DT::dataTableOutput("features_dt"))
   )
)

################################### Server ####################################
server <- function(input, output) {
   
   # icon
   output$icon <- renderImage(list(src = "../hexagons/expr.png",
                                   height = "90px", width = "85px"), 
                              deleteFile = F)
   
   # sample info table
   output$samples_dt <- DT::renderDataTable({
      si_summary <- reactvals$si_summary
      DT::datatable(si_summary, 
                    escape = F, 
                    rownames = F,
                    colnames = Hmisc::capitalize(gsub("_", " ", colnames(si_summary))),
                    selection = list(target = "row+column", 
                                     selected = list(rows = c(1:nrow(si_summary)), columns = 6)),
                    options = list(
                       dom = 'frtlip',
                       columnDefs = list(list(
                        targets = 7,
                        render = JS(
                          "function(data, type, row, meta) {",
                          "return type === 'display' && data.length > 15 ?",
                          "'<span title=\"' + data + '\">' + data.substr(0, 15) + '...</span>' : data;",
                          "}")))
                    ),
                    callback = JS('table.page(3).draw(false);'),
                    caption = "Select samples.")
   })
   samples_dt_proxy <- DT::dataTableProxy("samples_dt")
   observeEvent(input$samples_dt_sel, {
      if (isTRUE(input$samples_dt_sel)) {
         DT::selectRows(samples_dt_proxy, input$samples_dt_rows_all)
      } else {
         DT::selectRows(samples_dt_proxy, NULL)
      }
   })
   observeEvent(input$samples_dt_rows_selected, {
      sel_si <- get_sel_si()
      output$n_samples_sel <- renderText(paste0("Total selected samples: ",
                                                nrow(sel_si)))
   })
   observeEvent(input$samples_dt_columns_selected, {
      si_summary <- reactvals$si_summary
      idx <- input$samples_dt_columns_selected+1
      reactvals$grouping_col <- colnames(si_summary)[idx]
   })
   
   # return selected sample info
   get_sel_si <- reactive({
      si_summary <- reactvals$si_summary
      si_summary <- si_summary[input$samples_dt_rows_selected,]
      reactvals$group_levels <- unique(sample_info$group_name)
      apply(si_summary, 1, function(x) {
         sample_info[which(sample_info$tissue_type == x["tissue_type"] &
                              sample_info$disease_state == x["disease_state"] &
                              sample_info$cell_type == x["cell_type"] & 
                              sample_info$cell_subtype == x["cell_subtype"] &
                              sample_info$sample_type == x["sample_type"] & 
                              sample_info$group_name == x["group_name"]),]
      }) %>% bind_rows()
   })
   
   # summarise feature info
   get_fi_summary <- reactive({
      feature_info %>%
         group_by(feature_id) %>%
         summarize(platforms = paste0(unique(na.omit(platform)), collapse = ";"),
                   gene_symbols = paste0(unique(na.omit(gene_symbol)), collapse = ";"),
                   gene_aliases = paste0(unique(na.omit(gene_alias)), collapse = ";"),
                   entrez_IDs = paste0(unique(na.omit(entrez_id)), collapse = ";"),
                   ensembl_IDs = paste0(unique(na.omit(ensembl_id)), collapse = ";")) %>%
         as.data.frame()
   })
   
   # feature info table
   output$features_dt <- DT::renderDataTable({
      fi_summary <- get_fi_summary()
      DT::datatable(fi_summary,
                    escape = F, 
                    rownames = F,
                    selection = list(target = "row", selected = 1),
                    colnames = Hmisc::capitalize(gsub("_", " ", colnames(fi_summary))),
                    options = list(
                       pagelength = 5,
                       lengthMenu = list(c(5, 10, 20),
                                         c('5','10','20'))),
                    caption = "Select features.")
   })
   
   # return selected features
   get_sel_fi <- reactive({
      fi_summary <- get_fi_summary()
      fi_summary[input$features_dt_rows_selected,]
   })
   
   # get expr info
   get_expr_data <- reactive({
      sel_si <- get_sel_si()
      sel_fi <- get_sel_fi()
      if (nrow(sel_si) < 1 | nrow(sel_fi) < 1) return(NULL)
      db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
      query <- paste0('SELECT * FROM expr_data ',
                      "INNER JOIN sample_info ",
                      "ON expr_data.sample_accession = sample_info.sample_accession ",
                      'WHERE expr_data.feature_id IN ("',
                      paste0(sel_fi$feature_id, collapse = '", "'),
                      '") AND expr_data.sample_accession IN ("',
                      paste0(sel_si$sample_accession, collapse = '", "'),
                      '");')
      queryres <- dbSendQuery(db, query)
      exprdat <- dbFetch(queryres)
      dbClearResult(queryres)
      dbDisconnect(db)
      exprdat$gene_symbol <- sel_fi[match(exprdat$feature_id, sel_fi$feature_id),]$gene_symbol
      if (is_empty(reactvals$grouping_col)) {
         exprdat$group_label <- exprdat$group_name
      } else{
         exprdat$group_label <- exprdat[,reactvals$grouping_col]
      }
      reactvals$exprdat <- exprdat[, c(16:17,1:3,6:11)]
      return(NULL)
   })
   
   # generate expr summary
   get_expr_summary <- reactive({
      get_expr_data()
      exprdat <- reactvals$exprdat
      if (is_empty(exprdat)) return(NULL)
      exprdat %>%
         group_by(gene_symbol, group_label) %>%
         summarise(n = n(),
                   mean = round(mean(expr_val),2),
                   SD = round(sd(expr_val),2),
                   cell_type = paste0(unique(na.omit(cell_type)), collapse = ","),
                   subtype = paste0(unique(na.omit(cell_subtype)), collapse = ","),
                   tissue_type = paste0(unique(na.omit(tissue_type)), collapse = ","),
                   disease_state = paste0(unique(na.omit(disease_state)), collapse = ","),
                   sample_type = paste0(unique(na.omit(sample_type)), collapse = ","),
                   group_name = paste0(unique(na.omit(group_name)), collapse = ","))
   })
   
   # output expr summary dt
   output$expr_summary_dt <- DT::renderDataTable({
      exprsum <- get_expr_summary()
      if (is_empty(exprsum)) return(NULL)
      DT::datatable(exprsum, 
                    rownames = F,
                    colnames = Hmisc::capitalize(gsub("_", " ", colnames(exprsum))),
                    editable = 'column',
                    selection = list(target = 'row', mode = 'multiple'),
                    options = list(
                       dom = 'frtlip',
                       columnDefs = list(list(
                       targets = c(7,8),
                       render = JS(
                          "function(data, type, row, meta) {",
                          "return type === 'display' && data.length > 15 ?",
                          "'<span title=\"' + data + '\">' + data.substr(0, 15) + '...</span>' : data;",
                          "}")))
                    ),
                    callback = JS('table.page(3).draw(false);'),
                    caption = "Select columns or edit cells to change grouping and labels")
   })
  
   # group ranking
   output$group_ranks <- renderUI({
     if ((input$rank_toggle %% 2) != 0) {
       rank_list(text = "Drag groups to order",
                input_id = "group_rankings",
                labels = unique(reactvals$group_levels))
     } else return(NULL)
   })
   observeEvent(input$group_rankings, {
     reactvals$group_levels <- as.character(unique(input$group_rankings))
   })

   # update grouping column and names
   # expr_summary_dt_proxy <- DT::dataTableProxy("expr_summary_dt")
   # observeEvent(input$expr_summary_dt_cell_edit, {
   #    exprdat <- get_expr_data()
   #    original <- get_expr_summary()
   #    updated <- editData(original, input$expr_summary_dt_cell_edit,
   #                        'expr_summary_dt')
   #    idx <- match(exprdat$group_label, original$group_label)
   #    exprdat$group_label <- updated$group_label[idx]
   # })
   observeEvent(input$expr_summary_dt_rows_selected, {
      exprsum <- get_expr_summary()[input$expr_summary_dt_rows_selected,]
      reactvals$sel_groups <- unique(exprsum$group_label)
   })

   # generate distro plot
   plot_distribution <- reactive({
      exprdat <- reactvals$exprdat
      if (is_empty(exprdat)) return(NULL)
      sel_fi <- get_sel_fi()
      sel_groups <- reactvals$sel_groups
      plot_types <- input$plot_type
      if (is_empty(sel_groups)) {
         exprdat$col <- exprdat$group_label
      } else {
         exprdat$col <- ifelse(exprdat$group_label %in% sel_groups, "B","A")
      }
      p <- ggplot(exprdat, aes(x = group_label, y = expr_val)) +
         theme_bw(base_size = 18) + 
         ylab(paste0(unique(sel_fi$gene_symbol)," Expression [AU]")) +
         theme(panel.grid = element_blank(),
               legend.position = "none",
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.y = element_text(color = "black"),
               axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
      if ("dotplot" %in% plot_types) {
         p <- p + geom_point(aes(color = col),
                             alpha = .5,
                             position = position_dodge2(width = .5))
      }
      if ("violin" %in% plot_types) {
         p <- p + geom_violin(aes(fill = col))
      }
      if ("boxplot" %in% plot_types) {
         if ("dotplot" %in% plot_types & "violin" %in% plot_types) {
            p <- p + geom_boxplot(aes(fill = col), width = 0.2, outlier.shape = NA)
         } else if ("dotplot" %in% plot_types) {
            p <- p + geom_boxplot(aes(fill = col), width = 0.8, outlier.shape = NA)
         } else if ("violin" %in% plot_types) {
            p <- p + geom_boxplot(aes(fill = col), width = 0.2)
         } else {
            p <- p + geom_boxplot(aes(fill = col), width = 0.8)
         }
      }
      return(p)
   })
   
   # generate bivariate plot
   plot_bivariate <- reactive({
      plotdat <- reactvals$exprdat
      if (is_empty(plotdat)) return(NULL)
   })
   
   # output distro plot
   sel_expr_plot <- reactive({
      sel_fi <- get_sel_fi()
      if (is_empty(sel_fi)) return(NULL)
      if (length(unique(sel_fi$gene_symbol)) > 1) {
         plot_bivariate()
      } else {
         plot_distribution()
      }
   })
   
   # output expr plot
   output$expr_plot <- renderPlot({
      sel_expr_plot()
   })
   get_dims <- reactive({
      exprdat <- reactvals$exprdat
      ncond <- length(unique(exprdat$group_label))
      width <- min(70*ncond, 1000)
      list(w = width, h = 350)
   })
   output$expr_plot_ui <- renderUI({
      dims <- get_dims()
      plotOutput("expr_plot", height = dims$h, width = dims$w)
   })
   
   # download data as xlsx
   output$dl_expr_dat_xls <- downloadHandler(
      filename = function() {
         exprsum <- get_expr_summary()
         genes <- unique(exprsum$gene_name)
         paste0("Expression_level_", paste(genes, collapse="_"), ".xlsx")
      },
      content = function(file) {
         writexl::write_xlsx(reactvals$exprdat, path=file)
   })
   
   # download expr plot as png
   output$dl_expr_plot_png <- downloadHandler(
      filename = function() {
         exprsum <- get_expr_summary()
         if (length(unique(exprsum$gene_name)) > 2) {
            genes <- "multiple_genes"
         } else {
            genes <- paste0(unique(exprsum$gene_name), collapse = "_")
         }
         paste0("Expression_level_", genes, ".png")
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
}

# Run the application 
shinyApp(ui = ui, server = server)
