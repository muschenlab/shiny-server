library(shiny)
library(tidyverse)
library(DT)
library(rvg)
library(officer)
library(shinythemes)

## tmp
# di <- data.frame(id = c("JL_IFITM3_KO_BCRABL_RNAseq",
#                         "JL_IFITM3_KO_NRAS_RNAseq"),
#                  description = c("Differential expression analysis of Ifitm3-KO vs Ifitm3-WT cells in BCR-ABL1 transformed pre-B cells",
#                                  "Differential expression analysis of Ifitm3-KO vs Ifitm3-WT cells in NRAS-G12D transformed pre-B cells"),
#                  organism = c("mouse",
#                               "mouse"),
#                  celltype = c("pre-B",
#                               "pre-B"),
#                  assay = c("RNAseq",
#                            "RNAseq"),
#                  performed_by = c("Jaewoong Lee",
#                                   "Jaewoong Lee"),
#                  year = c("2020",
#                           "2020"))
# saveRDS(di, "functional_data/dataset_info.rds")

# load data
di <- readRDS("functional_data/dataset_info.rds")

# set stored values
reactvals <- reactiveValues(dat = NULL)

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

################### ui ########################
ui <- fluidPage(
   
  # title
  title = "Functional Studies", 
  theme = shinytheme("cosmo"),
  titlePanel(tags$h2(tags$a(
    imageOutput("icon", inline = TRUE),
    href="http://137.184.200.69:3838/"), "Functional Studies")),
  
  # dataset info table
  tags$h3("Datasets:"),
  tags$h5("Select a row from the table below to further explore the dataset"),
  div(DT::dataTableOutput("di_dt"),
      style = "font-size:90%"),
  
  # plots
  br(),
  tags$h3("Differential plots:"),
  column(width = 6, align = 'center',
         plotOutput('volcanoplot', 
                    click = clickOpts(id = "volcano_click"),
                    hover = hoverOpts(id = "volcano_hover"),
                    brush = brushOpts(id = "volcano_brush")),
         div(style="display:inline-block; padding-left:50px; ",
             downloadButton("dl_volcano_ppt", label = "PPT", style = "padding:4px; font-size:70%"),
             downloadButton("dl_volcano_png", label = "PNG", style = "padding:4px; font-size:70%"))
  ),
  column(width = 6, align = 'center',
         plotOutput('maplot', 
                    click = clickOpts(id = "ma_click"),
                    hover = hoverOpts(id = "ma_hover"),
                    brush = brushOpts(id = "ma_brush")),
         div(style="display:inline-block; padding-left:50px; ",
             downloadButton("dl_ma_ppt", label = "PPT", style = "padding:4px; font-size:70%"),
             downloadButton("dl_ma_png", label = "PNG", style = "padding:4px; font-size:70%"))
  )#,
  # column(width = 12, align="center",
  #        br(), 
  #        actionButton("reset_bttn", "Reset Selection", 
  #                     style = 'padding:4px; color: #fff; background-color: #337ab7; border-color: #2e6da4'),
  #        br(), br(), 
  #        div(style="display:inline-block", 
  #            textInput('id_input', 'Select protein by name:', 
  #                      value = "", 
  #                      placeholder = "Enter list of gene or protein names")),
  #        div(style="display:inline-block", 
  #            actionButton('input_go_bttn', 'Go!', style = 'font-size:90%')),
  #        textOutput("missing_prot"),
  #        br(), br(), 
  #        h4(strong("Selected Proteins")),
  #        DT::DTOutput('protein_info_dt')
  #  )
)

################### server ########################
server <- function(input, output, session) {
   
   # icon
   output$icon <- renderImage(list(src = "../assets/images/hex/functional_studies.png",
                                   height = "50px", width = "50px"), 
                              deleteFile = F)
   
   # dataset info table
   output$di_dt <- DT::renderDataTable({
     print("di_dt")
     DT::datatable(
       data = di,
       rownames = F,
       colnames = gsub("_", " ", colnames(di)),
       selection = list(mode = 'single', target = "row", selected = 1),
       options = list(columnDefs = list(list(
         targets = 0:(ncol(di)-1),
         render = JS(
           "function(data, type, row, meta) {",
           "return type === 'display' && data.length > 30 ?",
           "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
           "}")
       ))), callback = JS('table.page(3).draw(false);')
     )
   })
   
   # handle table selection
   get_dataset <- reactive({
     print("get_dataset")
     seldi <- di %>%
       dplyr::slice(input$di_dt_rows_selected)
     if(nrow(seldi) < 1) return(NULL)
     print(seldi)
     print(paste0("functional_data/",seldi$id,".rds"))
     dat <- readRDS(paste0("functional_data/",seldi$id,".rds"))
     if(is.null(dat)) return(NULL)
     reactvals$dat <- dat
     return(dat)
   })
   
   # plot scatter
   gen_volcano_plot <- reactive({
     dat <- get_dataset()
     if(is.null(dat)) { return(NULL) }
     top <- dat %>% mutate(rank_metric=-log10(padj+1e-100)*(L2FC/10)) %>% top_n(10, wt=rank_metric)
     bttm <- dat %>% mutate(rank_metric=-log10(padj+1e-100)*(-L2FC/10)) %>% top_n(10, wt=rank_metric)
     dat %>% 
       mutate(label = ifelse(gene_symbol %in% c(top$gene_symbol, bttm$gene_symbol),
                             gene_symbol, NA)) %>%
       mutate(color = ifelse(padj > 0.05, "unaltered", 
                             ifelse(L2FC > 0, "upregulated", "downregulated"))) %>%
       mutate(color = factor(color, levels=c("upregulated","downregulated","unaltered"))) %>%
       ggplot(aes(x = L2FC, y = -log10(padj+1e-200), label = label)) + 
       geom_point(aes(color=color), alpha = .5, size = 1) +
       ggrepel::geom_text_repel(aes(label = label), size = 5, arrow = arrow(angle=0)) +
       scale_color_manual(values = c("firebrick4","steelblue4","grey60")) +
       xlab("Fold Change [log2]") +
       ylab("Significance [-log10 adjusted p-value]") +
       theme_bw(base_size = 16) +
       theme(panel.grid = element_blank(),
             legend.background = element_blank(),
             legend.title = element_blank(), 
             legend.key = element_blank(),
             legend.position = "none") 
   })
   output$volcanoplot <- renderPlot({
      gen_volcano_plot()
   })
   
   # plot MA
   gen_ma_plot <- reactive({
     dat <- get_dataset()
     if(is.null(dat)) { return(NULL) }
     top <- dat %>% mutate(rank_metric=-log10(padj+1e-100)*(L2FC/10)) %>% top_n(10, wt=rank_metric)
     bttm <- dat %>% mutate(rank_metric=-log10(padj+1e-100)*(-L2FC/10)) %>% top_n(10, wt=rank_metric)
     dat %>% 
       mutate(label = ifelse(gene_symbol %in% c(top$gene_symbol, bttm$gene_symbol),
                             gene_symbol, NA)) %>%
       mutate(color = ifelse(padj > 0.05, "unaltered", 
                             ifelse(L2FC > 0, "upregulated", "downregulated"))) %>%
       mutate(color = factor(color, levels=c("upregulated","downregulated","unaltered"))) %>%
       ggplot(aes(x = av_expr, y = L2FC, label = label)) +
       geom_point(aes(color=color), alpha = .5, size = 1) +
       ggrepel::geom_text_repel(aes(label = label), size = 5, arrow = arrow(angle=0)) +
       scale_color_manual(values = c("firebrick4","steelblue4","grey60")) +
       xlab("Average expression [log normalized]") +
       ylab("Fold change [log2]") +
       theme_bw(base_size = 16) +
       theme(panel.grid = element_blank(),
             legend.background = element_blank(),
             legend.title = element_blank(), 
             legend.key = element_blank(),
             legend.position = "none") 
   })
   output$maplot <- renderPlot({
      gen_ma_plot()
   })
   
   # # update prot list from clicks/brushes
   # observeEvent(input$scatter_brush, {
   #    brushed <- row.names(brushedPoints(dat, input$scatter_brush))
   #    values$selprot <- unique(c(brushed, values$selprot))
   # })
   # observeEvent(input$scatter_click, {
   #    clicked <- row.names(nearPoints(dat, input$scatter_click))
   #    values$selprot <- unique(c(clicked, values$selprot))
   # })
   # observeEvent(input$volcano_brush, {
   #    brushed <- row.names(brushedPoints(dat, input$volcano_brush))
   #    values$selprot <- unique(c(brushed, values$selprot))
   # })
   # observeEvent(input$volcano_click, {
   #    clicked <- row.names(nearPoints(dat, input$volcano_click))
   #    values$selprot <- unique(c(clicked, values$selprot))
   # })
   # observeEvent(input$ma_brush, {
   #    brushed <- row.names(brushedPoints(dat, input$ma_brush))
   #    values$selprot <- unique(c(brushed, values$selprot))
   # })
   # observeEvent(input$ma_click, {
   #    clicked <- row.names(nearPoints(dat, input$ma_click))
   #    values$selprot <- unique(c(clicked, values$selprot))
   # })

   # # update prot list from input text
   # observeEvent(input$input_go_bttn, {
   #    id_list <- strsplit(input$id_input, split = ",| |\n")[[1]]
   #    id_list <- id_list[lapply(id_list, nchar) > 0]
   #    prot <- id_list[id_list %in% dat$protein_name]
   #    unipro <- id_list[id_list %in% dat$uniprot_id]
   #    if (length(unipro) > 0) {
   #       prot <- c(prot, dat[dat$uniprot_id %in% unipro,]$protein_name)
   #    }
   #    genelist <- sapply(dat$gene_symbols, function(x) {
   #       tmp <- strsplit(x, split = " |;")[[1]]
   #       tmp[!tmp==""] })
   #    genes <- id_list[id_list %in% unlist(genelist)]
   #    if (length(genes) > 0) {
   #       idx <- unlist(lapply(genelist, function(x) any(x %in% genes)))
   #       prot <- c(prot, dat[idx,]$protein_name)
   #    }
   #    missing <- id_list[!id_list %in% c(dat$protein_name, dat$uniprot_id, unlist(genelist))]
   #    if (length(missing) > 0) {
   #       print(missing)
   #       output$missing_prot <- renderText(paste0("Following names could not be found: ",
   #                                                paste(missing, collapse = ", ")))
   #    } else {
   #       output$missing_prot <- NULL
   #    }
   #    if (length(prot) > 0) {
   #      values$selprot <- unique(c(values$selprot, rownames(dat[dat$protein_name %in% prot,])))
   #    }
   # })
   
   # reset protein selection
   observeEvent(input$reset_bttn, {
      values$selprot <- vector()
   })
   
   # # datatable of info on selected proteins
   # output$protein_info_dt <- DT::renderDT({
   #    sel_prot <- values$selprot
   #    if (is_empty(sel_prot)) {
   #       seldat <- dat[, c(1:2,15,17:20,3:14)]
   #    } else {
   #       seldat <- dat[sel_prot, c(1:2,15,17:20,3:14)]
   #    }
   #    colnames(seldat) <- c("Protein Name", "UniProt ID", "Gene Symbols", 
   #                          "Specificity","TM Domain", "# TM Domain", "Membrane Type",
   #                          "L2FC", "p-value", "q-value",
   #                          "Av Expr", "Av Normal", "Av leukemia", 
   #                          "normal r1", "normal r2", "normal r3",
   #                          "leukemia r1", "leukemia r2", "leukemia r3")
   #    seldat[,sapply(seldat, is.numeric)] <- apply(seldat[,sapply(seldat, is.numeric)],
   #                                                 2, function(x) signif(x, 2))
   #    DT::datatable(seldat, rownames = F,
   #                  extensions = 'Buttons',
   #                  options = list(dom = 'Bfrtip',
   #                                 buttons = list(list(extend = 'copy',
   #                                                     title = NULL),
   #                                                list(extend = 'print',
   #                                                     title = "Surface Proteomics - Selected Proteins"),
   #                                                list(extend = 'csv',
   #                                                     filename = paste0("surfprot_selected_proteins_info.", Sys.Date()),
   #                                                     title = NULL),
   #                                                list(extend = 'excel',
   #                                                     filename = paste0("surfprot_selected_proteins_info.", Sys.Date()),
   #                                                     title = NULL),
   #                                                list(extend = 'pdf',
   #                                                     filename = paste0("surfprot_selected_proteins_info.", Sys.Date()),
   #                                                     title = "Surface Proteomics - Selected Protein")),
   #                                 columnDefs = list(list(
   #                                    targets = c(3,6),
   #                                    render = JS(
   #                                       "function(data, type, row, meta) {",
   #                                       "return type === 'display' && data.length > 10 ?",
   #                                       "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
   #                                       "}")))),
   #                  callback = JS('table.page(3).draw(false);'))
   # })
      
   # download scatter plot
   output$dl_scatter_png <- downloadHandler(
      filename = function() {
         paste0("Surf_proteomics_BCR-ABL1_scatterplot.", Sys.Date(), ".png")
      },
      content = function(file) {
         ggsave(plot = gen_scatter_plot(), filename = file, height = 5, width = 10)
      }
   )
   output$dl_scatter_ppt <- downloadHandler(
      filename = function() {
         paste0("Surf_proteomics_BCR-ABL1_scatterplot.", Sys.Date(), ".pptx")
      },
      content = function(file) {
         file_pptx <- tempfile(fileext = ".pptx")
         static_plot <- gen_scatter_plot() + theme(text = element_text(color = "white"),
                                                   axis.text = element_text(color = "white"))
         layer_length <- length(static_plot$layers)
         static_plot$layers[[layer_length]] <- NULL # to remove text
         vector_plot <- gen_scatter_plot() + theme(panel.background = element_blank(),
                                                   plot.background = element_blank())
         vector_plot$layers[[1]] <- NULL  # to remove points
         gen_pptx(static_plot, vector_plot,file_pptx)
         file.rename(from = file_pptx, to = file)
      }
   )
   
   # download volcano plot
   output$dl_volcano_png <- downloadHandler(
      filename = function() {
         paste0("Surf_proteomics_BCR-ABL1_volcanoplot.", Sys.Date(), ".png")
      },
      content = function(file) {
         ggsave(plot = gen_volcano_plot(), filename = file, height = 5, width = 10)
      }
   )
   output$dl_volcano_ppt <- downloadHandler(
      filename = function() {
         paste0("Surf_proteomics_BCR-ABL1_volcanoplot.", Sys.Date(), ".pptx")
      },
      content = function(file) {
         file_pptx <- tempfile(fileext = ".pptx")
         static_plot <- gen_volcano_plot() + theme(text = element_text(color = "white"),
                                                   axis.text = element_text(color = "white"))
         layer_length <- length(static_plot$layers)
         static_plot$layers[[layer_length]] <- NULL # to remove text
         vector_plot <- gen_volcano_plot() + theme(panel.background = element_blank(),
                                                   plot.background = element_blank())
         vector_plot$layers[[1]] <- NULL  # to remove points
         gen_pptx(static_plot, vector_plot,file_pptx)
         file.rename(from = file_pptx, to = file)
      }
   )
   
   # download MA plot
   output$dl_ma_png <- downloadHandler(
      filename = function() {
         paste0("Surf_proteomics_BCR-ABL1_MAplot.", Sys.Date(), ".png")
      },
      content = function(file) {
         ggsave(plot = gen_ma_plot(), filename = file, height = 5, width = 6.5)
      }
   )
   output$dl_ma_ppt <- downloadHandler(
      filename = function() {
         paste0("Surf_proteomics_BCR-ABL1_MAplot.", Sys.Date(), ".pptx")
      },
      content = function(file) {
         file_pptx <- tempfile(fileext = ".pptx")
         static_plot <- gen_ma_plot() + theme(text = element_text(color = "white"),
                                                   axis.text = element_text(color = "white"))
         layer_length <- length(static_plot$layers)
         static_plot$layers[[layer_length]] <- NULL # to remove text
         vector_plot <- gen_ma_plot() + theme(rect = element_rect(fill = "transparent",colour = NA),
                                              plot.background = element_rect(fill = "transparent",colour = NA))
         vector_plot$layers[[1]] <- NULL # to remove points
         gen_pptx(static_plot, vector_plot, file_pptx)
         file.rename(from = file_pptx, to = file)
      }
   )
}

# Run the application 
shinyApp(ui = ui, server = server)

