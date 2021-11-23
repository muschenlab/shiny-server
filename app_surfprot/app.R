library(shiny)
library(tidyverse)
library(DT)
library(rvg)
library(officer)
options(stringsAsFactors = F)

# test list: YLAT1 LAT3 S38A2  YLAT2 LAT1 S38A1 CTR1  4F2  Slc1a5 Slc1a4 Epha2 Epha5

# load data
dat <- readRDS("/poolio/internal_data/surf_proteomics/GD_surfprot_BCR_ABL1_4h.2019-11-08.rds")
dat <- as.data.frame(dat)

# set stored values
values <- reactiveValues(selprot = vector())

# function to generate pptx
gen_pptx <- function(static_plot, vector_plot, file) {
   read_pptx() %>% 
      add_slide(layout = "Title and Content", master = "Office Theme") %>% 
      ph_with_gg(value = static_plot, type = 'body', height = 5, width = 6.5) %>%
      ph_with_vg(ggobj = vector_plot, type = 'body', height = 5, width = 6.5) %>% 
      print(target = file)
}

################### ui ########################
ui <- fluidPage(
   
   title = "SurfProt", 
   theme = "bootstrap.flatly.css",
   
   titlePanel(tags$h3(tags$a(
      imageOutput("icon", height = "50px", width = "50px", inline = TRUE),
      href="http://10.197.211.94:80"), "Surface Proteomics")),
   
   column(width = 4, align = 'center',
          plotOutput('scatterplot', 
                     click = clickOpts(id = "scatter_click"),
                     hover = hoverOpts(id = "scatter_hover"),
                     brush = brushOpts(id = "scatter_brush")),
          div(style="display:inline-block; padding-left:50px; ",
              downloadButton("dl_scatter_ppt", label = "PPT", style = "padding:4px; font-size:70%"),
              downloadButton("dl_scatter_png", label = "PNG", style = "padding:4px; font-size:70%"))
   ),
   column(width = 4, align = 'center',
          plotOutput('volcanoplot', 
                     click = clickOpts(id = "volcano_click"),
                     hover = hoverOpts(id = "volcano_hover"),
                     brush = brushOpts(id = "volcano_brush")),
          div(style="display:inline-block; padding-left:50px; ",
              downloadButton("dl_volcano_ppt", label = "PPT", style = "padding:4px; font-size:70%"),
              downloadButton("dl_volcano_png", label = "PNG", style = "padding:4px; font-size:70%"))
   ),
   column(width = 4, align = 'center',
          plotOutput('maplot', 
                     click = clickOpts(id = "ma_click"),
                     hover = hoverOpts(id = "ma_hover"),
                     brush = brushOpts(id = "ma_brush")),
          div(style="display:inline-block; padding-left:50px; ",
              downloadButton("dl_ma_ppt", label = "PPT", style = "padding:4px; font-size:70%"),
              downloadButton("dl_ma_png", label = "PNG", style = "padding:4px; font-size:70%"))
   ),
   column(width = 12, align="center",
          br(), 
          actionButton("reset_bttn", "Reset Selection", 
                       style = 'padding:4px; color: #fff; background-color: #337ab7; border-color: #2e6da4'),
          br(), br(), 
          div(style="display:inline-block", 
              textInput('id_input', 'Select protein by name:', 
                        value = "", 
                        placeholder = "Enter list of gene or protein names")),
          div(style="display:inline-block", 
              actionButton('input_go_bttn', 'Go!', style = 'font-size:90%')),
          textOutput("missing_prot"),
          br(), br(), 
          h4(strong("Selected Proteins")),
          DT::DTOutput('protein_info_dt')
   )
)

################### server ########################
server <- function(input, output) {
   
   # icon
   output$icon <- renderImage(list(src = "www/moose.png",
                                   height = "50px", width = "50px"), 
                              deleteFile = F)
   
   # plot scatter
   gen_scatter_plot <- reactive({
      ggplot(dat, aes(x = av_norm, y = av_leuk, label = protein_name)) + 
         geom_point(aes(col = pval > -log10(0.05)), alpha = .5, size = 1) +
         geom_point(data = subset(dat, TM_domain_presence == "Y"),
                    col = "black", fill = NA, shape = 21, size = 1.5) +
         geom_point(data = subset(dat, rownames(dat) %in% values$selprot),
                    col = "purple", fill = NA, shape = 21, size = 2.5) +
         geom_text(data = subset(dat, rownames(dat) %in% values$selprot),
                   aes(label = protein_name, col = pval > -log10(0.05)), size = 5, 
                   position = position_nudge(y=.5)) +
         scale_color_manual(values = c("grey20","firebrick"),
                            labels = c("non-significant","p < 0.05", " ")) +
         xlab("Abundance - Normal [LFQ]") +
         ylab("Abundance - Leukemia [LFQ]") +
         scale_x_continuous(limits = c(min(dat$av_norm, na.rm = T),
                                       max(dat$av_norm, na.rm = T))) +
         scale_y_continuous(limits = c(min(dat$av_leuk, na.rm = T),
                                       max(dat$av_leuk, na.rm = T))) +
         theme_bw(base_size = 16) +
         theme(panel.grid = element_blank(),
               legend.background = element_blank(),
               legend.title = element_blank(), 
               legend.key = element_blank(),
               legend.position = "none") 
   })
   output$scatterplot <- renderPlot({
      gen_scatter_plot()
   })
   
   # plot scatter
   gen_volcano_plot <- reactive({
      ggplot(dat, aes(x = l2fc, y = pval, label = protein_name)) + 
         geom_point(aes(col = pval > -log10(0.05)), alpha = .5, size = 1) +
         geom_point(data = subset(dat, TM_domain_presence == "Y"),
                    col = "black", fill = NA, shape = 21, size = 1.5) +
         geom_point(data = subset(dat, rownames(dat) %in% values$selprot),
                    col = "purple", fill = NA, shape = 21, size = 2.5) +
         geom_text(data = subset(dat, rownames(dat) %in% values$selprot),
                   aes(label = protein_name, col = pval > -log10(0.05)), size = 5, 
                   position = position_nudge(y=.05)) +
         scale_color_manual(values = c("grey20","firebrick"),
                            labels = c("non-significant","p < 0.05", " ")) +
         scale_x_continuous(limits = c(-15,15)) +
         scale_y_continuous(limits = c(min(dat$pval, na.rm = T),
                                       max(dat$pval, na.rm = T))) +
         xlab("Fold Change - Leuk / Healthy [log2]") +
         ylab("Significance [-log10 p-value]") +
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
      ggplot(dat, aes(x = av_expr, y = l2fc, label = protein_name)) + 
         geom_point(aes(col = pval > -log10(0.05)), alpha = .5, size = 1) +
         geom_point(data = subset(dat, TM_domain_presence == "Y"),
                    col = "black", fill = NA, shape = 21, size = 1.5) +
         geom_point(data = subset(dat, rownames(dat) %in% values$selprot),
                    col = "purple", fill = NA, shape = 21, size = 2.5) +
         geom_text(data = subset(dat, rownames(dat) %in% values$selprot),
                   aes(label = protein_name, col = pval > -log10(0.05)), size = 5, 
                   position = position_nudge(x=1, y=1)) +
         scale_color_manual(values = c("grey20","firebrick"),
                            labels = c("non-significant","p < 0.05", " ")) +
         scale_y_continuous(limits = c(min(dat$l2fc, na.rm = T),
                                       max(dat$l2fc, na.rm = T))) +
         scale_x_continuous(limits = c(min(dat$av_expr, na.rm = T),
                                       max(dat$av_expr, na.rm = T))) +
         xlab("Abundance - Average [LFQ]") +
         ylab("Fold Change - Leuk / Healthy [log2]") +
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
   
   # update prot list from clicks/brushes
   observeEvent(input$scatter_brush, {
      brushed <- row.names(brushedPoints(dat, input$scatter_brush))
      values$selprot <- unique(c(brushed, values$selprot))
   })
   observeEvent(input$scatter_click, {
      clicked <- row.names(nearPoints(dat, input$scatter_click))
      values$selprot <- unique(c(clicked, values$selprot))
   })
   observeEvent(input$volcano_brush, {
      brushed <- row.names(brushedPoints(dat, input$volcano_brush))
      values$selprot <- unique(c(brushed, values$selprot))
   })
   observeEvent(input$volcano_click, {
      clicked <- row.names(nearPoints(dat, input$volcano_click))
      values$selprot <- unique(c(clicked, values$selprot))
   })
   observeEvent(input$ma_brush, {
      brushed <- row.names(brushedPoints(dat, input$ma_brush))
      values$selprot <- unique(c(brushed, values$selprot))
   })
   observeEvent(input$ma_click, {
      clicked <- row.names(nearPoints(dat, input$ma_click))
      values$selprot <- unique(c(clicked, values$selprot))
   })

   # update prot list from input text
   observeEvent(input$input_go_bttn, {
      id_list <- strsplit(input$id_input, split = ",| |\n")[[1]]
      id_list <- id_list[lapply(id_list, nchar) > 0]
      prot <- id_list[id_list %in% dat$protein_name]
      unipro <- id_list[id_list %in% dat$uniprot_id]
      if (length(unipro) > 0) {
         prot <- c(prot, dat[dat$uniprot_id %in% unipro,]$protein_name)
      }
      genelist <- sapply(dat$gene_symbols, function(x) {
         tmp <- strsplit(x, split = " |;")[[1]]
         tmp[!tmp==""] })
      genes <- id_list[id_list %in% unlist(genelist)]
      if (length(genes) > 0) {
         idx <- unlist(lapply(genelist, function(x) any(x %in% genes)))
         prot <- c(prot, dat[idx,]$protein_name)
      }
      missing <- id_list[!id_list %in% c(dat$protein_name, dat$uniprot_id, unlist(genelist))]
      if (length(missing) > 0) {
         print(missing)
         output$missing_prot <- renderText(paste0("Following names could not be found: ",
                                                  paste(missing, collapse = ", ")))
      } else {
         output$missing_prot <- NULL
      }
      if (length(prot) > 0) {
        values$selprot <- unique(c(values$selprot, rownames(dat[dat$protein_name %in% prot,])))
      }
   })
   
   # reset protein selection
   observeEvent(input$reset_bttn, {
      values$selprot <- vector()
   })
   
   # datatable of info on selected proteins
   output$protein_info_dt <- DT::renderDT({
      sel_prot <- values$selprot
      if (is_empty(sel_prot)) {
         seldat <- dat[, c(1:2,15,17:20,3:14)]
      } else {
         seldat <- dat[sel_prot, c(1:2,15,17:20,3:14)]
      }
      colnames(seldat) <- c("Protein Name", "UniProt ID", "Gene Symbols", 
                            "Specificity","TM Domain", "# TM Domain", "Membrane Type",
                            "L2FC", "p-value", "q-value",
                            "Av Expr", "Av Normal", "Av leukemia", 
                            "normal r1", "normal r2", "normal r3",
                            "leukemia r1", "leukemia r2", "leukemia r3")
      seldat[,sapply(seldat, is.numeric)] <- apply(seldat[,sapply(seldat, is.numeric)],
                                                   2, function(x) signif(x, 2))
      DT::datatable(seldat, rownames = F,
                    extensions = 'Buttons',
                    options = list(dom = 'Bfrtip',
                                   buttons = list(list(extend = 'copy',
                                                       title = NULL),
                                                  list(extend = 'print',
                                                       title = "Surface Proteomics - Selected Proteins"),
                                                  list(extend = 'csv',
                                                       filename = paste0("surfprot_selected_proteins_info.", Sys.Date()),
                                                       title = NULL),
                                                  list(extend = 'excel',
                                                       filename = paste0("surfprot_selected_proteins_info.", Sys.Date()),
                                                       title = NULL),
                                                  list(extend = 'pdf',
                                                       filename = paste0("surfprot_selected_proteins_info.", Sys.Date()),
                                                       title = "Surface Proteomics - Selected Protein")),
                                   columnDefs = list(list(
                                      targets = c(3,6),
                                      render = JS(
                                         "function(data, type, row, meta) {",
                                         "return type === 'display' && data.length > 10 ?",
                                         "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                                         "}")))),
                    callback = JS('table.page(3).draw(false);'))
   })
      
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

