library(shiny)
library(tidyverse)
library(shinyWidgets)
library(shinycssloaders)
library(shinythemes)
library(officer)
library(rvg)
library(pzfx)
library(ggpubr)
library(colourpicker)
library(ggnewscale)
library(readxl)
options(stringsAsFactors = F)

################################# Setup #######################################

# function to generate pptx
gen_pptx <- function(notxt_plot, plot_outline, file, height = 4, width = 4, left = 0, top = 0) {
    read_pptx() %>% 
        add_slide(layout = "Title and Content", master = "Office Theme") %>% 
        ph_with(value = notxt_plot,
                location = ph_location(height = height, width = width,
                                        left = left, top = top)) %>%
        ph_with(value = dml(ggobj = plot_outline),
                location = ph_location(height = height, width = width,
                                       left = left, top = top),
                bg = "transparent") %>%
        print(target = file)
}

# density function
get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}

# load internal datasets
internal_files <- list.files("/srv/shiny-server/app_scwb/scwb_data/", pattern = ".xlsx", full.names = T)
names(internal_files) <- sub("^.+/(.+)\\.xlsx", "\\1", internal_files)

#################################### UI #######################################
ui <- function(request) {
    fluidPage(
        
        # title
        title = "scWB", 
        theme = shinytheme("lumen"),
        titlePanel(tags$h3(tags$a(
            imageOutput("icon", height = "50px", width = "50px", inline = TRUE),
            href="http://10.197.211.94:80"), 
            "single cell Western Blot")),
        
        sidebarLayout(
            sidebarPanel(width = 3,
                         selectInput(inputId = "internal_file",
                                     label = "Choose an existing dataset:",
                                     choice = names(internal_files)),
                         fileInput(inputId = "input_file",
                                   label = "Upload an excel file:",
                                   multiple = F,
                                   accept = c(".xlsx")),
                         selectInput(inputId = "palette",
                                     label = "Select palette",
                                     selected = "OrRd",
                                     choices = c("Greys","Purples","Blues","Greens",
                                                 "Oranges","Reds","YlOrBr","YlOrRd",
                                                 "OrRd","PuRd","RdPu","PuBu","GnBu",
                                                 "PuBu","YlGnBu","PuBuGn","BuGn","YlGn",
                                                 "inferno","magma","plasma","viridis","cividis")),
                         div(style="display:inline-block;width:48%;", 
                             textInput("xlab_name", "X-axis target:")),
                         div(style="display:inline-block;width:48%;", 
                             textInput("ylab_name", "Y-axis target:")),
                         div(style="display:inline-block;width:48%;", 
                             sliderInput("xth", "X-axis theshold:", value = 2,
                                     min = 0, max = 20)),
                         div(style="display:inline-block;width:48%;", 
                             sliderInput("yth", "Y-axis theshold:", value = 2,
                                     min = 0, max = 20))
            ),
            mainPanel(width = 9,
                      uiOutput("plot_ui") %>% 
                          withSpinner(color="#0dc5c1"),
                      div(style="display:inline-block;",
                          downloadButton("dl_main_plot_ppt", label = "PPT",
                                         style = "font-size:12px;height:30px;padding:5px;"),
                          downloadButton("dl_main_plot_png", label = "PNG",
                                         style = "font-size:12px;height:30px;padding:5px;"),
                          bookmarkButton(label = "Save Session", 
                                         style = "font-size:13px;height:30px;padding:5px;")),
                      br(), br(),
                      tags$h4("Summary Stats:"),
                      DT::DTOutput('stats_dt')
            )
        )
        
    )
}

################################# Server ######################################
server <- function(input, output) {
    
    # icon
    output$icon <- renderImage(list(src = "../hexagons/scwb.png",
                                    height = "90px", width = "85px"), 
                               deleteFile = F)
    
    # load data
    load_data <- reactive({
        if (!is_empty(input$input_file)) {
            path <- input$input_file$datapath
        } else {
            path <- internal_files[[input$internal_file]]
        }
        sheets <- excel_sheets(path = path)
        df <- lapply(sheets, function(sheet) {
            tmp <- read_excel(path = path, sheet = sheet)
            tmp$condition <- sheet
            return(tmp)
        }) %>% bind_rows() %>% na.omit() %>% as.data.frame() 
        if (input$xlab_name == "") {
            colnames(df)[1] <- gsub("[ \\-]", "_", colnames(df)[1])
        } else {
            colnames(df)[1] <- input$xlab_name
        }
        if (input$ylab_name == "") {
            colnames(df)[2] <- gsub("[ \\-]", "_", colnames(df)[2])
        } else {
            colnames(df)[2] <- input$ylab_name
        }
        df[!is.finite(df[,1]), 1] <- 0
        df[!is.finite(df[,2]), 2] <- 0
        df[,1] <- df[,1] + runif(nrow(df), 0, 2)
        df[,2] <- df[,2] + runif(nrow(df), 0, 2)
        df[,1] <- ifelse(df[,1] < 1, df[,1], log2(df[,1]))
        df[,2] <- ifelse(df[,2] < 1, df[,2], log2(df[,2]))
        print(head(df))
        return(df)
    })
    
    # get name from file
    get_name <- reactive({
        if (!is_empty(input$input_file)) {
            gsub(" ", "_", sub("^(.+)\\.(xlsx)", "\\1", input$input_file$name))
        } else {
            gsub(" ", "_", sub("^(.+)\\.(xlsx)", "\\1", input$internal_file))
        }
    })
    
    # summary stats
    get_stats <- reactive({
        df <- load_data()
        xth <- input$xth
        yth <- input$yth
        xlab <- sym(colnames(df)[1])
        ylab <- sym(colnames(df)[2])
        if (is_empty(df)) return(NULL)
        stat_res <- df %>% group_by(condition) %>%
            summarise(mean_x = round(mean(!!xlab), 2),
                      mean_y = round(mean(!!ylab), 2),
                      percent_x = round(sum(!!xlab > xth)/length(!!xlab), 2),
                      percent_y = round(sum(!!ylab > yth)/length(!!ylab), 2),
                      percent_dp = round(sum(!!xlab > xth & !!ylab > yth) /
                                             length(!!xlab), 2) * 100) %>%
            as.data.frame()
        colnames(stat_res)[2:6] <- c(paste0(xlab," mean"), paste0(ylab," mean"),
                                     paste0(xlab, " percent"), paste0(ylab, " percent"),
                                     "percent DP")
        return(stat_res)
    })
    
    # show summary stats
    output$stats_dt <- DT::renderDT(
        DT::datatable(get_stats(), 
                      escape = F, 
                      rownames = T, 
                      selection = list(target = 'column', mode = 'single', selected = 6),
                      extensions = 'Buttons',
                      options = list(dom = 'frtipB',
                                     buttons = list('pageLength',
                                                    list(extend = 'copy',
                                                         title = NULL),
                                                    list(extend = 'csv',
                                                         filename = paste0(get_name(), "_stat_res.csv"),
                                                         title = NULL),
                                                    list(extend = 'excel',
                                                         filename = paste0(get_name(), "_stat_res.xlsx"),
                                                         title = NULL))),
                      caption = "Select rows to plot.")
    )
    
    # generate dotplot
    gen_plots <- reactive({
        df <- load_data()
        stats <- get_stats()
        selcond <-  stats[input$stats_dt_rows_selected,]$condition
        if (length(selcond) > 0) {
            df <- df[df$condition %in% selcond,]
            stats <- stats[stats$condition %in% selcond,]
        }
        df$condition <- factor(df$condition, levels = unique(df$condition))
        xlab <- sym(colnames(df)[1])
        ylab <- sym(colnames(df)[2])
        xth <- input$xth
        yth <- input$yth
        if (is_empty(df)) return(NULL)
        df$density <- get_density(df[,1], df[,2], n = 100)
        stats$x <- 21; stats$y <- 22
        baseplot <- ggplot(df, aes(x = !!xlab, y = !!ylab)) +
            scale_x_continuous(limits = c(-2,23), name = paste0(xlab, " intensity [AU]")) +
            scale_y_continuous(limits = c(-2,23), name = paste0(ylab, " intensity [AU]")) +
            facet_wrap(~condition, ncol = 4) +
            theme(panel.grid = element_blank(),
                  axis.title = element_text(size = 14, color = "grey20", family = "Arial"),
                  axis.text = element_text(size = 11, color = "grey20", family = "Arial"),
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, color = "black"),
                  strip.background = element_blank(),
                  strip.text = element_text(size = 14, color = "black", family = "Arial"),
                  legend.position = "none")
        if (input$palette %in% c("inferno","magma","plasma","viridis","cividis")) {
            baseplot <- baseplot + scale_fill_viridis_c(option = input$palette) +
                scale_color_viridis_c(option = input$palette)
        } else {
            baseplot <- baseplot + scale_fill_distiller(palette = input$palette, direction=-1) +
                scale_color_distiller(palette = input$palette, direction=-1)
        }
        if (!is_empty(input$stats_dt_columns_selected)) {
            stats$perc <- stats[,(input$stats_dt_columns_selected)]
            outline <- baseplot +
                geom_text(data = stats, aes(x = x, y = y, label = perc), size = 5)
        } else {
            outline <- baseplot
        }
        full <- baseplot + 
            geom_hline(yintercept = yth, lty = 2, color = "grey", alpha = .8) +
            geom_vline(xintercept = xth, lty = 2, color = "grey", alpha = .8) +
            geom_point(aes(color = density), alpha = .4, size = .3) +
            stat_density_2d(aes(fill = ..level..), geom = "polygon", h = 2, alpha = .6) 
        if (!is_empty(input$stats_dt_columns_selected)) {
            stats$perc <- stats[,(input$stats_dt_columns_selected)]
            full <- full +
                geom_text(data = stats, aes(x = x, y = y, label = perc), size = 5)
        } 
        notxt <- baseplot +
            geom_hline(yintercept = yth, lty = 2, color = "grey", alpha = .8) +
            geom_vline(xintercept = xth, lty = 2, color = "grey", alpha = .8) +
            geom_point(aes(color = density), alpha = .4, size = .3) +
            stat_density_2d(aes(fill = ..level..), geom = "polygon", h = 2, alpha = .6) +
            theme(text = element_text(color = "white"),
                  axis.text = element_text(color = "white"),
                  strip.text = element_text(color = "white"),
                  axis.title = element_text(color = "white"))
        return(list(full = full, outline = outline, notxt = notxt))
    })
    
    # output main plot
    output$main_plot <- renderPlot({
        plot <- gen_plots()$full
        validate(need(!is_empty(plot), "Upload a file."))
        return(plot)
    })
    get_dims <- reactive({
        stats <- get_stats()
        selcond <- stats[input$stats_dt_rows_selected,]$condition
        ncond <- ifelse(length(selcond) > 0, length(selcond), length(stats$condition))
        width <- ifelse(ncond < 4, 250*ncond, 800)
        height <- ifelse(ncond < 4, 270*ceiling(ncond/4), 220*ceiling(ncond/4))
        list(w = width, h = height)
    })
    output$plot_ui <- renderUI({
        dims <- get_dims()
        plotOutput("main_plot", height = dims$h, width = dims$w)
    })
    
    # download main plot as png
    output$dl_main_plot_png <- downloadHandler(
        filename = function() {
            name <- get_name()
            paste0(name, "_scWB_contourplots.png")
        },
        content = function(file) {
            plot <- gen_plots()$full
            dims <- get_dims()
            ggsave(plot = plot, filename = file, units = "mm",
                   height = dims$h/3, width = dims$w/3)
        }
    )
    
    # download main plot as ppt
    output$dl_main_plot_ppt <- downloadHandler(
        filename = function() {
            name <- get_name()
            paste0(name, "_scWB_contourplots.pptx")
        },
        content = function(file) {
            plots <- gen_plots()
            dims <- get_dims()
            file_pptx <- tempfile(fileext = ".pptx")
            gen_pptx(plots$notxt, plots$outline, file_pptx, 
                     height = (dims$h/3)*0.039,
                     width = (dims$w/3)*0.039)
            file.rename(from = file_pptx, to = file)
        }
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server, enableBookmarking = "server")
