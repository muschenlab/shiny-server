library(shiny)
library(tidyverse)
library(shinyWidgets)
library(shinythemes)
library(colourpicker)
library(org.Hs.eg.db)
library(survival)
library(survminer)
library(officer)
library(rvg)
library(DT)
options(stringsAsFactors = F)

# geneinfo <- AnnotationDbi::select(org.Hs.eg.db,
#                                   keys =  keys(org.Hs.eg.db,
#                                                keytype = "SYMBOL"),
#                                   keytype = "SYMBOL",
#                                   columns = c("ALIAS","ENTREZID",
#                                               "ENSEMBLID","GENENAME"))
# saveRDS(aliases, "/poolio/gene_info_hs.rds")
genes <- readRDS("/poolio/gene_info_hs.rds")

# tmp until i sort out db
# datasets <- list.files("/poolio/public_data/GEO_expr/rds", pattern = '.rds', full.names = T)
# names(datasets) <- sub("^.+/(.+)\\.[0-9\\-]+\\.rds", "\\1", datasets)
# datasets <- lapply(datasets, function(x) readRDS(x))
# saveRDS(datasets, "/poolio/public_data/GEO_expr/rds/combined_datsets.2019-10-24.rds")
datasets <- readRDS("/poolio/public_data/GEO_expr/rds/comb_temp.rds")
pmids <- sapply(datasets, function(x) {
  tmp <- experimentData(x)@pubMedIds
  tmp <- unlist(strsplit(tmp, split = "\n"))
  if (is_null(tmp)) {
    return("NA") 
  } else {
    paste0(sapply(tmp, function(y) {
      paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed/?term=",y,"'>",y,"</a>") 
    }), collapse = "; ") 
  }
})
surv_dat <- sapply(datasets, function(x) sum(grepl("censor_death", colnames(pData(x)))))
surv_dat <- names(datasets)[which(surv_dat == 1)]
relap_dat <- sapply(datasets, function(x) sum(grepl("censor_relapse", colnames(pData(x)))))
relap_dat <- names(datasets)[which(relap_dat == 1)]
dataset_info <- data.frame(row.names = NULL,
                           "Name" = c("ECOG B-ALL","MMMLNP DLBCL","COG B-ALL","BCCA DLBCL","TCGA AML"),
                           "Cancer Type" = c("Adult B-ALL","DLBCL","Childhood B-ALL","DLBCL","AML"),
                           "Patient #" = sapply(datasets, function(x) ncol(x)),
                           "Survival data" = ifelse(names(datasets) %in% surv_dat, "Y", "N"),
                           "Relapse data" = ifelse(names(datasets) %in% relap_dat, "Y", "N"),
                           "PMIDs" = pmids,
                           check.names = F)
names(datasets) <- dataset_info$Name

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

#################################### UI ####################################### 
ui <- function(request) {
  fluidPage(
    
    # title
    title = "Prognosis", 
    theme = shinytheme("flatly"),
    titlePanel(tags$h3(tags$a(
      imageOutput("icon", height = "50px", width = "50px", inline = TRUE),
      href="http://10.197.211.94:80"), 
      "Prognosis data")),
    
    sidebarLayout(
      
      #### sidebar ####
      sidebarPanel(width = 3,
        div(style="display:inline-block", textInput("gene", "Gene", value = "TP53")),
        div(style="display:inline-block", actionButton("go_button", "GO")),
        htmlOutput("gene_check"),
        br(),
        uiOutput("choose_probes"),
        tags$small("If multiple probes selected, average will be used"),
        br(), br(),
        checkboxGroupButtons("percentile", "Percentile", c("25%","50%"),
                             justified = T, selected = "50%"),
        br(),
        div(style="display:inline-block", 
            colourInput("col_expr_low","Low expr", value = "#3498dbff",
                        allowTransparent = TRUE)),
        div(style="display:inline-block", 
            colourInput("col_expr_high","High expr", value = "#C23C3C",
                        allowTransparent = TRUE)),
        br(),
        uiOutput("col_datasets")
      ),

      #### main ####
      mainPanel(width = 9,
        tabsetPanel(
          tabPanel("Study Info",
                   DT::DTOutput('dataset_dt')),
          tabPanel("Expression Level",
                   plotOutput("exprhist"),
                   br(),
                   tableOutput('expr_summary_dt'),
                   br(),
                   DT::DTOutput('expr_values_dt')),
          tabPanel("Survival",
                   plotOutput("survplot"),
                   div(style="display:inline-block", 
                       downloadButton("dl_surv_ppt", label = "Download PPT")),
                   div(style="display:inline-block", 
                       downloadButton("dl_surv_png", label = "Download PNG"))),
          tabPanel("Relapse",
                   plotOutput("relapseplot"),
                   div(style="display:inline-block", 
                       downloadButton("dl_relapse_ppt", label = "Download PPT")),
                   div(style="display:inline-block", 
                       downloadButton("dl_relapse_png", label = "Download PNG")))
        )
      )
    )
  )
}

################################# Server ##################################### 
server <- function(input, output) {

  # icon
  output$icon <- renderImage(list(src = "../hexagons/prognosis.png",
                                  height = "84px", width = "72px"), 
                             deleteFile = F)
  
  # get gene
  get_gene <- reactive({
    input$go_button
    input$gene
    return(isolate(input$gene))
  })

  # check gene
  output$gene_check <- renderUI({
    gene <- get_gene()
    if (is.null(gene) | gene == "") {
      return(tags$small("Gene field is empty - enter a gene name"))
    } else {
      if (gene %in% genes$SYMBOL) {
        return(NULL)
      } else {
        palias <- genes[which(genes$ALIAS == gene),]$SYMBOL
        pgenes <- genes[grep(gene, genes$SYMBOL, ignore.case = T),]$SYMBOL
        pdist <- adist(gene, pgenes)
        pgenes <- pgenes[which(pdist == min(pdist))]
        pmatch <- unique(c(palias, pgenes))
        pmatch <- ifelse(is_empty(pmatch), "none found", paste0(pmatch, collapse = ", "))
        return(tags$small(paste0("Gene symbol is not valid, potential matches: ", pmatch)))
      }
    }
  })

  # dataset info table
  output$dataset_dt <- DT::renderDT(
    DT::datatable(dataset_info, editable = "cell", escape = F, rownames = T,
                  selection = list(target = 'row', selected = c(1,3,5)))
  )
  
  # update dataset info
  observeEvent(input$dataset_dt_cell_edit, {
    dataset_info <<- editData(dataset_info, input$dataset_dt_cell_edit, 'dataset_dt')
    names(datasets) <- dataset_info$Name
  })

  # get dataset
  get_dat <- reactive({
    return(datasets[input$dataset_dt_rows_selected])
  })

  # probe choice ui
  output$choose_probes <- renderUI({
    dat <- get_dat()
    gene <- get_gene()
    if (is_empty(dat) | is_empty(gene)) {
      return()
    } else {
      lapply(seq_along(dat), function(x) {
        name <- names(dat)[x]
        probes <- fData(dat[[name]])[which(fData(dat[[name]])$gene_symbol == gene),]
        validate(need(nrow(probes) > 0, "Gene does not exist in selected dataset"))
        return(checkboxGroupButtons(paste0("probes_", x), paste0(name, " probes"),
                                    rownames(probes), justified = T,
                                    status = "info", selected = rownames(probes)[1],
                                    checkIcon = list(yes = icon("ok", lib = "glyphicon"),
                                                     no = icon("remove", lib = "glyphicon"))))
      })
    }
  })
  
  # dataset colour choice UI
  output$col_datasets <- renderUI({
    dat <- get_dat()
    col_presets <- c("grey80", "#22BD9E", "#3498DB", "#C23C3C")
    if (is_empty(dat)) {
      return()
    } else {
      lapply(seq_along(dat), function(x) {
        name <- names(dat)[x]
        precol <- col_presets[x]
        colourInput(paste0("col_dat_", x),
                    paste0(name, " col"),
                    value = precol,
                    allowTransparent = TRUE)
      })
    }
  })
  
  # get expr levels
  get_expr <- reactive({
    dat <- get_dat()
    lapply(seq_along(dat), function(x) {
      probes <- input[[paste0("probes_", x)]]
      if (is_empty(probes)) {
        return(NULL)
      }
      if (length(probes) > 1) {
        colMeans(exprs(dat[[x]])[rownames(dat[[x]]) %in% probes,], na.rm = T)
      } else {
        exprs(dat[[x]])[rownames(dat[[x]]) %in% probes,]
      }
    })
  })

  # get exprcats
  get_exprcats <- reactive({
    percs <- input$percentile
    expr <- get_expr()
    cats <- list()
    if (is.null(expr)) {
      return(NULL)
    } else {
      lapply(expr, function(x) {
        if ("50%" %in% percs) {
          th <- quantile(x, .5)
          cats$Median <- ifelse(x > th, "Upper", "Lower")
        }
        if ("25%" %in% percs) {
          lwth <- quantile(x, .25)
          upth <- quantile(x, .75)
          cats$Quartiles <- ifelse(x > upth, "Upper", NA)
          cats$Quartiles[x < lwth] <- "Lower"
        }
        return(cats)
      })
    }
  })

  # organise survival data
  get_osdat <- reactive({
    dat <- get_dat()
    expr <- get_expr()
    cats <- get_exprcats()
    if (is_empty(cats) | is_empty(dat)) {
      return(NULL)
    }
    lapply(seq_along(dat), function(x) {
      setname <- names(dat)[x]
      lapply(seq_along(cats[[x]]), function(y) {
        osdat <- data.frame(time = pData(dat[[x]])$follow_up_time,
                            events = pData(dat[[x]])$censor_death,
                            exprlevel = expr[[x]],
                            exprcats = cats[[x]][[y]],
                            exprset = names(cats[[x]])[y],
                            dataset = setname)
        osdat$exprcats <- factor(osdat$exprcats, levels = c("Lower", "Upper"))
        return(osdat)
      }) 
    })
  })

  # plot expr histogram
  output$exprhist <- renderPlot({
    dat <- get_dat()
    exprdat <- get_osdat()
    plotdat <- lapply(exprdat, function(x) { x %>% bind_rows() }) %>% bind_rows()
    plotdat <- plotdat[plotdat$exprset == unique(plotdat$exprset)[1],]
    plotdat$dataset <- factor(plotdat$dataset, levels = names(dat))
    cutoffs <- plotdat %>% group_by(dataset) %>% 
      summarise(median = median(exprlevel),
                lower = quantile(exprlevel, .25),
                upper = quantile(exprlevel, .75))
    collist <- list()
    for (x in seq_along(exprdat)) {
      collist[names(dat)[x]] <- input[[paste0("col_dat_", x)]] 
    }
    ggplot(plotdat, aes(x = exprlevel, fill = dataset)) + 
      geom_histogram(col = "black", binwidth = .05, lwd = .2) +
      geom_vline(data = cutoffs, aes(xintercept = median)) +
      geom_vline(data = cutoffs, aes(xintercept = upper), linetype = "dashed") + 
      geom_vline(data = cutoffs, aes(xintercept = lower), linetype = "dashed") +
      facet_wrap(~dataset, ncol = 1) +
      xlab("Expression level") + ylab("Frequency") +
      theme_bw(base_size = 16, base_family = "Arial") +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      scale_fill_manual(values = collist) +
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            axis.line = element_line(color = "black"),
            axis.text = element_text(color = "black"))
  })
  
  # render expr summary table
  output$expr_summary_dt <- renderTable({
    exprdat <- get_osdat()
    exprdat <- lapply(exprdat, function(x) { x %>% bind_rows() }) %>% bind_rows()
    exprdat <- exprdat[exprdat$exprset == unique(exprdat$exprset)[1],]
    exprdat <- group_by(exprdat, dataset) %>% 
      summarise("#" = length(exprlevel[!is.na(exprlevel)]),
                "Median" = round(median(exprlevel), 2),
                "Lower quartile" = round(quantile(exprlevel, .25), 2),
                "Upper quartile" = round(quantile(exprlevel, .75), 2))
    colnames(exprdat)[1] <- "Dataset"
    exprdat
  })

  # render expr full table
  output$expr_values_dt <- DT::renderDT({
    exprdat <- get_osdat()
    gene <- get_gene()
    tabdat <- lapply(exprdat, function(x) { 
      tab <- x[[1]][,c(6,1:3)]
      tmp <- lapply(x, function(y) {
        colnames(y)[4] <- paste0("Expr category: ", unique(y[,5]))
        return(y[, 4, drop = F])
      }) %>% bind_cols()
      cbind(tab, tmp)
    }) %>% bind_rows()
    tabdat$time <- round(tabdat$time,2)
    tabdat$exprlevel <- round(tabdat$exprlevel,2)
    tabdat$events <- factor(tabdat$events, levels = c(1,0), labels = c("Dead","Alive"))
    colnames(tabdat)[c(1:4)] <- c("Dataset", "Follow up time",
                                  "Status", paste0(gene," expr level"))
    DT::datatable(tabdat, escape = F, rownames = F, selection = "none")
  })
  
  # fit survival
  get_fits <- reactive({
    osdat <- get_osdat()
    if (is.null(osdat)) {
      return(NULL)
    }
    lapply(osdat, function(set) {
      lapply(set, function(cat) {
        fit <- survfit(Surv(time, events) ~ exprcats, data = cat)
        fit$pval_txt <- surv_pvalue(fit, data = cat)$pval.txt
        return(fit)
      })
    })
  })
  
  # generate KM plot
  plot_surv <- reactive({
    fits <- get_fits()
    cats <- get_exprcats()
    dat <- get_dat()
    if (is.null(fits)) return()
    plotdat <- lapply(seq_along(fits), function(x) {
      setname <- names(dat)[x]
      lapply(seq_along(fits[[x]]), function(y) {
        catname <- names(cats[[x]])[y]
        df <- data.frame(time = fits[[x]][[y]]$time,
                         surv = fits[[x]][[y]]$surv,
                         group = rep(names(fits[[x]][[y]]$strata),
                                     fits[[x]][[y]]$strata),
                         exprset = catname,
                         dataset = setname)
        df <- df[rep(seq_len(nrow(df)), each = 2), ]
        df$surv <- c(1,df$surv[1:(length(df$surv)-1)])
        switchidx <- (fits[[x]][[y]]$strata[[1]] + 1) * 2
        df <- df[c(1,1:switchidx,switchidx:nrow(df)),]
        df[switchidx:(switchidx+1),]$surv <- 1
        df[c(1,switchidx),]$time <- 0
        return(df)
      }) %>% bind_rows()
    }) %>% bind_rows()
    pvals <- lapply(seq_along(dat), function(x) {
      setname <- names(dat)[x]
      lapply(seq_along(fits[[x]]), function(y) {
        catname <- names(cats[[x]])[y]
        data.frame(dataset = setname,
                   exprset = catname,
                   pval = fits[[x]][[y]]$pval_txt,
                   nhigh = fits[[x]][[y]]$n[1],
                   nlow = fits[[x]][[y]]$n[2])
      }) %>% bind_rows()
    }) %>% bind_rows()
    pidx <- match(paste0(plotdat$dataset, "_", plotdat$exprset),
                  paste0(pvals$dataset, "_", pvals$exprset))
    plotdat$set_title <- paste0("(", str_trim(pvals[pidx,]$pval), 
                                "; #high = ", pvals[pidx,]$nhigh, 
                                "; #low = ", pvals[pidx,]$nlow, ")")
    ggplot(plotdat, aes(x = time, y = surv)) +
      geom_line(aes(col = group), lwd = 1.1, alpha = .6) +
      facet_wrap(dataset ~ set_title, scales = "free") +
      theme_classic(base_size = 16) +
      scale_color_manual(values = c(input$col_expr_low, input$col_expr_high), 
                         name = "Expression level", 
                         labels = c("Low","High")) +
      scale_y_continuous(limits = c(0,1), expand = c(0,0), name = "OS probability") +
      scale_x_continuous(expand = c(0,0), name = "Time [years]") +
      theme(strip.background = element_blank(),
            legend.position = "bottom", 
            panel.spacing = unit(1, "line"))
  })
  
  # output KM
  output$survplot <- renderPlot({
    if (is_empty(plot_surv())) return(NULL)
    plot_surv()
  })
  
  # download KM as ppt
  output$dl_surv_ppt <- downloadHandler(
    filename = function() {
      paste0(get_gene(), "_survplot_", Sys.Date(), ".pptx")
    },
    content = function(file) {
      fits <- get_fits()
      n <- length(flatten(fits))
      file_pptx <- tempfile(fileext = ".pptx")
      gen_pptx(plot_surv(), file_pptx, width = 3.8*n, height = 4.4)
      file.rename(from = file_pptx, to = file)
    }
  )
  
  # download KM as png
  output$dl_surv_png <- downloadHandler(
    filename = function() {
      paste0(get_gene(), "_survplot_", Sys.Date(), ".png")
    },
    content = function(file) {
      fits <- get_fits()
      n <- length(flatten(fits))
      plot <- plot_surv()
      ggsave(plot = plot, filename = file, width = 3.8*n, height = 4.4)
    }
  )
  
  # organise relapse data
  get_rldat <- reactive({
    dat <- get_dat()
    expr <- get_expr()
    cats <- get_exprcats()
    if (is.null(cats)) {
      return(NULL)
    }
    lapply(seq_along(dat), function(x) {
      setname <- names(dat)[x]
      relapsets <- dataset_info[dataset_info[,"Relapse data"] == "Y",]$Name
      if (!setname %in% relapsets) return(NULL)
      lapply(seq_along(cats[[x]]), function(y) {
        rldat <- data.frame(time = pData(dat[[x]])$years_to_relapse,
                            relapse_status = pData(dat[[x]])$censor_relapse,
                            exprlevel = expr[[x]],
                            exprcats = cats[[x]][[y]],
                            exprset = names(cats[[x]])[y],
                            dataset = setname)
        rldat$exprcats <- factor(rldat$exprcats, levels = c("Lower", "Upper"))
        rldat$relapse_status <- factor(rldat$relapse_status, levels = c(0,1),
                                       labels = c("censored","relapsed"))
        return(rldat)
      }) %>% bind_rows()
    }) %>% bind_rows()
  })
  
  # plot relapse
  gen_relapse_plot <- reactive({
    plotdat <- get_rldat()
    if (is_null(plotdat)) return(NULL)
    ggplot(plotdat, aes(x = relapse_status, y = exprlevel)) +
      geom_dotplot(binaxis = "y", stackdir = "center") +
      facet_wrap(dataset ~ exprset, scales = "free") +
      theme_classic(base_size = 16) +
      theme(strip.background = element_blank(),
            legend.position = "bottom", 
            panel.spacing = unit(1, "line"))
  })
  
  # output KM
  output$relapseplot <- renderPlot({
    gen_relapse_plot()
  })
}

############### run ############### 
shinyApp(ui = ui, server = server)
