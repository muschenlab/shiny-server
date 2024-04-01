library(tidyverse)
library(rstatix)
library(coin)
library(RMariaDB)

# ## gene info
# library(org.Hs.eg.db)
# orgdb <- org.Hs.eg.db
# keys <- AnnotationDbi::keys(orgdb, keytype='ENTREZID')
# gene_info <- AnnotationDbi::select(orgdb, keys = keys,
#                                    columns = c("SYMBOL", "GENENAME")) %>%
#   dplyr::rename(entrez_id = ENTREZID,
#          gene_symbol = SYMBOL,
#          gene_name = GENENAME)
# saveRDS(gene_info, "data/geneinfo.2023-10-31.rds")
gene_info <- readRDS("data/geneinfo.2023-10-31.rds")

### database config ###
dbname <- "DDDB"
cnf <- list.files(path = "data", pattern = paste0(dbname, ".cnf$"), full.names = T)

### load drug info ###
db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
query <- paste0('SELECT * FROM drug_info;')
queryres <- dbSendQuery(db, query)
di <- dbFetch(queryres)
dbClearResult(queryres)
dbDisconnect(db)

### load sample info ###
db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
query <- "SELECT * FROM sample_info;"
queryres <- dbSendQuery(db, query)
si <- dbFetch(queryres)
dbClearResult(queryres)
dbDisconnect(db)

### filter non-cancerous/unknown cell types 
si <- si %>%
  filter(!ds_type %in% c("NOS", "Unknown", "Non-cancerous", "Engineered", "Other"))

### set disease groupings ###
silist <- list("Solid tumor" = si[which(!si$ds_type %in% c("B-cell",
                                                           "T-cell",
                                                           "Myeloid",
                                                           "Plasma cell")),],
               "Haematopoietic" = si[which(si$ds_type %in% c("B-cell",
                                                             "T-cell",
                                                             "Myeloid",
                                                             "Plasma cell")),],
               "B-cell leukemia" = si[grep("B-ALL", si$ds_subtype),],
               "T-cell leukemia" = si[grep("T-ALL", si$ds_subtype),],
               "B-cell lymphoma" = si[which(si$ds_type == "B-cell" &
                                              grepl("lymphoma", si$ds_subtype)),],
               "Glioma" = si[which(si$ds_subtype %in% c("Glioblastoma",
                                                        "Glioma",
                                                        "Astrocytoma")),],
               "CML" = si[grep("CML", si$ds_subtype),])
cts <- unique(si$ds_type); names(cts) <- cts
tmp <- lapply(cts, function(ct) si[which(si$ds_type == ct),])
silist <- c(silist, tmp)
# remove low n
silist <- silist[!names(silist) %in% c("Testicular", "Embryonal", "Eye",
                                       "Gallbladder", "Skin carcinoma")]

### set comaprisons ###
comp <- data.frame(a = c("CML","B-cell","T-cell", "T-cell leukemia", names(silist)[-1]),
                   b = c("Solid tumor","Myeloid","B-cell", "B-cell leukemia", rep("Solid tumor", length(silist)-1)))

#################### calc stats for each disease group ###########################
dir.create("data/ctres3/")
lapply(1:nrow(comp), function(compidx) {
  
  ct1 <- comp[compidx, ]$a
  ct2 <- comp[compidx, ]$b
  print(paste0(ct1, " vs ", ct2))
  if (file.exists(paste0("data/ctres3/", ct1, "_vs_", ct2, "_diffsens.rds"))) return(NULL)
  
  ######### sample info #########
  ct1_si <- silist[[ct1]] %>%
    mutate(ct = ct1)
  ct2_si <- silist[[ct2]] %>%
    filter(!sample_id %in% ct1_si$sample_id) %>%
    mutate(ct = ct2)
  sisub <- rbind(ct1_si, ct2_si)
  
  ################### Gene Dep ######################
  
  # get dependency data
  db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
  query <- paste0('SELECT * FROM dep_data ',
                  'WHERE sample_id IN ("',
                  paste0(sisub$sample_id, collapse = '", "'),
                  '");')
  queryres <- dbSendQuery(db, query)
  depdata <- dbFetch(queryres)
  dbClearResult(queryres)
  dbDisconnect(db)
  
  ######### RNAi #########
  # rnai stats
  rnaidat <- depdata %>%
    filter(assay_type == "RNAi") %>%
    mutate(ct = sisub[match(sample_id, sisub$sample_id),]$ct) %>%
    mutate(ct = factor(ct, levels=c(ct1, ct2)))
  av <- rnaidat %>%
    group_by(ct, entrez_id) %>%
    summarise(avRNAi = mean(score, na.rm=T),
              nRNAi = n()) %>%
    pivot_wider(names_from = ct,
                values_from = c(avRNAi, nRNAi))
  nfilt <- filter_at(av, vars(starts_with("nRNAi")), all_vars(. > 3))
  if (nrow(nfilt) < 1) {
    rnai_summary <- av %>%
      mutate(dRNAi = NA,
             rRNAi = NA,
             pRNAi = NA,
             RNAi_prcntl = NA)
  } else {
    stats <- rnaidat %>%
      filter(entrez_id %in% nfilt$entrez_id) %>%
      group_by(entrez_id) %>% 
      rstatix::wilcox_test(score ~ ct, ref.group = ct1, detailed = T)
    eff <- rnaidat %>%
      filter(entrez_id %in% nfilt$entrez_id) %>%
      group_by(entrez_id) %>%
      rstatix::wilcox_effsize(score ~ ct, ref.group = ct1)
    stats <- merge(stats, eff[,4:5], by = c("entrez_id")) %>%
      dplyr::select(entrez_id, estimate, effsize, p) %>%
      mutate(effsize = ifelse(estimate < 0, -effsize, effsize)) %>%
      dplyr::rename(dRNAi = estimate,
             rRNAi = effsize,
             pRNAi = p) %>%
      mutate(RNAi_prcntl = rank(-rRNAi)/length(rRNAi)) 
    rnai_summary <- merge(av, stats, by = "entrez_id", all = T) %>%
      arrange(-RNAi_prcntl)
  }
  
  ########## CRISPR ##########
  crisprdat <- depdata %>%
    filter(assay_type == "CRISPR") %>%
    mutate(ct = sisub[match(sample_id, sisub$sample_id),]$ct) %>%
    mutate(ct = factor(ct, levels=c(ct1, ct2)))
  av <- crisprdat %>%
    group_by(ct, entrez_id) %>%
    summarise(avCRISPR = mean(score, na.rm=T),
              nCRISPR = n()) %>%
    pivot_wider(names_from = ct,
                values_from = c(avCRISPR, nCRISPR))
  nfilt <- filter_at(av, vars(starts_with("nCRISPR")), all_vars(. > 3))
  if (nrow(nfilt) < 1) {
    crispr_summary <- av %>%
      mutate(dCRISPR = NA,
             rCRISPR = NA,
             pCRISPR = NA,
             CRISPR_prcntl = NA)
  } else {
    stats <- crisprdat %>%
      filter(entrez_id %in% nfilt$entrez_id) %>%
      group_by(entrez_id) %>% 
      rstatix::wilcox_test(score ~ ct, ref.group = ct1, detailed = T)
    eff <- crisprdat %>%
      filter(entrez_id %in% nfilt$entrez_id) %>%
      group_by(entrez_id) %>%
      rstatix::wilcox_effsize(score ~ ct, ref.group = ct1)
    stats <- merge(stats, eff[,4:5], by = c("entrez_id")) %>%
      dplyr::select(entrez_id, estimate, effsize, p) %>%
      mutate(effsize = ifelse(estimate < 0, -effsize, effsize)) %>%
      dplyr::rename(dCRISPR = estimate,
             rCRISPR = effsize,
             pCRISPR = p) %>%
      mutate(CRISPR_prcntl = rank(-rCRISPR)/length(rCRISPR)) 
    crispr_summary <- merge(av, stats, by = "entrez_id", all = T) %>%
      arrange(-CRISPR_prcntl)
  }
  
  ####### combine CRISPR + RNAi #######
  dep_summary <- merge(crispr_summary, rnai_summary,
                       by="entrez_id", all=T)
  dep_summary <- merge(dep_summary, gene_info, by = "entrez_id") %>%
    arrange(-CRISPR_prcntl)
  
  ################### Compounds ######################
  
  # get compound data
  db <- dbConnect(RMariaDB::MariaDB(), default.file = cnf, group = dbname)
  query <- paste0('SELECT * FROM drug_metrics ',
                  'WHERE sample_id IN ("',
                  paste0(sisub$sample_id, collapse = '", "'),
                  '");')
  queryres <- dbSendQuery(db, query)
  drug_metrics <- dbFetch(queryres)
  dbClearResult(queryres)
  dbDisconnect(db)
  
  ######## CTD ########
  ctd_metrics <- drug_metrics %>%
    filter(dataset == "CTD") %>%
    mutate(ct = sisub[match(sample_id, sisub$sample_id),]$ct) %>%
    mutate(ct = factor(ct, levels=c(ct1, ct2))) %>%
    dplyr::rename(AAC = aac) %>%
    mutate(limEC50 = ifelse(logEC50 > logmaxc, log10(2*10^logmaxc), logEC50))
  ctd_summary <- lapply(c("limEC50", "AAC", "DSS1", "DSS2", "DSS3", "DSS4"), function(metric) {
    submetrics <- ctd_metrics %>%
      mutate(selmetric = .[, grep(metric, colnames(.))]) %>%
      dplyr::filter(!is.na(selmetric))
    av <- submetrics %>%
      group_by(ct, treatment_id) %>%
      summarise(n = n(),
                av = mean(selmetric, na.rm=T)) %>%
      pivot_wider(names_from = ct,
                  values_from = c(n, av))
    nfilt <- filter_at(av, vars(starts_with("n")), all_vars(. > 3))
    if (nrow(nfilt) < 1) {
      metric_summary <- av %>%
        mutate(d = NA,
               r = NA,
               p = NA,
               prcntl = NA)
    } else {
      tmp <- submetrics %>%
        filter(treatment_id %in% nfilt$treatment_id) %>%
        mutate(z = scale(selmetric)[,1]) %>%
        group_by(treatment_id) %>% 
        filter(!all(selmetric==0))
      stats <- tmp %>%
        rstatix::wilcox_test(z ~ ct, ref.group = ct1, detailed = T)
      eff <- tmp %>%
        rstatix::wilcox_effsize(z ~ ct, ref.group = ct1)
      stats <- merge(stats, eff[,4:5], by = c("treatment_id")) %>%
        dplyr::select(treatment_id, estimate, effsize, p) %>%
        mutate(effsize = ifelse(estimate < 0, -effsize, effsize)) %>%
        dplyr::rename(d = estimate,
                      r = effsize,
                      p = p) %>%
        mutate(prcntl = rank(-r)/length(r)) 
      metric_summary <- merge(av, stats, by = "treatment_id", all = T) %>%
        mutate(metric = metric)
    }
    return(metric_summary)
  }) %>% bind_rows()
  
  ######## GDSC ########
  gdsc_metrics <- drug_metrics %>%
    filter(dataset %in% c("GDSC1", "GDSC2")) %>%
    mutate(ct = sisub[match(sample_id, sisub$sample_id),]$ct) %>%
    mutate(ct = factor(ct, levels=c(ct1, ct2))) %>%
    dplyr::rename(AAC = aac) %>%
    mutate(limEC50 = ifelse(logEC50 > logmaxc, log10(2*10^logmaxc), logEC50))
  gdsc_summary <- lapply(c("limEC50", "AAC", "DSS1", "DSS2", "DSS3", "DSS4"), function(metric) {
    submetrics <- gdsc_metrics %>%
      mutate(selmetric = .[, grep(metric, colnames(.))]) %>%
      dplyr::filter(!is.na(selmetric))
    av <- submetrics %>%
      group_by(ct, treatment_id) %>%
      summarise(n = n(),
                av = mean(selmetric, na.rm=T)) %>%
      pivot_wider(names_from = ct,
                  values_from = c(n, av))
    nfilt <- filter_at(av, vars(starts_with("n")), all_vars(. > 3))
    if (nrow(nfilt) < 1) {
      metric_summary <- av %>%
        mutate(d = NA,
               r = NA,
               p = NA,
               prcntl = NA)
    } else {
      tmp <- submetrics %>%
        filter(treatment_id %in% nfilt$treatment_id) %>%
        mutate(z = scale(selmetric)[,1]) %>%
        group_by(treatment_id) %>% 
        filter(!all(selmetric==0))
      stats <- tmp %>%
        rstatix::wilcox_test(z ~ ct, ref.group = ct1, detailed = T)
      eff <- tmp %>%
        rstatix::wilcox_effsize(z ~ ct, ref.group = ct1)
      stats <- merge(stats, eff[,4:5], by = c("treatment_id")) %>%
        dplyr::select(treatment_id, estimate, effsize, p) %>%
        mutate(effsize = ifelse(estimate < 0, -effsize, effsize)) %>%
        dplyr::rename(d = estimate,
                      r = effsize,
                      p = p) %>%
        mutate(prcntl = rank(-r)/length(r)) 
      metric_summary <- merge(av, stats, by = "treatment_id", all = T) %>%
        mutate(metric = metric)
    }
    return(metric_summary)
  }) %>% bind_rows()
  
  #### merge GDSC + CTD + add drug info
  tmp <- drug_metrics %>%
    dplyr::select(treatment_id, cpd_name) %>%
    distinct()
  ctd_summary$dataset <- "CTD"
  gdsc_summary$dataset <- "GDSC"
  drug_summary <- rbind(ctd_summary, gdsc_summary) %>%
    mutate(cpd_name = tmp[match(treatment_id, tmp$treatment_id),]$cpd_name)
  drug_summary <- merge(drug_summary, di, by = "cpd_name")
  
  ### get fraction of drc within effective range for each treatment id
  inrange <- drug_metrics %>%
    mutate(ct = sisub[match(sample_id, sisub$sample_id),]$ct) %>%
    mutate(in_range = if_all(minr_logpval:EC50_logpval, ~ .x > 1.3)) %>%
    group_by(treatment_id, ct) %>%
    summarise(frac = sum(in_range, na.rm=T) / n()) %>%
    group_by(treatment_id) %>%
    summarise(frac = max(frac)) %>%
    mutate(frac = ifelse(frac==0, 0.001, frac))
  
  ### weighted average per cpd
  cpd_level <- drug_summary %>%
    mutate(frac = inrange[match(treatment_id, inrange$treatment_id),]$frac) %>%
    group_by(cpd_name, metric) %>%
    summarise(across(starts_with(c("av")), ~weighted.mean(., w = frac, na.rm = T)),
              across(d:p, ~weighted.mean(., w = frac, na.rm = T)),
              across(starts_with("n_"), ~sum(.)),
              datasets = paste0(unique(dataset), collapse=";")) %>%
    group_by(metric) %>%
    filter(!is.na(r)) %>%
    mutate(prcntl = rank(r, na.last = F)/length(r), .before=1) %>%
    arrange(-prcntl)
  cpd_level <- merge(cpd_level, di, by="cpd_name")
  
  ########### combine drug screen + gene dep ############
  tmp <- cpd_level %>%
    separate_rows(target_genes, sep=", ") %>%
    mutate(gene_symbol = sub("^ (.+)$", "\\1", target_genes)) %>%
    distinct()
  combined <- merge(dep_summary, tmp, by ="gene_symbol", all=T)
  
  #### gene lvl ####
  tmp1 <- combined %>% 
    filter(!is.na(entrez_id)) %>%
    group_by(entrez_id, metric) %>%
    summarise(gene_symbol = paste(unique(gene_symbol), collapse=";"),
              gene_name = paste(unique(gene_name), collapse=";"),
              compounds = paste(cpd_name, collapse=";"),
              CRISPR_score = unique(CRISPR_prcntl),
              RNAi_score = unique(RNAi_prcntl),
              cpd_score = mean(prcntl, na.rm=T),
              rCRISPR = unique(rCRISPR),
              dCRISPR = unique(dCRISPR),
              pCRISPR = unique(pCRISPR),
              rRNAi = unique(rRNAi),
              dRNAi = unique(dRNAi),
              pRNAi = unique(pRNAi)) %>%
    rowwise() %>% 
    mutate(gene_score = mean(c(CRISPR_score, cpd_score), na.rm=T), .before=6) %>%
    filter(!is.na(CRISPR_score) | !is.na(RNAi_score))
  tmp2 <- combined %>% 
    group_by(entrez_id, metric) %>%
    summarise(across(d:p, ~mean(., na.rm=T)),
              across(starts_with("av_", ignore.case = F), ~mean(., na.rm=T))) 
  genelvl <- merge(tmp1, tmp2, by=c("entrez_id", "metric")) %>%
    arrange(-gene_score)
  
  #### cpd lvl ####
  tmp1 <- combined %>% 
    dplyr::filter(!is.na(cpd_name)) %>%
    group_by(cpd_name, metric) %>%
    summarise(pubchem_cid = unique(pubchem_cid),
              target_genes = paste(unique(gene_symbol), collapse=";"),
              pathways = paste(unique(pathways), collapse=";"),
              datasets = paste(unique(datasets), collapse=";"),
              metric_score = unique(prcntl),
              CRISPR_score = mean(CRISPR_prcntl, na.rm=T),
              RNAi_score = mean(RNAi_prcntl, na.rm=T)) %>%
    rowwise() %>% 
    mutate(cpd_score = mean(c(metric_score, CRISPR_score), na.rm=T), .before=3) %>%
    filter(!is.na(cpd_score))
  tmp2 <- combined %>% 
    group_by(cpd_name, metric) %>%
    summarise(across(d:p, ~mean(., na.rm=T))) 
  cpdlvl <- merge(tmp1, tmp2, by=c("cpd_name","metric")) %>%
    arrange(-cpd_score)
  
  # save
  ctres <- list(si = sisub,
                ct1 = ct1,
                ct2 = ct2,
                genelvl = genelvl,
                cpdlvl = cpdlvl)
  saveRDS(ctres, paste0("data/ctres2/", ct1, "_vs_", ct2, "_diffsens.rds"))
})

