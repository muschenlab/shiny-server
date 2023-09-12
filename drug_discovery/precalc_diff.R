library(org.Hs.eg.db)
library(tidyverse)
library(rstatix)
library(coin)
library(RMariaDB)

## gene info
orgdb <- org.Hs.eg.db
keys <- AnnotationDbi::keys(orgdb, keytype='ENTREZID')
gene_info <- AnnotationDbi::select(orgdb, keys = keys, 
                                   columns = c("SYMBOL", "GENENAME")) %>%
  dplyr::rename(entrez_id = ENTREZID,
         gene_symbol = SYMBOL,
         gene_name = GENENAME)

# # run t-test only if sufficient observations
# filt_ttest <- function(x, metric, stat = "statistic", grouping_var = "ct") {
#   grpn <- table(x[[grouping_var]])
#   if(any(grpn < 2)) return(NA)
#   tres <- t.test(reformulate(grouping_var, metric), x)
#   return(tres[[stat]])
# }
# filt_wilcox_p <- function(x, metric, grouping_var = "ct") {
#   grpn <- table(x[[grouping_var]])
#   if(any(grpn < 2)) return(NA)
#   statres <- rstatix::wilcox_test(x, reformulate(grouping_var, metric))
#   return(statres$p)
# }
# filt_wilcox_eff <- function(x, metric, grouping_var = "ct") {
#   grpn <- table(x[[grouping_var]])
#   if(any(grpn < 2)) return(NA)
#   effres <- rstatix::wilcox_effsize(x, reformulate(grouping_var, metric))
#   return(effres$effsize)
# }
filt_wilcox <- function(x, metric, grouping_var = "ct", id_var = "entrez_id") {
  grpn <- table(x[[grouping_var]])
  if(any(grpn < 2)) return(NULL)
  statres <- rstatix::wilcox_test(x, reformulate(grouping_var, metric), 
                                  ref.group = levels(x[[grouping_var]])[[2]])
  effres <- rstatix::wilcox_effsize(x, reformulate(grouping_var, metric), 
                                    ref.group = levels(x[[grouping_var]])[[2]])
  statres <- merge(statres, effres, by = c(".y.","group1","group2","n1","n2")) %>%
    mutate(id = unique(x[[id_var]]), .before=1)
  return(statres)
}

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
                                                        "Astrocytoma")),])
cts <- unique(si$ds_type); names(cts) <- cts
tmp <- lapply(cts, function(ct) si[which(si$ds_type == ct),])
silist <- c(silist, tmp)
# remove low n
silist <- silist[!names(silist) %in% c("Testicular", "Embryonal", "Eye",
                                       "Gallbladder", "Skin carcinoma")]

### set comaprisons ###
comp <- data.frame(a = c("B-cell","T-cell", "T-cell leukemia", names(silist)[-1]),
                   b = c("Myeloid","B-cell", "B-cell leukemia", rep("Solid tumor", length(silist)-1)))

#################### calc stats for each disease group ###########################
lapply(1:nrow(comp), function(compidx) {
  
  ct1 <- comp[compidx, ]$a
  ct2 <- comp[compidx, ]$b
  print(paste0(ct1, " vs ", ct2))
  if (file.exists(paste0("data/ctres/", ct1, "_vs_", ct2, "_diffsens.rds"))) return(NULL)
  
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
    filter(!is.na(DSS4)) %>%
    mutate(limEC50 = ifelse(logEC50 > logmaxc, log10(2*10^logmaxc), logEC50))
  av <- ctd_metrics %>%
    group_by(ct, treatment_id) %>%
    summarise(n = n(),
              avEC50 = mean(limEC50, na.rm=T),
              avDSS4 = mean(DSS4, na.rm=T)) %>%
    pivot_wider(names_from = ct,
                values_from = c(avEC50, avDSS4, n))
  nfilt <- filter_at(av, vars(starts_with("n")), all_vars(. > 3))
  if (nrow(nfilt) < 1) {
    ctd_summary <- av %>%
      mutate(dDSS4 = NA,
             rDSS4 = NA,
             pDSS4 = NA,
             DSS4_prcntl = NA)
  } else {
    stats <- ctd_metrics %>%
      filter(treatment_id %in% nfilt$treatment_id) %>%
      mutate(zDSS4 = scale(DSS4)[,1]) %>%
      group_by(treatment_id) %>% 
      filter(!all(DSS4==0)) %>%
      rstatix::wilcox_test(zDSS4 ~ ct, ref.group = ct1, detailed = T)
    eff <- ctd_metrics %>%
      filter(treatment_id %in% nfilt$treatment_id) %>%
      mutate(zDSS4 = scale(DSS4)[,1]) %>%
      group_by(treatment_id) %>% 
      filter(!all(DSS4==0)) %>%
      rstatix::wilcox_effsize(zDSS4 ~ ct, ref.group = ct1)
    stats <- merge(stats, eff[,4:5], by = c("treatment_id")) %>%
      dplyr::select(treatment_id, estimate, effsize, p) %>%
      mutate(effsize = ifelse(estimate < 0, -effsize, effsize)) %>%
      dplyr::rename(dDSS4 = estimate,
             rDSS4 = effsize,
             pDSS4 = p) %>%
      mutate(DSS4_prcntl = rank(-rDSS4)/length(rDSS4)) 
    ctd_summary <- merge(av, stats, by = "treatment_id", all = T) %>%
      arrange(-DSS4_prcntl) 
  }
  
  ######## GDSC ########
  gdsc_metrics <- drug_metrics %>%
    filter(dataset %in% c("GDSC1", "GDSC2")) %>%
    mutate(ct = sisub[match(sample_id, sisub$sample_id),]$ct) %>%
    mutate(ct = factor(ct, levels=c(ct1, ct2))) %>%
    filter(!is.na(DSS4)) %>%
    mutate(limEC50 = ifelse(logEC50 > logmaxc, log10(2*10^logmaxc), logEC50))
  av <- gdsc_metrics %>%
    group_by(ct, treatment_id) %>%
    summarise(n = n(),
              avEC50 = mean(limEC50, na.rm=T),
              avDSS4 = mean(DSS4, na.rm=T)) %>%
    pivot_wider(names_from = ct,
                values_from = c(avEC50, avDSS4, n))
  nfilt <- filter_at(av, vars(starts_with("n")), all_vars(. > 3))
  if (nrow(nfilt) < 1) {
    gdsc_summary <- av %>%
      mutate(dDSS4 = NA,
             rDSS4 = NA,
             pDSS4 = NA,
             DSS4_prcntl = NA)
  } else {
    stats <- gdsc_metrics %>%
      filter(treatment_id %in% nfilt$treatment_id) %>%
      mutate(zDSS4 = scale(DSS4)[,1]) %>%
      group_by(treatment_id) %>% 
      filter(!all(DSS4==0)) %>%
      rstatix::wilcox_test(zDSS4 ~ ct, ref.group = ct1, detailed = T)
    eff <- gdsc_metrics %>%
      filter(treatment_id %in% nfilt$treatment_id) %>%
      mutate(zDSS4 = scale(DSS4)[,1]) %>%
      group_by(treatment_id) %>% 
      filter(!all(DSS4==0)) %>%
      rstatix::wilcox_effsize(zDSS4 ~ ct, ref.group = ct1)
    stats <- merge(stats, eff[,4:5], by = c("treatment_id")) %>%
      dplyr::select(treatment_id, estimate, effsize, p) %>%
      mutate(effsize = ifelse(estimate < 0, -effsize, effsize)) %>%
      dplyr::rename(dDSS4 = estimate,
             rDSS4 = effsize,
             pDSS4 = p) %>%
      mutate(DSS4_prcntl = rank(-rDSS4)/length(rDSS4)) 
    gdsc_summary <- merge(av, stats, by = "treatment_id", all = T) %>%
      arrange(-DSS4_prcntl) 
  }
  
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
    group_by(cpd_name) %>%
    summarise(across(starts_with(c("av","dD","rD","pD")), 
                     ~weighted.mean(., w = frac, na.rm = T)),
              across(starts_with("n_"), ~sum(.)),
              datasets = paste0(dataset, collapse=";")) %>%
    ungroup() %>%
    filter(!is.na(rDSS4)) %>%
    mutate(DSS4_prcntl = rank(rDSS4, na.last = F)/length(rDSS4), .before=1) %>%
    arrange(-DSS4_prcntl)
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
    group_by(entrez_id) %>%
    summarise(gene_symbol = paste(unique(gene_symbol), collapse=";"),
              gene_name = paste(unique(gene_name), collapse=";"),
              compounds = paste(cpd_name, collapse=";"),
              CRISPR_score = unique(CRISPR_prcntl),
              RNAi_score = unique(RNAi_prcntl),
              DSS4_score = mean(DSS4_prcntl, na.rm=T)) %>%
    mutate(gene_score = rowMeans(.[,c(5,7)], na.rm=T), .before=4) %>%
    filter(!is.na(CRISPR_score) | !is.na(RNAi_score))
  tmp2 <- combined %>% 
    group_by(entrez_id) %>%
    summarise(across(starts_with("d", ignore.case = F), ~mean(., na.rm=T)),
              across(starts_with("r", ignore.case = F), ~mean(., na.rm=T)),
              across(starts_with(c("pC","pR","pD"), ignore.case = F), ~mean(., na.rm=T)),
              across(starts_with("av", ignore.case = F), ~mean(., na.rm=T))) 
  genelvl <- merge(tmp1, tmp2, by=c("entrez_id")) %>%
    arrange(-gene_score)
  
  #### cpd lvl ####
  tmp1 <- combined %>% 
    group_by(cpd_name) %>%
    summarise(pubchem_cid = unique(pubchem_cid),
              target_genes = paste(gene_symbol, collapse=";"),
              pathways = unique(pathways),
              datasets = paste(unique(datasets), collapse=";"),
              DSS4_score = unique(DSS4_prcntl),
              CRISPR_score = mean(CRISPR_prcntl, na.rm=T),
              RNAi_score = mean(RNAi_prcntl, na.rm=T)) %>%
    rowwise() %>% 
    mutate(cpd_score = mean(c(DSS4_score, CRISPR_score), na.rm=T), .before=3) %>%
    filter(!is.na(cpd_score))
  tmp2 <- combined %>% 
    group_by(cpd_name) %>%
    summarise(across(starts_with(c("dD","dC","dR"), ignore.case = F), ~mean(., na.rm=T)),
              across(starts_with(c("rD","rC","rR"), ignore.case = F), ~mean(., na.rm=T)),
              across(starts_with(c("pD","pC","pR"), ignore.case = F), ~mean(., na.rm=T)),
              across(starts_with(c("avD","avC","avR"), ignore.case = F), ~mean(., na.rm=T))) 
  cpdlvl <- merge(tmp1, tmp2, by=c("cpd_name")) %>%
    arrange(-cpd_score)
  
  # save
  ctres <- list(si = sisub,
                ct1 = ct1,
                ct2 = ct2,
                genelvl = genelvl,
                cpdlvl = cpdlvl)
  saveRDS(ctres, paste0("data/ctres/", ct1, "_vs_", ct2, "_diffsens.rds"))
})

