# load data
load("data/depdata.rda")

# run t-test only if sufficient observations
filt_ttest <- function(x, metric, stat = "statistic", grouping_var = "ct") {
  x <- x[, c(grouping_var, metric)] %>% na.omit()
  if(any(table(x[,grouping_var]) < 2) | length(unique(x[,metric,drop=T])) < 3) return(NA)
  tres <- t.test(reformulate(grouping_var, metric), x)
  return(tres[[stat]])
}

### tmp - quick fix for data issues
data$si <- data$si %>% 
  dplyr::rename(COSMIC_ID = "COSMIC.identifier") %>%
  filter(!ds_type %in% c("Testicular", "Embryonal","NOS","Unknown",
                         "Engineered","Other","Eye","Gallbladder"))
data$gdsc_metrics <- data$gdsc_metrics %>%
  mutate(TARGET = ifelse(is.na(TARGET), "", TARGET))

# for testing
reactvals <- list()
reactvals$ct1_si <- data$si[grep("MCL", data$si$ds_subtype),]
reactvals$ct2_si <- data$si[which(data$si$ds_type == "Colorectal"),]
reactvals$ct1 <- "Mantle cell lymphoma (MCL)"
reactvals$ct2 <- "Colorectal"

# CTD metrics
ct1 <- data$ctd_metrics %>%
  dplyr::filter(!is.na(sampleid) &
                  sampleid %in% reactvals$ct1_si$sampleid) %>%
  mutate(ct = reactvals$ct1)
ct2 <- data$ctd_metrics %>%
  dplyr::filter(!is.na(sampleid) & 
                  sampleid %in% reactvals$ct2_si$sampleid) %>%
  mutate(ct = reactvals$ct2)
print(head(ct1)); print(head(ct2))
if (nrow(ct1) < 3 | nrow(ct2) < 3) return(NULL)
ctd_metrics <- rbind(ct1, ct2)
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





ctd_metrics %>%
  filter(treatmentid %in% c("AZD-7545","1S,3R-RSL-3")) %>%
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
         pDSS3 = map_dbl(data, ~filt_ttest(.x, metric="zDSS3", stat="p.value")))

ctd_metrics %>%
  filter(treatmentid %in% c("AZD-7545","1S,3R-RSL-3")) %>%
  select(treatmentid, genesymbol, ct, DSS1, DSS2, DSS3, EC50) %>%
  group_by(ct) %>%
  mutate(across(DSS1:DSS3, ~scale(.)[,1], .names = "z{.col}")) %>%
  ungroup() %>%
  mutate(ct = factor(ct, levels=c(reactvals$ct1, reactvals$ct2))) %>%
  nest(data = -c(treatmentid, genesymbol)) %>%
  mutate(dDSS3 = map_dbl(data, ~t.test(.x)),
         pDSS3 = map_dbl(data, ~filt_ttest(.x, metric="zDSS3", stat="p.value")))

select_if(is.numeric) %>%
  map_df(~ broom::tidy(t.test(. ~ grp)), .id = 'var')  