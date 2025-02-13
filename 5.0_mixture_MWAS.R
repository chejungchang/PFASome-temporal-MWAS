# Project PFAS-ome pilot, metabolomics
# MWAS 

# Library
library(tidyverse) # data cleaning and wrangling 
library(qgcomp)
library(tictoc)
library(foreach)
library(doParallel)


mix <- read_csv("clean_data/PFAS_ome_upload/sourcedata_figure5.csv")
ft_imp <- read_csv("clean_data/PFAS_ome_upload/ft_untar_tarPFAS_SS_imp.csv")
PFAS_list <- readRDS("clean_data/PFAS_ome_upload/List_PFASome_untar_ftID.RData")

nalist <- ft_imp %>% filter(is.na(tar.8.2FTS)) %>% select(PSID) %>% pull()
ft_imp1 <- ft_imp %>% filter(!PSID %in% nalist) %>% 
  select(starts_with("ft"), starts_with("tar"),PSID, timepoint, SS_race, SS_edu, SS_age_T0, SS_age_T1) %>%
  pivot_longer(cols = starts_with("SS_age_"), names_to = "age_type", values_to = "SS_age") %>%
  mutate(age_type = if_else(age_type == "SS_age_T0", "T0", "T1")) %>%
  filter(timepoint == age_type) 

grp1 <- mix %>% filter(membership == 1) %>% select(chemical) %>% pull()
# grp2 <- mix %>% filter(membership == 2) %>% select(chemical) %>% pull()
# grp3 <- mix %>% filter(membership == 3) %>% select(chemical) %>% pull()
# grp4 <- mix %>% filter(membership == 4) %>% select(chemical) %>% pull()
# grp5 <- mix %>% filter(membership == 5) %>% select(chemical) %>% pull()
# grp6 <- mix %>% filter(membership == 6) %>% select(chemical) %>% pull()
# grp7 <- mix %>% filter(membership == 7) %>% select(chemical) %>% pull()
# grp8 <- mix %>% filter(membership == 8) %>% select(chemical) %>% pull()
# grp9 <- mix %>% filter(membership == 9) %>% select(chemical) %>% pull()
# grp10 <- mix %>% filter(membership == 10) %>% select(chemical) %>% pull()
# grp11 <- mix %>% filter(membership == 11) %>% select(chemical) %>% pull()
# grp12 <- mix %>% filter(membership == 12) %>% select(chemical) %>% pull()
# grp13 <- mix %>% filter(membership == 13) %>% select(chemical) %>% pull()

grp1.str <- paste(grp1, collapse=" + ")
# grp2.str <- paste(grp2, collapse=" + ")
# grp3.str <- paste(grp3, collapse=" + ")
# grp4.str <- paste(grp4, collapse=" + ")
# grp5.str <- paste(grp5, collapse=" + ")
# grp6.str <- paste(grp6, collapse=" + ")
# grp7.str <- paste(grp7, collapse=" + ")
# grp8.str <- paste(grp8, collapse=" + ")
# grp9.str <- paste(grp9, collapse=" + ")
# grp10.str <- paste(grp10, collapse=" + ")
# grp11.str <- paste(grp11, collapse=" + ")
# grp12.str <- paste(grp12, collapse=" + ")
# grp13.str <- paste(grp13, collapse=" + ")


meta_C18_pos <- ft_imp %>% select(!any_of(PFAS_list) & matches("C18pos")) %>% colnames() 
meta_C18_neg <- ft_imp %>% select(!any_of(PFAS_list) & matches("C18neg")) %>% colnames() 
meta_HILIC_pos <- ft_imp %>% select(!any_of(PFAS_list) & matches("HILICpos")) %>% colnames() 
meta_HILIC_neg <- ft_imp %>% select(!any_of(PFAS_list) & matches("HILICneg")) %>% colnames() 


# C18 pos ----
res_C18_pos <- foreach(i = 1:length(meta_C18_pos), .packages = c("qgcomp", "tidyverse"), .combine = bind_rows, .errorhandling = "remove") %dopar% { # do is sequential execution and dopar is parallel 
  foru <- as.formula(paste0("log(", meta_C18_pos[i],") ~", grp3.str, "+ factor(SS_race) + factor(SS_edu) + SS_age"))
  qgcomp <- qgcomp.glm.ee(f = foru, expnms = grp3, q = 4, id = "PSID", data = ft_imp1, family = gaussian())
  mat <- tibble(ft = meta_C18_pos[i], 
                mz = gsub("^ft\\.(.*?)_.*$", "\\1", meta_C18_pos[i]),
                rt = gsub("^.*?_([^_]+)__.*$", "\\1", meta_C18_pos[i]),
                est = qgcomp$psi,
                pv = qgcomp$pval[2], 
                stat = qgcomp$tstat[2])
  return(mat)
} 

write_csv(res_C18_pos, "clean_data/PFAS_ome_upload/res_meta_C18_pos_grp1.csv" )


# C18 neg ----

res_C18_neg <- foreach(i = 1:length(meta_C18_neg), .packages = c("qgcomp", "tidyverse"), .combine = bind_rows, .errorhandling = "remove") %dopar% { # do is sequential execution and dopar is parallel 
  foru <- as.formula(paste0("log(", meta_C18_neg[i],") ~", grp3.str, "+ factor(SS_race) + factor(SS_edu) + SS_age"))
  qgcomp <- qgcomp.glm.ee(f = foru, expnms = grp3, q = 4, id = "PSID", data = ft_imp1, family = gaussian())
  mat <- tibble(ft = meta_C18_neg[i], 
                mz = gsub("^ft\\.(.*?)_.*$", "\\1", meta_C18_neg[i]),
                rt = gsub("^.*?_([^_]+)__.*$", "\\1", meta_C18_neg[i]),
                est = qgcomp$psi,
                pv = qgcomp$pval[2], 
                stat = qgcomp$tstat[2])
  return(mat)
} 

write_csv(res_C18_neg, "clean_data/PFAS_ome_upload/res_meta_C18_neg_grp1.csv" )


# HILIC pos ----
res_HILIC_pos <- foreach(i = 1:length(meta_HILIC_pos), .packages = c("qgcomp", "tidyverse"), .combine = bind_rows, .errorhandling = "remove") %dopar% { # do is sequential execution and dopar is parallel 
  foru <- as.formula(paste0("log(", meta_HILIC_pos[i],") ~", grp3.str, "+ factor(SS_race) + factor(SS_edu) + SS_age"))
  qgcomp <- qgcomp.glm.ee(f = foru, expnms = grp3, q = 4, id = "PSID", data = ft_imp1, family = gaussian())
  mat <- tibble(ft = meta_HILIC_pos[i], 
                mz = gsub("^ft\\.(.*?)_.*$", "\\1", meta_HILIC_pos[i]),
                rt = gsub("^.*?_([^_]+)__.*$", "\\1", meta_HILIC_pos[i]),
                est = qgcomp$psi,
                pv = qgcomp$pval[2], 
                stat = qgcomp$tstat[2])
  return(mat)
} 

write_csv(res_HILIC_pos, "clean_data/PFAS_ome_upload/res_meta_HILIC_pos_grp1.csv" )


# HILIC neg ----

res_HILIC_neg <- foreach(i = 1:length(meta_HILIC_neg), .packages = c("qgcomp", "tidyverse"), .combine = bind_rows, .errorhandling = "remove") %dopar% { # do is sequential execution and dopar is parallel 
  foru <- as.formula(paste0("log(", meta_HILIC_neg[i],") ~", grp3.str, "+ factor(SS_race) + factor(SS_edu) + SS_age"))
  qgcomp <- qgcomp.glm.ee(f = foru, expnms = grp3, q = 4, id = "PSID", data = ft_imp1, family = gaussian())
  mat <- tibble(ft = meta_HILIC_neg[i], 
                mz = gsub("^ft\\.(.*?)_.*$", "\\1", meta_HILIC_neg[i]),
                rt = gsub("^.*?_([^_]+)__.*$", "\\1", meta_HILIC_neg[i]),
                est = qgcomp$psi,
                pv = qgcomp$pval[2], 
                stat = qgcomp$tstat[2])
  return(mat)
} 

write_csv(res_HILIC_neg, "clean_data/PFAS_ome_upload/res_meta_HILIC_neg_grp1.csv" )


