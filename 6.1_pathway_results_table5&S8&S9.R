library(tidyverse)
library(metapone)

# Load data ---- 
for (i in 1:13) {
  file_name <- paste0("HPC_data/pathway_grp", i, ".rds")
  assign(paste0("grp", i), readRDS(file_name))
}

# Combine data

for (i in 1:13) {
 grp_object <- get(paste0("grp", i))
  
   tmp_list <- ptable(grp_object)[
    which(
      ptable(grp_object)[,"p_value"] < 0.05 & 
        ptable(grp_object)[,"n_significant_metabolites"] >= 2
    ),
  ] %>%
    mutate(mixture = paste0("M", i))
  
  assign(paste0("list", i), tmp_list)
  
  tmp_data <- data.frame(pathway = rownames(tmp_list))
  tmp_data[[paste0("M", i)]] <- 1
  
  assign(paste0("data", i), tmp_data)
}


# Combine all the significant pathway
library(purrr)
data_list <- mget(paste0("data", 1:13))
data_combine <- reduce(data_list, full_join, by = join_by(pathway))

write_csv(data_combine, "pathway_cat.csv")

list1_2<- ptable(grp1)
data_list <- mget(paste0("list", 1:13))
combined_df <- bind_rows(data_list)
combined_df2 <- combined_df %>% mutate(pathway = rep(rownames(list1_2), 2))
combined_df3 <- combined_df2 %>% select(pathway, p_value, mixture, n_significant_metabolites)
write_csv(combined_df3, "ind_pathway_cat.csv")

