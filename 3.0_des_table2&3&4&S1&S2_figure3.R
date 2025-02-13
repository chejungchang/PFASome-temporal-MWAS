# Project PFAS-ome pilot, metabolomics
# Concentrations and detection frequency


library(tidyverse)
library(dplyr)
library(xlsx)
library(Hmisc)
library(irr)
library(ggplot2)
library(ggpubr)

# Detection frequency for metabolomics -----

dat1 <- read_csv("clean_data/PFAS_ome_upload/ft_untar_tarPFAS_SS.csv") 
nalist <- dat1  %>% filter(is.na(tar.8.2FTS)) %>% select(PSID) %>% pull()


df1 <- dat1 |> 
  filter(!PSID %in% nalist) |> 
  filter(timepoint == "T0") |>
  summarise(
    across(
      .cols = everything(),
      .fns = ~sum(. > 0, na.rm = TRUE) / n() * 100,
     #.fns = ~sum(. != "N/F" & . > 0, na.rm = TRUE) / n() * 100, # those N/F or < 0 are considered non-detected
      .names = "{col}_detection"
    )
  ) 


df1_list <- df1[which(df1 >= 60)]
df1_list2 <- df1_list %>% colnames()
saveRDS(df1_list2, "clean_data/PFAS_ome_upload/List_over60dfT0_ftID.RData") 

df2 <- dat1 |> 
  filter(!PSID %in% nalist) |> 
  filter(timepoint == "T1") |>
  summarise(
    across(
      .cols = everything(),
      .fns = ~sum(. > 0, na.rm = TRUE) / n() * 100,
      #.fns = ~sum(. != "N/F" & . > 0, na.rm = TRUE) / n() * 100, # those N/F or < 0 are considered non-detected
      .names = "{col}_detection"
    )
  ) 

hist(t(df2))

df2_list <- df2[which(df2 >= 60)]
df2_list2 <- df2_list %>% colnames()

saveRDS(df2_list2, "clean_data/PFAS_ome_upload/List_over60dfT1_ftID.RData") 

# df1_1 <- as.data.frame(t(df1)) %>% rownames_to_column()
# df2_1 <- as.data.frame(t(df2)) %>% rownames_to_column()
# 
# write_csv(df1_1, "results/PFAS_ome/meta_detection_frequency_T0.csv")
# write_csv(df2_1, "results/PFAS_ome/meta_detection_frequency_T1.csv")

# Descriptive stats for PFAS -----


dat1 <- read_csv("clean_data/PFAS_ome_upload/PFASome_SS.csv") 
datimp1 <- read_csv("clean_data/PFAS_ome_upload/PFASome_SS_imp.csv") 
nalist <- dat1  %>% filter(is.na(tar.8.2FTS)) %>% select(PSID) %>% pull()


dat1 <- dat1 |> filter(!PSID %in% nalist) 
datimp1 <- datimp1 |> filter(!PSID %in% nalist) 


## Targeted PFAS concentrations ----

# detection frequency 

df <- dat1 %>% 
  group_by(timepoint) %>% 
  select(starts_with("tar")) %>%
  rename_all(~str_sub(., start = 5)) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = ~sum(. > 0, na.rm = TRUE) / n() * 100,
      .names = "{col}_detection"
    )
  ) 



# median, P25, P75, min, max

res <- datimp1 %>% 
  group_by(timepoint) %>% 
  select(starts_with("tar")) %>%
  rename_all(~str_sub(., start = 5)) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = list(min = ~min(., na.rm = TRUE), 
                  max = ~max(., na.rm = TRUE), 
                  median = ~median(., na.rm = TRUE), 
                  P25 = ~quantile(., 0.25, na.rm = TRUE), 
                  P75 = ~quantile(., 0.75, na.rm = TRUE), 
                  IQR = ~IQR(., na.rm = TRUE)), 
      .names = "{col}_{fn}"
    )
  )


sumT0 <- cbind(df[1,], res[1,]) # combine the detection frequency and other summary statistics 
sumT1 <- cbind(df[2,], res[2,]) # combine the detection frequency and other summary statistics 

sumT0 <- sumT0 %>% 
  select(-matches("point")) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = c("chemical", ".value"),
    names_sep = "_"
  ) %>% arrange(factor(chemical, levels = c("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnA", "PFDoA", "PFTrDA", "PFTeDA",
                                            "PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS", "PFNS", "PFDS", 
                                            "PFOSA", "4.2FTS", "6.2FTS", "8.2FTS", "N.MeFOSAA", "N.EtFOSAA")))

sumT1  <- sumT1 %>% 
  select(-matches("point")) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = c("chemical", ".value"),
    names_sep = "_"
  ) %>% arrange(factor(chemical, levels = c("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnA", "PFDoA", "PFTrDA", "PFTeDA",
                                            "PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS", "PFNS", "PFDS", 
                                            "PFOSA", "4.2FTS", "6.2FTS", "8.2FTS", "N.MeFOSAA", "N.EtFOSAA")))

write.xlsx(sumT0, "clean_data/PFAS_ome_upload/Res_des_tarPFAS.xlsx", sheetName = "T0", append = FALSE)
write.xlsx(sumT1, "clean_data/PFAS_ome_upload/Res_des_tarPFAS.xlsx", sheetName = "T1", append = TRUE)

## Untargeted PFAS concentrations ----

# detection frequency 

df <- dat1 %>% 
  group_by(timepoint) %>% 
  select(starts_with("ft")) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = ~sum(. > 0, na.rm = TRUE) / n() * 100, 
      .names = "{col}$detection"
    )
  )


# median, P25, P75, min, max

res <- datimp1 %>% 
  group_by(timepoint) %>% 
  select(starts_with("ft")) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = list(min = ~min(., na.rm = TRUE), 
                  max = ~max(., na.rm = TRUE), 
                  median = ~median(., na.rm = TRUE), 
                  P25 = ~quantile(., 0.25, na.rm = TRUE), 
                  P75 = ~quantile(., 0.75, na.rm = TRUE), 
                  IQR = ~IQR(., na.rm = TRUE)), 
      .names = "{col}${fn}"
    )
  )


sumT0 <- cbind(df[1,], res[1,])
sumT1 <- cbind(df[2,], res[2,])

sumT0 <- sumT0 %>% 
  select(-matches("point")) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = c("chemical", ".value"),
    names_sep = "\\$"
  ) %>% arrange(chemical)


sumT1  <- sumT1 %>% 
  select(-matches("point")) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = c("chemical", ".value"),
    names_sep = "\\$"
  ) %>% arrange(chemical)

write.xlsx(sumT0, "clean_data/PFAS_ome_upload/Res_des_untargetPFAS.xlsx", sheetName = "T0", append = FALSE)
write.xlsx(sumT1, "clean_data/PFAS_ome_upload/Res_des_untargetPFAS.xlsx", sheetName = "T1", append = TRUE)


## Detection frequency for PFASome ----

all_df <- dat1 %>% 
  select(starts_with("tar")|starts_with("ft")) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = ~sum(. != "N/F" & . > 0, na.rm = TRUE) / n() * 100, # those N/F or < 0 are considered non-detected
      .names = "{col}"
    )
  ) 

all_df_list <- all_df[which(all_df >= 60)]
all_df_list2 <- all_df_list %>% colnames()
saveRDS(all_df_list2, "clean_data/PFAS_ome_upload/List_over60df_PFASome.RData") 


# Based on only T0 
T0_df <- dat1 %>%
  filter(timepoint == "T0") %>%
  select(starts_with("tar")|starts_with("ft")) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = ~sum(. != "N/F" & . > 0, na.rm = TRUE) / n() * 100, # those N/F or < 0 are considered non-detected
      .names = "{col}"
    )
  ) 


T0_df_list <- T0_df[which(T0_df >= 60)]
T0_df_list2 <- T0_df_list %>% colnames()

saveRDS(T0_df_list2, "clean_data/PFAS_ome_upload/List_over60dfT0_PFASome.RData") 


# Based on only T1 
T1_df <- dat1 %>%
  filter(timepoint == "T1") %>%
  select(starts_with("tar")|starts_with("ft")) %>%
  summarise(
    across(
      .cols = everything(),
      .fns = ~sum(. != "N/F" & . > 0, na.rm = TRUE) / n() * 100, # those N/F or < 0 are considered non-detected
      .names = "{col}"
    )
  ) 


T1_df_list <- T1_df[which(T1_df >= 60)]
T1_df_list2 <- T1_df_list %>% colnames()

saveRDS(T1_df_list2, "clean_data/PFAS_ome_upload/List_over60dfT1_PFASome.RData") 


# Correlations and ICCs ----
datimp2 <- datimp1 %>% 
  select(PSID, starts_with("ft"), starts_with("tar"), timepoint) %>% 
  pivot_wider(names_from = timepoint,
              values_from = c(starts_with("ft"), starts_with("tar")),
              names_sep = "___") %>%
  mutate(across(-PSID, as.numeric))

cormat <- rcorr(as.matrix(datimp2[,-1]),  type="spearman")

r <- cormat[[1]] %>% data.frame() %>% 
  select(ends_with("T0")) %>%
  rownames_to_column("rowname") %>%
  filter(str_detect(rowname, "___T1$")) 

p <- cormat[[3]] %>% data.frame() %>% 
  select(ends_with("T0")) %>%
  rownames_to_column("rowname") %>%
  filter(str_detect(rowname, "___T1$")) 

cor_res <- rbind(chemical = colnames(r)[-1], r = diag(as.matrix(r[,-1])), p = diag(as.matrix(p[,-1]))) %>% t() %>% as_tibble
cor_res <- cor_res %>% mutate(chemical = str_sub(chemical, start =1, end=-6))

untar_cor <- cor_res %>% filter(str_detect(chemical,"ft")) %>% arrange(chemical)
tar_cor <- cor_res %>% filter(str_detect(chemical,"tar")) %>% mutate(chemical = str_sub(chemical, start = 5)) %>% arrange(factor(chemical, levels = c("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnA", "PFDoA", "PFTrDA", "PFTeDA",
                                                                                                                                                      "PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS", "PFNS", "PFDS", 
                                                                                                                                                      "PFOSA", "4.2FTS", "6.2FTS", "8.2FTS", "N.MeFOSAA", "N.EtFOSAA")))
library(lme4)
library(performance)

# ICC and standardized difference 
icc_res <- datimp2 %>% 
  pivot_longer(cols = !PSID, 
               names_to = c("chemical", "timepoint"),
               names_sep = "___",
               values_to = "value") %>%
  pivot_wider(names_from = timepoint, 
              values_from = value) %>%
  filter(!is.na(T0)) %>% 
  mutate(diff = T1-T0) %>% 
  group_by(chemical) %>% 
  summarise(median_T0 = median(T0, na.rm = TRUE), 
            median_T1 = median(T1, na.rm = TRUE),
            diff_IQR_std = round(median(diff / IQR(T0, na.rm = TRUE), na.rm = TRUE), digits = 4)) 


untar_icc <- icc_res2 %>% 
  filter(str_detect(chemical,"ft")) %>% 
  arrange(chemical)

tar_icc <- icc_res2 %>% 
  filter(str_detect(chemical,"tar")) %>% 
  mutate(chemical = str_sub(chemical, start = 5)) %>% 
  arrange(factor(chemical, levels = c("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnA", "PFDoA", "PFTrDA", "PFTeDA",
                                      "PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS", "PFNS", "PFDS", "PFOSA", "4.2FTS", "6.2FTS", "8.2FTS", 
                                      "N.MeFOSAA", "N.EtFOSAA")))


untar <- cbind(untar_cor, untar_icc) 
tar <- cbind(tar_cor, tar_icc)

write.xlsx(untar, "clean_data/PFAS_ome_upload/Res_des_untargetPFAS.xlsx", sheetName = "correlation_icc", append = TRUE)
write.xlsx(tar, "clean_data/PFAS_ome_upload/Res_des_tarPFAS.xlsx", sheetName = "correlation_icc", append = TRUE)

# Table 3 -----  
tar %>% select(-4) %>% summarise(r1 = sum(ifelse(abs(as.numeric(r)) >= 0.6, 1, 0)), 
                                 r2 = sum(ifelse(abs(as.numeric(r)) < 0.6 & abs(as.numeric(r)) >= 0.4, 1, 0)),
                                 r3 = sum(ifelse(abs(as.numeric(r)) < 0.4, 1, 0)), 
                                 icc1 = sum(ifelse(abs(as.numeric(ICC_adjusted)) >= 0.4, 1, 0), na.rm = TRUE), 
                                 icc2 = sum(ifelse(abs(as.numeric(ICC_adjusted)) < 0.4 | is.na(ICC_adjusted), 1, 0), na.rm = TRUE)) 
untar %>% select(-4) %>% summarise(r1 = sum(ifelse(abs(as.numeric(r)) >= 0.6, 1, 0)), 
                                   r2 = sum(ifelse(abs(as.numeric(r)) < 0.6 & abs(as.numeric(r)) >= 0.4, 1, 0)),
                                   r3 = sum(ifelse(abs(as.numeric(r)) < 0.4, 1, 0)), 
                                   icc1 = sum(ifelse(abs(as.numeric(ICC_adjusted)) >= 0.4, 1, 0), na.rm = TRUE), 
                                   icc2 = sum(ifelse(abs(as.numeric(ICC_adjusted))< 0.4|is.na(ICC_adjusted) , 1, 0), na.rm = TRUE)) 
# Figure 3 ----
# Plot only those with detection frequency >= 60% 
library(ggExtra)


T0_df <- readRDS("clean_data/PFAS_ome_upload/List_over60dfT0_PFASome.RData")
T1_df <- readRDS("clean_data/PFAS_ome_upload/List_over60dfT1_PFASome.RData")
figuredat <- read_csv("clean_data/PFAS_ome_upload/sourcedata_table3.csv") 

pdf("clean_data/PFAS_ome_upload/Fig3.pdf")

a <- figuredat  %>% 
  filter(chemical %in% T0_df & chemical %in% T1_df) %>% 
  mutate(approach = case_when(str_detect(chemical, "tar") ~ "Targeted", 
                              TRUE ~ "Untargeted")) %>%
  mutate(r = as.numeric(r)) %>%
  ggplot(aes(y=r,x=mean_dif/SD_dif, color = approach, shape = approach)) +
  geom_point()+
  scale_color_manual(values = c("Targeted" = "red", "Untargeted" = "black")) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(x = "Standardized T1-T0", y = "Spearman correlation coefficient")

a1 <- ggMarginal(a, type = "densigram",
                 groupColour = TRUE,
                 groupFill = TRUE)

b <- figuredat %>% 
  filter(chemical %in% T0_df & chemical %in% T1_df) %>% 
  mutate(approach = case_when(str_detect(chemical, "tar") ~ "Targeted", 
                              TRUE ~ "Untargeted")) %>%
  mutate(r = as.numeric(ICC)) %>%
  ggplot(aes(y=r,x=mean_dif/SD_dif, color = approach, shape = approach)) +
  geom_point()+
  scale_color_manual(values = c("Targeted" = "red", "Untargeted" = "black")) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  labs(x = "Standardized T1-T0", y = "Intraclass correlation coefficient")


b1 <- ggMarginal(b, type = "densigram",
                 groupColour = TRUE,
                 groupFill = TRUE)
ggarrange(a1, b1,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

dev.off()


