# Project PFAS-ome pilot, metabolomics
# Correlations between target and untarget metabolomics

library(tidyverse)
library(Hmisc)
library(irr)
library(ggplot2)
library(xlsx)

dat1 <- read_csv("clean_data/PFAS_ome_upload/PFASome_SS.csv") 
datimp1 <- read_csv("clean_data/PFAS_ome_upload/PFASome_SS_imp.csv") 
overlap <- readRDS("clean_data/PFAS_ome_upload/List_overlapPFAS_ftID.RData")
pfas <- read_csv("clean_data/PFAS_ome_upload/Annotation_PFASome.csv")
tarID <- read_csv("clean_data/PFAS_ome_upload/tarPFAS_InChiKey.csv")
T0_df <- readRDS("clean_data/PFAS_ome_upload/List_over60dfT0_PFASome.RData")
T1_df <- readRDS("clean_data/PFAS_ome_upload/List_over60dfT1_PFASome.RData")


dir <- pfas %>% filter(ft_ID %in% overlap) %>% select(ft_ID, InChiKey) %>% inner_join(tarID, join_by(InChiKey))
untar <- datimp1 %>% select(dir$ft_ID) %>% t() %>% as.data.frame %>% rownames_to_column()
tar <- datimp1 %>% select(dir$Abbre) %>% t() %>% as.data.frame %>% rownames_to_column()

v1 <- dir %>% inner_join(untar, join_by(ft_ID == rowname)) %>% mutate(type = "untarget")
v2 <- dir %>% inner_join(tar, join_by(Abbre == rowname)) %>% mutate(type = "target")
v3 <- rbind(v1, v2)

corr <- function(x){
  obj <- rcorr(as.matrix(t(x)),  type="spearman")
  r <- obj$r[1,2]
}


tar.untar.corr <- v3 %>% 
  select(-c(1:3), -type) %>% 
  group_by(Abbre) %>% nest() %>% mutate(corr = map(data, ~corr(.))) %>%
  unnest(corr) %>%
  inner_join(v3[,1:4], join_by("Abbre")) %>%
  inner_join(pfas[,c("mz", "time", "Mode", "ft_ID")], join_by("ft_ID"))%>%
  select(-data, -InChiKey) %>% 
  distinct() %>%
  mutate(Abbre = str_sub(Abbre, start = 5)) 

tar.untar.corr %>% select(Abbre, corr) %>% mutate(corr = round(corr, 2)) %>% print(n=25)


# How are those correlation related to confidence of annotation? 
list <- tar.untar.corr$ft_ID
pfas2 <- pfas %>% filter(ft_ID %in% list)



