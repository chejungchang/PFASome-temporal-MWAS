# Project PFAS-ome pilot, metabolomics
# pathway analysis

# Load packages  
library(metapone)
library(tidyverse)
library(ggplot2)

# This is the pathway data with only human metabolic pathways 
pa <- readRDS("clean_data/pa_human.rds") %>% filter(flag==1) %>% filter(KEGG.ID!="") 
data(hmdbCompMZ)
pos.adductsubset = c("M+","M+H", "M+Na", "M-H2O+H","2M+H","M+2H","M+ACN+H","M+ACN+2H","M+ACN+Na") # all adducts in RefMet/PubChem/Norman annotations. (same as xMSannotator list)
neg.adductsubset = c("M-H", "M+Cl", "M-2H+Na", "2M-H", "M+Hac-H", "M-H2O-H") # all adducts in RefMet/PubChem/Norman annotations 

# Load data ----
C18_pos_grp1 <- read_csv("HPC_data/res_meta_C18_pos_grp1.csv")
C18_neg_grp1 <- read_csv("HPC_data/res_meta_C18_neg_grp1.csv")
HILIC_pos_grp1 <- read_csv("HPC_data/res_meta_HILIC_pos_grp1.csv")
HILIC_neg_grp1 <- read_csv("HPC_data/res_meta_HILIC_neg_grp1.csv")

# Only save those with detection frequency >=60 in the MWAS analyses 
dfT0 <- readRDS("clean_data/PFAS_ome_upload/List_over60dfT0_ftID.RData")
dfT1 <- readRDS("clean_data/PFAS_ome_upload/List_over60dfT1_ftID.RData")

meta_df60 <- intersect(dfT0, dfT1)
meta_df60 <- str_sub(meta_df60, end=-11)


# Format the dataset to metapone ---- 
C18_pos_grp1 <- C18_pos_grp1 %>% filter(ft %in% meta_df60) %>% mutate(adj.pv=p.adjust(pv, method="BH")) %>% select(-1, -est, -pv) %>% select("mz", "rt", "adj.pv", "stat") %>% as.data.frame()
C18_neg_grp1 <- C18_neg_grp1 %>% filter(ft %in% meta_df60) %>% mutate(adj.pv=p.adjust(pv, method="BH")) %>% select(-1, -est, -pv) %>% select("mz", "rt", "adj.pv", "stat") %>% as.data.frame()
HILIC_pos_grp1 <- HILIC_pos_grp1 %>% filter(ft %in% meta_df60) %>% mutate(adj.pv=p.adjust(pv, method="BH")) %>% select(-1, -est, -pv) %>% select("mz", "rt", "adj.pv", "stat") %>% as.data.frame()
HILIC_neg_grp1 <- HILIC_neg_grp1 %>%  filter(ft %in% meta_df60) %>% mutate(adj.pv=p.adjust(pv, method="BH")) %>% select(-1, -est, -pv) %>% select("mz", "rt", "adj.pv", "stat") %>% as.data.frame()

colnames(C18_pos_grp1) <- c("mz", "time", "p-value", "F-statistics")
colnames(C18_neg_grp1) <- c("mz", "time", "p-value", "F-statistics")
colnames(HILIC_pos_grp1) <- c("mz", "time", "p-value", "F-statistics")
colnames(HILIC_neg_grp1) <- c("mz", "time", "p-value", "F-statistics")


dat <- list(C18_pos_grp1, HILIC_pos_grp1, C18_neg_grp1, HILIC_neg_grp1)
type <- list("pos", "pos","neg", "neg")

# permutation test

# permutation test
r2 <- metapone(dat, type, pa, hmdbCompMZ=hmdbCompMZ, pos.adductlist=pos.adductsubset, neg.adductlist=neg.adductsubset,  
               p.threshold=0.05,n.permu=2000, fractional.count.power=0.5, max.match.count=10, match.tol.ppm=5, use.fgsea=FALSE, use.meta=FALSE) 
saveRDS(r2, file = "clean_data/PFAS_ome_upload/pathway_grp1.RDS") 


