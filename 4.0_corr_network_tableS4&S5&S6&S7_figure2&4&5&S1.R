# Project PFAS-ome pilot, metabolomics
# Correlations, network analysis, community detection

# Library 
library(igraph) # visualize network analysis 
library(Hmisc) # correlation
library(tidyverse) # data cleaning and wrangling 
library(readxl) # save excel 
library(RColorBrewer)
library(ggplot2)

# Preparation -----
# load data
dat <- read_csv("clean_data/PFAS_ome_upload/PFASome_SS_imp.csv") 
pfas <- read_csv("clean_data/PFAS_ome_upload/Annotation_PFASome.csv")
overlap <- readRDS("clean_data/PFAS_ome_upload/List_overlapPFAS_ftID.RData")
dfT0 <- readRDS("clean_data/PFAS_ome_upload/List_over60dfT0_PFASome.RData")
dfT1 <- readRDS("clean_data/PFAS_ome_upload/List_over60dfT1_PFASome.RData")

df60 <- intersect(dfT0, dfT1)
keepPFAS <- df60[!df60 %in% LCPFAS_overlap_tar_untar]

nalist <- dat %>% filter(is.na(tar.8.2FTS)) %>% select(PSID) %>% pull()
dat1 <- dat %>% filter(!PSID %in% nalist) %>% select(all_of(keepPFAS), timepoint, PSID) 
dat2 <- dat %>% filter(!PSID %in% nalist) %>% arrange(PSID) %>% select(PSID)

# Network analysis (T0) -----
# use T0 concentrations to create network; because T0 network will be the mixtures that we investigate the in the main study 
T0 <- dat1 %>% filter(timepoint == "T0") %>% arrange(PSID) %>% select(c(starts_with("tar"),  starts_with("ft.")))
T1 <- dat1 %>% filter(timepoint == "T1") %>% arrange(PSID) %>% select(c(starts_with("tar"),  starts_with("ft.")))
cor <- rcorr(as.matrix(T0),  type="spearman") 

# extract the correlation coefficients and p-value; keep those r > 0.6 and p < 0.05  
r <- cor[[1]] %>% data.frame() %>% mutate_all(~ifelse(. >= 0.6, 1, 0))
p <- cor[[3]] %>% data.frame() %>% mutate_all(~ifelse(. < 0.05, 1, 0))

# identify those r > 0.6 AND p < 0.05 
com <- r+p ; 
com <- com %>% mutate_all(~ifelse(. == 2, 1, 0))

# rm PFAS without any correlations 
remove_col <- colnames(com)[colSums(com, na.rm = TRUE) == 0]
remove_row <- which(row.names(com) %in% remove_col)
com1 <- com %>% select(!all_of(remove_col)) %>% slice(-remove_row)

remove_col <- which(!colnames(r1) %in% colnames(com1))
remove_row <- which(!row.names(r1) %in% colnames(com1))

# run network analysis using directed model 
com1 <- as.matrix(com1)
# write_csv(data.frame(com1) , "clean_data/PFAS_ome_upload/sourcedata_figure2.csv")
network1 <- graph_from_adjacency_matrix(adjmatrix=com1, mode = "undirected")
summary(network1)


# identify communities
community1 <- cluster_louvain(network1)
n_cluster <- length(unique(community1$membership))
membershiptable <- data.frame(table(community1$membership))
colnames(membershiptable) <- c("old_mem", "num_PFAS")
membershiptable2 <- membershiptable %>% 
  mutate(num_PFAS = as.numeric(num_PFAS), 
         new_mem = rank(-num_PFAS, ties.method = "random"))

# plot the network 
layout <- layout_nicely(network1)
community1$names
ft_type <- c(rep(1, 15), rep(2, 379-15)) 
newmembership <- data.frame(old_mem = as.character(community1$membership)) %>% right_join(membershiptable2, by="old_mem")
membership <- paste0("M", newmembership$new_mem)
membership2 <- newmembership$new_mem
membership_list <- membership2

# extract information from the network
# compute node
nodes <- V(network1)
nodes <- attr(nodes, "names")

set.seed(4)
plot.igraph(network1, 
     layout = layout, 
     vertex.label = NA,
     mark.groups=communities(community1),
     mark.border=NA,
     margin=0,
     edge.color="darkgray",
     vertex.color=c("gray", "#28282B")[ft_type],
     vertex.label=nodes,
     vertex.color=rainbow(n_cluster, alpha=0.6)[community1$membership],
     vertex.shape=c("square", "circle")[ft_type],
     vertex.size = c(3,2)[ft_type],
     vertex.frame.color="black") # plot without labels 

# with labels
plot.igraph(network1, 
            layout = layout, 
            mark.groups=communities(community1),
            mark.border=NA,
            margin=0.05,
            edge.color="darkgray",
            vertex.color=c("gray", "#28282B")[ft_type],
            vertex.color=rainbow(n_cluster, alpha=0.6)[community1$membership],
            vertex.shape=c("square", "circle")[ft_type],
            vertex.size = c(3,2)[ft_type],
            vertex.frame.color="black",
            vertex.label=membership_list) 


# extract information from the network
# compute node
nodes <- V(network1)
nodes <- attr(nodes, "names")

# compute the degree centrality 
degr_cent <- centr_degree(network1, mode = 'all')
degr_cent <- degr_cent$res

# compute the eigenvector centrality 
eign_cent <- eigen_centrality(network1)
eign_cent <- eign_cent$vector

# compute the closeness centraility
clos_cent <- closeness(network1)

# compute betweeness centrality
betw_cent <- betweenness(network1)

# create data frame storing all of the measures of centrality
net_data1 <- data.frame(vertex = nodes,
                   membership = membership2, 
                   degree = degr_cent, 
                   eigen = eign_cent, 
                   closeness = clos_cent, 
                   betweeness = betw_cent)

net_data2 <- data.frame(membership = membership2, 
                        vertex = nodes,
                        degree = degr_cent, 
                        eigen = eign_cent, 
                        closeness = clos_cent, 
                        betweeness = betw_cent) %>%
  left_join(pfas, join_by(vertex == ft_ID)) %>%
  select(membership, vertex, mz, time, Mode, Level, `Lv1@Metabolite_Name`, `Lv34@Name`, degree, eigen, closeness, betweeness)


corricc <- read_csv("clean_data/PFAS_ome_upload/sourcedata_table3.csv") 

net_summary <- 
  net_data1 %>% 
  inner_join(corricc, join_by(vertex == chemical)) %>%
  group_by(membership) %>%
  summarise(n_PFAS = n(), 
            n_target = sum(str_detect(vertex, "tar.")),
            min_r = min(r),
            max_r = max(r), 
            median_r = median(r),
            min_icc = min(ICC),
            max_icc = max(ICC), 
            median_icc = median(ICC)
            )


large_community <- net_summary %>%  filter(n_PFAS >= 5) %>% select(membership) %>% pull() 
large_community_order <- net_summary %>%  filter(n_PFAS >= 5) %>% arrange(desc(n_PFAS)) %>% select(membership) %>% pull() 


savedat <- net_data1 %>% inner_join(corricc, join_by(vertex == chemical)) %>% filter(membership %in% large_community)
# write_csv(savedat, "clean_data/PFAS_ome_upload/sourcedata_figure5.csv")


a <- net_data1 %>% 
  inner_join(corricc, join_by(vertex == chemical)) %>%
  filter(membership %in% large_community) %>%
  mutate(membership = factor(membership, levels = large_community_order)) %>%
  ggplot() + 
  geom_boxplot(aes(x=as.factor(membership), y=r, fill=as.factor(membership))) +
  theme_bw() +
  labs(x="", y="Spearman correlation coefficient")

b <- net_data1 %>% 
  inner_join(corricc, join_by(vertex == chemical)) %>%
  filter(membership %in% large_community) %>%
  mutate(membership = factor(membership, levels = large_community_order)) %>%
  ggplot() + 
  geom_boxplot(aes(x=as.factor(membership), y=ICC, fill=as.factor(membership))) +
  theme_bw() +
  labs(x="", y="Intraclass correlation coefficient")

c <- net_data1 %>% 
  inner_join(corricc, join_by(vertex == chemical)) %>%
  filter(membership %in% large_community) %>%
  mutate(membership = factor(membership, levels = large_community_order)) %>%
  ggplot() + 
  geom_boxplot(aes(x=as.factor(membership), y=mean_dif/SD_dif, fill=as.factor(membership))) +
  theme_bw() +
  labs(x="", y="Standarized T1-T0")


d <- net_summary %>% 
  filter(membership %in% large_community) %>%
  mutate(membership = factor(membership, levels = large_community_order)) %>%
  ggplot(aes(x=as.factor(membership), y=n_PFAS)) + 
  geom_col(aes(fill=as.factor(membership))) +
  geom_text(aes(label=n_PFAS), color="black")+
  theme_bw() +
  labs(x="Mixture", y="Number of PFAS")
  
library(ggpubr)
ggarrange(a, b, c, d, nrow=4, labels = c("A", "B", "C", "D"), legend =FALSE) 
ggsave("clean_data/PFAS_ome_upload/Fig_mixtures.jpg", height = 12, width = 4.5, dpi = 500, bg = 'white')


# Network analysis (T1)  -----
cor <- rcorr(as.matrix(T1),  type="spearman") 
r <- cor[[1]] %>% data.frame() %>% mutate_all(~ifelse(. >= 0.6, 1, 0))
p <- cor[[3]] %>% data.frame() %>% mutate_all(~ifelse(. < 0.05, 1, 0))

# identify those r > 0.6 AND p < 0.05 
com <- r+p ; 
com <- com %>% mutate_all(~ifelse(. == 2, 1, 0))

# rm PFAS without any correlations 
remove_col <- colnames(com)[colSums(com, na.rm = TRUE) == 0]
remove_row <- which(row.names(com) %in% remove_col)
com1 <- com %>% select(!all_of(remove_col)) %>% slice(-remove_row)

remove_col <- which(!colnames(r1) %in% colnames(com1))
remove_row <- which(!row.names(r1) %in% colnames(com1))

# run network analysis using directed model 
com1 <- as.matrix(com1)
# write_csv(data.frame(com1) , "clean_data/PFAS_ome_upload/sourcedata_figureS1.csv")
network2 <- graph_from_adjacency_matrix(adjmatrix=com1, mode = "undirected")
summary(network2)

community1 <- cluster_louvain(network1)
communitya <- cluster_louvain(network2)

n_cluster <- length(unique(communitya$membership))
membershiptable <- data.frame(table(communitya$membership))
colnames(membershiptable) <- c("old_mem", "num_PFAS")
membershiptableb <- membershiptable %>% 
  mutate(num_PFAS = as.numeric(num_PFAS), 
         new_mem = rank(-num_PFAS, ties.method = "random"))

# plot the network 
layout <- layout_nicely(network2)
communitya$names
ft_type <- c(rep(1, 16), rep(2, 372-15)) 
newmembership <- data.frame(old_mem = as.character(communitya$membership)) %>% right_join(membershiptableb, by="old_mem")
membership <- paste0("M", newmembership$new_mem)
membership2 <- newmembership$new_mem
library(RColorBrewer)
membership_list <- membership2

# extract information from the network
# compute node
nodes <- V(network2)
nodes <- attr(nodes, "names")

set.seed(5)
plot.igraph(network2, 
            layout = layout, 
            vertex.label = NA,
            mark.groups=communities(communitya),
            mark.border=NA,
            margin=0,
            edge.color="darkgray",
            vertex.color=c("gray", "#28282B")[ft_type],
            vertex.label=nodes,
            vertex.color=rainbow(n_cluster, alpha=0.6)[communitya$membership],
            vertex.shape=c("square", "circle")[ft_type],
            vertex.size = c(3,2)[ft_type],
            vertex.frame.color="black") # plot without labels 

# with labels
plot.igraph(network2, 
            layout = layout, 
            mark.groups=communities(communitya),
            mark.border=NA,
            margin=0.05,
            edge.color="darkgray",
            vertex.color=c("gray", "#28282B")[ft_type],
            vertex.color=rainbow(n_cluster, alpha=0.6)[communitya$membership],
            vertex.shape=c("square", "circle")[ft_type],
            vertex.size = c(3,2)[ft_type],
            vertex.frame.color="black",
            vertex.label=membership_list) # plot without labels 


# compute the degree centrality 
degr_cent <- centr_degree(network2, mode = 'all')
degr_cent <- degr_cent$res

# compute the eigenvector centrality 
eign_cent <- eigen_centrality(network2)
eign_cent <- eign_cent$vector

# compute the closeness centraility
clos_cent <- closeness(network2)

# compute betweeness centrality
betw_cent <- betweenness(network2)

# create data frame storing all of the measures of centrality
net_data1 <- data.frame(vertex = nodes,
                        membership = membership2, 
                        degree = degr_cent, 
                        eigen = eign_cent, 
                        closeness = clos_cent, 
                        betweeness = betw_cent)


net_data2 <- data.frame(membership = membership2, 
                        vertex = nodes,
                        degree = degr_cent, 
                        eigen = eign_cent, 
                        closeness = clos_cent, 
                        betweeness = betw_cent) %>%
  left_join(pfas, join_by(vertex == ft_ID)) %>%
  select(membership, vertex, mz, time, Mode, Level, `Lv1@Metabolite_Name`, `Lv34@Name`, degree, eigen, closeness, betweeness)


net_summary <- 
  net_data1 %>% 
  inner_join(corricc, join_by(vertex == chemical)) %>%
  group_by(membership) %>%
  summarise(n_PFAS = n(), 
            n_target = sum(str_detect(vertex, "tar.")),
            min_r = min(r),
            max_r = max(r), 
            median_r = median(r),
            min_icc = min(ICC_adjusted),
            max_icc = max(ICC_adjusted), 
            median_icc = median(ICC_adjusted)
  )


twoT_member <- read_csv("clean_data/PFAS_ome_upload/networkT0T1alignment.csv")
twoT_member1 <- twoT_member %>% 
  mutate(new_T1 = ifelse(is.na(new_T1), "Other", new_T1)) %>%
  group_by(Mode, T0, new_T1) %>%
  summarise(count = n()) %>%
  mutate(T0 = as.character(T0)) %>% 
  filter(T0 %in% c(1:13))

twoT_member1$T0 <- factor(twoT_member1$T0, levels = c("1","2","3","4","5", "6", "7", "8","9","10", "11","12","13"))
twoT_member1$new_T1 <- factor(twoT_member1$new_T1, levels = c("1","2","3","4","5", "6", "7", "8","9","10", "11","12","13", "48", "49", "Other"))
colnames(twoT_member1)[3] <- "T1"
# write_csv(twoT_member1, "clean_data/PFAS_ome_upload/sourcedata_figure4.csv")

library(ggalluvial)
(memberchange <- ggplot(twoT_member1,
       aes(y = count, axis1 = T0, axis2 = T1)) +
  geom_alluvium(aes(fill = Mode), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_y_continuous(name = "Cumulative Number of PFAS Across Mixtures",
                     breaks = seq(0, 300, by = 50),
                     labels = function(x) 300 - x) + # Flip labels
  labs(title = "", x = "", y = "Cumulative Number of PFAS Across Mixtures") +
  theme_minimal()) +
  theme(axis.text.x = element_blank())

ggsave("clean_data/PFAS_ome_upload/memberchange.jpeg", plot = memberchange, width = 8, height = 9, dpi = 500)


twoT_member2 <- twoT_member %>% 
  mutate(new_T1 = ifelse(is.na(new_T1), "Other", new_T1)) %>%
  group_by(T1, new_T1)%>%
  summarise(count = n()) %>%
  mutate(T1 = as.character(T1))

twoT_member3 <- twoT_member %>% 
  mutate(new_T1 = ifelse(is.na(new_T1), "Other", new_T1)) %>%
  mutate(T1 = as.character(T1))

net_summary2 <- net_summary %>% 
  mutate(membership = as.character(membership)) %>%
  right_join(twoT_member2, by = c("membership"= "T1"))

net_data3 <- net_data2 %>% mutate(membership = as.character(membership)) %>%
  right_join(twoT_member3, by = c("vertex"= "ft"), relationship = "many-to-many")
  

