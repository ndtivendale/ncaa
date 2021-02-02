#####set working directory and load packages#####
setwd("C:/temp/noncanonical_amino_acids_r_processing/protein_and_peptide")
#PCA
library(devtools)
library(ggbiplot)
library(dplyr)
library(tidyr)
library(tibble)

#####HPG enriched vs. unenriched#####
hpg_unenriched <- read.csv("hpg_unenriched_for_pca.csv", stringsAsFactors = F) %>%
  select(c("firstID", "Experiment", "riBAQ")) %>% 
  filter(riBAQ > 10e-5)
hpg_enriched <- read.csv("hpg_enriched_for_pca.csv", stringsAsFactors = F) %>%
  select(c("firstID", "Experiment", "riBAQ")) %>% 
  filter(riBAQ > 10e-5)
hpg_combined <- rbind(hpg_unenriched, hpg_enriched) %>%
  mutate(riBAQ = as.numeric(riBAQ)) %>% 
  spread(key = firstID, value = riBAQ) %>% 
  filter(Experiment == "NT.19.47.10" | Experiment == "NT.19.47.13" | Experiment == "NT.19.47.16" |
           Experiment == "NT.20.21.25A" | Experiment == "NT.20.21.26A" | Experiment == "NT.20.21.27A") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames(., var = "Experiment") 
hpg_pca <- prcomp(hpg_combined, center = T, scale. = T)
hpg_treatment <- c(rep("Bulk", 3), rep("Enriched", 3))
png("Figure 4C hpg_enri_v_unenri.png", width = 500, height = 500)
ggbiplot(hpg_pca, scale = 1, var.axes = F, ellipse = T, groups = hpg_treatment, labels.size = 20) +
  ggtitle("C") +
  theme_minimal(base_size = 18) +
  theme(axis.title.y = element_text(size = 18, hjust = .4),
        axis.title.x = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  labs(x = "PC1 (54 %)", y = "PC2 (21 %)", colour = "Groups")
dev.off()
#####AHA enriched vs. unenriched#####
aha_unenriched <- read.csv("aha_unenriched_for_pca.csv", stringsAsFactors = F) %>%
  mutate(treatment = "unenriched") %>% 
  select(c("firstID", "Experiment", "riBAQ")) %>% 
  filter(riBAQ > 10e-5)
aha_enriched <- read.csv("aha_enriched_for_pca.csv", stringsAsFactors = F) %>%
  mutate(treatment = "enriched") %>% 
  select(c("firstID", "Experiment", "riBAQ")) %>% 
  filter(riBAQ > 10e-5)
aha_combined <- rbind(aha_unenriched, aha_enriched)%>%
  mutate(riBAQ = as.numeric(riBAQ)) %>% 
  spread(key = firstID, value = riBAQ) %>% 
  filter(Experiment == "NT.19.47.1" | Experiment == "NT.19.47.4" | Experiment == "NT.19.47.7" |
           Experiment == "NT.20.21.23A" | Experiment == "NT.20.21.24A" | Experiment == "NT.20.92.2A") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames(var = "Experiment") 
aha_pca <- prcomp(aha_combined, center = T, scale. = T)
aha_treatment <- c(rep("Bulk", 3), rep("Enriched", 3))
png("Figure 4B aha_enri_v_unenri.png", width = 500, height = 500)
ggbiplot(aha_pca, scale = 1, var.axes = F, ellipse = T, groups = aha_treatment) +
  ggtitle("B") +
  theme_minimal(base_size = 18) +
  theme(axis.title.y = element_text(size = 18, hjust = .4),
        axis.title.x = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  labs(x = "PC1 (55 %)", y = "PC2 (23 %)", colour = "Groups")
dev.off()

#####comparing the enriched AHA and HPG datasets######
enri_comb <- rbind(hpg_enriched, aha_enriched) %>%
  mutate(riBAQ = as.numeric(riBAQ)) %>% 
  spread(key = firstID, value = riBAQ) %>% 
  filter(Experiment == "NT.20.21.25A" | Experiment == "NT.20.21.26A" | Experiment == "NT.20.21.27A" |
           Experiment == "NT.20.21.23A" | Experiment == "NT.20.21.24A" | Experiment == "NT.20.92.2A") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames(var = "Experiment") 
enri_pca <- prcomp(enri_comb, center = T, scale. = T)
comb_treatment <- c(rep("AHA", 2), rep("HPG", 3), "AHA")
png("Figure 4D aha_v_hpg_enri.png", width = 500, height = 500)
ggbiplot(enri_pca, scale = 1, var.axes = F, ellipse = T, groups = comb_treatment) +
  ggtitle("D")+
  theme_minimal(base_size = 18) +
  theme(axis.title.y = element_text(size = 18, hjust = .4),
        axis.title.x = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  labs(x = "PC1 (55 %)", y = "PC2 (22 %)", colour = "Groups")
dev.off()

untagged <- read.csv("protGroups_untreated_clean.csv", stringsAsFactors = F) %>%
  mutate(treatment = "unenriched") %>% 
  select(c("firstID", "Experiment", "riBAQ")) %>% 
  filter(riBAQ > 10e-5)

#####iBAQ numbers compared between AHA, HPG and untagged datasets######
unenri_comb <- rbind(hpg_unenriched, aha_unenriched, untagged) %>%
  mutate(riBAQ = as.numeric(riBAQ)) %>% 
  spread(key = firstID, value = riBAQ) %>% 
  filter(Experiment == "NT.19.47.1" | Experiment == "NT.19.47.4" | Experiment == "NT.19.47.7" |
           Experiment == "NT.19.47.10" | Experiment == "NT.19.47.13" | Experiment == "NT.19.47.16" |
           Experiment == "NT.19.45.01" | Experiment == "NT.19.45.04" | Experiment == "NT.19.45.07") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames(var = "Experiment") 
unenri_pca <- prcomp(unenri_comb, center = T, scale. = T)
comb_treatment2 <- c(rep("untagged", 3), "AHA", rep("HPG", 3), rep("AHA", 2))
png("Figure 4A untagged_aha_hpg.png", width = 500, height = 500)
ggbiplot(unenri_pca, scale = 1, var.axes = F, ellipse = T, groups = comb_treatment2) +
  ggtitle("A")+
  theme_minimal(base_size = 18) +
  theme(axis.title.y = element_text(size = 18, hjust = .4),
        axis.title.x = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12,)) +
  labs(x = "PC1 (19 %)", y = "PC2 (18 %)", colour = "Treatment")
dev.off()
