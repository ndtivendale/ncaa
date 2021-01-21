#####set working directory and load packages#####
setwd("C:/temp/noncanonical_amino_acids_r_processing/protein_and_peptide")
library(VennDiagram)
library(ggplot2)
library(pacman)
library(scales)
library(dplyr)
library(stringr)
library(extrafont)
font_import()
loadfonts(device = "win")

#####import combined 15N, AHA, HPG data#####
ahTagUnenri <- read.csv("unenriTag.csv", stringsAsFactors = F)
ahTagEnri <- read.csv("enriTag.csv", stringsAsFactors = F)
data_full <- read.csv("15N_AHA_HPG_full.csv", stringsAsFactors = F) %>% mutate(n_tagged_pep_AHA = NA)
highConf <- read.csv("high_confidence_proteins_grouped.csv", stringsAsFactors = F)
wf <- windowsFonts()

#####find the overlaps#####
nahN <-  filter(data_full, is.na(LPF) == F)
n_N <- length(unique(nahN$prot_group_15N))
nahA <- filter(data_full, is.na(n_tagged_pep_AHA) == F)
n_nahA <- length(unique(nahA$prot))
nahH <- filter(data_full, tagged_unenriched_HPG == "Y")
n_nahH <- length(unique(nahH$prot))
nahNA <- filter(data_full, is.na(LPF) == F & is.na(n_tagged_pep_AHA) == F)
n_nahNA <- length(unique(nahNA$prot))
nahNH <- filter(data_full, is.na(LPF) == F & is.na(n_tagged_pep_HPG) == F)
n_nahNH <- length(unique(nahNH$prot))
nahAH <- filter(data_full, is.na(n_tagged_pep_AHA) == F & is.na(n_tagged_pep_HPG) == F)
n_nahAH <- length(unique(nahAH$prot))
nahNAH <- filter(data_full, is.na(LPF) == F,
                 is.na(n_tagged_pep_AHA) == F,
                 is.na(n_tagged_pep_HPG) == F)
n_nahNAH <- length(unique(nahNAH$prot))

naheH <-  filter(data_full, is.na(meanriBAQ_enriched_HPG) == F)
n_He <- length(unique(naheH$prot))
naheA <- filter(data_full, is.na(meanriBAQ_enriched_AHA) == F)
n_Ae <- length(unique(naheA$prot))
naheNH <- filter(naheH, is.na(LPF) == F)
n_NHe <- length(unique(naheNH$prot))
naheNA <- filter(naheA, is.na(LPF) == F)
n_NAe <- length(unique(naheNA$prot))
naheHA <- filter(naheH, is.na(meanriBAQ_enriched_AHA) == F)
n_HAe <- length(unique(naheHA$prot))
naheHAN <- filter(naheHA, is.na(LPF) == F)
n_HANe <- length(unique(naheHAN$prot))

naheH2 <-  filter(data_full, meanFoldEnri_HPG > 0)
n_He2 <- length(unique(naheH2$prot))
naheA2 <- filter(data_full, meanFoldEnri_AHA > 0)
n_Ae2 <- length(unique(naheA2$prot))
naheNH2 <- filter(naheH2, is.na(LPF) == F)
n_NHe2 <- length(unique(naheNH2$prot))
naheNA2 <- filter(naheA2, is.na(LPF) == F)
n_NAe2 <- length(unique(naheNA2$prot))
naheHA2 <- filter(naheH2, meanFoldEnri_AHA > 0)
n_HAe2 <- length(unique(naheHA2$prot))
naheHAN2 <- filter(naheHA2, is.na(LPF) == F)
n_HANe2 <- length(unique(naheHAN2$prot))

#####make venn diagrams showing the overlap of the datasets#####
#Global overlap of the three sets
png("Figure 3 venn_15n_ncaa_tag_overlap.png", width = 500, height = 500)
draw.triple.venn(area1 = n_N,
                 area2 = n_nahA,
                 area3 = n_nahH,
                 n12 = n_nahNA,
                 n23 = n_nahAH,
                 n13 = n_nahNH,
                 n123 = n_nahNAH,
                 fill = c("blue", "yellow", "red"),
                 alpha = 0.5,
                 cex = 2,
                 cat.cex = 2,
                 lty = "blank",
                 category = c(expression(""^{15}~"N"),"AHA", "HPG"),
                 fontfamily = wf$Arial,
                 cat.fontfamily = wf$Arial)
dev.off()
#tagged vs untagged for unenriched sets - AHA
grid.newpage()
draw.pairwise.venn(area1 = ahTagUnenri[1,2],
                   area2 = ahTagUnenri[1,3],
                   cross.area = ahTagUnenri[1,3],
                   fill = "yellow",
                   alpha = 0.5,
                   lty = "blank",
                   category = "AHA")

#tagged vs untagged for unenriched set - HPG
grid.newpage()
draw.pairwise.venn(area1 = ahTagUnenri[2,2],
                   area2 = ahTagUnenri[2,3],
                   cross.area = ahTagUnenri[2,3],
                   fill = "red",
                   alpha = 0.5,
                   lty = "blank",
                   category = "HPG")  

#tagged vs untagged for enriched sets - AHA
grid.newpage()
draw.pairwise.venn(area1 = ahTagEnri[1,2],
                   area2 = ahTagEnri[1,3],
                   cross.area = ahTagEnri[1,3],
                   fill = "yellow",
                   alpha = 0.5,
                   lty = "blank",
                   category = "AHA")

#tagged vs untagged for enriched sets - HPG
grid.newpage()
draw.pairwise.venn(area1 = ahTagEnri[2,2],
                   area2 = ahTagEnri[2,3],
                   cross.area = ahTagEnri[2,3],
                   fill = "red",
                   alpha = 0.5,
                   lty = "blank",
                   category = "HPG")          

#15N tagged compared to ncAA enriched
png("Figure 4E venn_15N_labelled_ncAA_enriched.png", width = 500, height = 500)
draw.triple.venn(area1 = n_N,
                   area2 = n_Ae,
                   area3 = n_He,
                   n12 = n_NAe,
                   n23 = n_HAe,
                   n13 = n_NHe,
                   n123 = n_HANe,
                   fill = c("blue", "yellow", "red"),
                   alpha = 0.5,
                   cex = 2,
                   cat.cex = 2,
                   lty = "blank",
                   category = c(expression(""^{15}~"N"),"AHA", "HPG"),
                   fontfamily = wf$Arial,
                   cat.fontfamily = wf$Arial)
dev.off()

#15N tagged compared to ncAA enriched > 0
png("Figure 4F venn_15N_labelled_ncAA_enriched_more_than_0.png", width = 500, height = 500)
draw.triple.venn(area1 = n_N,
                 area2 = n_Ae2,
                 area3 = n_He2,
                 n12 = n_NAe2,
                 n23 = n_HAe2,
                 n13 = n_NHe2,
                 n123 = n_HANe2,
                 fill = c("blue", "yellow", "red"),
                 alpha = 0.5,
                 cex = 2,
                 cat.cex = 2,
                 lty = "blank",
                 category = c(expression(""^{15}~"N"),"AHA", "HPG"),
                 fontfamily = wf$Arial,
                 cat.fontfamily = wf$Arial)
dev.off()
#####scatter plots for fold-enrichment vs.LPF*riBAQ#####
matching <- data_full
matching_HPG_15N <- matching %>% 
  filter(is.na(LPF.riBAQ) == F & is.na(meanFoldEnri_HPG) == F)
matching_HPG_15N2 <- matching %>% 
  filter(is.na(LPF) == F & is.na(meanFoldEnri_HPG) == F)
matching_AHA_15N <- matching %>% 
  filter(is.na(LPF.riBAQ) == F & is.na(meanFoldEnri_AHA) == F)
matching_AHA_15N2 <- matching %>% 
  filter(is.na(LPF) == F & is.na(meanFoldEnri_AHA) == F)

png("Supp fig S8 meanFoldEnri_HPG_v_LPFiBAQ.png", height = 500, width = 500)
ggplot(matching_HPG_15N, aes(x = LPF.riBAQ, y = meanFoldEnri_HPG, colour = tagged_unenriched_HPG)) +
  geom_point() + 
  theme_minimal(base_size = 18) +
  labs(x = "LPF*riBAQ", y = "Protein fold-enrichment", colour = "Tagged?")
dev.off()

png("Supp fig S7 meanFoldEnri_AHA_v_LPFiBAQ.png", height = 500, width = 500)
ggplot(matching_AHA_15N, aes(x = LPF.riBAQ, y = meanFoldEnri_AHA, colour = tagged_unenriched_AHA)) +
  geom_point()+
  theme_minimal(base_size = 18) +
  labs(x = "LPF*riBAQ", y = "Protein fold-enrichment", colour = "Tagged?")
dev.off()

HPG_N_plot2 <- ggplot(matching_HPG_15N2, aes(x = LPF, y = meanFoldEnri_HPG, colour = tagged_unenriched_HPG)) +
  geom_point()
print(HPG_N_plot2)

AHA_N_plot2 <- ggplot(matching_AHA_15N2, aes(x = LPF, y = meanFoldEnri_AHA, colour = tagged_unenriched_AHA)) +
  geom_point()
print(AHA_N_plot2)

#####generate bar graphs vs. AGI####
aEnri2 <- filter(data_full, meanFoldEnri_AHA > 0)
aEnriPlot <- ggplot(aEnri2, aes(x = factor(prot), y = as.numeric(meanFoldEnri_AHA), fill = tagged_unenriched_AHA)) +
  geom_bar(stat = "identity") +
  ylim(0,1)
print(aEnriPlot)

hEnri2 <- filter(data_full, meanFoldEnri_HPG > 0)
hEnriPlot <- ggplot(hEnri2, aes(x = factor(prot), y = meanFoldEnri_HPG, fill = tagged_unenriched_HPG)) +
  geom_bar(stat = "identity") +
  ylim(0,1)
print(hEnriPlot)

#####scatter and box-and-whisker plots for various protein attributes#####
full_data <- data_full
full_data2 <- filter(full_data, is.na(n_tagged_pep_HPG) == F)
full_data3 <- filter(full_data, is.na(n_tags_HPG) == F)
full_data4 <- filter(full_data, is.na(rel_HPG_pos_1) == F) %>%
  select(c("prot", "meanFoldEnri_HPG", "rel_HPG_pos_1", "rel_HPG_pos_2",
           "rel_HPG_pos_3", "rel_HPG_pos_4", "rel_HPG_pos_5",
           "rel_HPG_pos_6", "rel_HPG_pos_7", "rel_HPG_pos_8")) %>% 
  gather(key = "tag number", value = "relative_position", c("rel_HPG_pos_1", "rel_HPG_pos_2",
                                                  "rel_HPG_pos_3", "rel_HPG_pos_4",
                                                  "rel_HPG_pos_5", "rel_HPG_pos_6",
                                                  "rel_HPG_pos_7", "rel_HPG_pos_8")) %>% 
  mutate_all(funs(str_replace(., "rel_HPG_pos_", "")))
full_data4 <- full_data4[complete.cases(full_data4), ]

png("Supp Fig S4 fold_enri_HPG_v_length.png", height = 500, width = 500)
ggplot(full_data, aes(x = length, y = meanFoldEnri_HPG)) +
  geom_point(na.rm = T) +
  coord_cartesian(xlim = c(0, 1100), ylim = c(-1, 1)) +
  theme_minimal(base_size = 18) +
  labs(x = "Protein length (amino acids)", y = "Protein fold-enrichment")
dev.off()

png("Supp Fig S3 fold_enri_AHA_v_length.png", height = 500, width = 500)
ggplot(full_data, aes(x = length, y = meanFoldEnri_AHA)) +
  geom_point(na.rm = T) +
  coord_cartesian(xlim = c(0, 1100), ylim = c(-1, 1))  +
  theme_minimal(base_size = 18) +
  labs(x = "Protein length (amino acids)", y = "Protein fold-enrichment")
dev.off()

png("Supp Fig S2 fold_enri_HPG_v_percent_met.png", height = 500, width = 500)
ggplot(full_data, aes(x = percent_met, y = meanFoldEnri_HPG)) +
  geom_point() +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal(base_size = 18) +
  labs(x = "Percent Met", y = "Protein fold-enrichment")
dev.off()

png("Supp Fig S1 fold_enri_AHA_v_percent_met.png", height = 500, width = 500)
ggplot(full_data, aes(x = percent_met, y = meanFoldEnri_AHA)) +
  geom_point() +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal(base_size = 18) +
  labs(x = "Percent Met", y = "Protein fold-enrichment")
dev.off()

enri_v_tagPep <- ggplot(full_data2, aes(x = as.factor(n_tagged_pep_HPG), y = meanFoldEnri_HPG)) +
  geom_boxplot(outlier.colour = "indianred", outlier.shape = 16, outlier.size = 2, notch = F) +
  labs(y = "Protein fold-enrichment", x = "No. peptides tagged with HPG")
print(enri_v_tagPep)

png("Supp Fig S5 fold_enri_v_n_HPG_tags_per_protein.png", height = 500, width = 500)
ggplot(full_data3, aes(x = as.factor(n_tags_HPG), y = meanFoldEnri_HPG)) +
  geom_boxplot(outlier.colour = "indianred", outlier.shape = 16, outlier.size = 2, notch = F) +
  labs(y = "Protein fold-enrichment", x = "No. HPG tags per protein") +
  theme_minimal(base_size = 18)
dev.off()

png("Supp Fig S6 fold_enri_v_rel_pos_HPG_tag.png", height = 500, width = 500)
ggplot(full_data4, aes(x = as.numeric(relative_position), y = as.numeric(meanFoldEnri_HPG))) +
  geom_point(na.rm = T) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) + 
  labs(y = "Protein fold-enrichment", x = "Relative position of HPG tag") +
  theme_minimal(base_size = 18)
dev.off()

png("Supp Fig S13 tags_v_met_content.png", width = 500, height = 500)
ggplot(full_data, aes(x = percent_met, y = n_tags_HPG)) +
  geom_point(na.rm = T) +
  labs(y = "No. Met-to-HPG substitutions", x = "Percent Met") +
  theme_minimal(base_size = 18)
dev.off()
