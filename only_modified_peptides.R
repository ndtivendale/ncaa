#set working directory and load packages
setwd("C:/temp/noncanonical_amino_acids_r_processing/protein_and_peptide")
library(tidyverse)

#load and filter MQ output files
untagged <- read.delim("proteinGroupsUntagged.txt", stringsAsFactors = F)
aha_e <- read.delim("proteinGroupsEnrichedAHA.txt", stringsAsFactors = F)
aha_u <- read.delim("proteinGroupsunenrichedAHA.txt", stringsAsFactors = F)
hpg_e <- read.delim("proteinGroupsEnrichedHPG.txt", stringsAsFactors = F)
hpg_u <- read.delim("proteinGroupsUnenrichedHPG.txt", stringsAsFactors = F)

#clean and tidy data
untagged_low_conf <- untagged %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.19.45.01 > 0 & Intensity.NT.19.45.04 > 0) |
           (Intensity.NT.19.45.01 > 0 & Intensity.NT.19.45.07 > 0) |
           (Intensity.NT.19.45.04 > 0 & Intensity.NT.19.45.07 > 0)) %>% 
  filter(Score < 0) %>% 
  select(c("Protein.IDs", "Score", "Q.value", "Only.identified.by.site"))

aha_e_low_conf <- aha_e %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.20.92.2A > 0 & Intensity.NT.20.21.23A > 0) |
           (Intensity.NT.20.92.2A > 0 & Intensity.NT.20.21.24A > 0) |
           (Intensity.NT.20.21.23A > 0 & Intensity.NT.20.21.24A > 0)) %>% 
  filter(Score < 0) %>% 
  select(c("Protein.IDs", "Score", "Q.value", "Only.identified.by.site"))

aha_u_low_conf <- aha_u %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.19.47.1 > 0 & Intensity.NT.19.47.4 > 0) |
           (Intensity.NT.19.47.1 > 0 & Intensity.NT.19.47.7 > 0) |
           (Intensity.NT.19.47.4 > 0 & Intensity.NT.19.47.7 > 0)) %>% 
  filter(Score < 0) %>% 
  select(c("Protein.IDs", "Score", "Q.value", "Only.identified.by.site"))

hpg_e_low_conf <- hpg_e %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.20.21.25A > 0 & Intensity.NT.20.21.26A > 0) |
           (Intensity.NT.20.21.25A > 0 & Intensity.NT.20.21.27A > 0) |
           (Intensity.NT.20.21.26A > 0 & Intensity.NT.20.21.27A > 0)) %>% 
  filter(Score < 0) %>% 
  select(c("Protein.IDs", "Score", "Q.value", "Only.identified.by.site"))

hpg_u_low_conf <- hpg_u %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.19.47.10 > 0 & Intensity.NT.19.47.13 > 0) |
           (Intensity.NT.19.47.10 > 0 & Intensity.NT.19.47.16 > 0) |
           (Intensity.NT.19.47.13 > 0 & Intensity.NT.19.47.16 > 0)) %>% 
  filter(Score < 0) %>% 
  select(c("Protein.IDs", "Score", "Q.value", "Only.identified.by.site"))

#make one low confidence protein df and then turn it into a vector, splitting protein groups into idividual IDs
low_conf <- rbind(aha_e_low_conf, aha_u_low_conf, hpg_e_low_conf, hpg_u_low_conf, untagged_low_conf)
write_csv(low_conf, "low_confidence_proteins_combined.csv")

#the low confidence list is manually curated to find if the low confidence proteins were identified in any other analysis
#the resulting curated list is read back in as low_conf2 and used to make a vector which can be used to filter the relevant tables
low_conf2 <- read_csv("low_confidence_proteins_combined_2.csv")
l <- strsplit(as.character(low_conf2$Protein.IDs), ';')
low_conf_V <- unique(unlist(l))

#make sure the high confidence protein list doesn't contain any low confidence proteins
high_conf <- read_csv("high_confidence_proteins.csv") %>% 
  filter(prot %in% low_conf_V)

#make sure the HPG labelled protein list doesn't contain any low confidence proteins
hLabelled4 <- filter(hLabelled3, prot %in% low_conf_V)

#find the low confidence proteins in the supp tables so they can be lablled as such
S2 <- read_csv("S2.csv") %>% 
  filter(str_detect(`Protein group`, paste(low_conf_V, collapse = "|")))
S3 <- read_csv("S3.csv") %>% 
  filter(str_detect(`Protein group`, paste(low_conf_V, collapse = "|")))
S4 <- read_csv("S4.csv") %>% 
  filter(firstID %in% low_conf_V)
S5 <- read_csv("S5.csv") %>% 
  filter(firstID %in% low_conf_V)
S6 <- read_csv("S6.csv") %>% 
  filter(prot %in% low_conf_V)
S10 <- read_csv("S10.csv") %>% 
  filter(firstID %in% low_conf_V)
