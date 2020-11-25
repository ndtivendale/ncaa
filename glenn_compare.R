#This script compares the Glenn et al (2017) data to data from the current study
#it counts the number of proteins in the relevant protein groups files and then draws a venn diagram
#and writes the comparison to csv.
#It also counts the total number of peptides in Glenn's data and the total number of AHA-tagged peptides

setwd("C:/temp/noncanonical_amino_acids_r_processing/protein_and_peptide")
library(tidyverse)
library(VennDiagram)

hpg <- read.csv("protFoldEnriHPG.csv", stringsAsFactors = F) %>% 
  select(c("firstID", "meanriBAQ_enriched"))
aha <- read.csv("protFoldEnriAHA.csv", stringsAsFactors = F) %>% 
  select(c("firstID", "meanriBAQ_enriched"))
glenn_pep <- read.delim("peptidesGlenn.txt", stringsAsFactors = F) %>% 
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter(((Intensity.Wes_control1_1 > 0 & Intensity.Wes_control1_2 > 0) |
          (Intensity.Wes_control1_1 > 0 & Intensity.Wes_control1_3 > 0) |
          (Intensity.Wes_control1_2 > 0 & Intensity.Wes_control1_3 > 0)) |
         ((Intensity.Wes_control2_1 > 0 & Intensity.Wes_control2_2 > 0) |
          (Intensity.Wes_control2_1 > 0 & Intensity.Wes_control2_3 > 0) |
          (Intensity.Wes_control2_2 > 0 & Intensity.Wes_control2_3 > 0)) |
         ((Intensity.Wes_HS1 > 0 & Intensity.Wes_HS2 > 0) |
          (Intensity.Wes_HS1 > 0 & Intensity.Wes_HS3 > 0) |
          (Intensity.Wes_HS2 > 0 & Intensity.Wes_HS3 > 0)) |
         ((Intensity.Wes_Recovery1 > 0 & Intensity.Wes_Recovery2 > 0) |
          (Intensity.Wes_Recovery1 > 0 & Intensity.Wes_Recovery3 > 0) |
          (Intensity.Wes_Recovery2 > 0 & Intensity.Wes_Recovery3 > 0)))  
glenn_pep_count <- length(unique(glenn_pep$Sequence))
glenn_pep_AHA <- glenn_pep %>% 
  filter(Met..AHA.site.IDs != "")
glenn_pep_AHA_count <- length(unique(glenn_pep_AHA$Sequence))

glenn <- read.delim("proteinGroupsGlenn.txt", stringsAsFactors = F) %>% 
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter(((Intensity.Wes_control1_1 > 0 & Intensity.Wes_control1_2 > 0) |
            (Intensity.Wes_control1_1 > 0 & Intensity.Wes_control1_3 > 0) |
            (Intensity.Wes_control1_2 > 0 & Intensity.Wes_control1_3 > 0)) |
           ((Intensity.Wes_control2_1 > 0 & Intensity.Wes_control2_2 > 0) |
              (Intensity.Wes_control2_1 > 0 & Intensity.Wes_control2_3 > 0) |
              (Intensity.Wes_control2_2 > 0 & Intensity.Wes_control2_3 > 0)) |
           ((Intensity.Wes_HS1 > 0 & Intensity.Wes_HS2 > 0) |
              (Intensity.Wes_HS1 > 0 & Intensity.Wes_HS3 > 0) |
              (Intensity.Wes_HS2 > 0 & Intensity.Wes_HS3 > 0)) |
           ((Intensity.Wes_Recovery1 > 0 & Intensity.Wes_Recovery2 > 0) |
              (Intensity.Wes_Recovery1 > 0 & Intensity.Wes_Recovery3 > 0) |
              (Intensity.Wes_Recovery2 > 0 & Intensity.Wes_Recovery3 > 0))) %>% 
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  mutate(present_glenn = "Y") %>% 
  select(c("firstID", "present_glenn"))
glenn_all_prot_count <- length(unique(glenn$firstID))

yu <- read.csv("yu.csv", stringsAsFactors = F) %>% 
  dplyr::rename("firstID" = "Accession") %>% 
  select("firstID") %>% 
  mutate(present_yu = "Y")

glenn_control <- read.delim("proteinGroupsGlenn.txt", stringsAsFactors = F) %>% 
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter(((Intensity.Wes_control1_1 > 0 & Intensity.Wes_control1_2 > 0) |
            (Intensity.Wes_control1_1 > 0 & Intensity.Wes_control1_3 > 0) |
            (Intensity.Wes_control1_2 > 0 & Intensity.Wes_control1_3 > 0)) |
           ((Intensity.Wes_control2_1 > 0 & Intensity.Wes_control2_2 > 0) |
              (Intensity.Wes_control2_1 > 0 & Intensity.Wes_control2_3 > 0) |
              (Intensity.Wes_control2_2 > 0 & Intensity.Wes_control2_3 > 0))) %>% 
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  mutate(present_glenn = "Y") %>% 
  select(c("firstID", "present_glenn"))

#compare HPG and AHA enriched set to Glenn and Yu
#individual
comp <- full_join(hpg ,aha, by = "firstID", suffix = c("_HPG", "_AHA")) %>% 
  full_join(glenn, by = "firstID") %>% 
  full_join(yu, by = "firstID")
hpg_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_HPG) == F)
hpg_count <- length(unique(hpg_count$firstID))
aha_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_AHA) == F)
aha_count <- length(unique(aha_count$firstID))
glenn_count <- comp %>% 
  filter(present_glenn == "Y")
glenn_count <- length(unique(glenn_count$firstID))
yu_count <- comp %>% 
  filter(present_yu == "Y")
yu_count <- length(unique(yu_count$firstID))

#overlaps of two
hpg_aha_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_AHA) == F & is.na(meanriBAQ_enriched_HPG) == F)
hpg_aha_count <- length(unique(hpg_aha_count$firstID))
hpg_glenn_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_HPG) == F & present_glenn == "Y")
hpg_glenn_count <- length(unique(hpg_glenn_count$firstID))
hpg_yu_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_HPG) == F & present_yu == "Y")
hpg_yu_count <- length(unique(hpg_yu_count$firstID))
aha_glenn_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_AHA) == F & present_glenn == "Y")
aha_glenn_count <- length(unique(aha_glenn_count$firstID))
aha_yu_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_AHA) == F & present_yu == "Y")
aha_yu_count <- length(unique(aha_yu_count$firstID))
glenn_yu_count <- comp %>% 
  filter(present_glenn == "Y" & present_yu == "Y")
glenn_yu_count <- length(unique(glenn_yu_count$firstID))

#overlaps of three
hpg_aha_glenn_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_HPG) == F & is.na(meanriBAQ_enriched_AHA) == F & present_glenn == "Y")
hpg_aha_glenn_count <- length(unique(hpg_aha_glenn_count$firstID))
hpg_aha_yu_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_HPG) == F & is.na(meanriBAQ_enriched_AHA) == F & present_yu == "Y")
hpg_aha_yu_count <- length(unique(hpg_aha_yu_count$firstID))
hpg_glenn_yu_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_HPG) == F & present_glenn == "Y" & present_yu == "Y")
hpg_glenn_yu_count <- length(unique(hpg_glenn_yu_count$firstID))
aha_glenn_yu_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_AHA) == F & present_glenn == "Y" & present_yu == "Y")
aha_glenn_yu_count <- length(unique(aha_glenn_yu_count$firstID))

#overlap of four
hpg_aha_glenn_yu_count <- comp %>% 
  filter(is.na(meanriBAQ_enriched_HPG) == F & is.na(meanriBAQ_enriched_AHA) == F & present_glenn == "Y" & present_yu == "Y")
hpg_aha_glenn_yu_count <- length(unique(hpg_aha_glenn_yu_count$firstID))

venn_data <- data.frame(overlap = c("hpg", "aha", "glenn", "yu", "hpg_aha", "hpg_glenn", "hpg_yu", "aha_glenn", "aha_yu",
                                    "glenn_yu", "hpg_aha_glenn", "hpg_aha_yu", "hpg_glenn_yu", "aha_glenn_yu",
                                    "hpg_aha_glenn_yu"),
                        count = c(hpg_count, aha_count, glenn_count, yu_count, hpg_aha_count, hpg_glenn_count, hpg_yu_count,
                                  aha_glenn_count, aha_yu_count, glenn_yu_count, hpg_aha_glenn_count, hpg_aha_yu_count,
                                  hpg_glenn_yu_count, aha_glenn_yu_count, hpg_aha_glenn_yu_count))
write.csv(venn_data, "venn_data.csv", row.names = F)

grid.newpage()
png("Supp Fig S14 tivendale_glenn_yu_comparison.png", width = 1200, height = 1000)
draw.quad.venn(area1 = hpg_count,
                    area2 = aha_count,
                    area3 = glenn_count,
                    area4 = yu_count,  
                      n12 = hpg_aha_count,
                      n13 = hpg_glenn_count,
                      n14 = hpg_yu_count,
                      n23 = aha_glenn_count,
                      n24 = aha_yu_count,
                      n34 = glenn_yu_count,
                     n123 = hpg_aha_glenn_count,
                     n124 = hpg_aha_yu_count,
                     n134 = hpg_glenn_yu_count,
                     n234 = aha_glenn_yu_count,
                    n1234 = hpg_aha_glenn_yu_count, 
                      fill = c("blue", "yellow", "red", "green"),
                      alpha = 0.5,
                      cex = 1.8,
                      cat.cex = 1.8,
                      lty = "blank",
                      category = c("Tivendale et al., HPG", "Tivendale et al., AHA", "Glenn et al. (2017), AHA", "Yu et al. (2020), AHA"))
dev.off()

mapman <- read.delim("mapman.txt", stringsAsFactors = F, na.strings = "\''") %>% 
  select("IDENTIFIER", "BINCODE", "NAME", "DESCRIPTION") %>% 
  dplyr::rename("prot" = "IDENTIFIER") %>% 
  mutate_all(funs(str_replace(., "\'", ""))) %>% 
  mutate_all(funs(str_replace(., "\'", "")))
mapman <- setNames(mapman, tolower(names(mapman)))
mapman$prot <- toupper(mapman$prot)
mapman <- mapman[complete.cases(mapman), ]
comp <- full_join(comp, mapman, by = c("firstID" = "prot"))
write.csv(comp, "tivendale_glenn_yu_comparison.csv", row.names = F)

#for Glenn control only
comp_control <- full_join(hpg ,aha, by = "firstID", suffix = c("_HPG", "_AHA")) %>% 
  full_join(glenn_control, by = "firstID")
#single list
glenn_count_c <- comp_control %>% 
  filter(present_glenn == "Y")
glenn_count_c <- length(unique(glenn_count_c$firstID))
#overlaps of two
hpg_glenn_count_c <- comp_control %>% 
  filter(is.na(meanriBAQ_enriched_HPG) == F & present_glenn == "Y")
hpg_glenn_count_c <- length(unique(hpg_glenn_count_c$firstID))
aha_glenn_count_c <- comp_control %>% 
  filter(is.na(meanriBAQ_enriched_AHA) == F & present_glenn == "Y")
aha_glenn_count_c <- length(unique(aha_glenn_count_c$firstID))
#overlaps of three
hpg_aha_glenn_count_c <- comp_control %>% 
  filter(is.na(meanriBAQ_enriched_HPG) == F & is.na(meanriBAQ_enriched_AHA) == F & present_glenn == "Y")
hpg_aha_glenn_count_c <- length(unique(hpg_aha_glenn_count_c$firstID))

grid.newpage()
png("Supp Fig S15 tivendale_glenn_comare.png", width = 500, height = 500)
draw.triple.venn(area1 = hpg_count,
                    area2 = aha_count,
                    area3 = glenn_count_c,
                    n12 = hpg_aha_count,
                    n13 = hpg_glenn_count_c,
                    n23 = aha_glenn_count_c,
                    n123 = hpg_aha_glenn_count_c,
                    fill = c("blue", "yellow", "red"),
                    alpha = 0.5,
                    lty = "blank",
                    cex = 2,
                    cat.cex = 2,
                    category = c("HPG", "AHA", "Glenn et al. (2017), AHA"))
dev.off()
