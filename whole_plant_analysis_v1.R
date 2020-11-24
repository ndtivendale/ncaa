setwd("C:/temp/noncanonical_amino_acids_r_processing/protein_and_peptide")
library(tidyverse)
library(gridExtra)
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(VennDiagram)

#import and clean labelled datasets
#HPG
hLabelled <- fread("Met-_HPGSitesUnenriched.txt", na.strings = "")
hLabelled2 <- hLabelled %>%
  filter(is.na(Reverse) == T & is.na(`Potential contaminant`) == T) %>%
  filter(`Intensity NT-19-74-6` > 0) %>%  
  select(c("Protein", "Met->HPG Probabilities", "Intensity NT-19-74-6")) %>%
  dplyr::rename("seq" = "Met->HPG Probabilities", "prot" = "Protein") %>% #change names to match the 15N dataframe
  gather(key = "experiment", value = "Intensity", c("Intensity NT-19-74-6")) %>% #make a separate line for each experiment
  mutate_all(funs(str_replace(.,"Intensity NT", "NT"))) %>%
  select(-"Intensity")
hLabelled2$seq <- str_replace_all(hLabelled2$seq, "\\d", "")
hLabelled2$seq <- str_replace_all(hLabelled2$seq, "\\(", "")
hLabelled2$seq <- str_replace_all(hLabelled2$seq, "\\)", "")
data <- hLabelled2
protein_id <- unique(data$prot)
data_hold <- data.frame(prot = rep("NA", length(protein_id)), n_hits = NA, n_pep = NA)
n_outlier <- 0 # number of outliers
for(i in 1:length(protein_id)){
  sub_data <- filter(data, prot == protein_id[i])
  data_hold$prot[i] <- protein_id[i]
  data_hold$n_hits[i] <- nrow(sub_data) # number of spectrum hits
  data_hold$n_pep[i] <- length(unique(str_to_upper(sub_data$seq))) # transfer to all capital letters (number of unique peptides)
}
hLabelled3 <- data_hold

#AHA
aLabelled <- fread("Met-_AHASitesUnenriched.txt", na.strings = "")
aLabelled2 <- aLabelled %>%
  filter(is.na(Reverse) == T & is.na(`Potential contaminant`) == T) %>%
  filter(`Intensity NT-19-47-1` >0) %>%
  select(c("Protein", "Met->AHA Probabilities", "Intensity NT-19-74-2", "Fasta headers")) %>%
  dplyr::rename("seq" = "Met->AHA Probabilities", "prot" = "Protein") %>% #change names to match the 15N dataframe
  gather(key = "experiment", value = "Intensity", "Intensity NT-19-74-2") %>% #make a separate line for each experiment
  mutate_all(funs(str_replace(.,"Intensity NT", "NT"))) %>%
  select(-"Intensity") 
aLabelled2$seq <- str_replace_all(aLabelled2$seq, "\\d", "")
aLabelled2$seq <- str_replace_all(aLabelled2$seq, "\\(", "")
aLabelled2$seq <- str_replace_all(aLabelled2$seq, "\\)", "")
data <- aLabelled2
protein_id <- unique(data$prot)

data_hold <- data.frame(prot = rep("NA", length(protein_id)), n_hits = NA, n_pep = NA)
n_outlier <- 0 
for(i in 1:length(protein_id)){
  sub_data <- filter(data, prot == protein_id[i])
  data_hold$prot[i] <- protein_id[i]
  data_hold$n_hits[i] <- nrow(sub_data) # number of spectrum hits
  data_hold$n_pep[i] <- length(unique(str_to_upper(sub_data$seq))) # transfer to all capital letters (number of unique peptides)
}
aLabelled3 <- data_hold

#remove superfluous columns, clean and tidy data
list <- read.delim("proteinGroupsUnenrichedAHA.txt", stringsAsFactors = F)
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter(Intensity.NT.19.74.2 > 0)
list_lfq <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Intensity.NT.19.74.2)) %>%
  gather(key = "Experiment", value = "Intensity", "Intensity.NT.19.74.2") %>% 
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.19.74.2)) %>%
  gather(key = "Experiment", value = iBAQ, "iBAQ.NT.19.74.2")
iBAQ <- as.vector(list_iBAQ$iBAQ)

list_clean_complete <- list_lfq %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
  filter(Experiment == "NT.19.74.2") %>%
  group_by(Experiment) %>%
  mutate(meanExpiBAQ = mean(iBAQ, na.rm = T)) %>%
  mutate(sumExpiBAQ = sum(iBAQ, na.rm = T)) %>%
  mutate(riBAQ = iBAQ/sumExpiBAQ) %>%
  ungroup() %>%
  mutate(meanTotiBAQ = mean(sumExpiBAQ, na.rm = T)) %>%
  mutate(factor = meanTotiBAQ/meanExpiBAQ) %>%
  mutate(corriBAQ = iBAQ*factor) %>%
  group_by(Protein.IDs) %>%
  mutate(meanriBAQ = mean(riBAQ, na.rm = T)) %>%
  ungroup()

#use the aLabelled object to add a column (tagged) to the protU dataframe
protU <- list_clean_complete %>%
  group_by(Protein.IDs) %>%
  mutate(meanriBAQ = mean(as.numeric(riBAQ), na.rm = T)) %>%
  ungroup() %>%
  mutate(tagged = NA) %>%
  select(-"Experiment") %>% 
  dplyr::rename("prot" = "Protein.IDs")
av <- na.omit(unique(as.vector(aLabelled3$prot)))
protU$tagged <- na.omit(ifelse(str_detect(protU$prot, str_c(av, collapse = "|")), "Y", "N"))
aAll <- protU %>% 
  select(c("prot", "Fasta.headers", "riBAQ", "tagged"))

#HPG list
list <- read.delim("proteinGroupsUnenrichedHPG.txt", stringsAsFactors = F)

#remove superfluous columns, clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter(Intensity.NT.19.74.6 > 0)
list_lfq <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Intensity.NT.19.74.6)) %>%
  gather(key = "Experiment", value = "LFQIntensity", "Intensity.NT.19.74.6") %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_MSMS <- list_clean %>% 
  subset(select = c(MS.MS.count.NT.19.74.6)) %>%
  gather(key = "Experiment", value = "MS.MS.count", "MS.MS.count.NT.19.74.6")
MSMSCount <- as.vector(list_MSMS$MS.MS.count)

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.19.74.6)) %>%
  gather(key = "Experiment", value = iBAQ, "iBAQ.NT.19.74.6")
iBAQ <- as.vector(list_iBAQ$iBAQ)
list_clean_complete <- list_lfq %>%
  cbind(MSMSCount) %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
  filter(Experiment == "NT.19.74.6") %>%
  group_by(Experiment) %>%
  mutate(meanExpiBAQ = mean(iBAQ, na.rm = T)) %>%
  mutate(sumExpiBAQ = sum(iBAQ, na.rm = T)) %>%
  mutate(riBAQ = iBAQ/sumExpiBAQ) %>%
  ungroup() %>%
  mutate(meanTotiBAQ = mean(sumExpiBAQ, na.rm = T)) %>%
  mutate(factor = meanTotiBAQ/meanExpiBAQ) %>%
  mutate(corriBAQ = iBAQ*factor) %>%
  group_by(Protein.IDs) %>%
  mutate(meanriBAQ = mean(riBAQ, na.rm = T)) %>%
  ungroup()

#use the hLabelled object to add a column (tagged) to the protU dataframe
protU <- list_clean_complete %>%
  group_by(Protein.IDs) %>%
  mutate(meanriBAQ = mean(as.numeric(riBAQ), na.rm = T)) %>%
  ungroup() %>%
  mutate(tagged = NA) %>%
  select(-"Experiment") %>%
  dplyr::rename("prot" = "Protein.IDs")
hv <- na.omit(unique(as.vector(hLabelled3$prot)))
protU$tagged <- na.omit(ifelse(str_detect(protU$prot, str_c(hv, collapse = "|")), "Y", "N"))

hAll <- protU %>% 
  select(c("prot", "Fasta.headers", "riBAQ", "tagged"))

write.csv(aAll, "AHA_whole_plant.csv", row.names = F)
write.csv(hAll, "HPG_whole_plant.csv", row.names = F)

a <- length(unique(aAll$prot))
b <- length(unique(hAll$prot))
c <- filter(aAll, tagged == "Y")
c <- length(unique(c$prot))
d <- filter(hAll, tagged == "Y")
d <- length(unique(d$prot))

grid.newpage()
draw.pairwise.venn(area1 = a,
                   area2 = c,
                   cross.area = c)

grid.newpage()
draw.pairwise.venn(area1 = b,
                   area2 = d,
                   cross.area = d)
