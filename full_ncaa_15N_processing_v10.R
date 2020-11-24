#######################################################################################################################################
#This script takes 15N and ncAA dataframes and matches the protein groups between them. The approach uses a mixture of best ID (i.e.  #
#first ID matching and matching of individual proteins, regardless of where they are in the group.
#######################################################################################################################################

######set working directory and load packages######
setwd("C:/temp/noncanonical_amino_acids_r_processing/protein_and_peptide")
library(gridExtra)
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(dplyr)
library(utils)
library(tidyr)

######read in files#####
nLabelled <- fread("15NCellCulture24h.txt", na.strings = "") #15N dataset
aLabelled <- fread("Met-_AHASitesUnenriched.txt", na.strings = "") #AHA dataset
hLabelled <- fread("Met-_HPGSitesUnenriched.txt", na.strings = "") #HPG dataset

######do a data tidy on each file#####
#15N
nLabelled2 <- nLabelled %>%
  select(c("seq", "prot", "frac", "desc", "ria", "enrich")) %>%
  mutate(seq = toupper(seq)) %>% #some of the sequences had some lowercase letters; not sure why
  mutate_each(funs(str_replace_all(., "\"", ""))) %>%
  mutate(prot = str_replace(prot, "\\{", "")) %>% mutate(prot = str_replace(prot, "\\}", "")) %>% mutate(prot = str_replace_all(prot, "\\,", ";")) %>% #remove unneccessary characters from prot
  dplyr::rename("experiment" = "frac") %>% #change variable name to match ncAA datasets
  mutate(experiment = str_sub(experiment, end = 2))
nLabelled2$experiment <- paste("NT-19-59", nLabelled2$experiment, sep = "-") #add experiment number

#AHA
aLabelled2 <- aLabelled %>%
  filter(is.na(Reverse) == T & is.na(`Potential contaminant`) == T) %>%
  filter((`Intensity NT-19-47-1` > 0 & `Intensity NT-19-47-4` > 0) |
           (`Intensity NT-19-47-1` > 0 & `Intensity NT-19-47-7` > 0) |
           (`Intensity NT-19-47-4` > 0 & `Intensity NT-19-47-7` > 0)) %>%
  select(c("Protein", "Met->AHA Probabilities", "Intensity NT-19-47-1", "Intensity NT-19-47-4", "Intensity NT-19-47-7", "Fasta headers")) %>%
  dplyr::rename("seq" = "Met->AHA Probabilities", "prot" = "Protein") %>% #change names to match the 15N dataframe
  gather(key = "experiment", value = "Intensity", c("Intensity NT-19-47-1",
                                                    "Intensity NT-19-47-4", "Intensity NT-19-47-7")) %>% #make a separate line for each experiment
  mutate_all(funs(str_replace(.,"Intensity NT", "NT"))) %>%
  select(-"Intensity") #remove the intesity column

#this removes the Met->ncAA probability scores from the sequences
aLabelled2$seq <- str_replace_all(aLabelled2$seq, "\\d", "")
aLabelled2$seq <- str_replace_all(aLabelled2$seq, "\\(", "")
aLabelled2$seq <- str_replace_all(aLabelled2$seq, "\\)", "")

#HPG
hLabelled2 <- hLabelled %>%
  filter(is.na(Reverse) == T & is.na(`Potential contaminant`) == T) %>%
  filter((`Intensity NT-19-47-10` > 0 & `Intensity NT-19-47-13` > 0) |
           (`Intensity NT-19-47-10` > 0 & `Intensity NT-19-47-16` > 0) |
           (`Intensity NT-19-47-13` > 0 & `Intensity NT-19-47-16` > 0)) %>%  
  select(c("Protein", "Met->HPG Probabilities", "Intensity NT-19-47-10", "Intensity NT-19-47-13", "Intensity NT-19-47-16")) %>%
  dplyr::rename("seq" = "Met->HPG Probabilities", "prot" = "Protein") %>% #change names to match the 15N dataframe
  gather(key = "experiment", value = "Intensity", c("Intensity NT-19-47-10", "Intensity NT-19-47-13", "Intensity NT-19-47-16")) %>% #make a separate line for each experiment
  mutate_all(funs(str_replace(.,"Intensity NT", "NT"))) %>%
  select(-"Intensity")
#this removes the Met->ncAA probability scores from the sequences
hLabelled2$seq <- str_replace_all(hLabelled2$seq, "\\d", "")
hLabelled2$seq <- str_replace_all(hLabelled2$seq, "\\(", "")
hLabelled2$seq <- str_replace_all(hLabelled2$seq, "\\)", "")

#####package peptides into proteins#####
#15N
data <- nLabelled2 %>%
  mutate(ria = as.numeric(ria))
protein_id <- unique(data$prot) #create a list of unique protein groups (prot = protein group)
# main things calculation for-loop
data_hold <- data.frame(prot = rep("NA", length(protein_id)), LPF = NA, n_hits = NA, n_pep = NA, sd = NA, rsd = NA)
n_outlier <- 0 # number of outliers
for(i in 1:length(protein_id)){
  sub_data <- filter(data, prot == protein_id[i])
  data_hold$prot[i] <- protein_id[i]
  data_hold$LPF[i] <- mean(as.numeric(sub_data$ria), na.rm = T) # mean of LPF
  data_hold$n_hits[i] <- nrow(sub_data) # number of spectrum hits
  data_hold$n_pep[i] <- length(unique(str_to_upper(sub_data$seq))) # transfer to all capital letters (number of unique peptides)
  data_hold$sd[i] <- sd(sub_data$ria, na.rm = T) # standard deviation
  data_hold$rsd[i] <- round(data_hold$sd[i]/data_hold$LPF[i]*100, 2) # relative standard deviation by mean
  # if the actual abundance value is very small, use sd filter at mean(data$ria, na.rm = T)*0.25
  if(!is.na(data_hold$rsd[i]) & data_hold$rsd[i] > 25 & data_hold$sd[i] > mean(data$ria, na.rm = T)*0.25){ # remove outliers, which is rsd > 25.0% and sd > 0.15
    n_outlier <- n_outlier + 1
    data_hold$LPF[i] <- NA
    data_hold$sd[i] <- NA
    data_hold$rsd[i] <- NA
  }
}
nLabelled3 <- data_hold %>%
  filter(LPF >= 0.2 & rsd <= 30 & n_pep >= 2 & n_hits > 2) %>%
  mutate_each(funs(str_replace_all(., " ", "")))

#AHA
data <- aLabelled2
protein_id <- unique(data$prot)
#AHA has no labelled proteins, so it is just passed to the next variable unchanged
aLabelled3 <- aLabelled2

#HPG
data <- hLabelled2
protein_id <- unique(data$prot)
# main things calculation for-loop
data_hold <- data.frame(prot = rep("NA", length(protein_id)), n_hits = NA, n_pep = NA)
n_outlier <- 0 # number of outliers
for(i in 1:length(protein_id)){
  sub_data <- filter(data, prot == protein_id[i])
  data_hold$prot[i] <- protein_id[i]
  data_hold$n_hits[i] <- nrow(sub_data) # number of spectrum hits
  data_hold$n_pep[i] <- length(unique(str_to_upper(sub_data$seq))) # transfer to all capital letters (number of unique peptides)
}
hLabelled3 <- data_hold

######join the dataframes based upon each protein ID; since the ncAA tagged dataframes contain only one protein for each peptide
#they can be joined onto the 15N labelled dataframe without modification######
#this code chunk separates the protein IDs from the 15N labelled protein group into separate lines, but retains the protein group
#information
l <- strsplit(as.character(nLabelled3$prot), ';')
nLabelled4 <- data.frame(prot = unlist(l), prot_group_15N = rep(nLabelled3$prot, lengths(l)), LPF = rep(nLabelled3$LPF, lengths(l)),
                         n_hits = rep(nLabelled3$n_hits, lengths(l)), n_pep = rep(nLabelled3$n_pep, lengths(l)),
                         sd = rep(nLabelled3$sd, lengths(l)), rsd = rep(nLabelled3$rsd, lengths(l)))

cmplt <- full_join(nLabelled4, hLabelled3, by = "prot", suffix = c("_15N", "_HPG"))

######to get the iBAQ numbers, we use iBAQ numbers from untreated cells#####
#To get iBAQ numbers for untreated cells
list <- read.delim("proteinGroupsUntagged.txt", stringsAsFactors = F)
#clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.19.45.01 > 0 & Intensity.NT.19.45.04 > 0) |
           (Intensity.NT.19.45.01 > 0 & Intensity.NT.19.45.07 > 0) |
           (Intensity.NT.19.45.04 > 0 & Intensity.NT.19.45.07 > 0))
list_lfq <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Intensity.NT.19.45.01,
                    Intensity.NT.19.45.04, Intensity.NT.19.45.07)) %>% 
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  gather(key = "Experiment", value = "LFQIntensity", c("Intensity.NT.19.45.01", "Intensity.NT.19.45.04",
                                                       "Intensity.NT.19.45.07")) %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.19.45.01, iBAQ.NT.19.45.04, iBAQ.NT.19.45.07)) %>%
  gather(key = "Experiment", value = iBAQ, c("iBAQ.NT.19.45.01", "iBAQ.NT.19.45.04",
                                             "iBAQ.NT.19.45.07"))
iBAQ <- as.vector(list_iBAQ$iBAQ)
list_clean_complete <- list_lfq %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
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
  mutate(sdriBAQ = sd(riBAQ, na.rm = T)) %>% 
  select(c("firstID", "Protein.IDs", "Experiment",
           "meanriBAQ", "sdriBAQ", "corriBAQ", "riBAQ")) %>% 
  ungroup() %>% 
  distinct()
uAll <- list_clean_complete %>%
  dplyr::rename("prot_group_untagged" = "Protein.IDs")
write.csv(uAll, "protGroups_untreated_clean.csv")

#join cmplt and uAll by prot (i.e. by each protein ID in the untagged dataframe, matched to each protein ID in the 15N
#data frame)
l <- strsplit(as.character(uAll$prot_group_untagged), ';')
uAll2 <- data.frame(prot = unlist(l), prot_group_untagged = rep(uAll$prot_group_untagged, lengths(l)),
                    meanriBAQ_untagged = rep(uAll$meanriBAQ, lengths(l)), sdriBAQ_untagged = rep(uAll$sdriBAQ, lengths(l)))

cmplt2 <- left_join(cmplt, uAll2, by = "prot")

#calculate LPF*mean riBAQ for each protein ID
cmplt3 <- cmplt2 %>%
  mutate(LPF.riBAQ = as.numeric(LPF)*meanriBAQ_untagged)

#######generate a column indicating the position of each protein ID within the protein group#####
prot_groupsV <- na.omit(as.vector(unique(cmplt3$prot_group_15N)))
hold <- data.frame(prot_group_15N = NA, pos = NA)

for(i in 1:length(prot_groupsV)) {
  sub_data <- cmplt3 %>%
    filter(prot_group_15N == prot_groupsV[i]) %>%
    select(c("prot", "prot_group_15N")) %>%
    mutate(pos = NA) %>% 
    unique()
  sub_data$pos <- 1:nrow(sub_data)
  if(i == 1){
    hold <- sub_data
  }
  else{
    hold <- hold %>%
      rbind(sub_data) 
  }
}

#bind the position column to the main dataframe
cmplt4 <- cmplt3 %>%
  left_join(hold, by = c("prot", "prot_group_15N"))

#re-order and rename the variables
cmplt5 <- cmplt4 %>%
  select(c("prot", "pos", "LPF", "LPF.riBAQ", "n_pep_15N", "n_pep_HPG", "n_hits_15N", "n_hits_HPG",
           "prot_group_15N", "prot_group_untagged", "meanriBAQ_untagged", "sdriBAQ_untagged")) %>%
  dplyr::rename("pos_in_prot_group_15N" = "pos", "n_tagged_pep_HPG" = "n_pep_HPG", "n_tagged_hits_HPG" = "n_hits_HPG") %>% 
  distinct()

#####calculating enrichment values#####
#HPG unenriched dataset
list <- read.delim("proteinGroupsUnenrichedHPG.txt", stringsAsFactors = F)

#remove superfluous columns, clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.19.47.10 > 0 & Intensity.NT.19.47.13 > 0) |
           (Intensity.NT.19.47.10 > 0 & Intensity.NT.19.47.16 > 0) |
           (Intensity.NT.19.47.13 > 0 & Intensity.NT.19.47.16 > 0))
list_lfq <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Met..HPG.site.positions,Score, id, Intensity.NT.19.45.29,
                    Intensity.NT.19.45.32, Intensity.NT.19.45.35,
                    Intensity.NT.19.47.10, Intensity.NT.19.47.13, Intensity.NT.19.47.16)) %>%
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  gather(key = "Experiment", value = "LFQIntensity", c("Intensity.NT.19.45.29", "Intensity.NT.19.45.32",
                                                       "Intensity.NT.19.45.35", "Intensity.NT.19.47.10",
                                                       "Intensity.NT.19.47.13", "Intensity.NT.19.47.16")) %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_MSMS <- list_clean %>% 
  subset(select = c(MS.MS.count.NT.19.45.29, MS.MS.count.NT.19.45.32, MS.MS.count.NT.19.45.35,
                    MS.MS.count.NT.19.47.10, MS.MS.count.NT.19.47.13, MS.MS.count.NT.19.47.16)) %>%
  gather(key = "Experiment", value = "MS.MS.count", c("MS.MS.count.NT.19.45.29", "MS.MS.count.NT.19.45.32", "MS.MS.count.NT.19.45.35",
                                                      "MS.MS.count.NT.19.47.10", "MS.MS.count.NT.19.47.13", "MS.MS.count.NT.19.47.16"))
MSMSCount <- as.vector(list_MSMS$MS.MS.count)

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.19.45.29, iBAQ.NT.19.45.32, iBAQ.NT.19.45.35, iBAQ.NT.19.47.10,
                    iBAQ.NT.19.47.13, iBAQ.NT.19.47.16)) %>%
  gather(key = "Experiment", value = iBAQ, c("iBAQ.NT.19.45.29", "iBAQ.NT.19.45.32", "iBAQ.NT.19.45.35", "iBAQ.NT.19.47.10",
                                             "iBAQ.NT.19.47.13", "iBAQ.NT.19.47.16"))
iBAQ <- as.vector(list_iBAQ$iBAQ)
list_clean_complete <- list_lfq %>%
  cbind(MSMSCount) %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
  filter(Experiment == "NT.19.47.10" | Experiment == "NT.19.47.13" | Experiment == "NT.19.47.16") %>%
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
  mutate(sdriBAQ = sd(riBAQ, na.rm = T)) %>%
  select(c("firstID", "Protein.IDs", "meanriBAQ", "sdriBAQ", "Experiment", "corriBAQ", "Met..HPG.site.positions", "riBAQ")) %>% 
  ungroup()
hpg_unenriched <- list_clean_complete #create an object for 
write.csv(hpg_unenriched, "hpg_unenriched_for_pca.csv", row.names = F)
#use the hLabelled object to add a column (tagged) to the protU dataframe
protU <- list_clean_complete %>%
  mutate(tagged = NA) %>% 
  dplyr::rename("prot_group_HPG" = "Protein.IDs") %>% 
  select(-c("Experiment", "corriBAQ")) %>% 
  distinct()
hv <- na.omit(unique(as.vector(hLabelled3$prot)))
protU$tagged <- na.omit(ifelse(str_detect(protU$prot_group_HPG, str_c(hv, collapse = "|")), "Y", "N"))

hAll <- protU
hf <- filter(hAll, tagged == "Y")
n_hf <- na.omit(unique(hf$prot_group_HPG))

#HPG enriched dataset
list <- read.delim("proteinGroupsEnrichedHPG.txt", stringsAsFactors = F)

#remove superfluous columns, clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.20.21.25A > 0 & Intensity.NT.20.21.26A > 0) |
           (Intensity.NT.20.21.25A > 0 & Intensity.NT.20.21.27A > 0) |
           (Intensity.NT.20.21.26A > 0 & Intensity.NT.20.21.27A > 0))

list_instensity <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Intensity.NT.20.21.25A,
                    Intensity.NT.20.21.26A, Intensity.NT.20.21.27A, Met..HPG.site.positions)) %>%
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  gather(key = "Experiment", value = "LFQIntensity", c("Intensity.NT.20.21.25A", "Intensity.NT.20.21.26A",
                                                       "Intensity.NT.20.21.27A")) %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.20.21.25A, iBAQ.NT.20.21.26A, iBAQ.NT.20.21.27A)) %>%
  gather(key = "Experiment", value = iBAQ, c("iBAQ.NT.20.21.25A", "iBAQ.NT.20.21.26A", "iBAQ.NT.20.21.27A"))
iBAQ <- as.vector(list_iBAQ$iBAQ)

list_clean_complete <- list_instensity %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
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
  mutate(sdriBAQ = sd(riBAQ, na.rm = T)) %>% 
  select(c("firstID", "Protein.IDs", "meanriBAQ", "sdriBAQ", "Experiment", "corriBAQ", "Met..HPG.site.positions", "riBAQ")) %>% 
  ungroup()
hpg_enriched <- list_clean_complete #create an object for 
write.csv(hpg_enriched, "hpg_enriched_for_pca.csv", row.names = F)
protE <- list_clean_complete %>%
  mutate(tagged = NA) %>%
  select(-c("Experiment", "corriBAQ")) %>%
  dplyr::rename("prot_group_HPG" = "Protein.IDs")
hLabelledE <- fread("Met-_HPGSitesEnriched.txt", na.strings = "") %>%
  filter(is.na(Reverse) == T & is.na(`Potential contaminant`) == T)#HPG enriched dataset
hev <- na.omit(unique(as.vector(hLabelledE$Protein)))
protE$tagged <- ifelse(str_detect(protE$prot_group_HPG, str_c(hev, collapse = "|")), "Y", "N")

#match them by protein group and calculate the fold enrichment by dividing enriched riBAQ by unenriched riBAQ
protU <- hAll
prot <- protU %>%
  full_join(protE, by = "firstID", suffix = c("_unenriched", "_enriched")) %>%
  mutate(fold_enrichment = (meanriBAQ_enriched - meanriBAQ_unenriched)/(meanriBAQ_enriched + meanriBAQ_unenriched)) %>%
  mutate(sd_fold_enri = sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched - meanriBAQ_unenriched) +
                               sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched + meanriBAQ_unenriched)))) %>% 
  select(-c("sdriBAQ_unenriched", "sdriBAQ_enriched")) %>%
  select(c("firstID", "prot_group_HPG_unenriched", "prot_group_HPG_enriched", "meanriBAQ_unenriched", "meanriBAQ_enriched",
           "tagged_unenriched", "tagged_enriched", "fold_enrichment", "sd_fold_enri")) %>%
  distinct()

uf <- filter(prot, tagged_unenriched == "Y")
unique(uf$prot_group_HPG_unenriched)
write.csv(prot, "protFoldEnriHPG.csv", row.names = F)

#AHA unenriched dataset
list <- read.delim("proteinGroupsUnenrichedAHA.txt", stringsAsFactors = F)

#remove superfluous columns, clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.19.47.1 > 0 & Intensity.NT.19.47.4 > 0) |
           (Intensity.NT.19.47.1 > 0 & Intensity.NT.19.47.7 > 0) |
           (Intensity.NT.19.47.4 > 0 & Intensity.NT.19.47.7 > 0))
list_lfq <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Met..AHA.site.positions, Intensity.NT.19.45.20,
                    Intensity.NT.19.45.26, Intensity.NT.19.45.23,
                    Intensity.NT.19.47.1, Intensity.NT.19.47.4, Intensity.NT.19.47.7)) %>%
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  gather(key = "Experiment", value = "Intensity", c("Intensity.NT.19.45.20", "Intensity.NT.19.45.26",
                                                    "Intensity.NT.19.45.23", "Intensity.NT.19.47.1",
                                                    "Intensity.NT.19.47.4", "Intensity.NT.19.47.7")) %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.19.45.20, iBAQ.NT.19.45.26, iBAQ.NT.19.45.23, iBAQ.NT.19.47.1,
                    iBAQ.NT.19.47.4, iBAQ.NT.19.47.7)) %>%
  gather(key = "Experiment", value = iBAQ, c("iBAQ.NT.19.45.20", "iBAQ.NT.19.45.26", "iBAQ.NT.19.45.23", "iBAQ.NT.19.47.1",
                                             "iBAQ.NT.19.47.4", "iBAQ.NT.19.47.7"))
iBAQ <- as.vector(list_iBAQ$iBAQ)

list_clean_complete <- list_lfq %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
  filter(Experiment == "NT.19.47.1" | Experiment == "NT.19.47.4" | Experiment == "NT.19.47.7") %>%
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
  mutate(sdriBAQ = sd(riBAQ, na.rm = T)) %>% 
  select(c("firstID", "Protein.IDs", "meanriBAQ", "sdriBAQ", "Experiment", "corriBAQ", "Met..AHA.site.positions", "riBAQ")) %>% 
  ungroup()
aha_unenriched <- list_clean_complete #create an object for 
write.csv(aha_unenriched, "aha_unenriched_for_pca.csv", row.names = F)

#use the aLabelled object to add a column (tagged) to the protU dataframe
protU <- list_clean_complete %>%
  mutate(tagged = NA)
protU$tagged <- "N"
aAll <- protU %>% dplyr::rename("prot_group_AHA" = "Protein.IDs")

#AHA enriched dataset
list <- read.delim("proteinGroupsEnrichedAHA.txt", stringsAsFactors = F)

#remove superfluous columns, clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.20.92.2A > 0 & Intensity.NT.20.21.23A > 0) |
           (Intensity.NT.20.92.2A > 0 & Intensity.NT.20.21.24A > 0) |
           (Intensity.NT.20.21.23A > 0 & Intensity.NT.20.21.24A > 0))

list_instensity <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Intensity.NT.20.92.2A,
                    Intensity.NT.20.21.23A, Intensity.NT.20.21.24A)) %>%
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  gather(key = "Experiment", value = "Intensity", c("Intensity.NT.20.92.2A", "Intensity.NT.20.21.23A", "Intensity.NT.20.21.24A")) %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.20.92.2A, iBAQ.NT.20.21.23A, iBAQ.NT.20.21.24A)) %>%
  gather(key = "Experiment", value = iBAQ, c("iBAQ.NT.20.92.2A", "iBAQ.NT.20.21.23A", "iBAQ.NT.20.21.24A"))
iBAQ <- as.vector(list_iBAQ$iBAQ)

list_clean_complete <- list_instensity %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
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
  mutate(sdriBAQ = sd(riBAQ, na.rm = T)) %>% 
  select(c("firstID", "Protein.IDs", "meanriBAQ", "sdriBAQ", "Experiment", "corriBAQ", "riBAQ")) %>% 
  ungroup()
aha_enriched <- list_clean_complete #create an object for 
write.csv(aha_enriched, "aha_enriched_for_pca.csv")
protE <- list_clean_complete %>%
  mutate(tagged = NA) %>%
  select(-c("Experiment", "corriBAQ")) %>%
  dplyr::rename("prot_group_AHA" = "Protein.IDs")
aLabelledE <- fread("Met-_AHASitesEnriched.txt", na.strings = "") #AHA enriched dataset
ave <- na.omit(unique(as.vector(aLabelledE$Protein)))
protE$tagged <- ifelse(str_detect(protE$prot_group_AHA, str_c(ave, collapse = "|")), "Y", "N")

#match them by protein group and calculate the fold enrichment by dividing enriched riBAQ by unenriched riBAQ
protU <- aAll
prot <- protU %>%
  full_join(protE, by = "firstID", suffix = c("_unenriched", "_enriched")) %>%
  mutate(fold_enrichment = (meanriBAQ_enriched - meanriBAQ_unenriched)/(meanriBAQ_enriched + meanriBAQ_unenriched)) %>%
  mutate(sd_fold_enri = sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched - meanriBAQ_unenriched) +
                               sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched + meanriBAQ_unenriched)))) %>% 
  select(-c("sdriBAQ_unenriched", "sdriBAQ_enriched")) %>%
  select(c("firstID", "prot_group_AHA_unenriched", "prot_group_AHA_enriched", "meanriBAQ_unenriched", "meanriBAQ_enriched",
           "tagged_unenriched", "tagged_enriched", "fold_enrichment", "sd_fold_enri")) %>%
  distinct()

write.csv(prot, "protFoldEnriAHA.csv", row.names = F)

######This part of the script applies the same analysis steps to peptides that were cleaved from the resin
#HPG enriched dataset#####
list <- read.delim("proteinGroupsEnrichedHPG.txt", stringsAsFactors = F)

#remove superfluous columns, clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.20.21.25B > 0 & Intensity.NT.20.21.26B > 0) |
           (Intensity.NT.20.21.25B > 0 & Intensity.NT.20.21.27B > 0) |
           (Intensity.NT.20.21.26B > 0 & Intensity.NT.20.21.27B > 0))

list_instensity <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Intensity.NT.20.21.25B,
                    Intensity.NT.20.21.26B, Intensity.NT.20.21.27B, Met..HPG.site.positions)) %>%
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  gather(key = "Experiment", value = "LFQIntensity", c("Intensity.NT.20.21.25B", "Intensity.NT.20.21.26B",
                                                       "Intensity.NT.20.21.27B")) %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.20.21.25B, iBAQ.NT.20.21.26B, iBAQ.NT.20.21.27B)) %>%
  gather(key = "Experiment", value = iBAQ, c("iBAQ.NT.20.21.25B", "iBAQ.NT.20.21.26B", "iBAQ.NT.20.21.27B"))
iBAQ <- as.vector(list_iBAQ$iBAQ)

list_clean_complete <- list_instensity %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
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
  mutate(sdriBAQ = sd(riBAQ, na.rm = T)) %>% 
  select(c("firstID", "Protein.IDs", "meanriBAQ", "sdriBAQ", "Experiment", "corriBAQ", "Met..HPG.site.positions")) %>% 
  ungroup()
protE <- list_clean_complete %>%
  mutate(tagged = NA) %>%
  select(-c("Experiment", "corriBAQ")) %>%
  dplyr::rename("prot_group_HPG" = "Protein.IDs")
hLabelledE <- fread("Met-_HPGSitesEnriched.txt", na.strings = "") %>%
  filter(is.na(Reverse) == T & is.na(`Potential contaminant`) == T)#HPG enriched dataset
hev <- na.omit(unique(as.vector(hLabelledE$Protein)))
protE$tagged <- ifelse(str_detect(protE$prot_group_HPG, str_c(hev, collapse = "|")), "Y", "N")

#match them by protein group and calculate the fold enrichment by dividing enriched riBAQ by unenriched riBAQ
protU <- hAll
prot <- protU %>%
  full_join(protE, by = "firstID", suffix = c("_unenriched", "_enriched")) %>%
  mutate(fold_enrichment = (meanriBAQ_enriched - meanriBAQ_unenriched)/(meanriBAQ_enriched + meanriBAQ_unenriched)) %>%
  mutate(sd_fold_enri = sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched - meanriBAQ_unenriched) +
                               sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched + meanriBAQ_unenriched)))) %>% 
  select(-c("sdriBAQ_unenriched", "sdriBAQ_enriched")) %>%
  select(c("firstID", "prot_group_HPG_unenriched", "prot_group_HPG_enriched", "meanriBAQ_unenriched", "meanriBAQ_enriched",
           "tagged_unenriched", "tagged_enriched", "fold_enrichment", "sd_fold_enri")) %>%
  distinct()

uf <- filter(prot, tagged_unenriched == "Y")
unique(uf$prot_group_HPG_unenriched)

write.csv(prot, "protFoldEnriHPGCleaved.csv", row.names = F)

#AHA enriched dataset
list <- read.delim("proteinGroupsEnrichedAHA.txt", stringsAsFactors = F)

#remove superfluous columns, clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.20.21.22B.1 > 0 & Intensity.NT.20.21.23B > 0) |
           (Intensity.NT.20.21.22B.1 > 0 & Intensity.NT.20.21.24B > 0) |
           (Intensity.NT.20.21.23B > 0 & Intensity.NT.20.21.24B > 0))

list_instensity <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Intensity.NT.20.21.22B.1,
                    Intensity.NT.20.21.23B, Intensity.NT.20.21.24B)) %>%
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  gather(key = "Experiment", value = "Intensity", c("Intensity.NT.20.21.22B.1", "Intensity.NT.20.21.23B", "Intensity.NT.20.21.24B")) %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.20.21.22B.1, iBAQ.NT.20.21.23B, iBAQ.NT.20.21.24B)) %>%
  gather(key = "Experiment", value = iBAQ, c("iBAQ.NT.20.21.22B.1", "iBAQ.NT.20.21.23B", "iBAQ.NT.20.21.24B"))
iBAQ <- as.vector(list_iBAQ$iBAQ)

list_clean_complete <- list_instensity %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
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
  mutate(sdriBAQ = sd(riBAQ, na.rm = T)) %>% 
  select(c("firstID", "Protein.IDs", "meanriBAQ", "sdriBAQ", "Experiment", "corriBAQ")) %>% 
  ungroup()
protE <- list_clean_complete %>%
  mutate(tagged = NA) %>%
  select(-c("Experiment", "corriBAQ")) %>%
  dplyr::rename("prot_group_AHA" = "Protein.IDs")
aLabelledE <- fread("Met-_AHASitesEnriched.txt", na.strings = "") #AHA enriched dataset
ave <- na.omit(unique(as.vector(aLabelledE$Protein)))
protE$tagged <- ifelse(str_detect(protE$prot_group_AHA, str_c(ave, collapse = "|")), "Y", "N")

#match them by protein group and calculate the fold enrichment by dividing enriched riBAQ by unenriched riBAQ
protU <- aAll
prot <- protU %>%
  full_join(protE, by = "firstID", suffix = c("_unenriched", "_enriched")) %>%
  mutate(fold_enrichment = (meanriBAQ_enriched - meanriBAQ_unenriched)/(meanriBAQ_enriched + meanriBAQ_unenriched)) %>%
  mutate(sd_fold_enri = sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched - meanriBAQ_unenriched) +
                               sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched + meanriBAQ_unenriched)))) %>% 
  select(-c("sdriBAQ_unenriched", "sdriBAQ_enriched")) %>%
  select(c("firstID", "prot_group_AHA_unenriched", "prot_group_AHA_enriched", "meanriBAQ_unenriched", "meanriBAQ_enriched",
           "tagged_unenriched", "tagged_enriched", "fold_enrichment", "sd_fold_enri")) %>%
  distinct()

write.csv(prot, "protFoldEnriAHACleaved.csv", row.names = F)

######This part of the script applies the same analysis steps to peptides that were cleaved from the resin after tryptin digestion
#Import and clean enriched dataset#####
#HPG enriched dataset
list <- read.delim("proteinGroupsEnrichedHPG.txt", stringsAsFactors = F)

#remove superfluous columns, clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.20.21.25C > 0 & Intensity.NT.20.21.26C > 0) |
           (Intensity.NT.20.21.25C > 0 & Intensity.NT.20.21.27C > 0) |
           (Intensity.NT.20.21.26C > 0 & Intensity.NT.20.21.27C > 0))

list_instensity <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Intensity.NT.20.21.25C,
                    Intensity.NT.20.21.26C, Intensity.NT.20.21.27C, Met..HPG.site.positions)) %>%
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  gather(key = "Experiment", value = "LFQIntensity", c("Intensity.NT.20.21.25C", "Intensity.NT.20.21.26C",
                                                       "Intensity.NT.20.21.27C")) %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.20.21.25C, iBAQ.NT.20.21.26C, iBAQ.NT.20.21.27C)) %>%
  gather(key = "Experiment", value = iBAQ, c("iBAQ.NT.20.21.25C", "iBAQ.NT.20.21.26C", "iBAQ.NT.20.21.27C"))
iBAQ <- as.vector(list_iBAQ$iBAQ)

list_clean_complete <- list_instensity %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
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
  mutate(sdriBAQ = sd(riBAQ, na.rm = T)) %>% 
  select(c("firstID", "Protein.IDs", "meanriBAQ", "sdriBAQ", "Experiment", "corriBAQ", "Met..HPG.site.positions")) %>% 
  ungroup()
protE <- list_clean_complete %>%
  mutate(tagged = NA) %>%
  select(-c("Experiment", "corriBAQ")) %>%
  dplyr::rename("prot_group_HPG" = "Protein.IDs")
hLabelledE <- fread("Met-_HPGSitesEnriched.txt", na.strings = "") %>%
  filter(is.na(Reverse) == T & is.na(`Potential contaminant`) == T)#HPG enriched dataset
hev <- na.omit(unique(as.vector(hLabelledE$Protein)))
protE$tagged <- ifelse(str_detect(protE$prot_group_HPG, str_c(hev, collapse = "|")), "Y", "N")

#match them by protein group and calculate the fold enrichment by dividing enriched riBAQ by unenriched riBAQ
protU <- hAll
prot <- protU %>%
  full_join(protE, by = "firstID", suffix = c("_unenriched", "_enriched")) %>%
  mutate(fold_enrichment = (meanriBAQ_enriched - meanriBAQ_unenriched)/(meanriBAQ_enriched + meanriBAQ_unenriched)) %>%
  mutate(sd_fold_enri = sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched - meanriBAQ_unenriched) +
                               sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched + meanriBAQ_unenriched)))) %>% 
  select(-c("sdriBAQ_unenriched", "sdriBAQ_enriched")) %>%
  select(c("firstID", "prot_group_HPG_unenriched", "prot_group_HPG_enriched", "meanriBAQ_unenriched", "meanriBAQ_enriched",
           "tagged_unenriched", "tagged_enriched", "fold_enrichment", "sd_fold_enri")) %>%
  distinct()

uf <- filter(prot, tagged_unenriched == "Y")
unique(uf$prot_group_HPG_unenriched)

write.csv(prot, "protFoldEnriHPGCleavedPostDige.csv", row.names = F)

#AHA enriched dataset
list <- read.delim("proteinGroupsEnrichedAHA.txt", stringsAsFactors = F)

#remove superfluous columns, clean and tidy data
list_clean <- list %>%
  filter(Reverse != "+" & Potential.contaminant != "+") %>%
  filter((Intensity.NT.20.21.22C > 0 & Intensity.NT.20.21.23C > 0) |
           (Intensity.NT.20.21.22C > 0 & Intensity.NT.20.21.24C > 0) |
           (Intensity.NT.20.21.23C > 0 & Intensity.NT.20.21.24C > 0))

list_instensity <- list_clean %>%
  subset(select = c(Protein.IDs, Fasta.headers, Score, id, Intensity.NT.20.21.22C,
                    Intensity.NT.20.21.23C, Intensity.NT.20.21.24C)) %>%
  mutate(firstID = str_sub(Protein.IDs, end = 11)) %>% 
  gather(key = "Experiment", value = "Intensity", c("Intensity.NT.20.21.22C", "Intensity.NT.20.21.23C", "Intensity.NT.20.21.24C")) %>%
  mutate_all(funs(str_replace(.,"Intensity.NT", "NT"))) %>%
  select(-contains("Intensity."))

list_iBAQ <- list_clean %>%
  subset(select = c(iBAQ.NT.20.21.22C, iBAQ.NT.20.21.23C, iBAQ.NT.20.21.24C)) %>%
  gather(key = "Experiment", value = iBAQ, c("iBAQ.NT.20.21.22C", "iBAQ.NT.20.21.23C", "iBAQ.NT.20.21.24C"))
iBAQ <- as.vector(list_iBAQ$iBAQ)

list_clean_complete <- list_instensity %>%
  cbind(iBAQ)
list_clean_complete[list_clean_complete == 0] <- NA
list_clean_complete <- list_clean_complete %>%
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
  mutate(sdriBAQ = sd(riBAQ, na.rm = T)) %>% 
  select(c("firstID", "Protein.IDs", "meanriBAQ", "sdriBAQ", "Experiment", "corriBAQ")) %>% 
  ungroup()
protE <- list_clean_complete %>%
  mutate(tagged = NA) %>%
  select(-c("Experiment", "corriBAQ")) %>%
  dplyr::rename("prot_group_AHA" = "Protein.IDs")
aLabelledE <- fread("Met-_AHASitesEnriched.txt", na.strings = "") #AHA enriched dataset
ave <- na.omit(unique(as.vector(aLabelledE$Protein)))
protE$tagged <- ifelse(str_detect(protE$prot_group_AHA, str_c(ave, collapse = "|")), "Y", "N")

#match them by protein group and calculate the fold enrichment by dividing enriched riBAQ by unenriched riBAQ
protU <- aAll
prot <- protU %>%
  full_join(protE, by = "firstID", suffix = c("_unenriched", "_enriched")) %>%
  mutate(fold_enrichment = (meanriBAQ_enriched - meanriBAQ_unenriched)/(meanriBAQ_enriched + meanriBAQ_unenriched)) %>%
  mutate(sd_fold_enri = sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched - meanriBAQ_unenriched) +
                               sqrt((sqrt(sdriBAQ_enriched + sdriBAQ_unenriched)^2)/(meanriBAQ_enriched + meanriBAQ_unenriched)))) %>% 
  select(-c("sdriBAQ_unenriched", "sdriBAQ_enriched")) %>%
  select(c("firstID", "prot_group_AHA_unenriched", "prot_group_AHA_enriched", "meanriBAQ_unenriched", "meanriBAQ_enriched",
           "tagged_unenriched", "tagged_enriched", "fold_enrichment", "sd_fold_enri")) %>%
  distinct()

write.csv(prot, "protFoldEnriAHACleavedPostDige.csv", row.names = F)

######investigating enriched/unenriched data set#####
aEnri <- read.csv("protFoldEnriAHA_simplified.csv", stringsAsFactors = F)
aClea <- read.csv("protFoldEnriAHACleaved_simplified.csv", stringsAsFactors = F)
aCleaP <- read.csv("protFoldEnriAHACleavedPostDige_simplified.csv", stringsAsFactors = F)
hEnri <- read.csv("protFoldEnriHPG_simplified.csv", stringsAsFactors = F)
hClea <- read.csv("protFoldEnriHPGCleaved_simplified.csv", stringsAsFactors = F)
hCleaP <- read.csv("protFoldEnriHPGCleavedPostDige_simplified.csv", stringsAsFactors = F)
#digested on resin
#how many proteins in enriched set?
x <- filter(aEnri, is.na(meanriBAQ_enriched) == F)
x <- length(unique(x$firstID))
y <- filter(hEnri, is.na(meanriBAQ_enriched) == F)
y <- length(unique(y$firstID))

#of these, how many have a fold enrichment > 0?
i <- filter(aEnri, is.na(meanriBAQ_enriched) == F) %>%
  filter(fold_enrichment > 0)
i <- length(unique(i$firstID))
j <- filter(hEnri, is.na(meanriBAQ_enriched) == F) %>%
  filter(fold_enrichment > 0)
j <- length(unique(j$firstID))

#of those with enrichment > 1, how many are tagged in the unenriched set?
m <- filter(aEnri, fold_enrichment > 0) %>% 
  filter(tagged_unenriched == "Y")
m <- length(unique(m$firstID))
n <- filter(hEnri, fold_enrichment > 0) %>% 
  filter(tagged_unenriched == "Y")
n <- length(unique(n$firstID))

#how many of the enriched proteins are tagged in the uneneriched set?
mu <- filter(aEnri, tagged_unenriched == "Y") %>% 
  filter(is.na(meanriBAQ_enriched) == F)
mu <- length(unique(mu$firstID))
nu <- filter(hEnri, tagged_unenriched == "Y") %>% 
  filter(is.na(meanriBAQ_enriched) == F)
nu <- length(unique(nu$firstID))

#how many are tagged in the enriched set?
p <- filter(aEnri, tagged_enriched == "Y")
p <- length(unique(p$firstID))
q <- filter(hEnri, tagged_enriched == "Y")
q <- length(unique(q$firstID))

enri_table <- data.frame(tag = c("AHA", "HPG"), tagged_prot_unenri = rep("NA", 2), type = c("dige_on_resin", "dige_on_resin"),
                         n_prot_enri = c(x, y), n_tag_wo_enrich_thresh = c(mu, nu), n_above_0 = c(i, j),
                         n_tag_unenri_w_enri_above_0 = c(m, n), n_tag_enri = c(p, q))

#cleved from resin
#how many proteins in enriched set?
x <- filter(aClea, is.na(meanriBAQ_enriched) == F)
x <- length(unique(x$firstID))
y <- filter(hClea, is.na(meanriBAQ_enriched) == F)
y <- length(unique(y$firstID))

#of these, how many have a fold enrichment > 0?
i <- filter(aClea, is.na(meanriBAQ_enriched) == F) %>%
  filter(fold_enrichment > 0)
i <- length(unique(i$firstID))
j <- filter(hClea, is.na(meanriBAQ_enriched) == F) %>%
  filter(fold_enrichment > 0)
j <- length(unique(j$firstID))

#of those with enrichment > 0, how many are tagged in the unenriched set?
m <- filter(aClea, fold_enrichment > 0) %>% 
  filter(tagged_unenriched == "Y")
m <- length(unique(m$firstID))
n <- filter(hClea, fold_enrichment > 0) %>% 
  filter(tagged_unenriched == "Y")
n <- length(unique(n$firstID))

#how many of the enriched proteins are tagged in the uneneriched set?
mu <- filter(aClea, tagged_unenriched == "Y") %>% 
  filter(is.na(meanriBAQ_enriched) == F)
mu <- length(unique(mu$firstID))
nu <- filter(hClea, tagged_unenriched == "Y") %>% 
  filter(is.na(meanriBAQ_enriched) == F)
nu <- length(unique(nu$firstID))

#how many are tagged in the enriched set?
p <- filter(aClea, tagged_enriched == "Y")
p <- length(unique(p$firstID))
q <- filter(hClea, tagged_enriched == "Y")
q <- length(unique(q$firstID))

enri_table <- rbind(enri_table, data.frame(tag = c("AHA", "HPG"), tagged_prot_unenri = rep("NA", 2), type = c("cleaved", "cleaved"),
                                           n_prot_enri = c(x, y), n_tag_wo_enrich_thresh = c(mu, nu), n_above_0 = c(i, j),
                                           n_tag_unenri_w_enri_above_0 = c(m, n), n_tag_enri = c(p, q)))


#cleved from resin post digest
#how many proteins in enriched set?
x <- filter(aCleaP, is.na(meanriBAQ_enriched) == F)
x <- length(unique(x$firstID))
y <- filter(hCleaP, is.na(meanriBAQ_enriched) == F)
y <- length(unique(y$firstID))

#of these, how many have a fold enrichment > 0?
i <- filter(aCleaP, is.na(meanriBAQ_enriched) == F) %>%
  filter(fold_enrichment > 0)
i <- length(unique(i$firstID))
j <- filter(hCleaP, is.na(meanriBAQ_enriched) == F) %>%
  filter(fold_enrichment > 0)
j <- length(unique(j$firstID))

#of those with enrichment > 0, how many are tagged in the unenriched set?
m <- filter(aCleaP, fold_enrichment > 0) %>% 
  filter(tagged_unenriched == "Y")
m <- length(unique(m$firstID))
n <- filter(hCleaP, fold_enrichment > 0) %>% 
  filter(tagged_unenriched == "Y")
n <- length(unique(n$firstID))

#how many of the enriched proteins are tagged in the uneneriched set?
mu <- filter(aCleaP, tagged_unenriched == "Y") %>% 
  filter(is.na(meanriBAQ_enriched) == F)
mu <- length(unique(mu$firstID))
nu <- filter(hCleaP, tagged_unenriched == "Y") %>% 
  filter(is.na(meanriBAQ_enriched) == F)
nu <- length(unique(nu$firstID))

#how many are tagged in the enriched set?
p <- filter(aCleaP, tagged_enriched == "Y")
p <- length(unique(p$firstID))
q <- filter(hCleaP, tagged_enriched == "Y")
q <- length(unique(q$firstID))
enri_table <- rbind(enri_table, data.frame(tag = c("AHA", "HPG"), tagged_prot_unenri = rep("NA", 2), type = c("dige_post_cleave", "dige_post_cleave"),
                                           n_prot_enri = c(x, y), n_tag_wo_enrich_thresh = c(mu, nu), n_above_0 = c(i, j),
                                           n_tag_unenri_w_enri_above_0 = c(m, n), n_tag_enri = c(p, q)))

#what is the total number of tagged proteins in the unenriched sets?
a <- filter(aEnri, is.na(meanriBAQ_unenriched) == F)
a <- length(unique(a$firstID))
h <- filter(hEnri, is.na(meanriBAQ_unenriched) == F)
h <- length(unique(h$firstID))

#and how many of those are tagged?
b <- filter(aEnri, tagged_unenriched == "Y")
b <- length(unique(b$firstID))
i <- filter(hEnri, tagged_unenriched == "Y")
i <- length(unique(i$firstID))

enri_table$tagged_prot_unenri <- c(b, i)

ab <- (b/a)*100
hi <- (i/h)*100

unenriTag <- data.frame("treat" = c("AHA", "HPG"), no.prot = c(a, h), "no.tagged" = c(b, i), "percent.tagged" = c(ab, hi))
enriTag <- data.frame("treat" = c("AHA", "HPG"), "no.prot" = c(x, y), "no.tagged" = c(p, q), "precent.tagged" = c(p/x*100, q/y*100))
write.csv(unenriTag, "unenriTag.csv", row.names = F)
write.csv(enriTag, "enriTag.csv", row.names = F)
write.csv(enri_table, "enriAnalysis.csv", row.names = F)

######join the 15N/HPG/AHA tagging dataframe to the enrichment dataframe#####
hEnri2 <- hEnri %>% 
  select(c("firstID", "fold_enrichment", "sd_fold_enri", "tagged_unenriched", "tagged_enriched",
           "meanriBAQ_unenriched", "meanriBAQ_enriched", "prot_unenriched", "prot_enriched")) %>% 
  filter(is.na(meanriBAQ_enriched) == F | tagged_unenriched == "Y") %>% 
  dplyr::rename("prot_group_HPG_unenriched" = "prot_unenriched",
                "prot_group_HPG_enriched" = "prot_enriched")
aEnri2 <- aEnri %>% 
  select(c("firstID", "fold_enrichment", "sd_fold_enri", "tagged_unenriched", "tagged_enriched",
           "meanriBAQ_unenriched", "meanriBAQ_enriched", "prot_unenriched", "prot_enriched")) %>% 
  filter(is.na(meanriBAQ_enriched) == F | tagged_unenriched == "Y") %>% 
  dplyr::rename("prot_group_AHA_unenriched" = "prot_unenriched",
                "prot_group_AHA_enriched" = "prot_enriched")

#join hEnri and aEnri
cmplt6 <- cmplt5 %>% 
  full_join(hEnri2, by = c("prot" = "firstID"), suffix = c("_15N", "_HPG")) %>%
  full_join(aEnri2, by = c("prot" = "firstID"), suffix = c("_HPG", "_AHA"))
          
#remove surplus rows (i.e. rows containing protein matches when proteins have already been matched to IDs further up the protein group list)
#filtering is done based upon four successive filtering criteria:
#1. position in the 15N protein group (matches higher up the list are kept)
#2. if 1 leaves more than one observation, the number of peptides for that protein in the 15N dataset is used (the observation(s) with the 
#greates number of peptides is(are) kept)
#3. if 2 leaves more than one observation, the number of spectral hits for that protein in the 15N dataset is used (the observation(s) wtih
#the greatest number of hits is(are) kept)
#4. if 3 leaves more than one observation, the number of members in the protein group for the 15N dataset is used (the observation with the
#greater nuber of members is kept)
protV = unique(cmplt6$prot)
cmplt7 <- data.frame(prot = NA,
                     pos_in_prot_group_15N = NA,
                     LPF = NA,
                     LPF.riBAQ = NA,
                     n_pep_15N =NA,
                     n_tagged_pep_HPG =NA,
                     n_hits_15N = NA,
                     n_tagged_hits_HPG = NA,
                     prot_group_15N = NA,
                     prot_group_untagged = NA,
                     meanriBAQ_untagged = NA,
                     sdriBAQ_untagged = NA,
                     fold_enrichment_HPG = NA,
                     sd_fold_enri_HPG = NA,
                     tagged_unenriched_HPG = NA,
                     tagged_enriched_HPG = NA,
                     meanriBAQ_unenriched_HPG = NA,
                     meanriBAQ_enriched_HPG = NA,
                     prot_group_HPG_unenriched = NA,
                     prot_group_HPG_unenriched = NA,
                     fold_enrichment_AHA = NA,
                     sd_fold_enri_AHA = NA,
                     tagged_unenriched_AHA = NA,
                     tagged_enriched_AHA = NA,
                     meanriBAQ_unenriched_AHA = NA,
                     meanriBAQ_enriched_AHA = NA,
                     prot_group_AHA_unenriched = NA,
                     prot_group_AHA_unenriched = NA)
sub_data <- cmplt6 %>%
  filter(prot == protV[1])
if(all(is.na(sub_data$n_tagged_pep_HPG) == T) &
    all(is.na(sub_data$prot_group_HPG_unenriched) == T) &
    all(is.na(sub_data$prot_group_HPG_enriched) == T) &
    all(is.na(sub_data$prot_group_AHA_unenriched) == T) &
    all(is.na(sub_data$prot_group_AHA_unenriched) == T)) {
    sub_data <- filter(sub_data, pos_in_prot_group_15N == min(pos_in_prot_group_15N))
  } else {
    sub_data <- filter(sub_data, is.na(n_tagged_pep_HPG) == F | is.na(prot_group_HPG_unenriched) == F |
                         is.na(prot_group_HPG_enriched) == F |is.na(prot_group_AHA_unenriched) == F |
                         is.na(prot_group_AHA_enriched) == F)
  }
cmplt7 <- sub_data
for(i in 2:length(protV)){
  sub_data <- cmplt6 %>%
    filter(prot == protV[i])
if(all(is.na(sub_data$n_tagged_pep_HPG) == T) &
    all(is.na(sub_data$prot_group_HPG_unenriched) == T) &
    all(is.na(sub_data$prot_group_HPG_enriched) == T) &
    all(is.na(sub_data$prot_group_AHA_unenriched) == T) &
    all(is.na(sub_data$prot_group_AHA_unenriched) == T)) {
  sub_data <- filter(sub_data, pos_in_prot_group_15N == min(pos_in_prot_group_15N))
  } else {
    sub_data <- filter(sub_data, is.na(n_tagged_pep_HPG) == F | is.na(prot_group_HPG_unenriched) == F |
                         is.na(prot_group_HPG_enriched) == F |is.na(prot_group_AHA_unenriched) == F |
                         is.na(prot_group_AHA_enriched) == F)
  }
  cmplt7 <- cmplt7 %>%
    rbind(sub_data)
}
sub_data <- cmplt6 %>%
  filter(is.na(pos_in_prot_group_15N) == T)
cmplt7 <- cmplt7 %>%
  rbind(sub_data)

#filter out any remaining double matches between groups based upon protein groups
protV = unique(cmplt7$prot_group_15N)
cmplt8 <- data.frame(prot = NA,
                     pos_in_prot_group_15N = NA,
                     LPF = NA,
                     LPF.riBAQ = NA,
                     n_pep_15N =NA,
                     n_tagged_pep_HPG =NA,
                     n_hits_15N = NA,
                     n_tagged_hits_HPG = NA,
                     prot_group_15N = NA,
                     prot_group_untagged = NA,
                     meanriBAQ_untagged = NA,
                     sdriBAQ_untagged = NA,
                     fold_enrichment_HPG = NA,
                     sd_fold_enri_HPG = NA,
                     tagged_unenriched_HPG = NA,
                     tagged_enriched_HPG = NA,
                     meanriBAQ_unenriched_HPG = NA,
                     meanriBAQ_enriched_HPG = NA,
                     prot_group_HPG_unenriched = NA,
                     prot_group_HPG_unenriched = NA,
                     fold_enrichment_AHA = NA,
                     sd_fold_enri_AHA = NA,
                     tagged_unenriched_AHA = NA,
                     tagged_enriched_AHA = NA,
                     meanriBAQ_unenriched_AHA = NA,
                     meanriBAQ_enriched_AHA = NA,
                     prot_group_AHA_unenriched = NA,
                     prot_group_AHA_unenriched = NA)
sub_data <- cmplt7 %>%
filter(prot_group_15N == protV[1])
if(all(is.na(sub_data$n_tagged_pep_HPG) == T) &
   all(is.na(sub_data$prot_group_HPG_unenriched) == T) &
   all(is.na(sub_data$prot_group_HPG_enriched) == T) &
   all(is.na(sub_data$prot_group_AHA_unenriched) == T) &
   all(is.na(sub_data$prot_group_AHA_unenriched) == T)) {
  sub_data <- filter(sub_data, pos_in_prot_group_15N == min(pos_in_prot_group_15N))
 } else {
  sub_data <- filter(sub_data, is.na(n_tagged_pep_HPG) == F | is.na(prot_group_HPG_unenriched) == F |
                       is.na(prot_group_HPG_enriched) == F |is.na(prot_group_AHA_unenriched) == F |
                       is.na(prot_group_AHA_enriched) == F)
}
cmplt8 <- sub_data

for(i in 2:length(protV)){
  sub_data <- cmplt7 %>%
     filter(prot_group_15N == protV[i])
  if(all(is.na(sub_data$prot_group_HPG_unenriched) == T) &
     all(is.na(sub_data$n_tagged_pep_HPG) == T) &
     all(is.na(sub_data$prot_group_HPG_enriched) == T) &
     all(is.na(sub_data$prot_group_AHA_unenriched) == T) &
     all(is.na(sub_data$prot_group_AHA_unenriched) == T)) {
    sub_data <- filter(sub_data, pos_in_prot_group_15N == min(pos_in_prot_group_15N))
  } else {
      sub_data <- filter(sub_data, is.na(n_tagged_pep_HPG) == F | is.na(prot_group_HPG_unenriched) == F |
                           is.na(prot_group_HPG_enriched) == F |is.na(prot_group_AHA_unenriched) == F |
                           is.na(prot_group_AHA_enriched) == F)
  }
  cmplt8 <- cmplt8 %>%
    rbind(sub_data)
  }
  
sub_data <- cmplt7 %>%
  filter(is.na(prot_group_15N) == T)
cmplt8 <- cmplt8 %>%
    rbind(sub_data)

#group 15N data by turnover rate into medium, fast and slow.
cmplt8 <- cmplt8 %>%
  distinct()
cmplt8 <- cmplt8 %>% 
  mutate(turnoverGroup = rep("NA", NROW(cmplt8)))
for(i in 1:NROW(cmplt8)) {
  if(is.na(cmplt8$LPF[i] == T)) {
    cmplt8$turnoverGroup[i] <- "NA"
  } else {
  if(cmplt8$LPF[i] <= 0.5) {
   cmplt8$turnoverGroup[i] <- "slow"
  } else {
     if(cmplt8$LPF[i] < 0.7 & cmplt8$LPF[i] > 0.5) {
      cmplt8$turnoverGroup[i] <- "medium"
   }
    else {
      cmplt8$turnoverGroup[i] <- "fast"
      }
   }
  }
}


#####add the various extra columns for analysis#####
#add the Tair10 data based upon the first protein ID in prot_group, or by the ID that was matched between the 15N and HPG datasets
At_tair10 <- read.csv("tair10_fasta_dataframe.csv", stringsAsFactors = F)
cmplt9 <- left_join(cmplt8, At_tair10, by = "prot") %>%
  mutate(n_met = NA)

#count the number of methionines
for(i in 1:nrow(cmplt9)) {
  cmplt9$n_met[i] <- str_count(cmplt9$seq[i], "M")
}

#calculate the percent methionine in each protein
cmplt9 <- cmplt9 %>% 
  mutate(percent_met = (n_met/length)*100) %>%
  distinct()

formapman <- cmplt9 %>% 
  select("prot") %>% 
  distinct()
write.csv(formapman, "formapman.csv", row.names = F)

#add mapman data as new variables
mapman <- read.delim("mapman.txt", stringsAsFactors = F, na.strings = "\''") %>% 
  select("IDENTIFIER", "BINCODE", "NAME", "DESCRIPTION") %>% 
  dplyr::rename("prot" = "IDENTIFIER") %>% 
  mutate_all(funs(str_replace(., "\'", ""))) %>% 
  mutate_all(funs(str_replace(., "\'", "")))
mapman <- setnames(mapman, tolower(names(mapman)))
mapman$prot <- toupper(mapman$prot)
mapman <- mapman[complete.cases(mapman), ]
cmplt10 <- left_join(cmplt9, mapman, by = "prot")
aAll4 <- aAll %>% 
  select(c("firstID", "Met..AHA.site.positions"))
hAll4 <- hAll %>% 
  select(c("firstID", "Met..HPG.site.positions"))
cmplt11 <- cmplt10 %>% 
  left_join(aAll4, by = c("prot" = "firstID")) %>% 
  left_join(hAll4, by= c("prot" = "firstID")) %>% 
  distinct() %>% 
  mutate(n_tags_HPG = NA)

#The way MQ determines tagging for protein groups is a bit different to the way I've done it; this code chunk ensures that
#not spurious Met to HPG or Met to AHA site positions remain in the data
for(i in 1:nrow(cmplt11)) {
  if((is.na(cmplt11$tagged_unenriched_HPG[i]) == T | cmplt11$tagged_unenriched_HPG[i] == "N") &
     (is.na(cmplt11$tagged_enriched_HPG[i]) == T | cmplt11$tagged_enriched_HPG[i] == "N")) {
    cmplt11$Met..HPG.site.positions[i] <- NA
  }
  if((is.na(cmplt11$tagged_unenriched_AHA[i]) == T | cmplt11$tagged_unenriched_AHA[i] == "N") &
     (is.na(cmplt11$tagged_enriched_AHA[i]) == T | cmplt11$tagged_enriched_AHA[i] == "N")) {
    cmplt11$Met..AHA.site.positions[i] <- NA
  }
}

#count the number of tags per protein observation
for(i in 1:nrow(cmplt11)) {
  if(is.na(cmplt11$Met..HPG.site.positions[i]) == F) {
  cmplt11$n_tags_HPG[i] <- (length(unlist(strsplit(as.character(cmplt11$Met..HPG.site.positions[i]), ";"))))
  }
}

#calculate the relative position of each tag in the protein
cmplt12 <- cmplt11 %>%
  separate(Met..HPG.site.positions, into = c("HPG_pos_1", "HPG_pos_2", "HPG_pos_3", "HPG_pos_4", "HPG_pos_5",
                                             "HPG_pos_6", "HPG_pos_7", "HPG_pos_8"), sep = ";") %>%
  mutate(HPG_pos_1 = as.numeric(HPG_pos_1), HPG_pos_2 = as.numeric(HPG_pos_2), HPG_pos_3 = as.numeric(HPG_pos_3),
         HPG_pos_4 = as.numeric(HPG_pos_4), HPG_pos_5 = as.numeric(HPG_pos_5), HPG_pos_6 = as.numeric(HPG_pos_6),
         HPG_pos_7 = as.numeric(HPG_pos_7), HPG_pos_8 = as.numeric(HPG_pos_8)) %>% 
  mutate(rel_HPG_pos_1 = HPG_pos_1/length, rel_HPG_pos_2 = HPG_pos_2/length, rel_HPG_pos_3 = HPG_pos_3/length,
         rel_HPG_pos_4 = HPG_pos_4/length, rel_HPG_pos_5 = HPG_pos_5/length, rel_HPG_pos_6 = HPG_pos_6/length,
         rel_HPG_pos_7 = HPG_pos_7/length, rel_HPG_pos_8 = HPG_pos_8/length) %>%
  distinct()

#add GO and PO domains for each observation
go_domains <- read.delim("PO_GO_domains.txt", stringsAsFactors = F) %>%
  dplyr::rename("prot" = "Protein.stable.ID")

cmplt13 <- cmplt12 %>% 
  left_join(go_domains, by = "prot")

#write the final dataframe
write.csv(cmplt13,"15N_AHA_HPG_full.csv", row.names = F)
write.csv(cmplt9, "15N_AHA_HPG_wo_mapman.csv", row.names = F)

######generate a list of high confidence proteins#####
omega <- cmplt13 %>% 
  select(c("prot", "LPF", "LPF.riBAQ", "fold_enrichment_AHA", "fold_enrichment_HPG", "tagged_unenriched_AHA",
           "tagged_unenriched_HPG", "meanriBAQ_enriched_AHA", "meanriBAQ_enriched_HPG", "header",
           "bincode", "name", "description")) %>% 
  filter((is.na(LPF) == F & is.na(meanriBAQ_enriched_HPG) == F) | (is.na(LPF) == F & is.na(meanriBAQ_enriched_AHA) == F) |
           is.na(LPF) == F & tagged_unenriched_HPG == "Y") %>%
  distinct()
omega <- omega[order(omega$LPF.riBAQ, decreasing = T), ]
omega <- omega %>% 
  select(-c("LPF", "LPF.riBAQ")) %>% 
  distinct()
write.csv(omega, "high_confidence_proteins.csv", row.names = F)
count <- omega %>% 
  select("prot") %>% 
  distinct() %>% 
  NROW()

omega2 <- omega %>% 
  separate(name, into = c("group", "group2", "group3", "group4", "group5", "group6", "group7"), sep = "\\.") %>% 
  separate(bincode, into = c("bin", "bin2", "bin3", "bin4", "bin5", "bin6", "bin7"), sep = "\\.")
omega2 <- omega2 %>% 
  select(c("prot", "tagged_unenriched_HPG", "meanriBAQ_enriched_AHA", "meanriBAQ_enriched_HPG",
           "group", "group2", "group3", "group4", "group5", "group6", "group7",
           "bin", "bin2", "bin3", "bin4", "bin5", "bin6", "bin7")) %>% 
  distinct()
omega2 <- omega2[order(omega2$group),]
count2 <- omega2 %>% 
  select("prot") %>% 
  distinct() %>% 
  NROW()
#Total number of proteins in groups
a <- filter(omega2, group == "Amino acid metabolism")
a <- length(unique(a$prot))
b <- filter(omega2, group == "Carbohydrate metabolism")
b <- length(unique(b$prot))
c <- filter(omega2, group == "Lipid metabolism")
c <- length(unique(c$prot))
d <- filter(omega2, group == "Cellular respiration")
d <- length(unique(d$prot))
e <- filter(omega2, group == "Nucleotide metabolism" | group == "RNA processing")
e <- length(unique(e$prot))
f <- filter(omega2, group == "Coenzyme metabolism" | group == "Enzyme classification")
f <- length(unique(f$prot))
g <- filter(omega2, group == "Chromatin organisation" | group == "Cytoskeleton organisation" | group == "Cell wall organisation")
g <- length(unique(g$prot))
n <- filter(omega2, group == "not assigned")
n <- length(unique(n$prot))
m <- filter(omega2, group == "Nutrient uptake")
m <- length(unique(m$prot))
p <- filter(omega2, group == "Photosynthesis")
p <- length(unique(p$prot))
q <- filter(omega2, str_detect(group2, "ribosome")) 
q <- length(unique(q$prot))
quell <- filter(omega2, str_detect(group2, "translation") | str_detect(group3, "translation")) 
quell <- length(unique(quell$prot))
x <- filter(omega2, group2 == "proteolysis" | group2 == "ubiquitin-proteasome system")
x <- length(unique(x$prot))
xi <- filter(omega2, group2 == "protein quality control")
xi <- length(unique(xi$prot))
y <- filter(omega2, group == "Protein modification")
y <- length(unique(y$prot))
z <- filter(omega2, group == "Protein translocation")
z <- length(unique(z$prot))
r <- filter(omega2, group == "Redox homeostasis")
r <- length(unique(r$prot))
t <- filter(omega2, group == "Secondary metabolism")
t <- length(unique(t$prot))
u <- filter(omega2, group == "Solute transport")
u <- length(unique(u$prot))
v <- filter(omega2, group == "Vesicle trafficking")
v <- length(unique(v$prot))

#Number of proteins Enriched in AHA set
omega3 <- filter(omega2, is.na(meanriBAQ_enriched_AHA) == F)
a2 <- filter(omega3, group == "Amino acid metabolism")
a2 <- length(unique(a2$prot))
b2 <- filter(omega3, group == "Carbohydrate metabolism")
b2 <- length(unique(b2$prot))
c2 <- filter(omega3, group == "Lipid metabolism")
c2 <- length(unique(c2$prot))
d2 <- filter(omega3, group == "Cellular respiration")
d2 <- length(unique(d2$prot))
e2 <- filter(omega3, group == "Nucleotide metabolism" | group == "RNA processing")
e2 <- length(unique(e2$prot))
f2 <- filter(omega3, group == "Coenzyme metabolism" | group == "Enzyme classification")
f2 <- length(unique(f2$prot))
g2 <- filter(omega3, group == "Chromatin organisation" | group == "Cytoskeleton organisation" | group == "Cell wall organisation")
g2 <- length(unique(g2$prot))
n2 <- filter(omega3, group == "not assigned")
n2 <- length(unique(n2$prot))
m2 <- filter(omega3, group == "Nutrient uptake")
m2 <- length(unique(m2$prot))
p2 <- filter(omega3, group == "Photosynthesis")
p2 <- length(unique(p2$prot))
q2 <- filter(omega3, str_detect(group2, "ribosome"))
q2 <- length(unique(q2$prot))
quell2 <- filter(omega3, str_detect(group2, "translation") | str_detect(group3, "translation")) 
quell2 <- length(unique(quell2$prot))
x2 <- filter(omega3, group2 == "proteolysis" | group2 == "ubiquitin-proteasome system")
x2 <- length(unique(x2$prot))
xi2 <- filter(omega3, group2 == "protein quality control")
xi2 <- length(unique(xi2$prot))
y2 <- filter(omega3, group == "Protein modification")
y2 <- length(unique(y2$prot))
z2 <- filter(omega3, group == "Protein translocation")
z2 <- length(unique(z2$prot))
r2 <- filter(omega3, group == "Redox homeostasis")
r2 <- length(unique(r2$prot))
t2 <- filter(omega3, group == "Secondary metabolism")
t2 <- length(unique(t2$prot))
u2 <- filter(omega3, group == "Solute transport")
u2 <- length(unique(u2$prot))
v2 <- filter(omega3, group == "Vesicle trafficking")
v2 <- length(unique(v2$prot))

#Number of proteins Enriched in HPG set
omega4 <- filter(omega2, is.na(meanriBAQ_enriched_HPG) == F)
a3 <- filter(omega4, group == "Amino acid metabolism")
a3 <- length(unique(a3$prot))
b3 <- filter(omega4, group == "Carbohydrate metabolism")
b3 <- length(unique(b3$prot))
c3 <- filter(omega4, group == "Lipid metabolism")
c3 <- length(unique(c3$prot))
d3 <- filter(omega4, group == "Cellular respiration")
d3 <- length(unique(d3$prot))
e3 <- filter(omega4, group == "Nucleotide metabolism" | group == "RNA processing")
e3 <- length(unique(e3$prot))
f3 <- filter(omega4, group == "Coenzyme metabolism" | group == "Enzyme classification")
f3 <- length(unique(f3$prot))
g3 <- filter(omega4, group == "Chromatin organisation" | group == "Cytoskeleton organisation" | group == "Cell wall organisation")
g3 <- length(unique(g3$prot))
n3 <- filter(omega4, group == "not assigned")
n3 <- length(unique(n3$prot))
m3 <- filter(omega4, group == "Nutrient uptake")
m3 <- length(unique(m3$prot))
p3 <- filter(omega4, group == "Photosynthesis")
p3 <- length(unique(p3$prot))
q3 <- filter(omega4, str_detect(group2, "ribosome"))
q3 <- length(unique(q3$prot))
quell3 <- filter(omega4, str_detect(group2, "translation") | str_detect(group3, "translation")) 
quell3 <- length(unique(quell3$prot))
x3 <- filter(omega4, group2 == "proteolysis" | group2 == "ubiquitin-proteasome system")
x3 <- length(unique(x3$prot))
xi3 <- filter(omega4, group2 == "protein quality control")
xi3 <- length(unique(xi3$prot))
y3 <- filter(omega4, group == "Protein modification")
y3 <- length(unique(y3$prot))
z3 <- filter(omega4, group == "Protein translocation")
z3 <- length(unique(z3$prot))
r3 <- filter(omega4, group == "Redox homeostasis")
r3 <- length(unique(r3$prot))
t3 <- filter(omega4, group == "Secondary metabolism")
t3 <- length(unique(t3$prot))
u3 <- filter(omega4, group == "Solute transport")
u3 <- length(unique(u3$prot))
v3 <- filter(omega4, group == "Vesicle trafficking")
v3 <- length(unique(v3$prot))

#Number of proteins tagged with HPG
omega5 <- filter(omega2, tagged_unenriched_HPG == "Y")
a4 <- filter(omega5, group == "Amino acid metabolism")
a4 <- length(unique(a4$prot))
b4 <- filter(omega5, group == "Carbohydrate metabolism")
b4 <- length(unique(b4$prot))
c4 <- filter(omega5, group == "Lipid metabolism")
c4 <- length(unique(c4$prot))
d4 <- filter(omega5, group == "Cellular respiration")
d4 <- length(unique(d4$prot))
e4 <- filter(omega5, group == "Nucleotide metabolism" | group == "RNA processing")
e4 <- length(unique(e4$prot))
f4 <- filter(omega5, group == "Coenzyme metabolism" | group == "Enzyme classification")
f4 <- length(unique(f4$prot))
g4 <- filter(omega5, group == "Chromatin organisation" | group == "Cytoskeleton organisation" | group == "Cell wall organisation")
g4 <- length(unique(g4$prot))
n4 <- filter(omega5, group == "not assigned")
n4 <- length(unique(n4$prot))
m4 <- filter(omega5, group == "Nutrient uptake")
m4 <- length(unique(m4$prot))
p4 <- filter(omega5, group == "Photosynthesis")
p4 <- length(unique(p4$prot))
q4 <- filter(omega5, str_detect(group2, "ribosome"))
q4 <- length(unique(q4$prot))
quell4 <- filter(omega5, str_detect(group2, "translation") | str_detect(group3, "translation")) 
quell4 <- length(unique(quell4$prot))
x4 <- filter(omega5, group2 == "proteolysis" | group2 == "ubiquitin-proteasome system")
x4 <- length(unique(x4$prot))
xi4 <- filter(omega5, group2 == "protein quality control")
xi4 <- length(unique(xi4$prot))
y4 <- filter(omega5, group == "Protein modification")
y4 <- length(unique(y4$prot))
z4 <- filter(omega5, group == "Protein translocation")
z4 <- length(unique(z4$prot))
r4 <- filter(omega5, group == "Redox homeostasis")
r4 <- length(unique(r4$prot))
t4 <- filter(omega5, group == "Secondary metabolism")
t4 <- length(unique(t4$prot))
u4 <- filter(omega5, group == "Solute transport")
u4 <- length(unique(u4$prot))
v4 <- filter(omega5, group == "Vesicle trafficking")
v4 <- length(unique(v4$prot))


highConfCollapsed <- data.frame(Group = c("Amino acid metabolism", "Carbohydrate metabolism", "Cellular respiration",
                                          "Cellular organisation", "Photosynthesis", "Ribosomal subunits", "Translation",
                                          "Proteosome/proteolysis", "Chaperones", "Protein modification/translocation",
                                          "Other","Not assigned"),
                                # MapmanBin = c("4", "3", "2", "12, 20, 21", "1", "17", "17", "19", "19", "18,23",
                                #               "6-7, 9-10, 16, 22, 24-25, 50", "35"),
                                TotalNumberOfProteins = c(a, b, d, g, p, q, quell, x, xi, y + z, sum(c, e, f, m, r, t, u, v) , n),
                                NumberOfProteins_AHA_enriched = c(a2, b2, d2, g2, p2, q2, quell2, x2, xi2, y2 + z2, sum(c2, e2, f2, m2, r2, t2, u2, v2), n2),
                                NumberOfProteins_HPG_enriched = c(a3, b3, d3, g3, p3, q3, quell3, x3, xi3, y3 + z3, sum(c3, e3, f3, m3, r3, t3, u3, v3), n3),
                                NumberOfProteins_HPG_tagged = c(a4, b4, d4, g4, p4, q4, quell4, x4, xi, y4 + z4, sum(c4, e4, f4, m4, r4, t4, u4, v4), n4))
write.csv(highConfCollapsed, "high_confidence_proteins_collapsed.csv", row.names = F)

######does enrichment correlate with particular domains?######
hEnri3 <- hEnri %>% 
  dplyr::rename("prot" = "firstID") %>% 
  left_join(go_domains, by = "prot")
hEnri4 <- hEnri3 %>%
  filter(is.na(PO.domain) == F, is.na(GO.domain) == F, is.na(meanriBAQ_enriched) == F)
hEnri_PO_anatomy <- filter(hEnri4, PO.domain == "plant_anatomy")
hEnri_PO_anatomy <- length(unique(hEnri_PO_anatomy$prot))
hEnri_PO_dev <- filter(hEnri4, PO.domain == "plant_structure_development_stage")
hEnri_PO_dev <- length(unique(hEnri_PO_dev$prot))
hEnri_GO_proc <- filter(hEnri4, GO.domain == "biological_process")
hEnri_GO_proc <- length(unique(hEnri_GO_proc$prot))
hEnri_GO_comp <- filter(hEnri4, GO.domain == "cellular_component")
hEnri_GO_comp <- length(unique(hEnri_GO_comp$prot))
hEnri_GO_go <- filter(hEnri4, GO.domain == "go")
hEnri_GO_go <- length(unique(hEnri_GO_go$prot))
hEnri_GO_func <- filter(hEnri4, GO.domain == "molecular_function")
hEnri_GO_func <- length(unique(hEnri_GO_func$prot))

aEnri3 <- aEnri %>% 
  dplyr::rename("prot" = "firstID") %>% 
  left_join(go_domains, by = "prot")
aEnri4 <- aEnri3 %>% 
  filter(is.na(PO.domain) == F, is.na(GO.domain) == F, is.na(meanriBAQ_enriched) == F)
aEnri_PO_anatomy <- filter(aEnri4, PO.domain == "plant_anatomy")
aEnri_PO_anatomy <- length(unique(aEnri_PO_anatomy$prot))
aEnri_PO_dev <- filter(aEnri4, PO.domain == "plant_structure_development_stage")
aEnri_PO_dev <- length(unique(aEnri_PO_dev$prot))
aEnri_GO_proc <- filter(aEnri4, GO.domain == "biological_process")
aEnri_GO_proc <- length(unique(aEnri_GO_proc$prot))
aEnri_GO_comp <- filter(aEnri4, GO.domain == "cellular_component")
aEnri_GO_comp <- length(unique(aEnri_GO_comp$prot))
aEnri_GO_go <- filter(aEnri4, GO.domain == "go")
aEnri_GO_go <- length(unique(aEnri_GO_go$prot))
aEnri_GO_func <- filter(aEnri4, GO.domain == "molecular_function")
aEnri_GO_func <- length(unique(aEnri_GO_func$prot))

hEnri5 <- hEnri3 %>% 
  filter(is.na(PO.domain) == F, is.na(GO.domain) == F)
hEnri_PO_anatomy_background <- filter(hEnri5, PO.domain == "plant_anatomy")
hEnri_PO_anatomy_background <- length(unique(hEnri_PO_anatomy_background$prot))
hEnri_PO_dev_background <- filter(hEnri5, PO.domain == "plant_structure_development_stage")
hEnri_PO_dev_background <- length(unique(hEnri_PO_dev_background$prot))
hEnri_GO_proc_background <- filter(hEnri5, GO.domain == "biological_process")
hEnri_GO_proc_background <- length(unique(hEnri_GO_proc_background$prot))
hEnri_GO_comp_background <- filter(hEnri5, GO.domain == "cellular_component")
hEnri_GO_comp_background <- length(unique(hEnri_GO_comp_background$prot))
hEnri_GO_go_background <- filter(hEnri5, GO.domain == "go")
hEnri_GO_go_background <- length(unique(hEnri_GO_go_background$prot))
hEnri_GO_func_background <- filter(hEnri5, GO.domain == "molecular_function")
hEnri_GO_func_background <- length(unique(hEnri_GO_func_background$prot))

aEnri5 <- aEnri3 %>% 
  filter(is.na(PO.domain) == F, is.na(GO.domain) == F)
aEnri_PO_anatomy_background <- filter(aEnri5, PO.domain == "plant_anatomy")
aEnri_PO_anatomy_background <- length(unique(aEnri_PO_anatomy_background$prot))
aEnri_PO_dev_background <- filter(aEnri5, PO.domain == "plant_structure_development_stage")
aEnri_PO_dev_background <- length(unique(aEnri_PO_dev_background$prot))
aEnri_GO_proc_background <- filter(aEnri5, GO.domain == "biological_process")
aEnri_GO_proc_background <- length(unique(aEnri_GO_proc_background$prot))
aEnri_GO_comp_background <- filter(aEnri5, GO.domain == "cellular_component")
aEnri_GO_comp_background <- length(unique(aEnri_GO_comp_background$prot))
aEnri_GO_go_background <- filter(aEnri5, GO.domain == "go")
aEnri_GO_go_background <- length(unique(aEnri_GO_go_background$prot))
aEnri_GO_func_background <- filter(aEnri5, GO.domain == "molecular_function")
aEnri_GO_func_background <- length(unique(aEnri_GO_func_background$prot))

#does tagging correlate with particular domains?
hEnri6 <- hEnri3 %>% 
  filter(is.na(PO.domain) == F, is.na(GO.domain) == F, tagged_unenriched == "Y")
hEnri_PO_anatomy2 <- filter(hEnri6, PO.domain == "plant_anatomy")
hEnri_PO_anatomy2 <- length(unique(hEnri_PO_anatomy2$prot))
hEnri_PO_dev2 <- filter(hEnri6, PO.domain == "plant_structure_development_stage")
hEnri_PO_dev2 <- length(unique(hEnri_PO_dev2$prot))
hEnri_GO_proc2 <- filter(hEnri6, GO.domain == "biological_process")
hEnri_GO_proc2 <- length(unique(hEnri_GO_proc2$prot))
hEnri_GO_comp2 <- filter(hEnri6, GO.domain == "cellular_component")
hEnri_GO_comp2 <- length(unique(hEnri_GO_comp2$prot))
hEnri_GO_go2 <- filter(hEnri6, GO.domain == "go")
hEnri_GO_go2 <- length(unique(hEnri_GO_go2$prot))
hEnri_GO_func2 <- filter(hEnri6, GO.domain == "molecular_function")
hEnri_GO_func2 <- length(unique(hEnri_GO_func2$prot))

#summary of GO-PO enrichment analysis
domains <- data.frame(domain = c("PO_anatomy", "PO_dev", "GO_proc", "GO_comp","GO_func"),
                      aha_foreground_enriched = c(aEnri_PO_anatomy, aEnri_PO_dev, aEnri_GO_proc, aEnri_GO_comp,
                                                  aEnri_GO_func),
                      aha_background = c(aEnri_PO_anatomy_background, aEnri_PO_dev_background,
                                         aEnri_GO_proc_background, aEnri_GO_comp_background,
                                         aEnri_GO_func_background),
                      hpg_foreground_enriched = c(hEnri_PO_anatomy, hEnri_PO_dev, hEnri_GO_proc, hEnri_GO_comp,
                                                  hEnri_GO_func),
                      hpg_foreground_tagged = c(hEnri_PO_anatomy2, hEnri_PO_dev2, hEnri_GO_proc2, hEnri_GO_comp2,
                                                hEnri_GO_func2),
                      hpg_background = c(hEnri_PO_anatomy_background, hEnri_PO_dev_background,
                                         hEnri_GO_proc_background, hEnri_GO_comp_background,
                                         hEnri_GO_func_background)) %>% 
  mutate(aha_foreground_enriched_norm = aha_foreground_enriched/min(aha_foreground_enriched), aha_background_norm = aha_background/min(aha_background),
         hpg_foreground_enriched_norm = hpg_foreground_enriched/min(hpg_foreground_enriched), hpg_foreground_tagged_norm = hpg_foreground_tagged/min(hpg_foreground_tagged),
         hpg_background_norm = hpg_background/min(hpg_background))
write.csv(domains, "GO-PO_analysis.csv", row.names = F)

######Top 40 enriched proteins#####
bilbo <- cmplt13 %>% 
  select(c("prot", "n_pep_15N", "LPF", "LPF.riBAQ", "tagged_unenriched_HPG", "fold_enrichment_HPG", "meanriBAQ_enriched_HPG",
           "tagged_unenriched_AHA", "fold_enrichment_AHA", "meanriBAQ_enriched_AHA")) %>% 
  distinct()
bilbo$LPF.riBAQ <- as.numeric(bilbo$LPF.riBAQ)
bilbo <- bilbo[order(bilbo$LPF.riBAQ, decreasing = T),]
frodo <- data.frame(prot = NA, LPF = NA, LPF.riBAQ = NA, tagged_unenriched_HPG = NA, fold_enrichment_HPG = NA,
                    meanriBAQ_enriched_HPG = NA, tagged_unenriched_AHA = NA, fold_enrichment_AHA = NA,
                    meanriBAQ_enriched_AHA = NA)
gollum <- unique(bilbo$prot)
for(i in 1:length(gollum)) {
  sub_data <- bilbo %>% 
    filter(prot == gollum[i]) %>% 
    filter(LPF.riBAQ == max(LPF.riBAQ))
  if (i == 1) {
    frodo <- sub_data
  } else {
    frodo <- rbind(frodo, sub_data)
  }
}
frodo <- left_join(frodo, mapman, by = "prot")
frodo$name <- sub("\\..*", "", frodo$name)
frodo$bincode <- sub("\\..*", "", frodo$bincode)
frodo <- frodo %>% 
  select(-"bincode") %>% 
  distinct()
write.csv(frodo, "top_40.csv", row.names = F)
