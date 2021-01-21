#Set working directory
setwd("C:/temp/noncanonical_amino_acids_r_processing/modification_searches")

#Load packages
library(gridExtra)
library(utils)
library(tidyverse)
library(RNHANES)
library(stringi)
options(sringsAsFactors = F, na.strings = c("", "NA"))
#Import data files####
ahae1 <- read.delim("proteinAHAEnriched_1.tsv", na.strings = c("", "NA"))
ahae2 <- read.delim("proteinAHAEnriched_2.tsv", na.strings = c("", "NA"))
ahae3 <- read.delim("proteinAHAEnriched_3.tsv", na.strings = c("", "NA"))
ahau1 <- read.delim("proteinAHAUnenriched_1.tsv", na.strings = c("", "NA"))
ahau2 <- read.delim("proteinAHAUnenriched_2.tsv", na.strings = c("", "NA"))
ahau3 <- read.delim("proteinAHAUnenriched_3.tsv", na.strings = c("", "NA"))
hpge1 <- read.delim("proteinHPGEnriched_1.tsv", na.strings = c("", "NA"))
hpge2 <- read.delim("proteinHPGEnriched_2.tsv", na.strings = c("", "NA"))
hpge3 <- read.delim("proteinHPGEnriched_3.tsv", na.strings = c("", "NA"))
hpgu1 <- read.delim("proteinHPGUnenriched_1.tsv", na.strings = c("", "NA"))
hpgu2 <- read.delim("proteinHPGUnenriched_2.tsv", na.strings = c("", "NA"))
hpgu3 <- read.delim("proteinHPGUnenriched_3.tsv", na.strings = c("", "NA"))

ahau1f1 <- read.delim("NT-19-47-1_f1.tsv", na.strings = c("", "NA"))
ahau1f2 <- read.delim("NT-19-47-1_f2.tsv", na.strings = c("", "NA"))
ahau1f3 <- read.delim("NT-19-47-1_f3.tsv", na.strings = c("", "NA"))
ahau1f4 <- read.delim("NT-19-47-1_f4.tsv", na.strings = c("", "NA"))
ahau1f5 <- read.delim("NT-19-47-1_f5.tsv", na.strings = c("", "NA"))
ahau1f6 <- read.delim("NT-19-47-1_f6.tsv", na.strings = c("", "NA"))
ahau1f7 <- read.delim("NT-19-47-1_f7.tsv", na.strings = c("", "NA"))
ahau1f8 <- read.delim("NT-19-47-1_f8.tsv", na.strings = c("", "NA"))
ahau1f9 <- read.delim("NT-19-47-1_f9.tsv", na.strings = c("", "NA"))
ahau1f10 <- read.delim("NT-19-47-1_f10.tsv", na.strings = c("", "NA"))
ahau1f11 <- read.delim("NT-19-47-1_f11.tsv", na.strings = c("", "NA"))
ahau1f12 <- read.delim("NT-19-47-1_f12.tsv", na.strings = c("", "NA"))

ahau2f1 <- read.delim("NT-19-47-4_f1.tsv", na.strings = c("", "NA"))
ahau2f2 <- read.delim("NT-19-47-4_f2.tsv", na.strings = c("", "NA"))
ahau2f3 <- read.delim("NT-19-47-4_f3.tsv", na.strings = c("", "NA"))
ahau2f4 <- read.delim("NT-19-47-4_f4.tsv", na.strings = c("", "NA"))
ahau2f5 <- read.delim("NT-19-47-4_f5.tsv", na.strings = c("", "NA"))
ahau2f6 <- read.delim("NT-19-47-4_f6.tsv", na.strings = c("", "NA"))
ahau2f7 <- read.delim("NT-19-47-4_f7.tsv", na.strings = c("", "NA"))
ahau2f8 <- read.delim("NT-19-47-4_f8.tsv", na.strings = c("", "NA"))
ahau2f9 <- read.delim("NT-19-47-4_f9.tsv", na.strings = c("", "NA"))
ahau2f10 <- read.delim("NT-19-47-4_f10.tsv", na.strings = c("", "NA"))
ahau2f11 <- read.delim("NT-19-47-4_f11.tsv", na.strings = c("", "NA"))
ahau2f12 <- read.delim("NT-19-47-4_f12.tsv", na.strings = c("", "NA"))

ahau3f1 <- read.delim("NT-19-47-7_f1.tsv", na.strings = c("", "NA"))
ahau3f2 <- read.delim("NT-19-47-7_f2.tsv", na.strings = c("", "NA"))
ahau3f3 <- read.delim("NT-19-47-7_f3.tsv", na.strings = c("", "NA"))
ahau3f4 <- read.delim("NT-19-47-7_f4.tsv", na.strings = c("", "NA"))
ahau3f5 <- read.delim("NT-19-47-7_f5.tsv", na.strings = c("", "NA"))
ahau3f6 <- read.delim("NT-19-47-7_f6.tsv", na.strings = c("", "NA"))
ahau3f7 <- read.delim("NT-19-47-7_f7.tsv", na.strings = c("", "NA"))
ahau3f8 <- read.delim("NT-19-47-7_f8.tsv", na.strings = c("", "NA"))
ahau3f9 <- read.delim("NT-19-47-7_f9.tsv", na.strings = c("", "NA"))
ahau3f10 <- read.delim("NT-19-47-7_f10.tsv", na.strings = c("", "NA"))
ahau3f11 <- read.delim("NT-19-47-7_f11.tsv", na.strings = c("", "NA"))
ahau3f12 <- read.delim("NT-19-47-7_f12.tsv", na.strings = c("", "NA"))

hpgu1f1 <- read.delim("NT-19-47-10_f1.tsv", na.strings = c("", "NA"))
hpgu1f2 <- read.delim("NT-19-47-10_f2.tsv", na.strings = c("", "NA"))
hpgu1f3 <- read.delim("NT-19-47-10_f3.tsv", na.strings = c("", "NA"))
hpgu1f4 <- read.delim("NT-19-47-10_f4.tsv", na.strings = c("", "NA"))
hpgu1f5 <- read.delim("NT-19-47-10_f5.tsv", na.strings = c("", "NA"))
hpgu1f6 <- read.delim("NT-19-47-10_f6.tsv", na.strings = c("", "NA"))
hpgu1f7 <- read.delim("NT-19-47-10_f7.tsv", na.strings = c("", "NA"))
hpgu1f8 <- read.delim("NT-19-47-1_f8.tsv", na.strings = c("", "NA"))
hpgu1f9 <- read.delim("NT-19-47-10_f9.tsv", na.strings = c("", "NA"))
hpgu1f10 <- read.delim("NT-19-47-10_f10.tsv", na.strings = c("", "NA"))
hpgu1f11 <- read.delim("NT-19-47-10_f11.tsv", na.strings = c("", "NA"))
hpgu1f12 <- read.delim("NT-19-47-10_f12.tsv", na.strings = c("", "NA"))

hpgu2f1 <- read.delim("NT-19-47-13_f1.tsv", na.strings = c("", "NA"))
hpgu2f2 <- read.delim("NT-19-47-13_f2.tsv", na.strings = c("", "NA"))
hpgu2f3 <- read.delim("NT-19-47-13_f3.tsv", na.strings = c("", "NA"))
hpgu2f4 <- read.delim("NT-19-47-13_f4.tsv", na.strings = c("", "NA"))
hpgu2f5 <- read.delim("NT-19-47-13_f5.tsv", na.strings = c("", "NA"))
hpgu2f6 <- read.delim("NT-19-47-13_f6.tsv", na.strings = c("", "NA"))
hpgu2f7 <- read.delim("NT-19-47-13_f7.tsv", na.strings = c("", "NA"))
hpgu2f8 <- read.delim("NT-19-47-13_f8.tsv", na.strings = c("", "NA"))
hpgu2f9 <- read.delim("NT-19-47-13_f9.tsv", na.strings = c("", "NA"))
hpgu2f10 <- read.delim("NT-19-47-13_f10.tsv", na.strings = c("", "NA"))
hpgu2f11 <- read.delim("NT-19-47-13_f11.tsv", na.strings = c("", "NA"))
hpgu2f12 <- read.delim("NT-19-47-13_f12.tsv", na.strings = c("", "NA"))

hpgu3f1 <- read.delim("NT-19-47-16_f1.tsv", na.strings = c("", "NA"))
hpgu3f2 <- read.delim("NT-19-47-16_f2.tsv", na.strings = c("", "NA"))
hpgu3f3 <- read.delim("NT-19-47-16_f3.tsv", na.strings = c("", "NA"))
hpgu3f4 <- read.delim("NT-19-47-16_f4.tsv", na.strings = c("", "NA"))
hpgu3f5 <- read.delim("NT-19-47-16_f5.tsv", na.strings = c("", "NA"))
hpgu3f6 <- read.delim("NT-19-47-16_f6.tsv", na.strings = c("", "NA"))
hpgu3f7 <- read.delim("NT-19-47-16_f7.tsv", na.strings = c("", "NA"))
hpgu3f8 <- read.delim("NT-19-47-16_f8.tsv", na.strings = c("", "NA"))
hpgu3f9 <- read.delim("NT-19-47-16_f9.tsv", na.strings = c("", "NA"))
hpgu3f10 <- read.delim("NT-19-47-16_f10.tsv", na.strings = c("", "NA"))
hpgu3f11 <- read.delim("NT-19-47-16_f11.tsv", na.strings = c("", "NA"))
hpgu3f12 <- read.delim("NT-19-47-16_f12.tsv", na.strings = c("", "NA"))

ahae1f1 <- read.delim("NT-20-21-23A.tsv", na.strings = c("", "NA"))
ahae2f1 <- read.delim("NT-20-21-24A.tsv", na.strings = c("", "NA"))
ahae3f1 <- read.delim("NT-20-92-2A.tsv", na.strings = c("", "NA"))

hpge1f1 <- read.delim("NT-20-21-25A.tsv", na.strings = c("", "NA"))
hpge2f1 <- read.delim("NT-20-21-26A.tsv", na.strings = c("", "NA"))
hpge3f1 <- read.delim("NT-20-21-27A.tsv", na.strings = c("", "NA"))

alpha <- rbind(ahau1f1, ahau1f2, ahau1f3, ahau1f4, ahau1f5, ahau1f6, ahau1f7, ahau1f8, ahau1f9, ahau1f10, ahau1f11, ahau1f12)
beta <- rbind(ahau2f1, ahau2f2, ahau2f3, ahau2f4, ahau2f5, ahau2f6, ahau2f7, ahau2f8, ahau2f9, ahau2f10, ahau2f11, ahau2f12)
gamma <- rbind(ahau3f1, ahau3f2, ahau3f3, ahau3f4, ahau3f5, ahau3f6, ahau3f7, ahau3f8, ahau3f9, ahau3f10, ahau3f11, ahau3f12)
delta <- rbind(hpgu1f1, hpgu1f2, hpgu1f3, hpgu1f4, hpgu1f5, hpgu1f6, hpgu1f7, hpgu1f8, hpgu1f9, hpgu1f10, hpgu1f11, hpgu1f12)
epsilon <- rbind(hpgu2f1, hpgu2f2, hpgu2f3, hpgu2f4, hpgu2f5, hpgu2f6, hpgu2f7, hpgu2f8, hpgu2f9, hpgu2f10, hpgu2f11, hpgu2f12)
zeta <- rbind(hpgu3f1, hpgu3f2, hpgu3f3, hpgu3f4, hpgu3f5, hpgu3f6, hpgu3f7, hpgu3f8, hpgu3f9, hpgu3f10, hpgu3f11, hpgu3f12)

#Functions
#clean selects only the relevant columns, filters out proteins not from Arabidopsis, adds an extra column (Present)
#and renames the mod column####
clean <- function(x){
  x <- x %>% 
    filter(Total.Peptide.Ions > 0) %>% 
    select("Protein") %>% 
    filter(str_detect(Protein, pattern = "AT")) %>% 
    mutate(Present = "Y") %>% 
    distinct()
}
#triple_join joins three dataframes together and renames variables from the third so that they have the correct suffix####
triple_join <- function(w, x, y){
  z <- full_join(w, x, by = "Protein", suffix = c("_1", "_2"))
  z <- full_join(z, y, by = "Protein") %>% 
    rename("Mod_3" = "Mod", "Present_3" = "Present", "percentMod_3" = "percentMod", "n_mod_hit_3" = "n_mod_hit",
           "n_hit_3" = "n_hit") %>% 
    distinct()
}
#present_in_two filters out observations that are present in one sample and then creates a new column noting if modifications
#are present in two out of three observations
present_in_two <- function(x){
  x <- x %>% 
    filter(Present_1 == "Y" & Present_2 == "Y" |
           Present_1 == "Y" & Present_3 == "Y" |
           Present_2 == "Y" & Present_3 == "Y") %>% 
    mutate(modInTwo = NA, Present = "Y") %>% 
    select(-c("Present_1", "Present_2", "Present_3"))
  for(i in 1:nrow(x)){
      if((is.na(x$Mod_1[i]) == F & is.na(x$Mod_2[i]) == F) |
       (is.na(x$Mod_1[i]) == F & is.na(x$Mod_3[i]) == F) |
       (is.na(x$Mod_2[i])== F & is.na(x$Mod_3[i])== F)) {
      x$modInTwo[i] <- "Y"
    } else {
      x$modInTwo[i] <- "N"
    }
  }
  x <- x %>% 
    unite("Mods", c(Mod_1, Mod_2, Mod_3), sep = " ", remove = T, na.rm = T) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = stri_replace_last(Mods, "", fixed = ",")) %>%
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ",")) %>% 
    mutate(Mods = str_replace_all(Mods, ",,", ","))
  }
#count_mods sums the total number of modifications for each protein across all three replicates, regardless of whether
#the modification was found in two out of three####
count_mods <- function(x){
  x <- x %>% 
    mutate(Total_mods_across_3_reps_U = NA, Total_mods_across_3_reps_E = NA)
  for(i in 1:nrow(x)) {
    if(is.na(x$Mods_u[i]) == F) {
      x$Total_mods_across_3_reps_U[i] <- (length(unlist(strsplit(as.character(x$Mods_u[i]), ","))))
    } else {
      x$Total_mods_across_3_reps_U[i] <- 0
    }
  }
  for(i in 1:nrow(x)) {
    if(is.na(x$Mods_e[i]) == F) {
      x$Total_mods_across_3_reps_E[i] <- (length(unlist(strsplit(as.character(x$Mods_e[i]), ","))))
    } else
      x$Total_mods_across_3_reps_E[i] <- 0
  }
  x
}
#percent_mod calculates the percentage of spectral hits for each protein that are modified####
percent_mod <- function(x, y){
  #Filter out contaminants/decoys, remove all but three columns, modify the protein name to just containe the AGI
  x <- x %>% 
    filter(startsWith(protein, "A")) %>% 
    select(c("protein", "peptide", "modification_info")) %>% 
    mutate(protein = str_sub(protein, end = 11)) %>% 
    mutate(modification_info = str_replace_all(modification_info, ", \\d+C\\(57.021465\\)", "")) %>%
    mutate(modification_info = str_replace_all(modification_info, "\\d+C\\(57.021465\\), ", "")) %>%
    mutate(modification_info = str_replace_all(modification_info, "\\d+C\\(57.021465\\)", "")) %>% 
    mutate(modification_info = str_replace_all(modification_info, ", \\d+M\\(-4.9863\\)", "")) %>%
    mutate(modification_info = str_replace_all(modification_info, "\\d+M\\(-4.9863\\), ", "")) %>%
    mutate(modification_info = str_replace_all(modification_info, "\\d+M\\(-4.9863\\)", "")) %>% 
    mutate(modification_info = str_replace_all(modification_info, ", \\d+M\\(-21.9877\\)", "")) %>%
    mutate(modification_info = str_replace_all(modification_info, "\\d+M\\(-21.9877\\), ", "")) %>%
    mutate(modification_info = str_replace_all(modification_info, "\\d+M\\(-21.9877\\)", ""))
  
  #Create a vector that contains all of the AGIs from the PeptideProphe/ProteinProphet checked file
  prot <- as.vector(unique(y$Protein))
  #Create two empty dataframes
  data_hold <- data.frame(protein = NA, n_mod_hit = NA, n_hit = NA, Mod = NA)
  data_hold2 <- data.frame(protein = NA, n_mod_hit = NA, n_hit = NA, Mod = NA)
  #For each AGI in prot filter the unchecked dataframe, count the number of rows to get the number of spectral hits,
  #filter the filtered dataframe to get only spectra that are modified and then count the number of rows to get the 
  #modified spectral hits. Put this data into one of the empty dataframes and rbind to the second new dataframe.
  for(i in 1:length(prot)){
    sub_data <- filter(x, protein == prot[i])
    n_spec <- NROW(sub_data)
    sub_sub_data <- filter(sub_data, is.na(modification_info) == F &
                             modification_info != "")
    n_mod_spec <- NROW(sub_sub_data)
    if(n_mod_spec == 0){
      mod_info <- NA
    } else {
      mod_info <- as.vector(stri_paste(sub_sub_data$modification_info, sep = ", ", collapse = ", ", ignore_null = T)) %>% 
        str_split(pattern = ",")
      mod_info <- unlist(mod_info)
      mod_info <- str_c(mod_info, sep = ",", collapse = ",") %>%
        str_replace_all("\\)", "\\), ") %>%
        str_replace_all(", , , , , ,", ",") %>%
        str_replace_all(", , , , ,", ",") %>%
        str_replace_all(", , , ,", ",") %>%
        str_replace_all(", , ,", ",") %>%
        str_replace_all(", ,", ",") %>%
        str_replace_all("\\)", "\\), ") %>%
        str_replace_all(", , ", ", ") %>%
        str_replace_all("\\) ", "\\), ") %>%
        str_replace_all("\\)", "\\),") %>%
        str_replace_all(",,", ",") %>%
        str_trim(side = "right") %>%
        str_trim(side = "left") %>%
        stri_replace_last("", fixed = ",") %>%
        stri_replace_first("", fixed = ",") %>%
        str_trim(side = "left") %>% 
        str_replace_all("\\) ", "\\), ") %>% 
        str_replace_all("\\)", "\\),")
    }
    if(i == 1) {
      data_hold$protein <- prot[i]
      data_hold$n_mod_hit <- n_mod_spec
      data_hold$n_hit <- n_spec
      data_hold$Mod <- mod_info
    } else {
      data_hold2$protein <- prot[i]
      data_hold2$n_mod_hit <- n_mod_spec
      data_hold2$n_hit <- n_spec
      data_hold2$Mod <- mod_info
      data_hold <- rbind(data_hold, data_hold2)
    }
  }
  #Replace NAs with 0's
  data_hold[c("n_mod_hit", "n_hit")][is.na(data_hold[c("n_mod_hit", "n_hit")])] <- 0
  data_hold$n_mod_hit <- as.numeric(data_hold$n_mod_hit)
  data_hold$n_hit <- as.numeric(data_hold$n_hit)
  #Calculate the percent modified hits
  data_hold <- data_hold %>% 
    mutate(percentMod = round((n_mod_hit/n_hit)*100, digits = 0))
  #Join the percent modif
  y <- left_join(y, data_hold, by = c("Protein" = "protein"))
  for(i in 1:NROW(y)){
    if(is.na(y$Mod[i]) == T & y$n_mod_hit[i] != 0) {
      y$percentMod[i] <- NA
    }
  }
  for(i in 1:NROW(y)){
    if(y$percentMod[i] == "NaN" & is.na(y$Mod[i]) == F){
      y$Mod[i] <- NA
    }
  }
  y
}
#if a mod is only found in two out of three replicates (not necessarily the same mod) it is not considered reliable
#and the percentMod, Mods and n_mod_hit variables are set to NA.
mod_in_two <- function(x){
  for(i in 1:nrow(x)){
  if(x$modInTwo_e[i] == "N" & x$n_mod_hit_e[i] != 0){
    x$percentMod_e[i] <- 0
    x$Mods_e[i] <- NA
    x$n_mod_hit_e[i] <- 0
  }
  if(x$modInTwo_u[i] == "N" & x$n_mod_hit_u[i] != 0){
    x$percentMod_u[i] <- 0
    x$Mods_u[i] <- NA
    x$    n_mod_hit_u[i] <- 0
  }
  }
  x
}

empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

empty_as_na2 <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!=" ", x, NA)
}
  empty_as_na3 <- function(x){
    if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
    ifelse(as.character(x)!="  ", x, NA)  
  }
#####
#For each data file, select only: Protein and Razor Assigned Modifications and remove any proteins that are
#not from Arabidopsis####
ahae1_2 <- clean(ahae1)
ahae1_2 <- percent_mod(ahae1f1, ahae1_2)

ahae2_2 <- clean(ahae2)
ahae2_2 <- percent_mod(ahae2f1, ahae2_2)

ahae3_2 <- clean(ahae3)
ahae3_2 <- percent_mod(ahae3f1, ahae3_2)

ahau1_2 <- clean(ahau1)
ahau1_2 <- percent_mod(alpha, ahau1_2)

ahau2_2 <- clean(ahau2)
ahau2_2 <- percent_mod(beta, ahau2_2)

ahau3_2 <- clean(ahau3)
ahau3_2 <- percent_mod(gamma, ahau3_2)

hpge1_2 <- clean(hpge1)
hpge1_2 <- percent_mod(hpge1f1, hpge1_2)

hpge2_2 <- clean(hpge2)
hpge2_2 <- percent_mod(hpge2f1, hpge2_2)

hpge3_2 <- clean(hpge3)
hpge3_2 <- percent_mod(hpge3f1, hpge3_2)

hpgu1_2 <- clean(hpgu1)
hpgu1_2 <- percent_mod(delta, hpgu1_2)

hpgu2_2 <- clean(hpgu2)
hpgu2_2 <- percent_mod(epsilon, hpgu2_2)

hpgu3_2 <- clean(hpgu3)
hpgu3_2 <- percent_mod(zeta, hpgu3_2)

#Combine the replicates for each treatment and filter out any proteins that appear in only one replicate####
ahae <- triple_join(ahae1_2, ahae2_2, ahae3_2)
ahae <- present_in_two(ahae)
ahau <- triple_join(ahau1_2, ahau2_2, ahau3_2)
ahau <- present_in_two(ahau)
hpge <- triple_join(hpge1_2, hpge2_2, hpge3_2)
hpge <- present_in_two(hpge)
hpgu <- triple_join(hpgu1_2, hpgu2_2, hpgu3_2)
hpgu <- present_in_two(hpgu)

#Full join the AHA enriched table to the AHA unenriched table and the HPG enriched table to the HPG unenriched table####
aha <- full_join(ahau, ahae, by = "Protein", suffix = c("_u", "_e"))
hpg <- full_join(hpgu, hpge, by = "Protein", suffix = c("_u", "_e"))

#Create an output file that shows the total number of proteins with any modification in each dataset and the number of####
#modified proteins as proportion of the total number of proteins in each treatment group####
au <- length(unique(ahau$Protein))
ae <- length(unique(ahae$Protein))
hu <- length(unique(hpgu$Protein))
he <- length(unique(hpge$Protein))

auMod <- filter(ahau, modInTwo == "Y")
auMod <- length(unique(auMod$Protein))
aeMod <- filter(ahae, modInTwo == "Y")
aeMod <- length(unique(aeMod$Protein))
huMod <- filter(hpgu, modInTwo == "Y")
huMod <- length(unique(huMod$Protein))
heMod <- filter(hpge, modInTwo == "Y")
heMod <- length(unique(heMod$Protein))

output <- data.frame(Treatment = c("AHA_unenriched", "AHA_enriched", "HPG_unenriched", "HPG_enriched"),
                     Total_Prot_Count = c(au, ae, hu, he),
                     Mod_Prot_Count = c(auMod, aeMod, huMod, heMod),
                     Percent_Mod = c(auMod/au*100, aeMod/ae*100, huMod/hu*100, heMod/he*100))
write.csv(output, "mod_counts.csv", row.names = F)

#Inner join the AHA enriched table to the AHA unenriched table and the HPG enriched table to the HPG unenriched table####
aha2 <- inner_join(ahau, ahae, by = "Protein", suffix = c("_u", "_e"))%>%
  filter(modInTwo_u == "Y" | modInTwo_e == "Y") %>% 
  mutate(percentMod_u = NA, percentMod_e = NA, n_mod_hit_u = NA, n_mod_hit_e = NA, n_hit_u = NA, n_hit_e = NA, pValue = NA,
         sig = NA)
for(i in 1:nrow(aha2)){
  if(is.nan(aha2$percentMod_1_u[i]) == T){
    aha2$percentMod_1_u[i] <- NA
  }
  if(is.nan(aha2$percentMod_2_u[i]) == T){
    aha2$percentMod_2_u[i] <- NA
  }
  if(is.nan(aha2$percentMod_3_u[i]) == T){
    aha2$percentMod_3_u[i] <- NA
  }
  if(is.nan(aha2$percentMod_1_e[i]) == T){
    aha2$percentMod_1_e[i] <- NA
  }
  if(is.nan(aha2$percentMod_2_e[i]) == T){
    aha2$percentMod_1_u[i] <- NA
  }
  if(is.nan(aha2$percentMod_3_e[i]) == T){
    aha2$percentMod_1_u[i] <- NA
  }
}
aha2$percentMod_u <- round(rowMeans(aha2[, c("percentMod_1_u", "percentMod_2_u", "percentMod_3_u")], na.rm = T), digits = 0)
aha2$percentMod_e <- round(rowMeans(aha2[, c("percentMod_1_e", "percentMod_2_e", "percentMod_3_e")], na.rm = T), digits = 0)
aha2$n_mod_hit_u <- round(rowMeans(aha2[, c("n_mod_hit_1_u", "n_mod_hit_2_u", "n_mod_hit_3_u")], na.rm = T), digits = 0)
aha2$n_mod_hit_e <- round(rowMeans(aha2[, c("n_mod_hit_1_e", "n_mod_hit_2_e", "n_mod_hit_3_e")], na.rm = T), digits = 0)
aha2$n_hit_u <- round(rowMeans(aha2[, c("n_hit_1_u", "n_hit_2_u", "n_hit_3_u")], na.rm = T), digits = 0)
aha2$n_hit_e <- round(rowMeans(aha2[, c("n_hit_1_e", "n_hit_2_e", "n_hit_3_e")], na.rm = T), digits = 0)
for(i in 1:nrow(aha2)) {
  if((is.na(aha2$percentMod_1_e[i]) == F & is.na(aha2$percentMod_2_e[i]) == F) |
     (is.na(aha2$percentMod_1_e[i]) == F & is.na(aha2$percentMod_3_e[i]) == F) |
     (is.na(aha2$percentMod_2_e[i]) == F & is.na(aha2$percentMod_3_e[i]) == F)){
  aha2$pValue[i] <- tryCatch({round(t.test(aha2[i, c("percentMod_1_u", "percentMod_2_u", "percentMod_3_u")],
                                           aha2[i, c("percentMod_1_e", "percentMod_2_e", "percentMod_3_e")])$p.value, digits = 2)
  }, error = function(e) {
    "data are essentially constant"
  })
  }
}
for(i in 1:nrow(aha2)) {
  if(is.na(aha2$pValue[i]) == F & is.na(aha2$n_hit_u[i]) == F & is.na(aha2$n_hit_e[i]) == F){
    if(aha2$pValue[i] < 0.05 & aha2$n_hit_u[i] >= 5 & aha2$n_hit_e[i] >= 5) {
    aha2$sig[i] <- "Y"
  } else {
    aha2$sig[i] <- "N"
  }
}
}
aha2 <- aha2 %>% 
  select(-c("percentMod_1_u", "percentMod_2_u", "percentMod_3_u",
            "percentMod_1_e", "percentMod_2_e", "percentMod_3_e",
            "n_mod_hit_1_u", "n_mod_hit_2_u", "n_mod_hit_3_u",
            "n_mod_hit_1_e", "n_mod_hit_2_e", "n_mod_hit_3_e",
            "n_hit_1_u", "n_hit_2_u", "n_hit_3_u",
            "n_hit_1_e", "n_hit_2_e", "n_hit_3_e"))
hpg2 <- inner_join(hpgu, hpge, by = "Protein", suffix = c("_u", "_e")) %>%
  filter(modInTwo_u == "Y" | modInTwo_e == "Y") %>% 
  mutate(percentMod_u = NA, percentMod_e = NA, n_mod_hit_u = NA, n_mod_hit_e = NA, n_hit_u = NA, n_hit_e = NA, pValue = NA,
         sig = NA)
for(i in 1:nrow(hpg2)){
  if(is.nan(hpg2$percentMod_1_u[i]) == T){
    hpg2$percentMod_1_u[i] <- NA
  }
  if(is.nan(hpg2$percentMod_2_u[i]) == T){
    hpg2$percentMod_2_u[i] <- NA
  }
  if(is.nan(hpg2$percentMod_3_u[i]) == T){
    hpg2$percentMod_3_u[i] <- NA
  }
  if(is.nan(hpg2$percentMod_1_e[i]) == T){
    hpg2$percentMod_1_e[i] <- NA
  }
  if(is.nan(hpg2$percentMod_2_e[i]) == T){
    hpg2$percentMod_1_u[i] <- NA
  }
  if(is.nan(hpg2$percentMod_3_e[i]) == T){
    hpg2$percentMod_1_u[i] <- NA
  }
}
hpg2$percentMod_u <- round(rowMeans(hpg2[, c("percentMod_1_u", "percentMod_2_u", "percentMod_3_u")], na.rm = T), digits = 0)
hpg2$percentMod_e <- round(rowMeans(hpg2[, c("percentMod_1_e", "percentMod_2_e", "percentMod_3_e")], na.rm = T), digits = 0)
hpg2$n_mod_hit_u <- round(rowMeans(hpg2[, c("n_mod_hit_1_u", "n_mod_hit_2_u", "n_mod_hit_3_u")], na.rm = T), digits = 0)
hpg2$n_mod_hit_e <- round(rowMeans(hpg2[, c("n_mod_hit_1_e", "n_mod_hit_2_e", "n_mod_hit_3_e")], na.rm = T), digits = 0)
hpg2$n_hit_u <- round(rowMeans(hpg2[, c("n_hit_1_u", "n_hit_2_u", "n_hit_3_u")], na.rm = T), digits = 0)
hpg2$n_hit_e <- round(rowMeans(hpg2[, c("n_hit_1_e", "n_hit_2_e", "n_hit_3_e")], na.rm = T), digits = 0)
for(i in 1:nrow(hpg2)) {
  if((is.na(hpg2$percentMod_1_e[i]) == F & is.na(hpg2$percentMod_2_e[i]) == F) |
     (is.na(hpg2$percentMod_1_e[i]) == F & is.na(hpg2$percentMod_3_e[i]) == F) |
     (is.na(hpg2$percentMod_2_e[i]) == F & is.na(hpg2$percentMod_3_e[i]) == F)){
    hpg2$pValue[i] <- tryCatch({round(t.test(hpg2[i, c("percentMod_1_u", "percentMod_2_u", "percentMod_3_u")],
                                   hpg2[i, c("percentMod_1_e", "percentMod_2_e", "percentMod_3_e")])$p.value, digits = 2)
    }, error = function(e) {
      "data are essentially constant"
    })
  }
}
for(i in 1:nrow(hpg2)) {
  if(is.na(hpg2$pValue[i]) == F & is.na(hpg2$n_hit_u[i]) == F & is.na(hpg2$n_hit_e[i]) == F){
    if(hpg2$pValue[i] < 0.05 & hpg2$n_hit_u[i] >= 5 & hpg2$n_hit_e[i] >= 5) {
      hpg2$sig[i] <- "Y"
    } else {
      hpg2$sig[i] <- "N"
    }
  }
}
hpg2 <- hpg2 %>% 
  select(-c("percentMod_1_u", "percentMod_2_u", "percentMod_3_u",
            "percentMod_1_e", "percentMod_2_e", "percentMod_3_e",
            "n_mod_hit_1_u", "n_mod_hit_2_u", "n_mod_hit_3_u",
            "n_mod_hit_1_e", "n_mod_hit_2_e", "n_mod_hit_3_e",
            "n_hit_1_u", "n_hit_2_u", "n_hit_3_u",
            "n_hit_1_e", "n_hit_2_e", "n_hit_3_e"))

aha2 <- mod_in_two(aha2)
hpg2 <- mod_in_two(hpg2)

#Count the number of modifications on each protein in each group
aha2 <- count_mods(aha2)
hpg2 <- count_mods(hpg2)

#Create and output file that shows the number of modifications for each protein for each treatment group.
#add mapman data as new variables
mapman <- read.delim("mapman.txt", stringsAsFactors = F, na.strings = "\''") %>% 
  select("IDENTIFIER", "BINCODE", "NAME", "DESCRIPTION") %>% 
  dplyr::rename("prot" = "IDENTIFIER") %>% 
  mutate_all(funs(str_replace(., "\'", ""))) %>% 
  mutate_all(funs(str_replace(., "\'", "")))
mapman <- set_names(mapman, tolower(names(mapman)))
mapman$prot <- toupper(mapman$prot)
mapman <- mapman[complete.cases(mapman), ]
tair10 <- ahae1 %>% 
  select(c("Protein", "Protein.ID"))
output2 <- full_join(aha2, hpg2, by = "Protein", suffix = c("_AHA", "_HPG")) %>% 
  mutate(u_mod_n_greater_AHA = NA, u_mod_n_greater_HPG = NA, u_mod_percent_greater_AHA = NA, u_mod_percent_greater_HPG = NA) %>% 
  left_join(tair10, by = "Protein")
for(i in 1:nrow(output2)) {
  output2$u_mod_n_greater_AHA[i] <- output2$Total_mods_across_3_reps_U_AHA[i] > output2$Total_mods_across_3_reps_E_AHA[i]
  output2$u_mod_n_greater_HPG[i] <- output2$Total_mods_across_3_reps_U_HPG[i] > output2$Total_mods_across_3_reps_E_HPG[i]
  output2$u_mod_percent_greater_AHA[i] <- output2$percentMod_u_AHA[i] > output2$percentMod_e_AHA[i]
  output2$u_mod_percent_greater_HPG[i] <- output2$percentMod_u_HPG[i] > output2$percentMod_e_HPG[i]
}
output2 <- mutate_each(output2, funs(empty_as_na))
output2 <- mutate_each(output2, funs(empty_as_na2))
output2 <- mutate_each(output2, funs(empty_as_na3))
output2 <- left_join(output2, mapman, by = c("Protein" ="prot"))
write.csv(output2, "mod_counts_for_indiviaul_proteins.csv", row.names = F)

####To calculate the appearance of individual mods####
#select only AHA-relevant columns
aha_mods <- output2 %>% 
  select(c("Protein", "Mods_u_AHA", "modInTwo_u_AHA"))

#filter so that only proteins that are modified in two in the unenriched set remain
aha_mods <- filter(aha_mods, modInTwo_u_AHA == "Y")

#select the unenriched mods column only, remove the parentheses
aha_mods <- aha_mods %>% 
  select("Mods_u_AHA")
aha_mods$Mods_u_AHA <- str_replace_all(aha_mods$Mods_u_AHA, "\\)", "")
aha_mods$Mods_u_AHA <- str_replace_all(aha_mods$Mods_u_AHA, "\\(", "")

#split in the modifications into seperate strings and place them all in one column of the dataframe
aha_modsV <- str_trim(unlist(as.vector(str_split(aha_mods$Mods_u_AHA, pattern = ","))), side = "both") %>% 
  str_replace(., "[:digit:]+M", "M") %>% 
  str_replace(., "[:digit:]+P", "P") %>% 
  str_replace(., "[:digit:]+K", "K") %>% 
  str_replace(., "[:digit:]+S", "S") %>% 
  str_replace(., "[:digit:]+T", "T") %>% 
  str_replace(., "[:digit:]+Y", "Y") %>% 
  str_replace(., "[:digit:]+E", "E") %>% 
  str_replace(., "[:digit:]+R", "R")

aha_mods <- data.frame(PTM = aha_modsV)
#use table to count the number of each modification
n_aha_mods_u <- as.data.frame(table(aha_mods)) %>%
  rename("aha_u_freq" = "Freq", "mods" = "aha_mods")

#select only AHA-relevant columns
aha_mods <- output2 %>% 
  select(c("Protein", "Mods_e_AHA", "modInTwo_e_AHA"))

#filter so that only proteins that are modified in two in the unenriched set remain
aha_mods <- filter(aha_mods, modInTwo_e_AHA == "Y")

#select the unenriched mods column only, remove the parentheses
aha_mods <- aha_mods %>% 
  select("Mods_e_AHA")
aha_mods$Mods_e_AHA <- str_replace_all(aha_mods$Mods_e_AHA, "\\)", "")
aha_mods$Mods_e_AHA <- str_replace_all(aha_mods$Mods_e_AHA, "\\(", "")

#split in the modifications into seperate strings and place them all in one column of the dataframe
aha_modsV <- str_trim(unlist(as.vector(str_split(aha_mods$Mods_e_AHA, pattern = ","))), side = "both") %>% 
  str_replace(., "[:digit:]+M", "M") %>% 
  str_replace(., "[:digit:]+P", "P") %>% 
  str_replace(., "[:digit:]+K", "K") %>% 
  str_replace(., "[:digit:]+S", "S") %>% 
  str_replace(., "[:digit:]+T", "T") %>% 
  str_replace(., "[:digit:]+Y", "Y") %>% 
  str_replace(., "[:digit:]+E", "E") %>% 
  str_replace(., "[:digit:]+R", "R")

aha_mods <- data.frame(PTM = aha_modsV)
#use table to count the number of each modification
n_aha_mods_e <- as.data.frame(table(aha_mods)) %>%
  rename("aha_e_freq" = "Freq", "mods" = "aha_mods")

#select only HPG-relevant columns
hpg_mods <- output2 %>% 
  select(c("Protein", "Mods_u_HPG", "modInTwo_u_HPG"))

#filter so that only proteins that are modified in two in the unenriched set remain
hpg_mods <- filter(hpg_mods, modInTwo_u_HPG == "Y")

#select the unenriched mods column only, remove the parentheses
hpg_mods <- hpg_mods %>% 
  select("Mods_u_HPG")
hpg_mods$Mods_u_HPG <- str_replace_all(hpg_mods$Mods_u_HPG, "\\)", "")
hpg_mods$Mods_u_HPG <- str_replace_all(hpg_mods$Mods_u_HPG, "\\(", "")

#split in the modifications into seperate strings and place them all in one column of the dataframe
hpg_modsV <- str_trim(unlist(as.vector(str_split(hpg_mods$Mods_u_HPG, pattern = ","))), side = "both") %>% 
  str_replace(., "[:digit:]+M", "M") %>% 
  str_replace(., "[:digit:]+P", "P") %>% 
  str_replace(., "[:digit:]+K", "K") %>% 
  str_replace(., "[:digit:]+S", "S") %>% 
  str_replace(., "[:digit:]+T", "T") %>% 
  str_replace(., "[:digit:]+Y", "Y") %>% 
  str_replace(., "[:digit:]+E", "E") %>% 
  str_replace(., "[:digit:]+R", "R")

hpg_mods <- data.frame(PTM = hpg_modsV)
#use table to count the number of each modification
n_hpg_mods_u <- as.data.frame(table(hpg_mods)) %>%
  rename("hpg_u_freq" = "Freq", "mods" = "hpg_mods")

#select only HPG-relevant columns
hpg_mods <- output2 %>% 
  select(c("Protein", "Mods_e_HPG", "modInTwo_e_HPG"))

#filter so that only proteins that are modified in two in the unenriched set remain
hpg_mods <- filter(hpg_mods, modInTwo_e_HPG == "Y")

#select the unenriched mods column only, remove the parentheses
hpg_mods <- hpg_mods %>% 
  select("Mods_e_HPG")
hpg_mods$Mods_e_HPG <- str_replace_all(hpg_mods$Mods_e_HPG, "\\)", "")
hpg_mods$Mods_e_HPG <- str_replace_all(hpg_mods$Mods_e_HPG, "\\(", "")

#split in the modifications into seperate strings and place them all in one column of the dataframe
hpg_modsV <- str_trim(unlist(as.vector(str_split(hpg_mods$Mods_e_HPG, pattern = ","))), side = "both") %>% 
  str_replace(., "[:digit:]+M", "M") %>% 
  str_replace(., "[:digit:]+P", "P") %>% 
  str_replace(., "[:digit:]+K", "K") %>% 
  str_replace(., "[:digit:]+S", "S") %>% 
  str_replace(., "[:digit:]+T", "T") %>% 
  str_replace(., "[:digit:]+Y", "Y") %>% 
  str_replace(., "[:digit:]+E", "E") %>% 
  str_replace(., "[:digit:]+R", "R")

hpg_mods <- data.frame(PTM = hpg_modsV)
#use table to count the number of each modification
n_hpg_mods_e <- as.data.frame(table(hpg_mods)) %>%
  rename("hpg_e_freq" = "Freq", "mods" = "hpg_mods")
mod_counts_joined <- n_aha_mods_u %>% 
  full_join(n_aha_mods_e, by = "mods") %>% 
  full_join(n_hpg_mods_u, by = "mods") %>% 
  full_join(n_hpg_mods_e, by = "mods") %>% 
  mutate(u_freq = aha_u_freq + hpg_u_freq) %>% 
  select(-c("aha_u_freq", "hpg_u_freq"))
mod_counts_joined$mods <- str_replace_all(mod_counts_joined$mods, "14.0157", " methylation")
mod_counts_joined$mods <- str_replace_all(mod_counts_joined$mods, "42.0106", " acetylation")
mod_counts_joined$mods <- str_replace_all(mod_counts_joined$mods, "15.9949", " oxidation")
mod_counts_joined$mods <- str_replace_all(mod_counts_joined$mods, "42.0106", " acetylation")
mod_counts_joined$mods <- str_replace_all(mod_counts_joined$mods, "79.96633", " phosphorylation")
mod_counts_joined$mods <- str_replace_all(mod_counts_joined$mods, "-21.9877", "-to-HPG substitution")
mod_counts_joined$mods <- str_replace_all(mod_counts_joined$mods, "-4.9863", "-to-AHA substitution")
mod_tot_u <- sum(mod_counts_joined$u_freq, na.rm = T)
mod_tot_aha_e <- sum(mod_counts_joined$aha_e_freq, na.rm = T)
mod_tot_hpg_e <- sum(mod_counts_joined$hpg_e_freq, na.rm = T)
mod_counts_joined <- mod_counts_joined %>% 
  mutate(u_freq = u_freq/mod_tot_u*100)%>% 
  mutate(aha_e_freq = aha_e_freq/mod_tot_aha_e*100)%>% 
  mutate(hpg_e_freq = hpg_e_freq/mod_tot_hpg_e*100) %>% 
  rename("u_%_occurrence" = "u_freq", "aha_e_%_occurrence" = "aha_e_freq", "hpg_e_%_occurrence" = "hpg_e_freq")
write.csv(mod_counts_joined, "mod_frequencies.csv", row.names = F)

####Make a figure of the percentages####
mod_counts_joined2 <- mod_counts_joined %>% 
  pivot_longer(cols = 2:4, names_to = "Group", values_to = "occurrence") %>% 
  mutate(Group = str_sub(Group, end = 5)) %>% 
  mutate(Group = str_replace(Group, "aha_e", "AHA nascent"),
         Group = str_replace(Group, "hpg_e", "HPG nascent"),
         Group = str_replace(Group, "u_%_o", "Bulk"))
png("Figure 6 mod counts.png", width =  500, height = 500)
ggplot(mod_counts_joined2, aes(x = as.factor(mods), y = occurrence, fill = Group)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Modification", y = "Frequency (%)") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_minimal(base_size = 18) +
  theme(axis.title.y = element_text(size = 18, hjust = .4),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "bottom")
dev.off()
