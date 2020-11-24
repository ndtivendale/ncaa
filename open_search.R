library(ggplot2)
library(tidyr)
setwd("C:/temp/noncanonical_amino_acids_r_processing/protein_and_peptide")

load("C:/temp/noncanonical_amino_acids_r_processing/protein_and_peptide/ahahpgopensearch.RData")
shist <- function(x, y, z){
  ggplot(x, aes(y)) +
    geom_histogram(binwidth = 0.1, breaks = seq(-300, 350, by = 0.1)) +
    ylim(0, 1000) + 
    labs(title = z, x = "massdiff", y = "count")
}

carb <- 57.0215
#Import AGI lists of AHA and HPG-enriched proteins and HPG-tagged proteins
ahaEnri <- read.csv("protFoldEnriAHA.csv")
ahaEnriV <- as.vector(ahaEnri$firstID)
hpgEnri <- read.csv("protFoldEnriHPG.csv")
hpgEnriV <- as.vector(hpgEnri$firstID)
hpgEnriTag <- hpgEnri %>% 
  filter(tagged_unenriched == "Y")
hpgEnriTagV <- as.vector(hpgEnriTag$firstID)
#select only peptides containing methionine and reduce the number of dps in massdiff to 4
aha2 <- aha %>% 
  filter(str_detect(peptide, "M"))
notag2 <- notag %>% 
  filter(str_detect(peptide, "M"))
hpg2 <- hpg %>% 
  filter(str_detect(peptide, "M"))

#filter each dataset against the HPG-tagged AGIs
aha3 <- aha2 %>% 
  filter(seqsource %in% hpgEnriTagV)
notag3 <- notag2 %>% 
  filter(seqsource %in% hpgEnriTagV)
hpg3 <- hpg2 %>% 
  filter(seqsource %in% hpgEnriTagV)

#filter each dataset against the AHA-enriched AGIs
aha4 <- aha2 %>% 
  filter(seqsource %in% ahaEnriV)
notag4 <- notag2 %>% 
  filter(seqsource %in% ahaEnriV)
hpg4 <- hpg2 %>% 
  filter(seqsource %in% ahaEnriV)

#filter each dataset against the HPG-enriched AGIs
aha5 <- aha2 %>% 
  filter(seqsource %in% hpgEnriV)
notag5 <- notag2 %>% 
  filter(seqsource %in% hpgEnriV)
hpg5 <- hpg2 %>% 
  filter(seqsource %in% hpgEnriV)

#make histogram plots for each filtered dataset
aha_plot <- shist(aha, aha$massdiff, "aha1")
notag_plot <- shist(notag, notag$massdiff, "notag1")
hpg_plot <- shist(hpg, hpg$massdiff, "hpg1")

aha_plot2 <- shist(aha2, aha2$massdiff, "aha2")
notag_plot2 <- shist(notag2, notag2$massdiff, "notag2")
hpg_plot2 <- shist(hpg2, hpg2$massdiff, "hpg2")

aha_plot3 <- shist(aha3, aha3$massdiff, "aha3")
notag_plot3 <- shist(notag3, notag3$massdiff, "notag3")
hpg_plot3 <- shist(hpg3, hpg3$massdiff, "hpg3")

aha_plot4 <- shist(aha4, aha4$massdiff, "aha4")
notag_plot4 <- shist(notag4, notag4$massdiff, "notag4")
hpg_plot4 <- shist(hpg4, hpg4$massdiff, "hpg4")

aha_plot5 <- shist(aha5, aha5$massdiff, "aha5")
notag_plot5 <- shist(notag5, notag5$massdiff, "notag5")
hpg_plot5 <- shist(hpg5, hpg5$massdiff, "hpg5")

png("aha_plot.png", width = 1000, height = 1000)
aha_plot
dev.off()

png("hpg_plot.png", width = 1000, height = 1000)
hpg_plot
dev.off()

png("notag_plot.png", width = 1000, height = 1000)
notag_plot
dev.off()

png("aha_plot2.png", width = 1000, height = 1000)
aha_plot2
dev.off()

png("hpg_plot2.png", width = 1000, height = 1000)
hpg_plot2
dev.off()

png("notag_plot2.png", width = 1000, height = 1000)
notag_plot2
dev.off()

png("aha_plot3.png", width = 1000, height = 1000)
aha_plot3
dev.off()

png("hpg_plot3.png", width = 1000, height = 1000)
hpg_plot3
dev.off()

png("notag_plot3.png", width = 1000, height = 1000)
notag_plot3
dev.off()

png("aha_plot4.png", width = 1000, height = 1000)
aha_plot4
dev.off()

png("hpg_plot4.png", width = 1000, height = 1000)
hpg_plot4
dev.off()

png("notag_plot4.png", width = 1000, height = 1000)
notag_plot4
dev.off()

png("aha_plot5.png", width = 1000, height = 1000)
aha_plot5
dev.off()

png("hpg_plot5.png", width = 1000, height = 1000)
hpg_plot5
dev.off()

png("notag_plot5.png", width = 1000, height = 1000)
notag_plot5
dev.off()

aha_rounded <- aha
aha_rounded$massdiff <- round(aha_rounded$massdiff, digits = 4)
aha_rounded2 <- aha_rounded %>% 
  filter(str_detect(peptide, "C")) %>% 
  filter(str_detect(peptide, "M")) %>% 
  filter(massdiff == carb-4.9863)

hpg_rounded <- hpg
hpg_rounded$massdiff <- round(hpg_rounded$massdiff, digits = 4)
hpg_rounded2 <- hpg_rounded %>% 
  filter(str_detect(peptide, "C")) %>% 
  filter(str_detect(peptide, "M")) %>% 
  filter(massdiff == carb-21.9877)
                       