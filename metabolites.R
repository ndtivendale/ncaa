#set working directory
setwd("C:/temp/noncanonical_amino_acids_r_processing/metabolite")

#import data file
metabolites <- read.csv("metMetabolitesrr.csv", stringsAsFactors = F, header = T)

#load packages
library(dplyr)
library(tidyr)
library(ggplot2)

#Rename poorly named variables, replace 'NA's with '0's, remove IS results, group by treatment, compound and time point, get the mean level
metabolites2 <- metabolites %>%
  rename("time" = "time.h.", "5-Me-S-adenosine" = "X5.Me.S.Ado", "Treatment" = "treatment") %>%
  gather(compound, level, SAM:Met) %>%
  group_by(Treatment, compound, time) %>%
  summarise(meanLevel = mean(level), stDevLevel = (sd(level))/sqrt(3))

hpg <- as.numeric(metabolites2[57,4])
aha <- as.numeric(metabolites2[9,4])
met_h <- as.numeric(metabolites2[60,4])
met_a <- as.numeric(metabolites2[18,4])

aha_ratio <- aha/met_a
hpg_ratio <- hpg/met_h

#Remove AHA and HPG
metabolites3 <- metabolites2 %>%
  filter(compound != "HPG" & compound != "AHA")

#Make a scatter plot object
png("Figure 2 metabolites plot.png", width =  500, height = 500)
ggplot(metabolites3, aes(x = time, y = meanLevel, color = Treatment)) +
  scale_y_log10() +
  facet_wrap(~ compound, nrow = 3, ncol = 2, scales = "free_y") +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = meanLevel - stDevLevel, ymax = meanLevel + stDevLevel, width = 1)) + 
  geom_errorbarh(aes(xmin = time - 80/120, xmax = time + 80/120, height = 0.2)) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
  labs(y = expression(paste("Mean abundance (", mu, "g/g (FW)")), x = "Time (h)") +
  theme_minimal(base_size = 18)  +
  theme(axis.title.y = element_text(size = 18, hjust = .4),
        axis.title.x = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.position = c(0.75, 0.13),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 12))
dev.off()

#pairwise ttest for significance
metabolites4 <- metabolites %>%
  filter(time.h. == 24 & treatment != "Control") %>%
  rename("MeSAdo" = "X5.Me.S.Ado")
pairwise_ttest <- pairwise.t.test(metabolites4$Met, metabolites4$treatment,
                                  p.adjust.method = "none",
                                  paired = TRUE,
                                  pool.sd = FALSE,
                                  alternative = "two.sided")
pairwise_ttest_df <- do.call(rbind.data.frame, pairwise_ttest)
write.csv(pairwise_ttest_df, "pairwise_ttest_Met.csv")