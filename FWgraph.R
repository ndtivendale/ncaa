#set working directory
setwd("C:/temp/noncanonical_amino_acids_r_processing/cell_culture_growth")

#import data file
fresh_weights <- read.csv("NT-20-18,19FWs.csv", stringsAsFactors = F, header = T)

#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

#Rename poorly named variables, replace 'NA's with '0's, remove IS results, group by treatment, compound and time point, get the mean level
fresh_weights2 <- fresh_weights %>%
  rename("time" = "Time..h.") %>%
  mutate(percentDensityChange = Density.change..relative.to.initial.density.*100) %>%
  group_by(Treatment, time) %>%
  summarise(meanDensityChange = mean(percentDensityChange), stDevDensityChange = (sd(percentDensityChange))/sqrt(3)) %>%
  ungroup() %>% 
  mutate(Treatment = factor(Treatment, levels = c("AHA", "HPG", "Met", "Nitrogen-15", "Control")))
fresh_weights2$Treatment = factor(fresh_weights2$Treatment, labels = c("AHA", "HPG", "Met", "{}^{15}~N", "Control"))
fresh_weights2 <- filter(fresh_weights2, Treatment != "Met")
#Make a scatter plot object
postscript("Figure 1D FW_growth_bar_plot_separate.ps")
ggplot(fresh_weights2, aes(x = factor(time), y = meanDensityChange, fill = Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~ Treatment, ncol = 4, labeller = label_parsed) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
  geom_errorbar(aes(ymin = meanDensityChange - stDevDensityChange, ymax = meanDensityChange + stDevDensityChange),
                width = .2, position=position_dodge()) +
  theme(axis.text = element_text(size = 18), strip.text = element_text(size = 18, face = "bold")) + 
  labs(y = expression(atop("Mean density change", paste("(% relative to initial density)"))), x = "Time (h)")
dev.off()

#Find statistically significant differences
fresh_weights8 <- fresh_weights %>%
  rename("time" = "Time..h.") %>%
  filter(time == 8.5) %>%
  group_by(Treatment, time)
pairwise_ttest8 <- pairwise.t.test(fresh_weights8$Density.change..relative.to.initial.density., fresh_weights8$Treatment,
                                  p.adjust.method = "none",
                                  paired = FALSE,
                                  pool.sd = FALSE,
                                  alternative = "two.sided")
pairwise_ttest_df8 <- do.call(rbind.data.frame, pairwise_ttest8) 

fresh_weights24 <- fresh_weights %>%
  rename("time" = "Time..h.") %>%
  filter(time == 24) %>%
  group_by(Treatment, time)
pairwise_ttest24 <- pairwise.t.test(fresh_weights24$Density.change..relative.to.initial.density., fresh_weights24$Treatment,
                                   p.adjust.method = "none",
                                   paired = FALSE,
                                   pool.sd = FALSE,
                                   alternative = "two.sided")
pairwise_ttest_df24 <- do.call(rbind.data.frame, pairwise_ttest24) 

write.csv(pairwise_ttest_df8, "pairwise_ttest_8h.csv")
write.csv(pairwise_ttest_df24, "pairwise_ttest_24h.csv")
