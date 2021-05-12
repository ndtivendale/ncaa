#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

#set working directory
setwd("C:/temp/noncanonical_amino_acids_r_processing/whole_plant_growth")

#import data file
a <- read.csv("combined_whole_plant_analysis.csv", stringsAsFactors = F, header = T) %>% 
  select("treatment", "time", "mean_area_change", "se") %>% 
  mutate(mean_area_change = mean_area_change*100, se = se*100) %>% 
  filter(time != 72)
a$treatment <- factor(a$treatment, levels = c("aha", "hpg", "met", "control", "no_treatment"))
a <- a[complete.cases(a),]
labels <- c(aha = "AHA", hpg = "HPG", met = "Met", control = "Mock", no_treatment = "No treatment")
a <- filter(a, treatment != "met")
#Make a barplot
postscript("Figure 1C (whole plant treatment effect).ps")
ggplot(a, aes(x = factor(time), y = mean_area_change, fill = treatment)) +
  geom_bar(stat="identity", position = position_dodge()) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#009E73")) +
  facet_wrap(~ treatment, labeller = labeller(treatment = labels), ncol = 5) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin = mean_area_change - se, ymax = mean_area_change + se),
                width = .2, position=position_dodge())+
  theme(axis.text = element_text(size = 18), strip.text = element_text(size = 18)) + 
  labs(y = expression(atop("Mean grean area change", paste("(% relative to initial density)"))), x = "Time (h)")
dev.off()

#T-tests
a24 <- read.csv("combined_whole_plant_analysis.csv", stringsAsFactors = F, header = T) %>% 
  select("treatment", "time", "green_area_change") %>% 
  filter(time == 24) %>%
  group_by(treatment, time)
pairwise_ttest24 <- pairwise.t.test(a24$green_area_change, a24$treatment,
                                   p.adjust.method = "none",
                                   paired = FALSE,
                                   pool.sd = FALSE,
                                   alternative = "two.sided")
pairwise_ttest_df24 <- do.call(rbind.data.frame, pairwise_ttest24) 

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

write_csv(pairwise_ttest_df24, "pairwise_ttest_24h.csv")
write_csv(pairwise_ttest_df48, "pairwise_ttest_48h.csv")