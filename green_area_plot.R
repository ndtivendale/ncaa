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
labels <- c(aha = "AHA", hpg = "HPG", met = "Met", control = "Control", no_treatment = "No treatment")
#Make a scatter plot object
png("Figure 1B (whole plant treatment effect).png", width = 750, height = 250)
ggplot(a, aes(x = factor(time), y = mean_area_change)) +
  geom_bar(stat="identity", fill ="#009E73", position = position_dodge()) +
  facet_wrap(~ treatment, labeller = labeller(treatment = labels), ncol = 5) +
  theme_minimal(base_size = 18) +
  geom_errorbar(aes(ymin = mean_area_change - se, ymax = mean_area_change + se),
                width = .2, position=position_dodge())+
  theme(axis.text = element_text(size = 18), strip.text = element_text(size = 18)) + 
  labs(y = expression(atop("Mean grean area change", paste("(% relative to initial density)"))), x = "Time (h)")
dev.off()