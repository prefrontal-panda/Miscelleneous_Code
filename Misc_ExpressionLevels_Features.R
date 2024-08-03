# This is to make plots to compare the expression levels of various features between groups

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Load library
library(ggplot2)
library(tidyverse)
library(ggrepel)

#----GAERS: Alpha-ketoglutarate----
dt <- read.csv("GAERS_Tha_Met_AllTimePt.csv", stringsAsFactors = F, header = T)
dt1 <- dt %>%
  filter(grepl("C00026", feature))
# Making plot
pdf("Tha_Alpha_ketoglutarate_ExpressionLevels.pdf")
# All on one plot
ggplot(dt1, aes(x=as.factor(week), y=value, fill=group)) +
  geom_violin(alpha=0.8, position = position_dodge(0.9)) +
  geom_boxplot(width=0.1, color="black", alpha=0.5, position = position_dodge(0.9)) +
  xlab("") +
  ylab("Alpha-ketoglutarate") +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Alpha-ketoglutarate expression levels in Tha")
# individual plots
for (i in as.factor(c("3","7","16"))) {
  #print(i)
  d <- dt1 %>%
    filter(grepl(i, week))
  #print(d)
  plt <- ggplot(d, aes(x=as.factor(group), y=value, fill=group)) +
    geom_violin(alpha=0.8) +
    geom_boxplot(width=0.1, color="black", alpha=0.5) +
    xlab("") +
    ylab("Alpha-ketoglutarate") +
    theme_bw() +
    theme(legend.position = "none") +
    ggtitle(paste0("Alpha-ketoglutarate expression levels in Tha of ", i, "wk animals"))
  print(plt)
}
dev.off()

#----PTE metabolites----
# Read in dataframe
dt <- read.csv("PTE_Metabolomics_Pre-processed_SampImp_PCA.csv", stringsAsFactors = F,
               header = T)

# Removing Naives and Sham
dt_PTE_TBI <- dt %>%
  filter(label == "PTE" | label == "TBI")
# Removing outliers
dt_noOut <- dt[-c(1,28),]

# Violin plots
lab <- c("X","label","Carbamoyl.phosphate", "LysoPC.16.0.")
dt_violin <- dt[,lab]
# Making long dataframe
dt_long <- dt_violin %>%
  pivot_longer(!c("X","label"), #using the metabolite and values are reference
               names_to = "Metabolite", # metabolites to column
               values_to = "value") %>% # normalised intensities to this column
  arrange(Metabolite)
# Separating metabolites
dt_CP <- dt_long %>%
  filter(Metabolite == "Carbamoyl.phosphate")
dt_LPC <- dt_long %>%
  filter(Metabolite == "LysoPC.16.0.")
# Plotting
# Carbamoyl Phosphate
dt_CP %>%
  ggplot(aes(x=label, y=value, fill = label, alpha = 0.5)) +
  geom_violin() + 
  geom_jitter(color="black", size=1, alpha=1, 
              position = position_jitter(seed = 1)) +
  geom_text_repel(aes(label = X), size = 3, 
                  position = position_jitter(seed = 1), max.overlaps = Inf) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("Carbamoyl Phosphate") +
  xlab("")
# LysoPC(16:0)
dt_LPC %>%
  ggplot(aes(x=label, y=value, fill = label, alpha = 0.5)) +
  geom_violin() + 
  geom_jitter(color="black", size=1, alpha=1, 
              position = position_jitter(seed = 1)) +
  geom_text_repel(aes(label = X), size = 3, 
                  position = position_jitter(seed = 1), max.overlaps = Inf) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(size = 11)
  ) +
  ggtitle("LysoPC(16:0)") +
  xlab("")