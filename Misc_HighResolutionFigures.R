# This is to make higher resolution .tiff figures for the thesis.

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

#----GAERS Volcano Plots----
# Load libraries
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(tidyverse)

# Reading in datafiles
scx <- read.csv("GAERS_Scx_Trans_16wks_DEAll_FC2.csv", stringsAsFactors = F, sep = ",", row.names = 1)
tha <- read.csv("GAERS_Tha_Trans_16wks_DEAll_FC2.csv", stringsAsFactors = F, sep = ",", row.names = 1)

# Making TRANSCRIPTOMICS Volcano Plots
top_n <- 10
# Scx
top_genes <- rbind(
  scx %>%
    filter(Sig.DE == "Up") %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top_n),
  scx %>%
    filter(Sig.DE == "Down") %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top_n)  
)
scx_vol <- ggplot(scx, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = Sig.DE), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"(P-Value)")) +
  scale_color_manual(labels = c("Down-regulated", "Not Significant", "Up-regulated"),
                     values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(data = top_genes,
                  aes(logFC, -log(PValue,10), 
                      label = gene_name),
                  size = 3,
                  max.overlaps = 50)
scx_vol <- ggpar(scx_vol, legend.title = "Differential Expression")
# Thalamus
top_tha <- rbind(
  tha %>%
    filter(Sig.DE == "Up") %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top_n),
  tha %>%
    filter(Sig.DE == "Down") %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(top_n)  
)
tha_vol <- ggplot(tha, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = Sig.DE), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"(P-Value)")) +
  scale_color_manual(labels = c("Down-regulated", "Not Significant", "Up-regulated"),
                     values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(data = top_tha,
                  aes(logFC, -log(PValue,10), 
                      label = gene_name),
                  size = 3,
                  max.overlaps = 50)
tha_vol <- ggpar(tha_vol, legend.title = "Differential Expression")

# Putting together on same plot
final <- ggarrange(scx_vol, tha_vol, 
                   labels = c("(a)", "(b)"),
                   font.label = list(size = 11),
                   ncol = 2, nrow = 1,
                   common.legend = T, legend = "bottom")
# Saving
ggsave(filename = "GAERS_TransVolPlot.tiff", plot = final, 
       width = 20, height = 16, units = "cm", device='tiff', dpi=300)

#----PTE Volcano Plot----
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(tidyverse)

# Reading in datafiles
trans <- read.csv("PTE_Trans_TBI+PTE vs TBI-PTE_DEAll_FC2.csv", stringsAsFactors = F, sep = ",", row.names = 1)
met <- read.csv("TBI+PTE vs TBI-PTE_FC1.5_NoOutliers.csv", stringsAsFactors = F, sep = ",")

# Making TRANSCRIPTOMICS Volcano Plot
gene <- trans %>%
  filter(gene_name == "Cd22")
trans_vol <- ggplot(trans, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = Sig.DE), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"(P-Value)")) +
  scale_color_manual(labels = c("Not Significant", "Up-regulated"),
                     values = c("gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(data = gene,
                  aes(logFC, -log(PValue,10), 
                      label = gene_name),
                  size = 3,
                  max.overlaps = 50) +
  theme_bw()
trans_vol <- ggpar(trans_vol, legend.title = "Differential Expression")

# METABOLOMICS Volcano Plot
# Changing names
met$X[1:3] <- c("LysoPC(16:0)", "Carbamoyl Phosphate", "Tagatose,1-6,bisphosphate")
metabolites <- met %>%
  filter(X %in% c("LysoPC(16:0)", "Carbamoyl Phosphate", "Tagatose,1-6,bisphosphate"))
met_vol <- ggplot(met, aes(logFC, -log(P.Value,10))) +
  geom_point(aes(color = Sig.DE), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"(P-Value)")) +
  scale_color_manual(labels = c("Down-regulated", "Not Significant", "Up-regulated"),
                     values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_text_repel(data = metabolites,
                  aes(logFC, -log(P.Value,10), 
                      label = X),
                  size = 3,
                  max.overlaps = 50) +
  theme_bw()
met_vol <- ggpar(met_vol, legend.title = "Differential Expression")

# Putting together on same plot
final <- ggarrange(trans_vol, met_vol, 
                   labels = c("(a)", "(b)"),
                   font.label = list(size = 11),
                   ncol = 2, nrow = 1,
                   common.legend = T, legend = "bottom")
# Saving
ggsave(filename = "PTE_VolPlot.tiff", plot = final, 
       width = 20, height = 16, units = "cm", device='tiff', dpi=300)

#---- MEFISTO Common Features----
# Loading libraries
library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)

# Reading in MEFISTO trained object
filepath <- file.path(getwd(), "MEFISTO_Input_nogroup_20230627.hdf5") 
trained_model <- load_model(filepath)

# Plotting factor plots
samp_metadata <- read.csv("Aim 1 - GAERS/Sample metadata.csv", stringsAsFactors = F, sep = ",") 
samp_metadata$sample <- gsub("s$", "t", samp_metadata$sample) 
metadata <- samp_metadata[samp_metadata$sample %in% samples_metadata(trained_model)$sample,] 
metadata <- metadata[order(match(metadata$sample, samples_metadata(trained_model)$sample)),] 
colnames(metadata) <- c("sample", "strain", "strainCode", "week", "3 week", "7 week", "16 week", "Time in Open Field Centre", "% Sucrose Preference") 
samples_metadata(trained_model) <- metadata
# Plot
Fac1 <- plot_factors_vs_cov(trained_model,
                            factors = 1,
                            covariates = "week",
                            color_by = "strain",
                            legend = T)
Fac2 <- plot_factors_vs_cov(trained_model,
                            factors = 2,
                            covariates = "week",
                            color_by = "strain",
                            legend = T)
# arranging on the same plot
fin_plt <- ggarrange(Fac1, Fac2,
                     labels = c("(a)", "(b)"),
                     font.label = list(size = 11),
                     ncol = 2, nrow = 1,
                     common.legend = T, legend = "bottom")
# Saving
ggsave(filename = "MEFISTO_Scx_FactorPlts.tiff", plot = fin_plt, 
       width = 20, height = 16, units = "cm", device='tiff', dpi=300)

# Plotting weigh plots
# Changing the labels of the molecules
feat_name <- features_names(trained_model)
# Proteomics
p_lab <- read.csv("Prot_ED1_Genes_ProteinID_ForMEFISTO.csv", stringsAsFactors = F)
p_trained <- feat_name[["Proteomic"]]
p_labMatch <- p_lab[p_lab$Protein.IDs %in% p_trained,] 
p_labMatch <- p_labMatch[order(match(p_labMatch$Protein.IDs, p_trained)),] 
# Metabolomics
m_lab <- read.csv("Met_Scx_KEGG_CompoundName.csv", stringsAsFactors = F, header = 1)
m_trained <- feat_name[["Metabolomic"]]
m_labMatch <- m_lab[m_lab$KEGG.ID %in% m_trained,]
m_labMatch <- m_labMatch[order(match(m_labMatch$KEGG.ID, m_trained)),]
m_labMatch$Compound <- gsub("Dimethylallylpyrophosphate","Isopentenyl pyrophosphate", m_labMatch$Compound)
m_labMatch$Compound <- gsub("_",",", m_labMatch$Compound)
# Making a list
feat_list <- list(Metabolomic = m_labMatch$Compound,
                  Proteomic = p_labMatch$Gene.Name)
features_names(trained_model) <- feat_list

# Colouring in the features that were found to be common in the single-omics
# Proteomics
p <- plot_top_weights(trained_model,
                      view = "Proteomic",
                      factors = 2,
                      nfeatures = 20,
                      scale = T)
# Read in the table with the common findings
#p_common <- read.csv("Prot_All_ForMEFISTO.csv", stringsAsFactors = F, row.names = 1)
p_common <- read.csv("Prot_Scx_CommonSigDE_ForMEFISTOPlot.csv", 
                     stringsAsFactors = F, row.names = 1)
p_match <- rev(unlist(p[["data"]][["feature_id"]]) %in% row.names(p_common))
p[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(p_match == T, "deeppink1", "black") #deeppink1, darkorchid2

# Metabolomics
m <- plot_top_weights(trained_model,
                      view = "Metabolomic",
                      factors = 2,
                      nfeatures = 20,
                      scale = T)
#m_common <- read.csv("Met_Scx_CommonSigDE_ForMEFISTOPlot.csv", stringsAsFactors = F, row.names = 1)
m_common <- read.csv("Met_All_ForMEFISTO.csv", stringsAsFactors = F, row.names = 1)
row.names(m_common) <- gsub("Dimethylallylpyrophosphate","Isopentenyl pyrophosphate", row.names(m_common))
row.names(m_common) <- gsub("_",",", row.names(m_common))
m_match <- rev(unlist(m[["data"]][["feature_id"]]) %in% row.names(m_common))
m[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(m_match == T, "deeppink1", "black")

# arranging on the same plot
final_wght <- ggarrange(p, m,
                        labels = c("(a)", "(b)"),
                        font.label = list(size = 11),
                        ncol = 2, nrow = 1,
                        common.legend = T, legend = "bottom",
                        widths = c(1,1.5))
# Saving
ggsave(filename = "MEFISTO_Scx_AllDE_WeightPlts.tiff", plot = final_wght, 
       width = 20, height = 16, units = "cm", device='tiff', dpi=300)

#----MOFA Common Features----
# Loading libraries
library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggpubr)

## Reading in the trained model
filepath <- file.path(getwd(), "GAERS_MOFA_Input_Tha_TransFilt_20240502.hdf5")
trained_model <- load_model(filepath)

## Adding metadata to model
samp_metadata <- read.csv("Aim 1 - GAERS/Sample metadata.csv", stringsAsFactors = F, sep = ",")
samp_metadata$sample <- gsub("s$", "t", samp_metadata$sample)
metadata <- samp_metadata[samp_metadata$sample %in% samples_metadata(trained_model)$sample,]
metadata <- metadata[order(match(metadata$sample, samples_metadata(trained_model)$sample)),]
metadata <- metadata[,-c(4:7)]
colnames(metadata) <- c("sample", "strain", "strainCode", "Time in Open Field Centre", "% Sucrose Preference")
samples_metadata(trained_model) <- metadata

## Plotting factor values
Factor <- plot_factor(trained_model,
                      factors = 1,
                      color_by = "strain",
                      legend = T)

## Plotting Feature weights
# Changing the labels of the molecules
feat_name <- features_names(trained_model)
# Proteomics
p_lab <- read.csv("Prot_ED1_Genes_ProteinID_ForMEFISTO.csv", stringsAsFactors = F)
p_trained <- feat_name[["Proteomic"]]
p_labMatch <- p_lab[p_lab$Protein.IDs %in% p_trained,]
p_labMatch <- p_labMatch[order(match(p_labMatch$Protein.IDs, p_trained)),]
# Metabolomics
m_lab <- read.csv("Met_Tha_KEGG_CompoundName.csv", stringsAsFactors = F, header = 1)
m_trained <- feat_name[["Metabolomic"]]
m_labMatch <- m_lab[m_lab$KEGG.ID %in% m_trained,]
m_labMatch <- m_labMatch[order(match(m_labMatch$KEGG.ID, m_trained)),]
m_labMatch$Compound <- gsub("Dimethylallylpyrophosphate","Isopentenyl pyrophosphate", m_labMatch$Compound)
m_labMatch$Compound <- gsub("_",",", m_labMatch$Compound)
# Making a list
feat_list <- list(Metabolomic = m_labMatch$Compound,
                  Proteomic = p_labMatch$Gene.Name,
                  Transcriptomic = feat_name$Transcriptomic)
# Changing feature name
features_names(trained_model) <- feat_list

# Colouring in the features of interest
# Proteomics
p <- plot_top_weights(trained_model,
                      view = "Proteomic",
                      factors = 1,
                      nfeatures = 20,
                      scale = T)
# Highlight proteins of interest
p_common <- c("Nit2", "Ak1", "Tmlhe")
p_match <- rev(unlist(p[["data"]][["feature_id"]]) %in% p_common)
p[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(p_match == T, "red", "black")

# Metabolomics
m <- plot_top_weights(trained_model,
                      view = "Metabolomic",
                      factors = 1,
                      nfeatures = 20,
                      scale = T)
# To highlight only seubset of interest
m_common <- c("2-Keto-glutaramic acid", "Isopentenyl pyrophosphate", "Cystathionine")
m_match <- rev(unlist(m[["data"]][["feature_id"]]) %in% m_common)
m[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(m_match == T, "red", "black")

# Transcriptomics
t <- plot_top_weights(trained_model,
                      view = "Transcriptomic",
                      factors = 1,
                      nfeatures = 20,
                      scale = T)
# Highlighting NIT2
t_common <- c("Nit2")
t_match <- rev(unlist(t[["data"]][["feature_id"]]) %in% t_common)
t[["theme"]][["axis.text.y"]][["colour"]] <- ifelse(t_match == T, "red", "black")

# arranging on the same plot
final_plts <- ggarrange(Factor, t, p, m,
                        labels = c("(a)", "(b)", "(c)", "(d)"),
                        font.label = list(size = 11),
                        ncol = 2, nrow = 2)
# Saving
ggsave(filename = "MOFA_AllAnalysis.tiff", plot = final_plts, 
       width = 30, height = 24, units = "cm", device='tiff', dpi=220)

#----Bubble Plots for Enrichment Analysis----
# Loading libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Reading in the files
df <- read.csv("PTE Enrichment Analysis/PTE_PTETBIOnly_ReactomeGSA_Protein.csv", 
               stringsAsFactors = F, row.names = 1)
df <- read.csv("MA_PTETBIOnly_Enrichment/msea_qea_result.csv", stringsAsFactors = F)

# Removing any terms with less than 10 hits
# Reactome
df1 <- df %>%
  filter(NGenes.Transcriptomics > 10) %>%
  arrange(FDR.Transcriptomics)
df1 <- df %>%
  filter(NGenes.Proteomics > 10) %>%
  arrange(FDR.Proteomics)
# MetaboAnalyst
df1 <- df %>%
  filter(Hits > 10) %>%
  arrange(Raw.p)

# Taking the top 20
df2 <- df1[1:20,]
# New lines for transcriptomics
df2$Name[5] <- "Nonsense Mediated Decay (NMD) independent \nof the Exon Junction Complex (EJC)"
df2$Name[7] <- "Energy dependent regulation of \nmTOR by LKB1-AMPK"  
df2$Name[8] <- "Activated PKN1 stimulates transcription of AR \n(androgen receptor) regulated genes KLK2 and KLK3"
df2$Name[16] <- "Response of EIF2AK4 (GCN2) \nto amino acid deficiency" 
df2$Name[20] <- "Nonsense Mediated Decay (NMD) enhanced by \nthe Exon Junction Complex (EJC)"
# Putting in new lines for proteomics
df2$Name[2] <- "HSP90 chaperone cycle for steroid hormone \nreceptors (SHR) in the presence of ligand"
df2$Name[15] <- "Factors involved in megakaryocyte \ndevelopment and platelet production"
df2$Name[17] <- "Regulation of actin dynamics for \nphagocytic cup formation" 
# Metabolomics
df2$X[13] <- "Arginine: Glycine Amidinotransferase \nDeficiency (AGAT Deficiency)" 
df2$X[14] <- "Creatine deficiency, guanidinoacetate \nmethyltransferase deficiency"
df2$X[15] <- "Guanidinoacetate Methyltransferase \nDeficiency (GAMT Deficiency)"   
df2$X[17] <- "Hyperornithinemia-hyperammonemia-homocitrullinuria \n[HHH-syndrome]" 

# Making the plots
trans_plt <- ggplot(df2, aes(x = FDR.Transcriptomics,
                             y = reorder(Name,FDR.Transcriptomics))) +
  geom_point(aes(size = NGenes.Transcriptomics,
                 color = FDR.Transcriptomics),
             alpha = 0.7) +
  geom_text(label = df2$NGenes.Transcriptomics, nudge_x = 0.002) +
  labs(size = "Number of Genes", color = "False Discovery Rate") +
  xlab("False Discovery Rate") +
  ylab("") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9))
# For Proteomics
prot_plt <- ggplot(df2, aes(x = FDR.Proteomics,
                            y = reorder(Name,FDR.Proteomics))) +
  geom_point(aes(size = NGenes.Proteomics,
                 color = FDR.Proteomics),
             alpha = 0.7) +
  geom_text_repel(label = df2$NGenes.Proteomics, nudge_x = 0.001) +
  labs(size = "Number of Genes", color = "False Discovery Rate") +
  xlab("False Discovery Rate") +
  ylab("") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9))
# Metabolomics
met_plt <- ggplot(df2, aes(x = -log10(Raw.p),
                           y = reorder(X,-log10(Raw.p)))) +
  geom_point(aes(size = Hits,
                 color = -log10(Raw.p)),
             alpha = 0.7) +
  geom_text(label = df2$Hits, nudge_x = 0.17) +
  labs(size = "Number of Genes", color = expression("-log"[10]*"(p-Value)")) +
  xlab(expression("-log"[10]*"(p-Value)")) +
  ylab("") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9))

# As required
# arranging on the same plot
grp_plts <- ggarrange(trans_plt, prot_plt,
                      labels = c("(a)", "(b)"),
                      font.label = list(size = 11),
                      ncol = 2, nrow = 1,
                      common.legend = T,
                      legend = "bottom")

# Saving
ggsave(filename = "TransProt_BubblePlt.tiff", plot = grp_plts, 
       width = 24, height = 18, units = "cm", device='tiff', dpi=300)

#----Alpha-ketoglutarate BoxPlots----
# Loading libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Reading in dataframe
wk_16 <- read.csv("Aim 1 - GAERS/Met_Scx_16wk_noOut_timpt.csv", stringsAsFactors = F)
wk_7 <- read.csv("Aim 1 - GAERS/Met_Tha_7wk_noOut_timpt.csv", stringsAsFactors = F)
wk_3 <- read.csv("Aim 1 - GAERS/Met_Tha_3wk_noOut_timpt.csv", stringsAsFactors = F)

# Selecting alpha-ketoglutarate (C00026)
lab <- c("X","C00026")
wk_3t <- wk_3 %>%
  select(lab) %>%
  mutate(strain = c(rep("GAERS",10), rep("NEC",10)), .before = 2) %>%
  # 16 week_scx: rep("GAERS",10), rep("NEC",12)
  # 7 week_scx: rep("GAERS",8), rep("NEC",11)
  #3 week_scx: rep("GAERS",11), rep("NEC",11)
  # 16 week_tha: rep("GAERS",10), rep("NEC",11)
  # 7 week_tha: rep("GAERS",8), rep("NEC",11)
  #3 week_tha: rep("GAERS",10), rep("NEC",10)
  pivot_longer(!c("X","strain"),
               names_to = "Metabolite",
               values_to = "value") %>%
  mutate(Timpt_BR = "3week_Tha", .before = 3) %>%
  arrange(Metabolite)

# Joining together in long dataframe (for faceting)
dt_plt <- rbind(wk_3s, wk_7s, wk_16s,
                wk_3t, wk_7t, wk_16t)
# Make boxplot
box_plt <- dt_plt %>%
  ggplot(aes(x=strain, y=value, fill = strain, alpha = 0.5)) +
  geom_boxplot() + #geom_violin for a violin plot
  facet_wrap(~factor(Timpt_BR, c("3week_Scx", "7week_Scx", "16week_Scx", 
                                 "3week_Tha", "7week_Tha", "16week_Tha"))) +
  geom_jitter(color="black", size=1, alpha=1, 
              position = position_jitter(seed = 1)) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_blank()
  ) +
  ylab("Normalised Concentration")

# arranging on the same plot (if making individual boxplots)
final_box <- ggarrange(wk3_plt, wk7_plt, wk16_plt,
                       labels = c("(a)", "(b)", "(c)"),
                       font.label = list(size = 11),
                       ncol = 2, nrow = 2)
# Saving
ggsave(filename = "Alpha-ketoglutarate_AllGrps.tiff", plot = box_plt, 
       width = 20, height = 20, units = "cm", device='tiff', dpi=300)
save.image(file = "GAERS_AlphaKetoLevels.RData")