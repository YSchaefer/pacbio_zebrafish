setwd("C:/Dump/")

library(ggplot2)
library(forcats)
library(patchwork)
library(dplyr)
library(RColorBrewer)


orthocluster_counts = read.table("orthocluster_counts.tsv", header = T)

lab_wild_reordered = orthocluster_counts %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "DP", "KG", "SN", "CHT"))

lab_wild_reordered$Exon <- factor(lab_wild_reordered$Exon,
                         levels = c("B30.2", "FISNACHT", "NLR_B30.2", "other_B30.2"))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/FISNACHT_cluster_counts.pdf", width = 15, height = 9)

ggplot(data = lab_wild_reordered[lab_wild_reordered$Exon == "FISNACHT",], 
       aes(x = pop_names, y = Orthoclusters, fill = pop_names)) +
  geom_boxplot(outlier.shape = NA, show.legend = F, size = 1) +
  ggbeeswarm::geom_beeswarm(size = 5, show.legend = F) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 580)) +
  xlab("Population") + ylab("NLR-C genes") + #ggtitle("FISNACHT") +
  scale_fill_manual(values = c("red","orange", "white", "royalblue3","deepskyblue","green")) +
  theme(
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    legend.background = element_rect(colour = "black"),
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", linewidth = 2)

dev.off()

summary(aov(Orthoclusters ~ Population, data = lab_wild_reordered))
TukeyHSD(aov(Orthoclusters ~ Population, data = lab_wild_reordered))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/B30.2_cluster_counts.pdf", width = 15, height = 9)

ggplot(data = lab_wild_reordered[lab_wild_reordered$Exon == "NLR_B30.2",], 
       aes(x = pop_names, y = Orthoclusters, fill = pop_names)) +
  geom_boxplot(outlier.shape = NA, show.legend = F, size = 1) +
  ggbeeswarm::geom_beeswarm(size = 5, show.legend = F) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 220)) +
  scale_fill_manual(values = c("red","orange", "white", "royalblue3","deepskyblue","green")) +
  xlab("Population") + ylab("Unique sequences") + #ggtitle("NLR-B30.2") +
  theme(
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    axis.title.y = element_blank(),
    #axis.title.x = element_blank(),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    legend.background = element_rect(colour = "black"),
    legend.title = element_text(size = 35, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 30)
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", linewidth = 2)

dev.off()


# Repeat without DP

counts_without_DP = orthocluster_counts[orthocluster_counts$Population != "DP",]

lab_wild_reordered = counts_without_DP %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "KG", "SN", "CHT"))

lab_wild_reordered$Exon <- factor(lab_wild_reordered$Exon,
                                  levels = c("B30.2", "FISNACHT", "NLR_B30.2", "other_B30.2"))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/FISNACHT_cluster_counts_noDP.pdf", width = 15, height = 9)

ggplot(data = lab_wild_reordered[lab_wild_reordered$Exon == "FISNACHT",], 
       aes(x = pop_names, y = Orthoclusters, fill = pop_names)) +
  geom_boxplot(outlier.shape = NA, show.legend = F, size = 1) +
  ggbeeswarm::geom_beeswarm(size = 5, show.legend = F) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 580)) +
  xlab("Population") + ylab("NLR-C genes") + #ggtitle("FISNACHT") +
  scale_fill_manual(values = c("red","orange", "royalblue3","deepskyblue","green")) +
  theme(
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    legend.background = element_rect(colour = "black"),
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", linewidth = 2)

dev.off()

summary(aov(Orthoclusters ~ Population, data = lab_wild_reordered[lab_wild_reordered$Exon == "FISNACHT",]))
TukeyHSD(aov(Orthoclusters ~ Population, data = lab_wild_reordered[lab_wild_reordered$Exon == "FISNACHT",]))

# All lab strains differ significantly from all wild populations


ggplot(data = lab_wild_reordered[lab_wild_reordered$Exon == "NLR_B30.2",], 
       aes(x = pop_names, y = Orthoclusters, fill = pop_names)) +
  geom_boxplot(outlier.shape = NA, show.legend = F, size = 1) +
  ggbeeswarm::geom_beeswarm(size = 5, show.legend = F) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 220)) +
  xlab("Population") + ylab("NLR-C genes") + #ggtitle("FISNACHT") +
  scale_fill_manual(values = c("red","orange", "royalblue3","deepskyblue","green")) +
  theme(
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    legend.background = element_rect(colour = "black"),
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", linewidth = 2)

summary(aov(Orthoclusters ~ Population, data = lab_wild_reordered[lab_wild_reordered$Exon == "NLR_B30.2",]))
TukeyHSD(aov(Orthoclusters ~ Population, data = lab_wild_reordered[lab_wild_reordered$Exon == "NLR_B30.2",]))
