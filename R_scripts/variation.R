# Variation statistics from vcftools

setwd("C:/Dump/")

library(dplyr)
library(forcats)
library(ggplot2)
library(scales)
library(tidyr)

B30.2_statistics = read.table("B30.2_variation_fig4.tsv", header = T)
FISNACHT_statistics = read.table("FISNACHT_variation_fig4.tsv", header = T)

B30.2_statistics$Class = rep("Unknown", times = 4351)

both_variation = rbind(FISNACHT_statistics, B30.2_statistics)

variation_without_DP = both_variation[both_variation$Population != "DP",]

lab_wild_reordered = both_variation[both_variation$No_Fish >= 4,] %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "DP", "KG", "SN", "CHT"))

# Variation parameters by exon

reshaped_data <- lab_wild_reordered %>%
  gather(Parameter, Value, Watterson, Pi)

pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/FISNACHT_variation_parameters.pdf", width = 15, height = 9)

ggplot(data = reshaped_data[reshaped_data$Group != "no" & reshaped_data$Value != 0 & reshaped_data$Exon == "FISNACHT",], aes(x = pop_names, y = Value, fill = Parameter)) +
  geom_boxplot(size = 1, show.legend = F, outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = c("yellow","blue3"), labels=c("Pi","Watterson")) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.0005, 0.028)) +
  xlab("Population") + ylab("Value") + ggtitle("FISNA-NACHT") + labs(fill = "Estimator") +
  theme(
    plot.title = element_text(size = 50, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 40, face = "bold"),
    legend.text = element_text(size = 35),
    legend.background = element_rect(colour = "black"),
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "white")
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2, colour = "black")

dev.off()

summary(aov(data = FISNACHT_statistics, Pi ~ Population))
TukeyHSD(aov(data = FISNACHT_statistics, Pi ~ Population))

summary(aov(data = FISNACHT_statistics, Watterson ~ Population))
TukeyHSD(aov(data = FISNACHT_statistics, Watterson ~ Population))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/B30.2_variation_parameters.pdf", width = 15, height = 9)

ggplot(data = reshaped_data[reshaped_data$Group != "no" & reshaped_data$Value != 0 & reshaped_data$Exon == "B30.2",], aes(x = pop_names, y = Value, fill = Parameter)) +
  geom_boxplot(size = 1, show.legend = F, outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = c("yellow","blue3"), labels=c(expression("Nucleotide Diversity \u03B8"["\u03C0"]), expression("Watterson's \u03B8"[w]))) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.0005, 0.028)) +
  xlab("Population") + ylab("Value") + ggtitle("NLR-B30.2") + labs(fill = "Estimator") +
  theme(
    plot.title = element_text(size = 50, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 35, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 30, hjust = 0),
    legend.background = element_rect(colour = "black"),
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "white")
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2, colour = "black")

dev.off()


summary(aov(data = B30.2_statistics, Pi ~ Population))
TukeyHSD(aov(data = B30.2_statistics, Pi ~ Population))

summary(aov(data = B30.2_statistics, Watterson ~ Population))
TukeyHSD(aov(data = B30.2_statistics, Watterson ~ Population))


# Fractions of exons without variation

nrow(lab_wild_reordered[lab_wild_reordered$Population == "SN" &
                     lab_wild_reordered$Group != "no" &
                     lab_wild_reordered$Exon == "FISNACHT" &
                     lab_wild_reordered$Watterson == 0 &
                     lab_wild_reordered$Pi == 0,])

zero_variation = c(108,171,295,417,463,423,37,66,75,148,116,131)
Population = c(rep(c("TU","CGN","DP","KG","SN","CHT"),2))
Exon = c(rep("FISNACHT",6),rep("B30.2",6))

nrow(lab_wild_reordered[lab_wild_reordered$Population == "TU" &
                        lab_wild_reordered$Exon == "FISNACHT" &
                        lab_wild_reordered$No_Fish >= 4,])

Total = c(192,267,566,788,850,714,180,244,328,599,628,573)

zero_variation_df = data.frame(Exon, Population, Total, Zeros = zero_variation)

zeros_df_reordered = zero_variation_df %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "DP", "KG", "SN", "CHT"))

pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Exons_without_variation.pdf", width = 15, height = 9)

ggplot(data = zeros_df_reordered, aes(x = pop_names, y = Zeros/Total, fill = Exon)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", show.legend = F, size = 1) +
  scale_fill_manual(values = c("white","darkgrey")) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.9)) +
  xlab("Population") + ylab("Exons without variation") +
  theme(
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 35, face = "bold"),
    legend.text = element_text(size = 30),
    legend.background = element_rect(colour = "black")
  ) +
  geom_vline(xintercept = 2.5, linewidth = 2, linetype = "dashed", colour = "black")

dev.off()


# Plot of the ratio of thetas

pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Theta_ratios.pdf", width = 15, height = 9)

ggplot(data = lab_wild_reordered[lab_wild_reordered$Group != "no" & lab_wild_reordered$Watterson != 0,], aes(x = pop_names, y = Ratio, fill = Exon)) +
  geom_boxplot(size = 1, show.legend = F, outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = c("white","darkgrey"), labels=c("NLR-B30.2","FISNACHT")) +
  scale_y_log10(limits = c(0.2, 1.8)) +
  xlab("Population") + ylab(expression("\u03B8"["\u03C0"]*"/\u03B8"[w])) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1) +
  theme(
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    axis.title.y = element_text(size = 40, face = "bold"),
    legend.title = element_text(size = 35, face = "bold"),
    legend.text = element_text(size = 30),
    legend.background = element_rect(colour = "black"),
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "white")
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2, colour = "black")

dev.off()

expression("\u03B8"["\u03C0"]*"/\u03B8"[w])

summary(aov(data = FISNACHT_statistics, Ratio ~ Population))
TukeyHSD(aov(data = FISNACHT_statistics, Ratio ~ Population))

summary(aov(data = B30.2_statistics, Ratio ~ Population))
TukeyHSD(aov(data = B30.2_statistics, Ratio ~ Population))


# repeat without DP

variation_without_DP = both_variation[both_variation$Population != "DP",]

lab_wild_reordered = variation_without_DP[variation_without_DP$No_Fish >= 4,] %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "KG", "SN", "CHT"))


# Variation parameters by exon

reshaped_data <- lab_wild_reordered %>%
  gather(Parameter, Value, Watterson, Pi)

pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/FISNACHT_variation_parameters_noDP.pdf", width = 15, height = 9)

ggplot(data = reshaped_data[reshaped_data$Group != "no" & reshaped_data$Value != 0 & reshaped_data$Exon == "FISNACHT",], aes(x = pop_names, y = Value, fill = Parameter)) +
  geom_boxplot(size = 1, show.legend = F, outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = c("yellow","blue3"), labels=c("Pi","Watterson")) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.0005, 0.028)) +
  xlab("Population") + ylab("Value") + ggtitle("FISNA-NACHT") + labs(fill = "Estimator") +
  theme(
    plot.title = element_text(size = 50, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 40, face = "bold"),
    legend.text = element_text(size = 35),
    legend.background = element_rect(colour = "black"),
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "white")
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2, colour = "black")

dev.off()

summary(aov(data = FISNACHT_statistics, Pi ~ Population))
TukeyHSD(aov(data = FISNACHT_statistics, Pi ~ Population))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/B30.2_variation_parameters_noDP.pdf", width = 15, height = 9)

ggplot(data = reshaped_data[reshaped_data$Group != "no" & reshaped_data$Value != 0 & reshaped_data$Exon == "B30.2",], aes(x = pop_names, y = Value, fill = Parameter)) +
  geom_boxplot(size = 1, show.legend = F, outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = c("yellow","blue3"), labels=c(expression("Nucleotide Diversity \u03B8"["\u03C0"]), expression("Watterson's \u03B8"[w]))) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.0005, 0.028)) +
  xlab("Population") + ylab("Value") + ggtitle("NLR-B30.2") + labs(fill = "Estimator") +
  theme(
    plot.title = element_text(size = 50, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 35, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 30, hjust = 0),
    legend.background = element_rect(colour = "black"),
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "white")
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2, colour = "black")

dev.off()


# Fractions of exons without variation

zero_variation = c(108,171,417,463,423,37,66,148,116,131)
Population = c(rep(c("TU","CGN","KG","SN","CHT"),2))
Exon = c(rep("FISNACHT",5),rep("B30.2",5))

Total = c(192,267,788,850,714,180,244,599,628,573)

zero_variation_df = data.frame(Exon, Population, Total, Zeros = zero_variation)

zeros_df_reordered = zero_variation_df %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "KG", "SN", "CHT"))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Exons_without_variation_noDP.pdf", width = 15, height = 9)

ggplot(data = zeros_df_reordered, aes(x = pop_names, y = Zeros/Total, fill = Exon)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", show.legend = F, size = 1) +
  scale_fill_manual(values = c("white","darkgrey")) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.9)) +
  xlab("Population") + ylab("Exons without variation") +
  theme(
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    legend.title = element_text(size = 40, face = "bold"),
    legend.text = element_text(size = 35),
    legend.background = element_rect(colour = "black")
  ) +
  geom_vline(xintercept = 2.5, linewidth = 2, linetype = "dashed", colour = "black")

dev.off()


# Plot of the ratio of thetas

pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Theta_ratios_noDP.pdf", width = 15, height = 9)

ggplot(data = lab_wild_reordered[lab_wild_reordered$Group != "no" & lab_wild_reordered$Watterson != 0,], aes(x = pop_names, y = Ratio, fill = Exon)) +
  geom_boxplot(size = 1, show.legend = F, outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = c("white","darkgrey"), labels=c("NLR-B30.2","FISNACHT")) +
  scale_y_log10(limits = c(0.2, 1.8)) +
  xlab("Population") + ylab(expression("\u03B8"["\u03C0"]*"/\u03B8"[w])) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1) +
  theme(
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
    axis.title.y = element_text(size = 40, face = "bold"),
    legend.title = element_text(size = 35, face = "bold"),
    legend.text = element_text(size = 30),
    legend.background = element_rect(colour = "black"),
    strip.text = element_text(size = 30, face = "bold"),
    strip.background = element_rect(fill = "white")
  ) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2, colour = "black")

dev.off()

expression("\u03B8"["\u03C0"]*"/\u03B8"[w])