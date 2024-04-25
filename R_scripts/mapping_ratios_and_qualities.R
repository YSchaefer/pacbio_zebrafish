setwd("C:/Dump/")

mapping_ratios = read.table("mapping_ratios_fig2S2.tsv", header = T)

NLR_mapping_stats = read.table("orthocluster_mapping_stats_fig2S2.tsv", header = T)

no_none = NLR_mapping_stats[NLR_mapping_stats$Mapping_Quality!="None",]

no_none$Mapping_Quality = as.numeric(no_none$Mapping_Quality)

library(ggplot2)
library(forcats)
library(dplyr)

lab_wild_reordered = mapping_ratios %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "DP", "KG", "SN", "CHT"))

exons_reordered = lab_wild_reordered %>%
  mutate(Exon_new = fct_relevel(Exon, "NLR_B30.2", "FISNACHT"))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Mapping_reatios.pdf", width = 15, height = 9)

ggplot(exons_reordered, aes(x = pop_names, y = Mapping_Ratio, fill = Exon_new)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", show.legend = T, size = 1) +
  theme_bw() +
  scale_fill_manual(values = c("white", "darkgrey"), labels = c("NLR-B30.2", "FISNACHT")) + labs(fill = "Exon") +
  xlab("Population") + ylab("Mapping ratio") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  theme(
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 45, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 40),
    legend.background = element_rect(colour = "black"),
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
  ) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2)

dev.off()


no_none_reordered = no_none %>%
  mutate(Exon_new = fct_relevel(Exon, "NLR-B30.2", "FISNACHT"))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Mapping_qualities.pdf", width = 15, height = 9)

ggplot(no_none_reordered, aes(x = Mapping_Quality, fill = Exon_new)) +
  geom_histogram(binwidth = 1, position = position_dodge(), colour = "black", show.legend = F, size = 1) +
  theme_bw() +
  scale_fill_manual(values = c("white", "darkgrey"), labels = c("FISNACHT", "NLR-B30.2")) +
  xlab("Mapping quality") + ylab("Unique sequences") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 61)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 820)) +
  theme(
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 45, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 40),
    legend.background = element_rect(colour = "black"),
    axis.title = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 40),
  )

dev.off()
