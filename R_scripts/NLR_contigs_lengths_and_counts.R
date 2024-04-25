setwd("C:/Dump")

NLR_contig_lengths = read.table("NLR_contig_lengths_fig2S1.tsv", header = T)

library(ggplot2)
library(RColorBrewer)
library(forcats)
library(dplyr)

NLR_contig_lengths_reordered = NLR_contig_lengths %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "DP", "KG", "SN", "CHT"))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Contig_lengths.pdf", width = 15, height = 9)

ggplot(NLR_contig_lengths_reordered, aes(x = pop_names, y = Contig_Length, fill = Domain)) +
  geom_boxplot(show.legend = F, outlier.shape = NA, size = 1) +
  theme_bw() +
  scale_fill_manual(values = c("white", "darkgrey"), labels = c("B30.2","FISNA-NACHT")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10200)) +
  xlab("Population") + ylab("Length [bp]") +
  theme(axis.text=element_text(size=40), 
        axis.title=element_text(size=45,face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.background = element_rect(colour = "black"),
        axis.title.y = element_blank()) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2)

dev.off()


# NLR Contig Counts

NLR_contig_counts = read.table("NLR_contig_counts_fig2S1.tsv", header = T)

NLR_contig_counts_reordered = NLR_contig_counts %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "DP", "KG", "SN", "CHT"))

NLR_contig_counts_reordered$Domain = gsub("FISNACHT","FISNA-NACHT",NLR_contig_counts_reordered$Domain)


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Contig_counts.pdf", width = 15, height = 9)

ggplot(NLR_contig_counts_reordered, aes(x = pop_names, y = Contig_Count, fill = Domain)) +
  geom_boxplot(show.legend = F, outlier.shape = NA, size = 1) +
  theme_bw() +
  ggtitle("Contigs of interest") +
  scale_fill_manual(values = c("white", "darkgrey")) +
  scale_y_continuous(expand = c(0, 10), limits = c(0, 900)) +
  xlab("Population") + ylab("Contigs") +
  theme(axis.text=element_text(size=40), 
        axis.title=element_text(size=45,face = "bold"),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 35, face = "bold"),
        legend.background = element_rect(colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 50, face = "bold", hjust = 0.5)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_point(aes(x = "TU", y = 378), shape = 24, colour = "black", fill = "darkgrey", size = 4, stroke = 1.5,
             show.legend = F, position = position_nudge(x = 0.2)) +
  geom_point(aes(x = "TU", y = 608), shape = 24, colour = "black", fill = "white", size = 4, stroke = 1.5,
             show.legend = F, position = position_nudge(x = -0.2)) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2)

dev.off()


# Counts of B30.2 contigs in NLRs and otherwise, informed by the special exon

B30.2_contig_counts = read.table("B30.2_contig_counts_fig2S3.tsv", header = T)

B30.2_contig_counts_reordered = B30.2_contig_counts %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "DP", "KG", "SN", "CHT"))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/B30.2_contig_counts.pdf", width = 15, height = 9)

ggplot(B30.2_contig_counts_reordered, aes(x = pop_names, y = Contig_Count, fill = Class)) +
  geom_boxplot(show.legend = T, outlier.shape = NA, size = 1) +
  theme_bw() +
  scale_fill_manual(values = c("white", "darkgrey")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 450)) +
  xlab("Population") + ylab("Contigs with B30.2") + labs(fill = "Gene") +
  theme(axis.text=element_text(size=40), 
        axis.title=element_text(size=45,face = "bold"),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 35, face = "bold", hjust = 0.5),
        legend.background = element_rect(colour = "black"),
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.title.y = element_blank(),
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 8, color = "black",
               position = position_dodge2(width = 0.75, preserve = "single"), show.legend = F) +
  geom_vline(xintercept = 2.5, linetype = "dashed", linewidth = 2)

dev.off()
