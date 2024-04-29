setwd("C:/Dump/")

read_lengths = read.table("NLR_read_lengths_fig2S1.tsv", header = T)

library(ggplot2)
library(RColorBrewer)
library(forcats)
library(dplyr)

read_lengths_reordered = read_lengths %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "DP", "KG", "SN", "CHT"))

min(read_lengths$Read_Length)
max(read_lengths$Read_Length)
mean(read_lengths[read_lengths$Population=="CHT","Read_Length"])
mean(read_lengths[read_lengths$Population=="DP","Read_Length"])


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Read_lengths.pdf", width = 15, height = 9)

ggplot(data = read_lengths_reordered, aes(x = pop_names, y = Read_Length, fill = pop_names)) +
  # geom_violin(show.legend = F) +
  geom_boxplot(show.legend = F, outlier.shape = NA, size = 1) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4700)) +
  scale_fill_manual(values = c("red","orange","white","royalblue3","deepskyblue","green")) +
  ylab("Length [bp]") + xlab("Population") +
  theme(axis.text=element_text(size=40), 
        axis.title=element_text(size=45,face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25, face = "bold")) +
  geom_point(stat = "summary", size = 8, shape = 18, fill = "black") +
  geom_vline(xintercept = 2.5, linewidth = 2, linetype = "dashed")

dev.off()
