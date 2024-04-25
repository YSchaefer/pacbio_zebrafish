setwd("C:/Dump/")

read_percentages = read.table("all_read_percentages_fig2S1.tsv", header = T)

library(ggplot2)
library(dplyr)
library(forcats)
library(ggbeeswarm)

lab_wild_reordered = read_percentages %>%
  mutate(pop_names = fct_relevel(Population, "TU", "CGN", "DP", "KG", "SN", "CHT"))


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/Read_counts.pdf", width = 15, height = 9)

ggplot(data = lab_wild_reordered, aes(x = pop_names, y = NLR_Reads, colour = pop_names)) +
  geom_beeswarm(size = 6, show.legend = F) +
  theme_bw() +
  scale_colour_manual(values = c("red","orange","black","royalblue3","deepskyblue","green")) +
  ylab("Count") + xlab("Population") +
  ggtitle("NLR reads with at least 3 passes") +
  scale_y_continuous(expand = c(0, 1000), limits = c(0, 103000)) +
  theme(axis.text=element_text(size=40), 
        axis.title=element_text(size=45,face = "bold"),
        plot.title = element_text(size = 50, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  geom_vline(xintercept = 2.5, linewidth = 2, linetype = "dashed")

dev.off()
