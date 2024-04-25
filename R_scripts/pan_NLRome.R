setwd("C:/Dump/")

library(ggplot2)

three_wilds = read.table("three_wild_pan_fig3.tsv")
four_wilds = read.table("four_wild_pan_fig3S2.tsv")


scale = c(0,0.2,0.4,0.6,0.8,1)


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/FISNACHT_three_pans.pdf", width = 15, height = 9)

ggplot(data = three_wilds, aes(x = V2, y = V1)) +
  geom_point(size = 3) +
  geom_smooth(linewidth = 2, colour = "black", method = "loess", level = 0.95) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0,90)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), labels = scale, breaks = scale) +
  xlab("Relative Amount of Fish") + ylab("NLR Genes") +
  geom_vline(xintercept = 0.2) + geom_vline(xintercept = 0.8) +
  theme(
    axis.title = element_text(size = 35, face = "bold"),
    axis.text = element_text(size = 30)
  )

dev.off()

pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/FISNACHT_four_pans.pdf", width = 15, height = 9)

ggplot(data = four_wilds, aes(x = V1, y = V2)) +
  geom_point(size = 3) +
  geom_smooth(linewidth = 2, colour = "black", method = "loess", level = 0.95) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0,85)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), labels = scale, breaks = scale) +
  xlab("Relative Amount of Fish") + ylab("NLR Genes") +
  geom_vline(xintercept = 0.2) + geom_vline(xintercept = 0.8) +
  theme(
    axis.title = element_text(size = 35, face = "bold"),
    axis.text = element_text(size = 30)
  )

dev.off()