# Heat maps with ggplot2 geom_tile

setwd("C:/Dump/")

library(ggplot2)
library(dplyr)
library(forcats)

reduced_FISNACHT_data = read.table("FISNACHT.heatmap_data_fig3.tsv", header = T)


pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/FISNACHT_heatmap.pdf", width = 15, height = 9)

ggplot(data = reduced_FISNACHT_data, aes(x = reorder(Fish, Order), y = reorder(Cluster, Relative_Depth), fill = Relative_Depth)) +
  geom_tile(show.legend = T) +
  theme_minimal() +
  scale_fill_gradient(low="white", high="black", limit = c(0,0.0001)) +
  labs(x ="Fish", y ="Cluster", fill = "Relative Depth") +
  theme(
    legend.title = element_text(size = 35, face = "bold"),
    legend.text = element_text(size = 30),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

dev.off()