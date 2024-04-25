## For B30.2

setwd("C:/Dump")


library(ggvenn)

zfdat <- read.table("NLR-B30.2.venn_diagram_fig2S1.tsv", header = T)
zfdat <- zfdat[rowSums(zfdat[,c(7:99)]) > 0,] # exclude contigs not supported by any HiFi reads
zfcounts <- zfdat[,c(1:6)] # create a new table for population counts, then fill it in the next steps
zfcounts$CGN <- rowSums(zfdat[,colnames(zfdat)[grepl("CGN[0-9]",colnames(zfdat)) == T]])
zfcounts$TU <- rowSums(zfdat[,colnames(zfdat)[grepl("TU[0-9]",colnames(zfdat)) == T]])
zfcounts$CHT <- rowSums(zfdat[,colnames(zfdat)[grepl("CHT[0-9]",colnames(zfdat)) == T]])
zfcounts$DP <- rowSums(zfdat[,colnames(zfdat)[grepl("DP[0-9]",colnames(zfdat)) == T]])
zfcounts$KG <- rowSums(zfdat[,colnames(zfdat)[grepl("KG[0-9]",colnames(zfdat)) == T]])
zfcounts$SN <- rowSums(zfdat[,colnames(zfdat)[grepl("SN[0-9]",colnames(zfdat)) == T]])
zfcounts$nlrsum <- rowSums(zfdat[,c(7:99)])

#plot lists as venn charts
lists.wild <- list("KG" = zfcounts[zfcounts$KG > 0, 1], "DP" = zfcounts[zfcounts$DP > 0, 1], "SN" = zfcounts[zfcounts$SN > 0, 1], "CHT" = zfcounts[zfcounts$CHT > 0, 1])

pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/B30.2_venn_four.pdf", width = 15, height = 9)

ggvenn(lists.wild, fill_color = c("lightblue", "khaki", "lightpink", "lightgrey"), text_size = 7, stroke_size = 2, set_name_size = 12)

dev.off()

## For FISNACHT

zfdat <- read.table("FISNACHT.venn_diagram_fig2S1.tsv", header = T)
zfdat <- zfdat[rowSums(zfdat[,c(8:100)]) > 0,] # exclude contigs not supported by any HiFi reads
zfcounts <- zfdat[,c(1:7)] # create a new table for population counts, then fill it in the next steps
zfcounts$CGN <- rowSums(zfdat[,colnames(zfdat)[grepl("CGN[0-9]",colnames(zfdat)) == T]])
zfcounts$TU <- rowSums(zfdat[,colnames(zfdat)[grepl("TU[0-9]",colnames(zfdat)) == T]])
zfcounts$CHT <- rowSums(zfdat[,colnames(zfdat)[grepl("CHT[0-9]",colnames(zfdat)) == T]])
zfcounts$DP <- rowSums(zfdat[,colnames(zfdat)[grepl("DP[0-9]",colnames(zfdat)) == T]])
zfcounts$KG <- rowSums(zfdat[,colnames(zfdat)[grepl("KG[0-9]",colnames(zfdat)) == T]])
zfcounts$SN <- rowSums(zfdat[,colnames(zfdat)[grepl("SN[0-9]",colnames(zfdat)) == T]])
zfcounts$nlrsum <- rowSums(zfdat[,c(8:100)])

#plot lists as venn charts

lists.wild <- list("KG" = zfcounts[zfcounts$KG > 0, 1], "SN" = zfcounts[zfcounts$SN > 0, 1], "DP" = zfcounts[zfcounts$DP > 0, 1], "CHT" = zfcounts[zfcounts$CHT > 0, 1])

pdf("C:/Uni/Papers/Our_papers/PacBio_Zebrafish/Figures/FISNACHT_venn_four.pdf", width = 15, height = 9)

ggvenn(lists.wild, fill_color = c("lightblue", "khaki", "lightpink", "lightgrey"), text_size = 7, stroke_size = 2, set_name_size = 12)

dev.off()