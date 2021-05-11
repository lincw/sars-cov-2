# tissue specificity assay of protein interactome
## date: 2021.04.08

library(openxlsx)
library(ggplot2)
library(reshape2)
library(ggsignif)
library(ggpubr)

interactome_raw <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/7interactome_statistics.xlsx", sheet = 1) # mainly generated from "coronavirus_tissue_specificity.Rmd"
i_raw_df <- interactome_raw[c(1:8, 11), c(1, 12, 13)]
names(i_raw_df) <- c("interactome", "specific", "common")
ird_melt <- melt(i_raw_df)
#1 ird_melt$interactome <- factor(ird_melt$interactome, levels = c("HuSCI", "HuRI", "Gordon", "Stukalov", "Li", "Nabeel", "BioPlex3.0"))
ird_melt$interactome <- factor(ird_melt$interactome, levels = c("HuSCI", "Gordon", "Stukalov", "Li", "Nabeel", "Laurent", "St_Germain", "Samavarchi", "HPA"))
interactome_p <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/7interactome_statistics.xlsx", sheet = 3)

gp_style <- theme_bw() +
    theme(axis.title = element_text(color = "black", size = 14, face = "bold"),
        axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ird_melt$pattern <- ifelse(ird_melt$interactome == "HPA", "Yes", "No")

gp_percent <- ggplot(ird_melt, aes(x = interactome, y = value)) +
    geom_bar(aes(fill = variable), stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
    labs(x = "", y = "% of gene", fill = "Specificity") +
    # when manual add p onto the stacked barplot, the pvalue position is crashed!
    # stat_pvalue_manual(interactome_p, y.position = 1.05, step.increase = 0.01, label = "p") +
    # p between HuSCI and Gordon
    annotate("text", x = 1.5, y = 1.05, label = "p < 0.00001", size = 2) +
    annotate("segment", x = 0.98343, xend = 2.01657, y = 1.03, yend = 1.03, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.01, yend = 1.03, color = "black") +
    annotate("segment", x = 2, xend = 2, y = 1.01, yend = 1.03, color = "black") +
    # p between HuSCI and Stukalov
    annotate("text", x = 2, y = 1.1, label = "p < 0.0001", size = 2) +
    annotate("segment", x = 0.98343, xend = 3.01657, y = 1.08, yend = 1.08, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.05, yend = 1.08, color = "black") +
    annotate("segment", x = 3, xend = 3, y = 1.01, yend = 1.08, color = "black") +
    # p between HuSCI and Li
    annotate("text", x = 2.5, y = 1.15, label = "p < 0.00001", size = 2) +
    annotate("segment", x = 0.98343, xend = 4.01657, y = 1.13, yend = 1.13, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.1, yend = 1.13, color = "black") +
    annotate("segment", x = 4, xend = 4, y = 1.01, yend = 1.13, color = "black") +
    # p between HuSCI and Nabeel
    annotate("text", x = 3, y = 1.2, label = "p < 0.00001", size = 2) +
    annotate("segment", x = 0.98343, xend = 5.01657, y = 1.18, yend = 1.18, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.15, yend = 1.18, color = "black") +
    annotate("segment", x = 5, xend = 5, y = 1.01, yend = 1.18, color = "black") +
    # p between HuSCI and Laurent
    annotate("text", x = 3.5, y = 1.25, label = "p < 0.00001", size = 2) +
    annotate("segment", x = 0.98343, xend = 6.01657, y = 1.23, yend = 1.23, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.2, yend = 1.23, color = "black") +
    annotate("segment", x = 6, xend = 6, y = 1.01, yend = 1.23, color = "black") +
    # p between HuSCI and St_Germain
    annotate("text", x = 4, y = 1.3, label = "p < 0.00001", size = 2) +
    annotate("segment", x = 0.98343, xend =7.01657, y = 1.28, yend = 1.28, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.25, yend = 1.28, color = "black") +
    annotate("segment", x = 7, xend = 7, y = 1.01, yend = 1.28, color = "black") +
    # p between HuSCI and Samavarchi
    annotate("text", x = 4.5, y = 1.35, label = "p < 0.00001", size = 2) +
    annotate("segment", x = 0.98343, xend = 8.01657, y = 1.32, yend = 1.32, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.3, yend = 1.32, color = "black") +
    annotate("segment", x = 8, xend = 8, y = 1.01, yend = 1.32, color = "black") +
    # p between HuSCI and HPA
    annotate("text", x = 5, y = 1.4, label = "p = 0.041", size = 2) +
    annotate("segment", x = 0.98343, xend = 9.01657, y = 1.38, yend = 1.38, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.35, yend = 1.38, color = "black") +
    annotate("segment", x = 9, xend = 9, y = 1.01, yend = 1.38, color = "black") +
    gp_style

ggsave(gp_percent, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/interactome_specificity_percent.pdf", width = 5, height = 5)


gp_count <- ggplot(ird_melt, aes(x = interactome, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "# of gene count", fill = "Specificity") +
    gp_style
ggsave(gp_count, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/7interactome_specificity_count_withBioID.pdf", width = 5, height = 5)

## statistics ----
pvalue <- c()
for (i in 2:9) {
    fp <- fisher.test(matrix(as.numeric(c(i_raw_df[1, c(2, 3)], i_raw_df[i, c(2, 3)])), ncol = 2))
    pvalue <- c(pvalue, fp$p.value)
}
pvalue_adj <- p.adjust(pvalue, method = "fdr")
