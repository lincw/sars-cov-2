# tissue specificity assay of protein interactome ----
## date: 2021.04.08

library(openxlsx)
library(ggplot2)
library(reshape2)
library(ggsignif)
library(ggpubr)

## loading data ----
interactome_raw <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/7interactome_statistics.xlsx", sheet = 1) # mainly generated from "coronavirus_tissue_specificity.Rmd"
i_raw_df <- interactome_raw[c(1:8, 11), c(1, 12, 13)]
names(i_raw_df) <- c("interactome", "specific", "common")
i_raw_df$interactome <- factor(i_raw_df$interactome, levels = c("HuSCI", "Gordon", "Stukalov", "Li", "Nabeel", "Laurent", "St_Germain", "Samavarchi", "HPA"))
ird_melt <- melt(i_raw_df)
ird_melt$variable <- factor(ird_melt$variable, levels = c("common", "specific"))
interactome_p <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/7interactome_statistics.xlsx", sheet = 3)
total <- c(19670, 171, 384, 875, 285, 277, 2128, 1010, 2242)

## graph style ----
gp_style <- theme_bw() +
    theme(axis.title = element_text(color = "black", size = 14, face = "bold"),
        axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ird_melt$pattern <- ifelse(ird_melt$interactome == "HPA", "Yes", "No")
colour <- c("deeppink4", "turquoise4", "turquoise4", "turquoise4", "turquoise4", "cyan4", "cyan4", "cyan4", "grey50", "deeppink", "turquoise", "turquoise", "turquoise", "turquoise", "cyan", "cyan", "cyan", "grey70")

gp_percent <- ggplot(ird_melt, aes(x = interactome, y = value)) +
    geom_bar(stat = "identity", position = "fill", fill = colour) +
    scale_linetype_manual(values = c("common" = "dashed", "specific" = "blank")) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0%", "20%", "40%", "60%", "80%", "100%"), expand = expansion(mult = c(0, 0.1))) +
    labs(x = "", y = "fraction") +
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
    annotate("segment", x = 0.98343, xend = 8.01657, y = 1.33, yend = 1.33, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.3, yend = 1.33, color = "black") +
    annotate("segment", x = 8, xend = 8, y = 1.01, yend = 1.33, color = "black") +
    # p between HuSCI and HPA
    annotate("text", x = 5, y = 1.4, label = "p = 0.041", size = 2) +
    annotate("segment", x = 0.98343, xend = 9.01657, y = 1.38, yend = 1.38, color = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.35, yend = 1.38, color = "black") +
    annotate("segment", x = 9, xend = 9, y = 1.01, yend = 1.38, color = "black") +
    gp_style

ggsave(gp_percent, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/interactome_specificity_percent_rev.pdf", width = 5, height = 5)

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

## butterfly chart ----
i_raw_df$common <- i_raw_df$common * -1
ird_melt2 <- melt(i_raw_df)
bg <- ggplot(ird_melt2, aes(x = interactome, y = value / total[c(2:9, 1)])) +
    geom_bar(position = "identity", stat = "identity", fill = rep(c("pink", "lightblue"), each = 9)) +
    scale_y_continuous("", breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4), labels = c("80%", "60%", "40%", "20%", "0", "20%", "40%")) +
    labs(x = "PPI dataset", y = "% of Gene", color = "Specificity") +
    coord_flip() + theme_bw() +
    theme(axis.title = element_text(color = "black", size = 14, face = "bold"),
        axis.text = element_text(color = "black", size = 14))
ggsave(bg, file = "~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/7interactome_specificity_count_butterfly.pdf", width = 5, height = 3)
