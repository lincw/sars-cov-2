# 1. cluster with BETA/SE score
pdf("/tmp/1.pdf", width = 10)
h1 <- Heatmap(t(plot_data), column_title = "FDR < 0.1", name = "BETA/SE", col = col_fun, column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"))
# 2. collapse everything between 1 and -1 as 0
h2 <- Heatmap(t(plot_data_2), column_title = paste0("FDR < 0.1 \n  collapse ", expression(beta), "/SE score between 1 and -1 to 0"), name = "BETA/SE", col = col_fun, column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal"))
ht_list <- h1 + h2
draw(ht_list,
    ht_gap = unit(1, "cm"),
    row_dend_side = "right",
    heatmap_legend_side = "bottom",
    # auto_adjust = FALSE
    )
# 1.1 & 2.1 flip 90 degree of the table
# NOT necessary

# h3 <- Heatmap(plot_data, column_title = "FDR < 0.1", name = "BETA/SE", col = col_fun, column_names_gp = gpar(fontsize = 8),
#     heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"))
# h4 <- Heatmap(plot_data_2, column_title = "cluster based on BETA/SE \n collapse score 1 to -1 as 0", name = "BETA/SE", col = col_fun, column_names_gp = gpar(fontsize = 8),
#     heatmap_legend_param = list(direction = "horizontal"))
# ht_list <- h3 + h4
# draw(ht_list,
#     ht_gap = unit(1, "cm"),
#     row_dend_side = "right",
#     heatmap_legend_side = "bottom",
#     # auto_adjust = FALSE
#     )

# 1.3 & 2.3 from 1, with cell gaps
h5 <- Heatmap(t(plot_data), column_title = "FDR < 0.1", name = "BETA/SE", col = col_fun, column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"),
    rect_gp = gpar(col = "white", lwd = 1))
h6 <- Heatmap(t(plot_data_2), column_title = "FDR < 0.1 \n collapse score 1 to -1 as 0", name = "BETA/SE", col = col_fun, column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal"),
    rect_gp = gpar(col = "white", lwd = 1))
ht_list <- h5 + h6
draw(ht_list,
    ht_gap = unit(1, "cm"),
    row_dend_side = "right",
    heatmap_legend_side = "bottom",
    # auto_adjust = FALSE
    )
# 3. clustering based on p-value
h7 <- Heatmap(t(plot_p), column_title = "clustering based on pvalue", name = "p", col = col_fun_p, column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal"),
    rect_gp = gpar(col = "white", lwd = 1))
h8 <- Heatmap(t(plot_p10), column_title = "clustering based on pvalue", name = "-log10(p)", col = col_fun_p10, column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal"),
    rect_gp = gpar(col = "white", lwd = 1))
draw(h7 + h8,
    ht_gap = unit(1, "cm"),
    row_dend_side = "right",
    heatmap_legend_side = "bottom"
    )
# 4. filter by FDR < 0.05
h9 <- Heatmap(t(plot_data[row_dim, col_dim]), column_title = "clustering based on pvalue", name = "BETA/SE", col = col_fun, column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal"),
    rect_gp = gpar(col = "white", lwd = 1))
h10 <- Heatmap(t(plot_p[row_dim, col_dim]), column_title = "clustering based on pvalue", name = "p", col = col_fun_p, column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal"),
    rect_gp = gpar(col = "white", lwd = 1))
draw(h9,
    ht_gap = unit(1, "cm"),
    row_dend_side = "right",
    heatmap_legend_side = "bottom"
    )
draw(h10,
    ht_gap = unit(1, "cm"),
    row_dend_side = "right",
    heatmap_legend_side = "bottom"
    )
dev.off()

# 4.1 using pheatmap
pdf("/tmp/2.0.pdf", width = 10)
pheatmap::pheatmap(t(plot_data[row_dim, col_dim]),
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, border_color = "white",
        display_numbers = ifelse(t(plot_fdr[row_dim, col_dim]) < 0.05, "*", ""),
        fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05)")
pheatmap::pheatmap(t(plot_p10[row_dim, col_dim]),
    color = colorRampPalette(c("white", "red"))(50), border_color = "white",
        display_numbers = ifelse(t(plot_fdr[row_dim, col_dim]) < 0.05, "*", ""),
        fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(0, 1, 2, 3, 4, 5, 5.85),
    legend_labels = c("0", "1", "2", "3", "4", "5", "-log10(p)"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05)")
dev.off()

######
# for candidates
pdf("~/Documents/INET-work/virus_network/figure_results/GWAS/cluster_heatmap.pdf", width = 8)
pheatmap::pheatmap(t(plot_data),
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, border_color = "white",
        display_numbers = ifelse(t(plot_fdr) < 0.05, "*", ""),
        fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.1) \n hierarchical cluster based on BETA/SE ratio")
pheatmap::pheatmap(t(plot_data_2),
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, border_color = "white",
        display_numbers = ifelse(t(plot_fdr) < 0.05, "*", ""),
        fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.1) \n hierarchical cluster based on BETA/SE ratio (collapsed score <= |1| as 0)")
pheatmap::pheatmap(t(plot_p10[row_dim, col_dim]),
    color = colorRampPalette(c("white", "red"))(50), border_color = "white",
        display_numbers = ifelse(t(plot_fdr[row_dim, col_dim]) < 0.05, "*", ""),
        fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(0, 1, 2, 3, 4, 5, 5.85),
    legend_labels = c("0", "1", "2", "3", "4", "5", "-log10(p)"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05) \n hierarchical cluster based on -log10(p)")
dev.off()
