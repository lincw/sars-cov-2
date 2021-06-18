phtmap <- pheatmap::pheatmap(t(plot_p10[comm_position, row_dim2]),
    cluster_rows = F,
    color = colorRampPalette(c("white", "red"))(50), border_color = "white",
    display_numbers = ifelse(t(plot_fdr[comm_position, row_dim2]) < 0.05, "*", ""),
    fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(0, 1, 2, 3, 4, 5, 5.85),
    legend_labels = c("0", "1", "2", "3", "4", "5", "-log10(p)"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05) \n hierarchical cluster based on -log10(p)")
col_dend <- phtmap[[2]]
# reorder columns
col_dend <- rotate(col_dend, order = c("3682", "2769", "731", "316", "4227", "1353", "4160", "2831", "1833", "2545", "692", "926", "2398", "525", "1882", "1652", "3941", "415", "377", "571"))
pdf("~/workplace/gwas_10.pdf")
pheatmap::pheatmap(t(plot_p10[comm_position, row_dim2]),
    cluster_cols = as.hclust(col_dend), cluster_rows = F,
    color = colorRampPalette(c("white", "red"))(50), border_color = "white",
    display_numbers = ifelse(t(plot_fdr[comm_position, row_dim2]) < 0.05, "*", ""),
    fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(0, 1, 2, 3, 4, 5, 5.85),
    legend_labels = c("0", "1", "2", "3", "4", "5", "-log10(p)"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05) \n hierarchical cluster based on -log10(p)")
pheatmap::pheatmap(t(plot_data[comm_position, row_dim2]),
    cluster_cols = as.hclust(col_dend), cluster_rows = F,
    color = colorRampPalette(c("blue", "white", "red"))(50), breaks = myBreaks, border_color = "white",
    display_numbers = ifelse(t(plot_fdr[comm_position, row_dim2]) < 0.05, "*", ""),
    fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 4.6),
    legend_labels = c("-3", "-2", "-1", "0", "1", "2", "3", "4", "BETA/SE\n"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05) \n hierarchical cluster based on -log10(p)")
dev.off()
# have dendgrogram structure from above heatmap
# have dendgrogram structure from above heatmap
row_dend <- phtmap[[1]]
# reorder rows
row_dend <- rotate(row_dend, order = c("NEUROT_UKB", "HRET", "RET", "PHF", "HC_UKBS", "SCZ_UKBS", "BMIA", "ADPN", "FAT_UKB", "HYPOTHY_UKBS", "IBD_UKBS", "OST_UKBS", "T2D_UKBS", "HEIGHT", "HIP"))
# replot
pheatmap::pheatmap(t(plot_p10[row_dim, col_dim]),
    cluster_cols = as.hclust(col_dend), cluster_rows = as.hclust(row_dend),
    color = colorRampPalette(c("white", "red"))(50), border_color = "white",
    display_numbers = ifelse(t(plot_fdr[row_dim, col_dim]) < 0.05, "*", ""),
    fontsize_number = 12, fontsize_col = 12, fontsize_row = 12, number_color = "black",
    legend_breaks = c(0, 1, 2, 3, 4, 5, 5.85),
    legend_labels = c("0", "1", "2", "3", "4", "5", "-log10(p)"),
    legend = TRUE,
    main = "Community GWAS hits (FDR < 0.05) \n hierarchical cluster based on -log10(p)")
