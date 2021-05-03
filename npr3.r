library(ComplexHeatmap)
library(circlize)

source("~/Documents/INET-work/scripts/2methods.r")
npr3a <- read.csv("~/Documents/INET-work/edge-score/data/npr3.a.list", header = F)
npr3b <- read.csv("~/Documents/INET-work/edge-score/data/npr3.b.list", header = F)
npr.iba <- as.vector(npr3a$V1)
npr.ibb <- as.vector(npr3b$V1)
npr.datalist <- log2edgeCal(merged.jasaaba.log2, npr.iba, npr.ibb)
npr.small_data <- do.call(rbind, npr.datalist)
npr.small_mean <- rowMeans(as.matrix(npr.small_data))
npr.small_std <- apply(npr.small_data, 1, sd)
npr.small_zscore <- (npr.small_data - npr.small_mean) / npr.small_std
npr.small_zscore <- npr.small_zscore[!grepl("*NA*", row.names(npr.small_zscore)), ]
colnames(npr.small_zscore) <- k_category_name$Description[match(colnames(npr.small_zscore), k_category_name$Abbreviation)]
row.names(npr.small_zscore) <- c("AT1G02450, NIMIN1", "AT3G25882, NIMIN2", "AT1G09415, NIMIN3", "AT5G65210, TGA1", "AT5G06950, TGA2", "AT1G22070, TGA3", "AT1G64280, NPR1", "AT5G45110, NPR3", "AT1G15750, TPL", "AT3G27890, quinone reductase")
pdf("/tmp/npr3_stress.pdf", width = 25, height = 12)
npr.map <- Heatmap(npr.small_zscore, cluster_rows = F, cluster_columns = F, column_names_max_height = max_text_width(colnames(npr.small_zscore), gp = gpar(fontsize = 4)), height = unit(8, "cm"),
    col = c("blue", "white", "red"),
    # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    width = unit(50, "cm"),column_title = "NPR3 interactions", row_title = "edge-score", heatmap_legend_param = list(title = "edge-score"),
  # cell_fun = function(j, i, x, y, width, height, fill) {
    # ifelse(npr.small_zscore[i, j] >= 2, grid.text(sprintf("%.1f", npr.small_zscore[i, j]), x, y, gp = gpar(fontsize = 7)), "")
  # }
    )
draw(npr.map, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()
