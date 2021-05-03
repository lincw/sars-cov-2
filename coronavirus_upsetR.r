library(UpSetR)
library(grid)
library(gridExtra)
library(ggplot2)
setwd("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified")
data <- read.csv("SARS2_multiple_resource.csv", header = T)
data_intersect <- fromList(data)
row.names(data_intersect) <- unique(unlist(data))
data_intersect <- data_intersect[row.names(data_intersect) != "",]
write.table(data_intersect, file = "SARS2_multiple_resource_target-intersect.csv", sep = ",")

data <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/SARS2_multiple_resource_target-intersect.csv", header = T, row.names = 1)
names(data) <- c("INET (Y2H)", "Roth (fBFG)", "Gordon (AP-MS)", "Li (AP-MS)", "Stukalov (AP-MS)", "Laurent (BioID)")

##################################################################
# move horizontal bar to right
# ref: https://stackoverflow.com/questions/58131712/upset-plot-with-set-size-bars-in-right-side
# library(UpSetR)
# library(grid)
# library(gridExtra)
# library(ggplot2)

# NoAttBasePlot <- function (legend, size_plot_height, Main_bar_plot, Matrix_plot,
#     hratios, Size_plot, query_legend, set_metadata, set_metadata_plots,
#     newpage) {
#     top <- 1
#     bottom <- 100
#     if ((!is.null(legend)) && (query_legend != tolower("none"))) {
#         if (query_legend == tolower("top")) {
#             top <- 3
#             bottom <- 102
#             legend_top <- 1
#             legend_bottom <- 3
#             size_plot_height <- (size_plot_height + 2)
#         }
#         else if (query_legend == tolower("bottom")) {
#             legend_top <- 101
#             legend_bottom <- 103
#         }
#     }
#     if (is.null(set_metadata)) {
#         matrix_and_mainbar_right <- 100
#         matrix_and_mainbar_left <- 21
#         size_bar_right <- 20
#         size_bar_left <- 1
#     }
#     else if (!is.null(set_metadata)) {
#         matrix_and_mainbar_right <- set_metadata$ncols + 100
#         matrix_and_mainbar_left <- set_metadata$ncols + 21
#         size_bar_right <- set_metadata$ncols + 20
#         size_bar_left <- set_metadata$ncols + 1
#         metadata_right <- set_metadata$ncols
#         metadata_left <- 1
#     }
#     if (newpage) {
#         grid::grid.newpage()
#     }
#     if ((!is.null(legend)) && (query_legend != tolower("none"))) {
#         if (query_legend == tolower("top")) {
#             pushViewport(viewport(layout = grid.layout(102, matrix_and_mainbar_right)))
#         }
#         else if (query_legend == tolower("bottom")) {
#             pushViewport(viewport(layout = grid.layout(103, matrix_and_mainbar_right)))
#         }
#     }
#     else if ((is.null(legend)) || (query_legend == tolower("none"))) {
#         pushViewport(viewport(layout = grid.layout(100, matrix_and_mainbar_right)))
#     }
#     # Modified
#     vp = UpSetR:::vplayout(top:bottom, 1:(matrix_and_mainbar_right-matrix_and_mainbar_left))
#     pushViewport(vp)
#     grid.draw(arrangeGrob(Main_bar_plot, Matrix_plot, heights = hratios))
#     popViewport()
#     # Modified
#     vp = UpSetR:::vplayout(size_plot_height:bottom, (matrix_and_mainbar_right-matrix_and_mainbar_left-1):96)
#     pushViewport(vp)
#     grid.draw(arrangeGrob(Size_plot))
#     popViewport()
#     if (!is.null(set_metadata)) {
#         for (i in 1:length(set_metadata_plots)) {
#             if (i != 1) {
#                 metadata_left <- 1 + metadata_right
#                 metadata_right <- metadata_right + set_metadata$plots[[i]]$assign
#             }
#             else {
#                 metadata_left <- 1
#                 metadata_right <- set_metadata$plots[[i]]$assign
#             }
#             vp = UpSetR:::vplayout(size_plot_height:bottom, metadata_left:metadata_right)
#             pushViewport(vp)
#             grid.draw(arrangeGrob(set_metadata_plots[[i]]))
#             popViewport()
#         }
#     }
#     if ((!is.null(legend)) && (query_legend != tolower("none"))) {
#         vp = UpSetR:::vplayout(legend_top:legend_bottom, matrix_and_mainbar_left:matrix_and_mainbar_right)
#         pushViewport(vp)
#         grid.draw(arrangeGrob(legend))
#         popViewport()
#     }
# }

# Make_size_plot <- function (Set_size_data, sbar_color, ratios, ylabel, scale_sets,
#     text_scale, set_size_angle, set_size.show, set_size.scale_max,
#     set_size.number_size) {
#     if (length(text_scale) > 1 && length(text_scale) <= 6) {
#         x_axis_title_scale <- text_scale[3]
#         x_axis_tick_label_scale <- text_scale[4]
#     }
#     else {
#         x_axis_title_scale <- text_scale
#         x_axis_tick_label_scale <- text_scale
#     }
#     if (ylabel == "Set Size" && scale_sets != "identity") {
#         ylabel <- paste("Set Size", paste0("( ",
#             scale_sets, " )"))
#         if (scale_sets == "log2") {
#             Set_size_data$y <- log2(Set_size_data$y)
#         }
#         if (scale_sets == "log10") {
#             Set_size_data$y <- log10(Set_size_data$y)
#         }
#     }
#     if (!is.null(set_size.number_size)) {
#         num.size <- (set_size.number_size/2.845276) * x_axis_tick_label_scale
#     }
#     else {
#         num.size <- (7/2.845276) * x_axis_tick_label_scale
#     }
#     Size_plot <- (ggplot(data = Set_size_data, aes_string(x = "x",
#         y = "y")) + geom_bar(stat = "identity", colour = sbar_color,
#         width = 0.4, fill = sbar_color, position = "identity") +
#         scale_x_continuous(limits = c(0.5, (nrow(Set_size_data) +
#             0.5)), breaks = c(0, max(Set_size_data)), expand = c(0,
#             0)) + theme(panel.background = element_rect(fill = "white"),
#         plot.margin = unit(c(-0.11, -1.3, 0.5, 0.5), "lines"),
#         axis.title.x = element_text(size = 8.3 * x_axis_title_scale),
#         axis.text.x = element_text(size = 7 * x_axis_tick_label_scale,
#             vjust = 1, hjust = 0.5), axis.line = element_line(colour = "gray0"),
#         axis.line.y = element_blank(), axis.line.x = element_line(colour = "gray0",
#             size = 0.3), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#         panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
#         xlab(NULL) + ylab(ylabel) + coord_flip())
#     if (set_size.show == TRUE) {
#         Size_plot <- (Size_plot + geom_text(aes(label = y, vjust = 0.5,
#             hjust = 1.2, angle = set_size_angle), size = num.size))
#     }
#     if (scale_sets == "log10") {
#         if (!is.null(set_size.scale_max)) {
#             Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max,
#                 0), trans = log10_reverse_trans()))
#         }
#         else {
#             Size_plot <- (Size_plot + scale_y_continuous(trans = log10_reverse_trans()))
#         }
#     }
#     else if (scale_sets == "log2") {
#         if (!is.null(set_size.scale_max)) {
#             Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max,
#                 0), trans = log2_reverse_trans()))
#         }
#         else {
#             Size_plot <- (Size_plot + scale_y_continuous(trans = log2_reverse_trans()))
#         }
#     }
#     else {
#         if (!is.null(set_size.scale_max)) {
#             Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max,
#                 0), trans = "reverse"))
#         }
#         else {
#             # Modified
#             #Size_plot <- (Size_plot + scale_y_continuous(trans = "reverse"))
#         }
#     }
#     Size_plot <- ggplot_gtable(ggplot_build(Size_plot))
#     return(Size_plot)
# }

# assignInNamespace(x="NoAttBasePlot", value=NoAttBasePlot, ns="UpSetR")
# assignInNamespace(x="Make_size_plot", value=Make_size_plot, ns="UpSetR")
##################################################################

pdf("SARS2_multiple_dataset_UpSetR.pdf", width = 12, height = 5)
par(mar = c(5, 15, 5, 2))
upset(data, order.by = "freq", nsets = 6, text.scale = 1.5, sets.bar.color = c("orange", "red", "red", "red", "blue", "blue"), sets.x.label = "interactor candidate (n)",
    set_size.show = T
    )
dev.off()

# subset of intersection
pdf("SARS2_subset_UpSetR.pdf", width = 8, height = 5)
upset(data, sets = c("INET (Y2H)", "Roth (fBFG)"), text.scale = 1.5, sets.x.label = "interactor candidate (n)")
grid.text("Human interaction candiates distribution (in Y2H screening)",x = 0.6, y = 0.95
    # gp = gpar(fontsize = 12)
    )

upset(data, sets = c("Gordon (AP-MS)", "Li (AP-MS)", "Stukalov (AP-MS)"), text.scale = 1.5, sets.x.label = "interactor candidate (n)")
grid.text("Human interaction candiates distribution (in AP-MS screening)",x = 0.6, y = 0.95)

upset(data, sets = c("INET (Y2H)", "Roth (fBFG)", "Laurent (BioID)"), text.scale = 1.5, sets.x.label = "interactor candidate (n)")
grid.text("Human interaction candiates distribution (Y2H and BioID)",x = 0.6, y = 0.95)

upset(data, sets = c("Gordon (AP-MS)", "Li (AP-MS)", "Stukalov (AP-MS)", "Laurent (BioID)"), text.scale = 1.5, sets.x.label = "interactor candidate (n)")
grid.text("Human interaction candiates distribution (AP-MS and BioID)",x = 0.6, y = 0.95)

upset(data, sets = c("INET (Y2H)", "Roth (fBFG)", "Gordon (AP-MS)", "Li (AP-MS)", "Stukalov (AP-MS)"), text.scale = 1.5, sets.x.label = "interactor candidate (n)")
grid.text("Human interaction candiates distribution (Y2H and AP-MS)",x = 0.6, y = 0.95)
dev.off()

##################################################################
# distribution of edge
edge <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/SARS2_edge_multiple_resource.csv", header = T)
names(edge) <- c("INET (Y2H)", "Roth (fBFG)", "Gordon (AP-MS)", "Li (AP-MS)", "Stukalov (AP-MS)", "Laurent (BioID)")
edge_intersect <- fromList(edge)
row.names(edge_intersect) <- unique(unlist(edge))
edge_intersect <- edge_intersect[row.names(edge_intersect) != "",]
write.table(edge_intersect, file = "~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/SARS2_edge_multiple_resource_intersect.csv", sep = ",")

edge <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/SARS2_edge_multiple_resource_intersect.csv", header = T, row.names = 1)
names(edge) <- c("INET (Y2H)", "Roth (fBFG)", "Gordon (AP-MS)", "Li (AP-MS)", "Stukalov (AP-MS)", "Laurent (BioID)")

## output node distribution as pdf
pdf("SARS2_interactome_UpSetR.pdf", width = 14, height = 6)
upset(data, sets = c("Laurent (BioID)", "Stukalov (AP-MS)", "Li (AP-MS)", "Gordon (AP-MS)", "Roth (fBFG)", "INET (Y2H)"), order.by = "degree", nsets = 6, text.scale = 1.5, point.size = 3.5, line.size = 2, sets.bar.color = c("orange", "red", "red", "red", "blue", "blue"), sets.x.label = "Candidator", keep.order = TRUE
    # set_size.show = T
    )
grid.text("Human interaction candidator distribution",x = 0.6, y = 0.95, gp = gpar(fontsize = 14))
## output edge distribution as pdf
upset(edge, sets = c("Laurent (BioID)", "Stukalov (AP-MS)", "Li (AP-MS)", "Gordon (AP-MS)", "Roth (fBFG)", "INET (Y2H)"), order.by = "degree", nsets = 6, text.scale = 1.5, point.size = 3.5, line.size = 2, sets.bar.color = c("orange", "red", "red", "red", "blue", "blue"), sets.x.label = "Protein interactions", keep.order = TRUE
    # set_size.show = T
    )
grid.text("Virhost interactome edge distribution",x = 0.6, y = 0.95, gp = gpar(fontsize = 14))
dev.off()
