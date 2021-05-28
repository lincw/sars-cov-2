library(reshape2)
library(ggplot2)
library(openxlsx)
library(msa)
# library(ggmsa)
# library(ape)
# library(ggtree)
############################################################
# not run now
# Get lower triangle of the correlation matrix
# get_lower_tri <- function(cormat){
#     cormat[upper.tri(cormat, diag = TRUE)] <- NA
#     return(cormat)
# }
# Get upper triangle of the correlation matrix
# get_upper_tri <- function(cormat){
#     cormat[lower.tri(cormat, diag = TRUE)] <- NA
#     return(cormat)
# }
############################################################

# label <- c(aa = "amino acid", nt = "nucleotide")
# for (now in 1:2) {
    sheetName <- getSheetNames(file = "~/Documents/INET-work/virus_network/statistic_results/align_score.xlsx")
    aa_list <- list()
    for (i in 1:length(sheetName)) {
        plot <- read.xlsx("~/Documents/INET-work/virus_network/statistic_results/align_score.xlsx", sheet = i, colNames = T, rowNames = F)
        plot_met <- melt(plot)
        plot_met$protein <- factor(plot_met$protein, levels = c("h229E", "HKU1", "NL63", "OC43", "MERS", "SARS2", "SARS"))
        plot_met$variable <- factor(plot_met$variable, levels = c("h229E", "HKU1", "NL63", "OC43", "MERS", "SARS2", "SARS"))
        aa_list[[i]] <-
        ggplot(plot_met, aes(x = protein, y = variable, fill = value)) +
            geom_tile(color = "white") +
            geom_text(aes(x = protein, y = variable, label = value), color = "black", size = 4) +
            scale_fill_gradient(low = "white", high = "blue", limit = c(0, 100), name = "percent identity") +
            labs(x = "coronavirus", y = "coronavirus", title = paste0(sheetName[i], " alignment")) +
            geom_label(aes(x = 6, y = 0.6, label = "amino acid (lower triangular)"), fill = "white") +
            geom_label(aes(x = 2, y = 7.4, label = "nucleotide (upper triangular)"), fill = "white") +
            theme_minimal()+
            coord_fixed()
        }
        pdf(file = "coronaviral_alignment.pdf")
        print(aa_list)
        dev.off()
# }

proteins <- data.frame(nsp1 = c(1, 22, 44, 64, 86, 110, 137),
                       nsp2	= c(2, 23, 45, 65, 87, 111, 138),
                       nsp3 = c(3, 24, 46, 66, 88, 112, 139),
                       nsp4 = c(4, 25, 47, 67, 89, 113, 140),
                       nsp5 = c(5, 26, 48, 68, 90, 114, 141),
                       nsp6 = c(6, 27, 49, 69, 91, 115, 142),
                       nsp7 = c(7, 28, 50, 70, 92, 116, 143),
                       nsp8 = c(8, 29, 51, 71, 93, 117, 144),
                       nsp9 = c(9, 30, 52, 72, 94, 118, 145),
                       nucleocapsid = c(21, 42, 63, 85, 108, 135, 162),
                       spike = c(16, 38, 59, 81, 101, 125, 152),
                       nsp10 = c(10, 31, 53, 73, 95, 119, 146),
                       nsp12 = c(11, 32, 54, 74, 96, 120, 147),
                       nsp13 = c(12, 33, 55, 75, 97, 121, 148),
                       nsp14 = c(13, 34, 56, 76, 98, 122, 149),
                       nsp15 = c(14, 35, 57, 77, 99, 123, 150),
                       nsp16 = c(15, 36, 58, 78, 100, 124, 151),
                       envelope = c(19, 40, 61, 83, 106, 128, 155),
                       membrane = c(20, 41, 62, 84, 107, 129, 156))

setwd("~/Documents/INET-work/virus_network/text_results/sequences/")
aa <- readAAStringSet("aa.fa")
nt <- readDNAStringSet("nt.fa")

## nsp2
align <- msa(aa[c("SARS-CoV-2-nsp2", "SARS-CoV-nsp2")], method = "ClustalW", order = "input")
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "nsp2_SARSs.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\feature{top}{1}{85..85}{fill:$\\downarrow$}{S.501Y.V2 T265I}"
        # "\\feature{top}{1}{501..501}{fill:$\\downarrow$}{S.501Y.V3 F681L}",
        # "\\feature{top}{1}{580..580}{fill:$\\downarrow$}{S.501Y.V3 I760T}"
    ))

## nsp3
align <- msa(aa[c("SARS-CoV-2-nsp3", "SARS-CoV-nsp3")], method = "ClustalW", order = "input")
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "nsp3_SARSs.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\feature{top}{1}{183..183}{fill:$\\downarrow$}{B.1.1.7 T1001I}",
        "\\feature{top}{1}{370..370}{fill:$\\downarrow$}{S.501Y.V3 S1188L}",
        "\\feature{top}{1}{837..837}{fill:$\\downarrow$}{S.501Y.V2 K1655N}",
        "\\feature{top}{1}{890..890}{fill:$\\downarrow$}{B.1.1.7 A1708D}",
        "\\feature{top}{1}{977..977}{fill:$\\downarrow$}{S.501Y.V3 K1795Q}",
        "\\feature{top}{1}{1412..1412}{fill:$\\downarrow$}{B.1.1.7 I2230T}"
    ))

## nsp5
align <- msa(aa[proteins[, "nsp5"]][c(6, 7, 5, 3, 4, 1, 2)], method = "ClustalW", order = "input")
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "nsp5_aa_full_align.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\frameblock{1}{90..90}{Red[1pt]}",
        "\\feature{top}{1}{90..90}{fill:$\\downarrow$}{S.501Y.V2 K3353R}"
    ))

## nsp6
align <- msa(aa[c("SARS-CoV-2-nsp6", "SARS-CoV-nsp6")], method = "ClustalW", order = "input")
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "nsp6_SARSs.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\frameblock{1}{SGF}{Red[1pt]}",
        "\\feature{top}{1}{106..108}{}{B.1.1.7 and S.501Y.V3 SGF 3675-3677 deletion}"
    ))

## nsp13
align <- msa(aa[c("SARS-CoV-2-nsp13", "SARS-CoV-nsp13")], method = "ClustalW", order = "input")
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "nsp13_SARSs.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\feature{top}{1}{341..341}{fill:$\\downarrow$}{S.501Y.V3 E5665D}"
    ))

## spike
# align <- msa(aa[proteins[, "spike"]][c(6, 7, 5, 3, 4, 1, 2)], method = "ClustalW", order = "input")
align <- msa(aa[c("SARS-CoV-2-Spike protein", "SARS-CoV-Spike protein")], method = "ClustalW", order = "input")
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "spike_aa_full_align.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\feature{top}{1}{18..18}{fill:$\\downarrow$}{S.501Y.V3 L18F}",
        "\\feature{bottom}{1}{20..20}{fill:$\\uparrow$}{S.501Y.V3 T20N}",

        # shift label position
        "\\feature{top}{1}{26..26}{fill:$\\downarrow$}{}",
        "\\feature{top}{1}{30..30}{}{S.501Y.V3 P26S}",

        "\\frameblock{1}{69..70}{Red[1pt]}",
        "\\feature{top}{1}{69..70}{}{B.1.1.7 HV69-70 deletion}",
        "\\feature{top}{1}{138..138}{fill:$\\downarrow$}{S.501Y.V3 D138Y}",
        "\\frameblock{1}{144..144}{Red[1pt]}",
        "\\feature{bottom}{1}{144..144}{}{B.1.1.7 Y144 deletion}",
        "\\feature{top}{1}{190..190}{fill:$\\downarrow$}{S.501Y.V3 R190S}",
        "\\frameblock{1}{417..417}{Red[1pt]}",
        "\\feature{top}{1}{417..417}{fill:$\\downarrow$}{S.501Y.V2 K417N}",
        "\\feature{bottom}{1}{484..484}{fill:$\\uparrow$}{S.501Y.V2 E484K}",
        "\\feature{top}{1}{501..501}{fill:$\\downarrow$}{B.1.1.7, S.501Y.V2 and S.501Y.V3 N501Y}",
        "\\feature{top}{1}{570..570}{fill:$\\downarrow$}{B.1.1.7 A570D}",
        "\\feature{top}{1}{614..614}{fill:$\\downarrow$}{B.1.1.7, S.501Y.V2 and S.501Y.V3 D614G}",
        "\\feature{top}{1}{655..655}{fill:$\\downarrow$}{S.501Y.V3 H655Y}",
        "\\feature{top}{1}{681..681}{fill:$\\downarrow$}{B.1.1.7 P681H}",
        "\\feature{top}{1}{701..701}{fill:$\\downarrow$}{S.501Y.V2 A701V}",
        "\\feature{top}{1}{761..761}{fill:$\\downarrow$}{B.1.1.7 T716I}",
        "\\feature{top}{1}{982..982}{fill:$\\downarrow$}{B.1.1.7 S982A}",
        "\\feature{top}{1}{1027..1027}{fill:$\\downarrow$}{S.501Y.V3 T1027I}",
        "\\feature{top}{1}{1118..1118}{fill:$\\downarrow$}{B.1.1.7 D1118H}"
    ))

## orf3a
align <- msa(aa[c(164, 165)], method = "ClustalW", order = "input")
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "orf3a_aa_full_align.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\frameblock{1}{57..57}{Red[1pt]}",
        "\\feature{top}{1}{57..57}{fill:$\\downarrow$}{S.501Y.V2 Q57H}",
        "\\frameblock{1}{171..171}{Red[1pt]}",
        "\\feature{top}{1}{171..171}{fill:$\\downarrow$}{S.501Y.V2 S171L}"
    ))

## orf8a
align <- msa(aa[c("SARS-CoV-2-ORF8", "SARS-CoV-ORF8a")], method = "ClustalW", order = "input")
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "orf8_SARSs.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\feature{top}{1}{92..92}{fill:$\\downarrow$}{S.501Y.V3 E92K}"
    ))

## nucleocapsid, all pan coronavirus
align <- msa(aa[proteins[, "nucleocapsid"]][c(6, 7, 5, 3, 4, 1, 2)], method = "ClustalW", order = "input")
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "nucleocapsid_aa_full_align.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\frameblock{1}{3..3}{Red[1pt]}",
        "\\feature{top}{1}{3..3}{fill:$\\downarrow$}{B.1.1.7 D3L}",
        "\\frameblock{1}{205..205}{Red[1pt]}",
        "\\feature{top}{1}{205..205}{fill:$\\downarrow$}{S.501Y.V2 T205I}{blue}",
        "\\frameblock{1}{235..235}{Red[1pt]}",
        "\\feature{top}{1}{235..235}{fill:$\\downarrow$}{B.1.1.7 S235F}"
    ))

## nucleocapsid, SARS1 and SARS2
align <- msa(aa[c("SARS-CoV-2-Nucleocapsid protein", "SARS-CoV-Nucleocapsid protein")])
msaPrettyPrint(align, shadingMode = "functional", shadingModeArg = "chemical", file = "SARSs_N.pdf", askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
    "\\feature{top}{1}{3..3}{fill:$\\downarrow$}{B.1.1.7 D3L}",
    "\\feature{top}{1}{80..80}{fill:$\\downarrow$}{B.501Y.V3 P80R}",
    "\\frameblock{1}{180..209}{Red[1pt]}",
    "\\feature{top}{1}{180..209}{}{SR-rich motif}",
    "\\feature{bottom}{1}{205..205}{fill:$\\uparrow$}{S.501Y.V2 T205I}",
    "\\feature{top}{1}{235..235}{fill:$\\downarrow$}{B.1.1.7 S235F}")
    )

for (i in names(proteins)) {
    tryCatch({
        align <- msa(nt[proteins[, i]])
#        msaPrettyPrint(align, output = "pdf", shadingMode = "similar", file = paste0(i, "_nt_fingerprint_align.pdf"), askForOverwrite = FALSE, code = c("\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\threshold[80]{50}"))
    }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    tryCatch({
        align <- msa(nt[proteins[, i]])
        msaPrettyPrint(align, output = "pdf", shadingMode = "similar", file = paste0(i, "_nt_full_align.pdf"), askForOverwrite = FALSE, code = c("\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\threshold[80]{50}", "\\showsequencelogo{top}"))
    }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n") })
}

#####
# for nsp3 alignment, which contains very long sequences. Latex cannot operate it now.

align <- msa(aa[proteins[, "nsp3"]][c(6, 7, 5, 3, 4, 1, 2)], method = "ClustalW", order = "input")

chunkSize <- 4000
for (start in seq(1, ncol(align), by = chunkSize)) {
    end <- min(start + chunkSize - 1, ncol(align))
    alnPart <- DNAMultipleAlignment(subseq(unmasked(align), start, end))
  msaPrettyPrint(alnPart, shadingMode = "functional", shadingModeArg = "chemical", file = paste0("nsp3", "_nt_full_align", start, "-", end, ".pdf") askForOverwrite = FALSE, code = c("\\shadingmode[chemical]{functional}", "\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\showcaption[top]{alignment with ClustalW, shading with chemical properties}",
        "\\frameblock{1}{SGF}{Red[1pt]}",
        "\\feature{top}{1}{106..108}{fill:$\\downarrow$}{B.1.1.7 SGF 3675-3677 deletion}"
    ))
    msaPrettyPrint(alnPart, shadingMode = "similar", showLegend = TRUE, file = paste0("nsp3", "_nt_full_align", start, "-", end, ".pdf"), askForOverwrite = FALSE, code = c("\\shownumbering{right}", "\\showruler{1}{top}", "\\showlegend", "\\threshold[80]{50}", "\\showsequencelogo{top}"))
   # msaPrettyPrint(x = alnPart, output = "pdf", subset = NULL, file = paste0("aln_", start, "-", end, ".pdf"))
}
