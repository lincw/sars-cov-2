# checking human TFs in HuSCI

library(rethinking)
library(openxlsx)
source("~/Documents/Programming/R_functions/random_sig.r")

huri <- read.csv("~/Documents/INET-work/references/HuRI_binaryPPI/HuRI_Tong_withSymbol.csv", header = T)
huri_uniq <- unique(c(huri$symbol_a, huri$symbol_b))

horf <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/ORF_searchSpace.xlsx", sheet = 1)
horf_sym <- unique(horf$ensembl_geneName)
horf <- unique(horf$orf_id)

fbfg <- read.xlsx("~/Documents/INET-work/references/HuRI_binaryPPI/ORF_searchSpace.xlsx", sheet = 2)
fbfg_sym <- unique(fbfg$ensembl_geneName)
fbfg <- unique(fbfg$ORF_id)
tot_space <- unique(c(horf, fbfg))

total_sym <- unique(c(horf_sym, fbfg_sym))

husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)
husci <- husci[husci$group == "human", "node"]
tf <- read.table("~/Documents/INET-work/virus_network/references/human_TFome.txt")
husci_tf <- as.numeric(table(husci %in% tf$V1)[2])

random_tf <- simula(10000, huri_uniq, husci, tf$V1)
random_tf_sym <- simula(10000, total_sym, husci, tf$V1)

sign_tf <- round(1 - (abs(sig(random_tf, husci_tf)) / 10000), 4)
sign_tf_sym <- round(1 - (abs(sig(random_tf_sym, husci_tf)) / 10000), 4)

dens <- hist(random_tf, breaks = 30, plot = FALSE)
dens_sym <- hist(random_tf_sym, breaks = 30, plot = FALSE)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/figures/random_TF.pdf", width = 3, height = 3)
par(mgp = c(2, 0.7, 0), ps = 8)
plot(dens_sym, xlim = c(0, 35), ylim = c(0, 0.12), col = rgb(0.75, 0.75, 0.75, 1/2), freq = FALSE, border = NA, las = 1, xlab = "Number of identified TFs", ylab = "Frequency density", main = NULL, xaxt = "n")
plot(dens, xlim = c(0, 35), ylim = c(0, 0.12), col = rgb(1, 0.65, 0, 1/2), freq = FALSE, border = NA, las = 1, xaxt = "n", add = T)
axis(side = 1, at = seq(1, 35, by = 5) - 0.5, labels = c(seq(1, 35, by = 5)))
box(col = "black")
arrows(husci_tf - 0.5, 200/10000, husci_tf - 0.5, 10/10000, col = "red", lwd = 2, length = 0.1)
text(11, 350/10000, label = paste0("p = ", sign_tf_sym, "(integrate search space)\np = ", sign_tf, " (HuRI)"), cex = 0.6, pos = 4)
legend(17, 0.13, inset = 0.02, legend = c("Integrated search space", "HuRI"), fill = c(rgb(0.75, 0.75, 0.75, 1/2), rgb(1, 0.65, 0, 1/2)), cex = 0.5, bty = "n")
dev.off()
