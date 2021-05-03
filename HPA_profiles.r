# shown the profiles of transcriptome data used in the Human Protein Atlas

library(dplyr)
library(ggplot2)

output <- function(file, name) {
    consens <- read.delim(file, header = T, sep = "\t")
    names(consens) <- c("Gene", "Gene_name", "Tissue", "TPM", "pTPM", "NX")
    df <- consens %>% group_by(Tissue) %>% count(Tissue)
    write.csv(df, file = paste0("~/workplace/database/Homo_sapiens/ProteinAtlas/statistics/", name, "_profile.csv"), row.names = F)
    gp_TPM <-
        ggplot(consens, aes(x = Tissue, y = log(TPM + 1, 2))) +
        geom_boxplot() +
        geom_hline(yintercept = 1, color = "red") +
        geom_text(aes(1, 1.5, label = 1), color = "red") +
        labs(x = "tissues", y = "log2(TPM + 1)", title = name) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    ggsave(gp_TPM, file = paste0("~/workplace/database/Homo_sapiens/ProteinAtlas/figures/", name, "_allBoxplot_log2TPM.pdf"))
    gp_NX <-
        ggplot(consens, aes(x = Tissue, y = log(NX + 1, 2))) +
        geom_boxplot() +
        geom_hline(yintercept = 1, color= "red") +
        geom_text(aes(1, 1.5, label = 1), color = "red") +
        labs(x = "tissues", y = "log2(NX + 1)", title = name) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    ggsave(gp_NX, file = paste0("~/workplace/database/Homo_sapiens/ProteinAtlas/figures/", name, "_allBoxplot_log2NX.pdf"))

    tumo <- consens[consens$Gene_name %in% tumorigenesis, ]
    gp_tumo <-
        ggplot(tumo, aes(x = Tissue, y = NX)) +
        geom_bar(aes(fill = Gene_name), stat = "identity", position = "dodge") +
        facet_wrap(~Gene_name, nrow = 1) +
        coord_flip() +
        labs(x = "tissues", y = "NX", title = name) +
        theme_bw() +
        theme(axis.text = element_text(color = "black", size = 10), legend.position = "none")
    ggsave(gp_tumo, file = paste0("~/workplace/database/Homo_sapiens/ProteinAtlas/figures/", name, "_tumorigenesis_NX.pdf"), width = 12, height = 8)
}

rna1 <- file.path("~/workplace/database/Homo_sapiens/ProteinAtlas/rna_consensus.tsv")
consens <- read.delim(rna1, header = T, sep = "\t")
name <- "consensus"
gp_NX <-
    ggplot(consens, aes(x = Tissue, y = log(NX + 1, 2))) +
    geom_boxplot() +
    geom_hline(yintercept = 1, color = "red") +
    geom_text(aes(1, 1.5, label = 1), color = "red") +
    labs(x = "tissues", y = "log2(NX + 1)", title = name) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(gp_NX, file = paste0("~/workplace/database/Homo_sapiens/ProteinAtlas/figures/", name, "_allBoxplot_log2NX.pdf"))

tumorigenesis <- c("BRCA1", "HTR1A", "PCSK9", "TP53", "ACE2")
tumo <- consens[consens$"Gene.name" %in% tumorigenesis, ]
gp_tumo <-
    ggplot(tumo, aes(x = Tissue, y = NX)) +
    geom_bar(aes(fill = Gene.name), stat = "identity", position = "dodge") +
    facet_wrap(~ Gene.name, nrow = 1) +
    coord_flip() +
    labs(x = "tissues", y = "NX", title = name) +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 10), legend.position = "none")
ggsave(gp_tumo, file = paste0("~/workplace/database/Homo_sapiens/ProteinAtlas/figures/", name, "_tumorigenesis_NX.pdf"), width = 12, height = 8)

rna2 <- file.path("~/workplace/database/Homo_sapiens/ProteinAtlas/rna_tissue_fantom.tsv")
output(rna2, "fantom")

rna3 <- file.path("~/workplace/database/Homo_sapiens/ProteinAtlas/rna_tissue_gtex.tsv")
output(rna3, "gtex")

rna4 <- file.path("~/workplace/database/Homo_sapiens/ProteinAtlas/rna_tissue_hpa.tsv")
output(rna4, "HPA")
