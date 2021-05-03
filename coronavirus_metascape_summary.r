---
title: "Output of functional enrichment result from metascape.org"
output: pdf_document
always_allow_html: true
---

## 21.10.2020 Frist result
```{r}
# devtools::install_github("hrbrmstr/ggalt")
library(ggplot2)
library(stringr)

library(ggalt)

theme_set(theme_classic())
meta_plot <- function(x, main) {
    names(x) <- c("term", "log10p", "observeRatio")
    x$term <- factor(x$term, levels = x$term[rev(order(x$log10p))])
    ggplot(x, aes(x = term, y = -log10p)) +
        geom_point(aes(, size = observeRatio * 100), stat = "identity") +
        geom_segment(aes(y = 0, x = term, yend = -log10p, xend = term), color = "black") +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 60)) +
        labs(x = "", y = "-Log10(p)", size = "Observed ratio (%)", title = main) +
        theme(axis.text = element_text(color = "black", size = 14, hjust = 0,),
            axis.title = element_text(size = 16, face = "bold")) +
        coord_flip()
}
meta_plot_dumbbel <- function(x, main) {
    names(x) <- c("term", "log10p", "Enrichment", "Background")
    x$term <- factor(x$term, levels = x$term[order(x$Enrichment)])
    ggplot(x, aes(x = Background * 100, xend = Enrichment * 100, y = term, group = term)) +
        geom_dumbbell(color = "#e3e2e1", colour_x = "red", colour_xend = "blue", size = 2, dot_guide = TRUE, dot_guide_size = 0.25) +
        labs(y = "Functional category", x = "Enrichment (%)", title = main) +
        theme_bw() +
        theme(axis.text = element_text(color = "black", size = 14, hjust = 0,),
            axis.title = element_text(size = 16, face = "bold"))
}
meta_inet_nonspecific <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/INET_nonspecific/metascape_summary.csv", header = T)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/inet_non_specific_meta.pdf", width = 10, height = 8)
meta_plot(meta_inet_nonspecific, "INET common targets")
dev.off()

meta_inet_specific <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/INET_specific/metascape_summary.csv", header = T)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/inet_specific_meta.pdf", width = 10, height = 8)
meta_plot(meta_inet_specific, "INET unique targets")
dev.off()

meta_y2h_nonspecific <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/Y2H_binary_common/metascape_summary.csv", header = T)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/Y2Hbinary_common_meta.pdf", width = 10, height = 8)
meta_plot(meta_y2h_nonspecific, "Y2H binary common targets")
dev.off()

meta_y2h_specific <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/Y2H_binary_unique/metascape_summary.csv", header = T)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/Y2Hbinary_unique_meta.pdf", width = 10, height = 8)
meta_plot(meta_y2h_specific, "Y2H binary unique targets")
dev.off()
```

## 25.10.2020 Second result
According to Fritz and Pascal, it makes less sense to separate unique and common proteins for functional enrichment assay.
***
The DisGeNET is used as another functional enrichment term for query.
> the DisGeNET terms with duplicated geneID were removed.

```{r}
meta_inet_124 <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/INET_all124_GO_reactome_KEGG/metascape_summary.csv", header = T)
disgenet_inet_124 <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/INET_all124_GO_reactome_KEGG/Enrichment_QC/GO_DisGeNET_summary.csv", header = T)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/INET_all124_GO_reactome_KEGG_meta.pdf", width = 10, height = 9)
meta_plot(meta_inet_124, "Viral targeted human interactors (INET, n = 124)")
meta_plot_dumbbel(meta_inet_124, "INET")
meta_plot(disgenet_inet_124, "Gene-disease enrichment (INET, n = 124)")
dev.off()

meta_roth_229 <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/Roth_all228_GO_reactome_KEGG/metascape_summary.csv", header = T)
disgenet_roth_229 <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/Roth_all228_GO_reactome_KEGG/Enrichment_QC/GO_DisGeNET_summary.csv", header = T)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/Roth_all228_GO_reactome_KEGG_meta.pdf", width = 10, height = 8)
meta_plot(meta_roth_229, "Viral targeted human interactors (Roth, n = 229)")
meta_plot_dumbbel(meta_roth_229, "Roth")
meta_plot(disgenet_roth_229, "Gene-disease enrichment (Roth, n = 229)")
dev.off()

meta_y2h_346 <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/2_Y2H_all346_GO_reactome_KEGG/metascape_summary.csv", header = T)
disgenet_y2h_346 <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/meta_enrichment/2_Y2H_all346_GO_reactome_KEGG/Enrichment_QC/GO_DisGeNET_summary.csv", header = T)
pdf("~/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/2_Y2H_all346_GO_reactome_KEGG_meta.pdf", width = 10, height = 9)
meta_plot(meta_y2h_346, "Viral targeted human interactors (binary, n = 346)")
meta_plot_dumbbel(meta_y2h_346, "Binary")
meta_plot(disgenet_y2h_346, "Gene-disease enrichment (binary, n = 346)")
dev.off()
```
