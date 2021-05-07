# Evaluate the gene-disease association (GDA) for HuSCI dataset
# date: 03.05.2021
# author: Lin Chung-wen

## environment ----
library(gprofiler2)
library(disgenet2r)
library(KEGGREST)

#!!!!!!!!!!!!!!!!!!!!!!!!
# Not used not
## custom GMT ----
# gda <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/DisGeNet_curated_GDA.gmt") # generated 03.05.2021
# gda <- upload_GMT_file(gmtfile = "~/workplace/database/Homo_sapiens/disgenet.curated.v7.symbols.gmt")

## HuSCI ----
husci <- read.csv("~/Documents/INET-work/virus_network/Y2H_screening/20201104_final/binary_node_1126.csv", header = TRUE)
husci <- husci[husci$group == "human", c("Ensembl_ID", "node", "category")]
husci_convert <- gconvert(query = husci$Ensembl_ID, 
                          organism = "hsapiens",
                          target = "ENTREZGENE_ACC",
                          mthreshold = Inf, 
                          filter_na = TRUE)

sars2_gda <- disease_enrichment(entities = binary_node, 
                  vocabulary = "HGNC",
                  database = "HPO")

## KEGG ----
# id mapping (ref: https://www.researchgate.net/post/I-have-a-list-of-ENSEMBLE-gene-id-and-want-to-map-these-genes-to-pathways)
library(biomaRt)
# Retrieve the ensembl info
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = biomaRtOrg)
# Filter the ensembl ID
ids <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", 'entrezgene_id', 'entrezgene_accession'), values = , mart = ensembl)

# ref: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html
pathway_list <- keggList("pathway", "hsa")
pathway.codes <- sub("path:", "", names(pathway_list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                               pw <- keggGet(pwid)
                               if (is.null(pw[[1]]$GENE)) return(NA)
                               pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                               pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                               return(pw2)
                           }
)
head(genes.by.pathway)

geneList <- husci_convert$target
names(geneList) <- husci_convert$target

pVals.by.pathway <- t(sapply(names(genes.by.pathway),
    function(pathway) {
        pathway.genes <- genes.by.pathway[[pathway]]
        list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
        list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
        scores.in.pathway <- geneList[list.genes.in.pathway]
        scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
        if (length(scores.in.pathway) > 0){
            p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
        } else {
            p.value <- NA
        }
        return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
        }
    ))
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathway_list[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),]
head(outdat)
