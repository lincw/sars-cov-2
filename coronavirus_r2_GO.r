################################################
# date: 29.09.2020
# output: failed! looks like too less query items

library(openxlsx)
library(topGO)
library(org.Hs.eg.db)

go_calculate <- function(x, GO) {
    if (missing(GO)) {GO = "BP"}
    myGOdata <- new("topGOdata",
                    description = "Human",
                    ontology = GO,
                    allGenes = x,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO
                    )
    resultFisher <- runTest(myGOdata,
                    algorithm = "weight01",
                    statistic = "fisher"
                    )
    go <- p.adjust(score(resultFisher), "fdr", n = length(score(resultFisher)))
    geneTable <- GenTable(myGOdata,
                 classicFisher = resultFisher,
                 topNodes = sum(score(resultFisher) < .01),
                 ranksOf = "classicFisher",
                 orderBy = "classicFisher",
                 numChar = 100
                 )
    go <- data.frame(fdr = go[names(go) %in% geneTable$GO.ID])
    go$GO.ID <- rownames(go)
    geneTable <- merge(geneTable, go, by = "GO.ID")
    geneTable$genes <- sapply(sapply(geneTable$GO.ID, simplify = T, function(x) {
            genes <- genesInTerm(myGOdata, x)
            genes[[1]][genes[[1]] %in% sigGenes(myGOdata)]
            }), paste, collapse = ",")
        return(geneTable[order(geneTable$fdr), ])
}

x <- org.Hs.egGO
mapped_genes <- mappedkeys(x)
horf <- read.table("~/Documents/INET-work/references/HuRI_binaryPPI/HORFeome_all.tsv", sep = "\t", header = T)
horf_clip <- horf[, c(1, 2, 5)]
orf_screen <- mapped_genes[mapped_genes %in% horf_clip$entrez_gene_id]
geneID2GO <- as.list(x[orf_screen])
geneUniverse <- names(geneID2GO)

inet <- read.xlsx("/Users/chung-wen.lin/Documents/INET-work/virus_network/Y2H_screening/20200923_verified/INET_All_Interactions_200923.xlsx", sheet = 1)
candidate <- names(table(inet$"Viral.Protein"))

for (i in 1:length(candidate)) {
    go_query <- inet[inet$"Viral.Protein" == candidate[i], "entrezID"]
    go_query <- factor(as.integer(geneUniverse %in% go_query))
    names(go_query) <- geneUniverse
    tryCatch({
        result <- go_calculate(go_query, GO = "BP")
        write.csv(result, file = paste0(candidate[i], "_GO_enrichment.csv"), row.names = FALSE)
    }, error = function(e) {
        cat("ERROR: ", conditionMessage(e), "\n")
    })
}
