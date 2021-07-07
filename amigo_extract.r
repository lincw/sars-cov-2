# Parse GO associated annotated genes
# Lin Chung-wen
# 07.07.2021

library(rvest)

go <- list(ubiquitination = "GO:0016567", # protein ubiquitination
		   cytoskeleton = "GO:0005856", # cytoskeleton
		   trafficking = "GO:0016192", # vesicle mediate trafficking
		   viral = "GO:0016032", # viral process
		   immune = "GO:0006955" # immune response
		)

url5 <- "http://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity_label,bioentity,bioentity_name,qualifier,annotation_class,assigned_by,taxon,isa_partof_closure,evidence_type,evidence_with,reference,date&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22annotation%22&fq=regulates_closure:%22"

url6 <- "%22&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&facet.field=aspect&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=evidence_subset_closure_label&facet.field=regulates_closure_label&facet.field=annotation_class_label&facet.field=qualifier&facet.field=annotation_extension_class_closure_label&facet.field=assigned_by&facet.field=panther_family_label&q=*:*"

# download GO annotated from AmiGO
# split html into text
go_query <- list()
for (i in 1:length(go)) {
	go_txt <- html_text(read_html(paste0(url5, go[[i]], url6)))
	go_query[[names(go[i])]] <- lapply(strsplit(go_txt, "\n"), function(x) strsplit(x, "\t"))
}
# convert list into data frame format for saving
go_summary <- list()
df1 <- lapply(go_query, function(x) data.frame(x[[1]]))
df2 <- lapply(df1, function(x) t(x))
write.xlsx(df2, file = "~/Documents/INET-work/virus_network/raw/GO_gene_annotation_AmiGO_raw.xlsx")
# extract unique gene list from above downloaded data
go_query_gene <- list()
for (i in 1:length(go)) {
	gene_id <- unique(unlist(lapply(go_query[[i]][[1]], function(x) 
		c(x[[1]], x[[2]]))))
	go_query_gene[[names(go[i])]] <- data.frame(matrix(gene_id, ncol = 2, byrow = T))
	names(go_query_gene[[names(go[i])]]) <- c("GeneName", "UniProt_ID")
}

write.xlsx(go_query_gene, file = "~/Documents/INET-work/virus_network/raw/go_gene_annotation.xlsx", overwrite = T)

### NOT accessable, different data length
# go_df_last <- data.frame(matrix(unlist(go_df), ncol = 15, byrow = T))
# names(go_df_last) <- c("source","bioentity_internal_id","bioentity_label","qualifier","annotation_class","reference","evidence_type","evidence_with","aspect","bioentity_name","synonym","type","taxon","date","assigned_by","annotation_extension_class")

# from QuickGO API
library(httr)
library(jsonlite)
library(xml2)
library(openxlsx)

goa_human <- read.delim("~/workplace/database/Homo_sapiens/GO/goa_human.gaf", header = F, skip = 41)
names(goa_human) <- c("db", "uniprot_id", "symbol", "qualifier", "GO_ID", "reference", "evidence_code", "WithorFrom", "aspect", "name", "synonym", "type", "taxon", "date", "assigned", "annotation_extension", "product_form_ID")

go_gene <- list()

for (i in 1:length(go)) {
	requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", go[[i]], "/children")
	r <- GET(requestURL, accept("application/json"))

	stop_for_status(r)

	json <- toJSON(content(r))
	df <- as.data.frame(unlist(fromJSON(json)))
	children_df <- df[grep("children.id*", row.names(df)), ]
	go_gene[[names(go[i])]] <- unique(goa_human[goa_human$"GO_ID" %in% children_df, c("uniprot_id", "symbol")])
}

write.xlsx(go_gene, file = "~/Documents/INET-work/virus_network/raw/go_gene_annotation.xlsx", overwrite = T)