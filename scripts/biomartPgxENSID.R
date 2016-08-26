library(biomaRt)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org",
                   path="/biomart/martservice",archive=FALSE, verbose=TRUE)
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

my_genes <- read.table("../data/all_pharmacogenes")
myGeneDetails <- getBM(attributes=c("chromosome_name", "hgnc_symbol", "ensembl_gene_id"), filters = "hgnc_symbol", values = my_genes , mart = ensembl)

# Remove non-numeric chromosome Ensembl gene IDs
autosomal_chr <- !is.na(as.numeric(as.character(myGeneDetails$chromosome_name)))
x_chr <- myGeneDetails$chromosome_name=="X"
myGeneDetails <- myGeneDetails[autosomal_chr|x_chr,]

write.table(myGeneDetails, "../data/all_pharmacogenes_details", quote=F, row.names = F) 
