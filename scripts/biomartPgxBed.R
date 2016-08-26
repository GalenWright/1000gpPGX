library(biomaRt)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org",path="/biomart/martservice",archive=FALSE, verbose=TRUE)
#building a query, requires filters, attributes and values

filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

my_ensemblID <- read.table("../pgx_analyses/gene_ens_details/final_pgx_ensID")

myBedDetails <- getBM(attributes=c("chromosome_name", "exon_chrom_start", "exon_chrom_end","ensembl_gene_id", "ensembl_transcript_id"), filters = "ensembl_gene_id", values = my_ensemblID , mart = ensembl)
geneID <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = my_ensemblID , mart = ensembl)
myBedDetails <- merge (myBedDetails, geneID)
write.table(myBedDetails, "../pgx_analyses/gene_ens_details/pgxMartExport.txt", quote=F, row.names = F)
