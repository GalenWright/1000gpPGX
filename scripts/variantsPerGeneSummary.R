library(biomaRt)
library(dplyr)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org",path="/biomart/martservice",archive=FALSE, verbose=TRUE)
#building a query, requires filters, attributes and values
#listFilters shows all filters
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

my_ensemblID <- read.table("../pgx_analyses/gene_ens_details/final_pgx_ensID")

maskGeneInfo <- read.table("../pgx_analyses/consequences/genes/maskGeneSummary.txt",header= TRUE)
variantsGeneInfo <- read.table("../pgx_analyses/consequences/genes/summaryVariantsPerGene",header= TRUE)

myAttributes <- getBM(attributes=c('ensembl_gene_id',"ensembl_transcript_id","cds_length"), filters = "ensembl_gene_id", values = my_ensemblID , mart = ensembl)
geneID <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = my_ensemblID , mart = ensembl)

geneSummaryMask <- merge(geneID, myAttributes)
names(geneSummaryMask)[2] <- "gene"

geneSummaryMask <- merge(geneSummaryMask, maskGeneInfo, by="gene")
geneSummaryMask <- merge(geneSummaryMask, variantsGeneInfo, by="gene")


# CDS length /3 - 1
# The CdS includes the stop-codon, if present. The reason for this is there are many annotated transcripts in the Ensembl database the length of whose transcribed exons are not divisible by 3. Hence we leave it to the user to decide how to deal with that, but mention here that determining the number of complete codons is trivial and you can slice the Cds so that it’s length is divisible by 3.
geneSummaryMask$protein_length <- geneSummaryMask$cds_length/3-1
geneSummaryMask <- na.omit(geneSummaryMask)
geneSummaryMask <- geneSummaryMask %>% group_by(gene) %>% filter(cds_length == max(cds_length))  
geneSummaryMask$missensePerCDS <- geneSummaryMask$missense_variant/geneSummaryMask$cds_length
geneSummaryMask$ensembl_transcript_id <- NULL
geneSummaryMask <- geneSummaryMask[!duplicated(geneSummaryMask),]
geneSummaryMask$inaccessOrSegmental <- geneSummaryMask$segmental_variants+geneSummaryMask$inaccessible_variants-geneSummaryMask$inaccess_segmental_variants
geneSummaryMask$percentageInaccess <- geneSummaryMask$inaccessible_variants/geneSummaryMask$total_variants
geneSummaryMask$percentageSegmental <- geneSummaryMask$segmental_variants/geneSummaryMask$total_variants
geneSummaryMask$percentageInaccessOrSegmental <- geneSummaryMask$inaccessOrSegmental/geneSummaryMask$total_variants
geneSummaryMask$variantsPerLength <- geneSummaryMask$total_variants/geneSummaryMask$geneBedLengthTotal
geneSummaryMask$missensePerCDSrank <- rank(-geneSummaryMask$missensePerCDS,ties.method= "first")
geneSummaryMask$variantsPerLengthRank <- rank(-geneSummaryMask$variantsPerLength,ties.method= "first")

write.csv(geneSummaryMask, "../pgx_analyses/consequences/genes/geneSummaryVariantsMask.csv", quote=F, row.names = F)

