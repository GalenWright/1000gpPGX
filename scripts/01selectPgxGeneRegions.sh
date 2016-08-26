#!/bin/bash
ROOT_DIR=..
DATA=$ROOT_DIR/data

mkdir -p $ROOT_DIR/pgx_analyses
PGX=$ROOT_DIR/pgx_analyses

# The following pharmacogenes were selected
# 1) VIP pharmaocogenes from PharmGKB (https://www.pharmgkb.org/downloads/)
# These are found in $DATA/vip_genes

# 2) Pharmacogenes that have SNPs with PharmGKB clinical evidence level 1 and 2
# These are found in $DATA/pharmgkb_1_2_clinical

# 3) Pharmacogenes mentioned the following recent reviews: PMID 24684508 PMID 25870137
# These are found in $DATA/emerging_pgx

# Create a unique list of pharmacogenes that meet any of these criteria
cat $DATA/vip_genes $DATA/pharmgkb_1_2_clinical $DATA/emerging_pgx |\
sort | uniq > $DATA/all_pharmacogenes

# Find Ensembl IDs and chromosomes for pharmacogenes
Rscript biomartPgxENSID.R

# Remove the HLA and UGT genes due to complexities of the loci
# Exclude IFNL4 as it is a pseudogene
# Remove and X chromosome genes
grep -v 'IFNL4\|UGT\|HLA\|^X' $DATA/all_pharmacogenes_details > $DATA/final_all_pharmacogenes_details

mkdir -p $PGX/gene_ens_details
ENS=$PGX/gene_ens_details

# Select Ensembl ID for extracting information from biomaRt
awk 'NR>1 {print $3}' $DATA/final_all_pharmacogenes_details > $ENS/final_pgx_ensID

# Use Ensembl archive site 
# http://grch37.ensembl.org/index.html

# Run R script to extract relevant information from BioMart using the relevant R package
Rscript biomartPgxBed.R

# Make *bed file
awk 'NR>1 {print $2"\t"$3"\t"$4"\t"$1","$5","$6}' $ENS/pgxMartExport.txt  > $ENS/pgxMartExport.bed

# Download the 1000GP exome target regions
curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets/20130108.exome.targets.bed
mv 20130108.exome.targets.bed $ENS

# Remove chr from 1000GP bed file 
sed 's/^chr//g' $ENS/20130108.exome.targets.bed > $ENS/nochr
mv $ENS/nochr $ENS/20130108.exome.targets.bed

# Intersect Ensembl exons with the 1000GP exome capture region, only keeping the regions common to 1000GP
bedtools intersect -a $ENS/pgxMartExport.bed -b $ENS/20130108.exome.targets.bed \
> $ENS/pgx1000GPexome.bed

# Add 25bp flanking sequences
awk -F"\t" '{print $1"\t"$2-25"\t"$3+25"\t"$4}' $ENS/pgx1000GPexome.bed > $ENS/pgx1000GPexomeSlop25.bed

# Make a unique region bed file and
# Calculate initial bp of sequence extracted
sort -k1,1 -k2,2n $ENS/pgx1000GPexomeSlop25.bed | bedtools merge -i stdin \
> $ENS/pgx1000GPexomeSlop25unique.bed

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM,"bp"}' $ENS/pgx1000GPexomeSlop25unique.bed \
> $ENS/pgx1000GPexomeSlop25uniqueSequence
