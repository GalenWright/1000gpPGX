#!/bin/sh
ROOT_DIR=..
PGX=$ROOT_DIR/pgx_analyses
VCF=$PGX/vcf
ENS=$PGX/gene_ens_details
CSQ=$PGX/consequences
VEP=~/Programs/ensembl-tools-release-83/scripts/variant_effect_predictor/
in_VCF=$VCF/pgx1000GP.vcf
out_VCF=$VCF/pgx1000GP.unfiltered.vcf
# The GRCh37 v76 databases are hosted on a different port on ensembldb.ensembl.org (3337)

perl $VEP/variant_effect_predictor.pl \
--port 3337 \
--input $in_VCF \
--output_file $out_VCF \
--format vcf --vcf \
--species homo_sapiens \
--cache --fork 4 \
--compress "gunzip -c" \
--per_gene --check_existing \
-force --symbol --regulatory \
--sift b --polyphen b \
--plugin LoF,human_ancestor_fa:~/Resources/human_ancestor.fa.rz \
--plugin ExAC,/Users/galenwright/Resources/ExAC/ExAC.r0.3.sites.vep.vcf.gz

# Only keep VIP genes (i.e. exclude non-PGx genes that overlap)
perl $VEP/filter_vep.pl -i $out_VCF \
--filter "Gene in $ENS/final_pgx_ensID" --force \
-o $VCF/pgx1000GP.VEP.vcf --only_matched

# Extract Ensembl VEP consequence field from vcf
mkdir -p $PGX/consequences
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --get-INFO AC \
--get-INFO AF --get-INFO CSQ --out $CSQ/pgx1000GP.VEP.CSQ

# Make a information file for annotating variants from R analyses
# Make header for file
echo "CHROM POS ID REF ALT AC AF FUNCTION GENE" > $CSQ/pgx1000GP.VEP.CSQ.annotation.INFO
# Select position, ID, and allele information information
grep -v ^# $VCF/pgx1000GP.VEP.vcf | awk '{print $1,$2,$3,$4,$5}' \
> $CSQ/pgx1000GP.VEP.CSQ.position.INFO
# Select allele count and allele frequency
awk 'NR>1 {print $5,$6}' $CSQ/pgx1000GP.VEP.CSQ.INFO > $CSQ/pgx1000GP.VEP.CSQ.allele.INFO
# Select gene name and functional consequence from INFO file
# Use underscore to deal with empty rsID columns
awk -F'|' 'NR>1 {print $2,$4}' $CSQ/pgx1000GP.VEP.CSQ.INFO > $CSQ/pgx1000GP.VEP.CSQ.geneFunction.INFO
# Combine the three files
paste $CSQ/pgx1000GP.VEP.CSQ.position.INFO $CSQ/pgx1000GP.VEP.CSQ.allele.INFO \
$CSQ/pgx1000GP.VEP.CSQ.geneFunction.INFO >> $CSQ/pgx1000GP.VEP.CSQ.annotation.INFO

rm $CSQ/pgx1000GP.VEP.CSQ.log $CSQ/pgx1000GP.VEP.CSQ.position.INFO \
$CSQ/pgx1000GP.VEP.CSQ.allele.INFO $CSQ/pgx1000GP.VEP.CSQ.geneFunction.INFO \
$VCF/pgx1000GP.unfiltered.vcf
