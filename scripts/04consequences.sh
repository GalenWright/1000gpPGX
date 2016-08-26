#!/bin/bash
ROOT_DIR=..
DATA=$ROOT_DIR/data
PGX=$ROOT_DIR/pgx_analyses
ENS=$PGX/gene_ens_details
VCF=$PGX/vcf
CSQ=$PGX/consequences
GENE=$CSQ/genes
VEP=~/Programs/ensembl-tools-release-83/scripts/variant_effect_predictor/

# Create different vcfs for different frequency classes
# Singletons
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --recode --max-mac 1 \
--recode-INFO-all --out $CSQ/pgx1000GP.VEP.singletons

# Non-singletons to MAF 0.01
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --recode --mac 2 --max-maf 0.01 \
--recode-INFO-all --out $CSQ/pgx1000GP.VEP.nonsing-0.01MAF

# MAF 0.01-0.05
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --recode --maf 0.01 --max-maf 0.05 \
--recode-INFO-all --out $CSQ/pgx1000GP.VEP.0.01-0.05MAF

# MAF > 0.05
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --recode --maf 0.05 \
--recode-INFO-all --out $CSQ/pgx1000GP.VEP.greater.0.05MAF

# Count number of appearances of each of the Sequence Ontology (SO) terms
# $DATA/SO.terms was modified from:
# http://www.ensembl.org/info/genome/variation/predicted_data.html

while read SO; do
	grep $SO $VCF/pgx1000GP.VEP.vcf | wc -l
done < $DATA/SO.terms > $CSQ/allCount

while read SO; do
        grep $SO $CSQ/pgx1000GP.VEP.singletons.recode.vcf | wc -l 
done < $DATA/SO.terms > $CSQ/singletonsCount

while read SO; do
        grep $SO $CSQ/pgx1000GP.VEP.nonsing-0.01MAF.recode.vcf | wc -l 
done < $DATA/SO.terms > $CSQ/nonsing-0.01MAFCount

while read SO; do
        grep $SO $CSQ/pgx1000GP.VEP.0.01-0.05MAF.recode.vcf | wc -l 
done < $DATA/SO.terms > $CSQ/0.01-0.05MAFCount

while read SO; do
        grep $SO $CSQ/pgx1000GP.VEP.greater.0.05MAF.recode.vcf | wc -l 
done < $DATA/SO.terms > $CSQ/greater0.05MAFCount

echo "terms 	countAll	countSingletons	countnonSing-0.01MAF	count0.01-0.05MAF	countGreater0.05" > $CSQ/SO.terms.MAF.summary 
paste $DATA/SO.terms $CSQ/allCount $CSQ/singletonsCount  $CSQ/nonsing-0.01MAFCount \
$CSQ/0.01-0.05MAFCount $CSQ/greater0.05MAFCount >> $CSQ/SO.terms.MAF.summary

cp $VCF/pgx1000GP.VEP.vcf $CSQ
rm $CSQ/*log $CSQ/*Count

# Create UniqueID for each of the vcf frequency classes
grep -v ^# $CSQ/pgx1000GP.VEP.singletons.recode.vcf | awk '{print $1":"$2":"$5}' \
>  $CSQ/pgx1000GP.VEP.singletons.uniqueID

grep -v ^# $CSQ/pgx1000GP.VEP.nonsing-0.01MAF.recode.vcf | awk '{print $1":"$2":"$5}' \
>  $CSQ/pgx1000GP.VEP.nonsing-0.01MAF.uniqueID

grep -v ^# $CSQ/pgx1000GP.VEP.0.01-0.05MAF.recode.vcf | awk '{print $1":"$2":"$5}' \
>  $CSQ/pgx1000GP.VEP.0.01-0.05MAF.uniqueID

grep -v ^# $CSQ/pgx1000GP.VEP.greater.0.05MAF.recode.vcf | awk '{print $1":"$2":"$5}' \
>  $CSQ/pgx1000GP.VEP.greater.0.05.uniqueID

# Create individual files for each coding SO term for all genes
mkdir -p $CSQ/genes

while read SO; do echo $SO > $GENE/$SO.gene.txt ; done < $DATA/SO.coding.terms

# Add heading to VIP genes gene symbol
echo gene > $GENE/symbol.txt
awk 'NR>1 {print $2}' $DATA/final_all_pharmacogenes_details >> $GENE/symbol.txt

# Analyse each of the genes in the study independently, appending counts

while read gene;
do
perl $VEP/filter_vep.pl -i $VCF/pgx1000GP.VEP.vcf \
--filter "Gene is $gene" --force \
-o $GENE/$gene.vcf --only_matched
    while read term; do grep $term $GENE/$gene.vcf | wc -l >> $GENE/$term.gene.txt; done < $DATA/SO.coding.terms
done < $ENS/final_pgx_ensID

paste $GENE/symbol.txt $GENE/*gene* > $GENE/summaryVariantsPerGene
rm $GENE/*txt rm $CSQ/pgx1000GP.VEP.vcf
