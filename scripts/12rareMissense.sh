#!/bin/bash
ROOT_DIR=..
PGX=$ROOT_DIR/pgx_analyses
VCF=$PGX/vcf
FRQ=$PGX/freq
POP=$PGX/popdiff
ENS=$PGX/gene_ens_details
CSQ=$PGX/consequences
GENE=$CSQ/genes
VEP=~/Programs/ensembl-tools-release-83/scripts/variant_effect_predictor/

# Analyse each of the genes in the study independently, appending counts
echo "rare_missense" > $GENE/temp.rare.missense.gene.txt
while read gene;
do
vcftools --vcf $GENE/$gene.vcf --max-maf 0.005 --recode-INFO-all --recode --stdout |\
grep missense_variant | wc -l >> $GENE/temp.rare.missense.gene.txt
done < $ENS/final_pgx_ensID
rm *log

paste $GENE/symbol.txt $GENE/temp.rare.missense.gene.txt > $GENE/rare.missense.gene.txt
rm $GENE/temp.rare.missense.gene.txt
