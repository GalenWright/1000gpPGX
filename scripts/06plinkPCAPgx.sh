#!/bin/bash

ROOT_DIR=..
PGX=$ROOT_DIR/pgx_analyses
ENS=$PGX/gene_ens_details
VCF=$PGX/vcf
PCA=$PGX/pca

mkdir -p $PGX/pca

# Remove low frequency variants
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --maf 0.01 \
--recode --out $PCA/pgx1000GP.VEP.maf.0.01

# Convert to tped first
vcftools --vcf $PCA/pgx1000GP.VEP.maf.0.01.recode.vcf --plink-tped \
--out $PCA/pgx1000GP.VEP.maf.0.01.recode.plink

# Make bed files
plink --tped $PCA/pgx1000GP.VEP.maf.0.01.recode.plink.tped \
--tfam $PCA/pgx1000GP.VEP.maf.0.01.recode.plink.tfam --make-bed \
--out $PCA/pgx1000GP.VEP.maf.0.01.recode.plink
rm $PCA/*tfam $PCA/*tped

# Prune for LD
plink --bfile $PCA/pgx1000GP.VEP.maf.0.01.recode.plink --indep-pairwise 50 5 0.2 \
--out $PCA/pgx1000GP.VEP.maf.0.01.recode.plink.LD
plink --bfile $PCA/pgx1000GP.VEP.maf.0.01.recode.plink \
--extract $PCA/pgx1000GP.VEP.maf.0.01.recode.plink.LD.prune.in \
--make-bed --out $PCA/pgx1000GP.VEP.maf.0.01.recode.plink.LD.pruned

# Run pca with parfile, which is located where data is stored
cp pgx1000GP.VEP.parfile.txt $PCA
cd $PCA
smartpca -p pgx1000GP.VEP.parfile.txt > pgx1000GP.VEP.PCA.logfile

# Clean up first column of evec file which has format SampleName:SampleName
# So that it only has the sample name once
awk -F":" 'NR>1 {print $2}' pgx1000GP.VEP.pruned.evec > pgx1000GP.VEP.samples.pruned.evec
