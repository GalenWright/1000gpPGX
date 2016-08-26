#!/bin/bash
ROOT_DIR=..
DATA=$ROOT_DIR/data
PGX=$ROOT_DIR/pgx_analyses
ENS=$PGX/gene_ens_details
VCF=$PGX/vcf

mkdir -p $PGX/vcf

FTP=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502

for chr in 22 {1..22}
do
	grep -w ^${chr} $ENS/pgx1000GPexomeSlop25unique.bed > $VCF/splitChr${chr}Regions.bed
	tabix -fh $FTP/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
	-R $VCF/splitChr${chr}Regions.bed  > $VCF/pgx1000GPchr${chr}.vcf
	sleep 30s
done
rm $VCF/split*
rm ALL.chr*

# Concatenate all VCFs into a single file
vcf-concat $VCF/pgx1000GPchr{1..22}.vcf > $VCF/pgx1000GP.vcf

# Remove individual chromosome files
rm $VCF/pgx1000GPchr*vcf

# Remove AC=0 variants (i.e. those that were in the individuals excluded from the 1000GP final list)
vcftools --vcf $VCF/pgx1000GP.vcf --recode-INFO-all \
--out $VCF/pgx1000GPnoAC0 --recode --mac 1

# Remove exome structural variants as they cannot be analysed adequately
grep -v esv $VCF/pgx1000GPnoAC0.recode.vcf > $VCF/pgx1000GP.vcf 

rm $VCF/*log $VCF/pgx1000GPnoAC0.recode.vcf
