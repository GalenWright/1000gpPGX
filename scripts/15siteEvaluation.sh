#!/bin/bash
ROOT_DIR=..
DATA=$ROOT_DIR/data
PGX=$ROOT_DIR/pgx_analyses
CSQ=$PGX/consequences
ENS=$PGX/gene_ens_details
SVM=$PGX/site_eval
mkdir -p $PGX/site_eval

for chr in {1..22}
do
curl -O ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/site_assessment/ALL.chr${chr}.unfiltered_union_sites_with_svm.20130502.biallelic_snps.sites.vcf.gz
curl -O ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/site_assessment/ALL.chr${chr}.unfiltered_union_sites_with_svm.20130502.biallelic_snps.sites.vcf.gz.tbi
mv *biallelic_snps.sites.vcf.gz* $SVM/
grep -w ^${chr} $ENS/pgx1000GPexomeSlop25unique.bed > $SVM/splitChr${chr}Regions.bed
tabix -fh $SVM/ALL.chr${chr}.unfiltered_union_sites_with_svm.20130502.biallelic_snps.sites.vcf.gz \
-R $SVM/splitChr${chr}Regions.bed > $SVM/pgx1000GPchr${chr}.unfiltered.svm.vcf
rm $SVM/splitChr${chr}Regions.bed
done

# Combine vcf files - need to edit the header as some data are missing
# Make header
grep ^# $SVM/pgx1000GPchr1.unfiltered.svm.vcf > $SVM/pgx1000GP.unfiltered.svm.vcf
# Add missing header information
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" \
>> $SVM/pgx1000GP.unfiltered.svm.vcf
# Add variant information
cat $SVM/pgx1000GPchr*.unfiltered.svm.vcf | grep -v ^# | sort -n -k1 -k2 >> $SVM/pgx1000GP.unfiltered.svm.vcf

rm $SVM/ALL.chr*.unfiltered_union_sites_with_svm.20130502.biallelic_snps.sites.vcf.gz* \
$SVM/pgx1000GPchr*.unfiltered.svm.vcf

#Annotate with VEP
VEP=~/Programs/ensembl-tools-release-83/scripts/variant_effect_predictor/
in_VCF=$SVM/pgx1000GP.unfiltered.svm.vcf
out_VCF=$SVM/pgx1000GP.svm.VEP.vcf
# The GRCh37 v76 databases are hosted on a different port on ensembldb.ensembl.org (3337)
perl $VEP/variant_effect_predictor.pl \
--port 3337 \
--input $in_VCF \
--output_file stdout \
--format vcf --vcf \
--species homo_sapiens \
--cache --fork 4 \
--compress "gunzip -c" \
--per_gene --check_existing \
-force --symbol --regulatory \
--sift b --polyphen b \
--plugin LoF,human_ancestor_fa:~/Resources/human_ancestor.fa.rz \
--plugin ExAC,/Users/galenwright/Resources/ExAC/ExAC.r0.3.sites.vep.vcf.gz |\
perl $VEP/filter_vep.pl \
--filter "Gene in $ENS/final_pgx_ensID" --force \
-o $out_VCF --only_matched
rm *summary.html

vcftools --vcf $SVM/pgx1000GP.svm.VEP.vcf --get-INFO DP --get-INFO AC --get-INFO AF \
--get-INFO SVM  --get-INFO CSQ --out $SVM/pgx1000GP.svm.VEP
rm $SVM/*log

# Calculate proportion of different types of variants per gene
# All variants, pass, marginal and fail
echo "all_variants" > $SVM/all_variants
echo "high_qc_variants" > $SVM/high_qc_variants
echo "marginal_variants" > $SVM/marginal_variants
echo "fail_variants" > $SVM/fail_variants

while read gene
do
grep -w $gene $SVM/pgx1000GP.svm.VEP.INFO | wc -l >> $SVM/all_variants
grep -w $gene $SVM/pgx1000GP.svm.VEP.INFO | awk '$8>0.3 {print}' | wc -l >> $SVM/high_qc_variants
grep -w $gene $SVM/pgx1000GP.svm.VEP.INFO | awk '$8>0 {print}' | awk '$8<0.3 {print}' | wc -l >> $SVM/marginal_variants
grep -w $gene $SVM/pgx1000GP.svm.VEP.INFO | awk '$8<0 {print}' | wc -l >> $SVM/fail_variants
done < $ENS/final_pgx_ensID

paste $CSQ/genes/symbol.txt $SVM/all_variants $SVM/high_qc_variants $SVM/marginal_variants $SVM/fail_variants \
> $SVM/pgx1000GP.svm.analyses.txt 
rm $SVM/all_variants $SVM/high_qc_variants $SVM/marginal_variants $SVM/fail_variants

# Generate list of marginal variants that may be artefacts for the supplementary
echo "CHROM	POS	REF	ALT	total_depth	AC	AF	SVM	annotation	gene	ID" \
> $SVM/pgx1000GP.svm.marginal.variants.txt
awk '$8>0 {print}' $SVM/pgx1000GP.svm.VEP.INFO | awk '$8<0.3 {print}' | awk -F"|" '{print $1,$2,$4,$18}' |\
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12}' >> $SVM/pgx1000GP.svm.marginal.variants.txt
