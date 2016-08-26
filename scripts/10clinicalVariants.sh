#!/bin/bash
ROOT_DIR=..
PGX=$ROOT_DIR/pgx_analyses
DATA=$ROOT_DIR/data
VCF=$PGX/vcf
FRQ=$PGX/freq
CLIN=$PGX/clinical

mkdir -p $PGX/clinical
# The clinicalAnnotations.csv was downloaded from PharmGKB on 26 August 2014
# This was used to find rsIDs for genes in this study (saved $DATA/pgxClinicalLevel1.txt)

# Select rsIDs from file and use this list to search annotated vcf
cut -f1 $DATA/pgxClinicalLevel1.txt > $CLIN/pgxClinicalrsID.txt

vcftools --vcf $VCF/pgx1000GP.VEP.vcf --snps $CLIN/pgxClinicalrsID.txt --recode \
--recode-INFO-all --out $CLIN/pgx1000GP.VEP.clinical

while read POPULATION; do
          vcf-subset -c $FRQ/$POPULATION.samples.list $CLIN/pgx1000GP.VEP.clinical.recode.vcf \
	  > $CLIN/$POPULATION.clinical.vcf
          vcftools --freq2 --vcf $CLIN/$POPULATION.clinical.vcf --out $CLIN/$POPULATION.clinical.ind
done < $FRQ/phase3populations.txt
rm $CLIN/*log

# Create an index from vcf file that can be used to map back rsIDs as freq only outputs coordinates
awk '!/^($|#)/ {print $3"\t"$1":"$2}' $CLIN/ACB.clinical.vcf > $CLIN/rsIndexClin.txt

# Increase number of open files allowed (temp)
ulimit -n 4000

# Calculate the number of non-ref alleles per sample
vcftools --vcf $CLIN/pgx1000GP.VEP.clinical.recode.vcf --012 --out $CLIN/ClinPerSample

# Return number of open files allowed to original
ulimit -n 256
