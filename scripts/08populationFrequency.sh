#!/bin/bash

ROOT_DIR=..
PGX=$ROOT_DIR/pgx_analyses
VCF=$PGX/vcf
FRQ=$PGX/freq

mkdir -p $PGX/freq

# Download 1000GP Phase 3 sample information files
curl -O ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
mv *panel $FRQ

# Make files with unique (super)population
awk 'NR>1 {print $2}' $FRQ/integrated_call_samples_v3.20130502.ALL.panel | sort | uniq > $FRQ/phase3populations.txt
awk 'NR>1 {print $3}' $FRQ/integrated_call_samples_v3.20130502.ALL.panel | sort | uniq > $FRQ/phase3superPopulations.txt

# Calculate bi-allelic SNPs
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --min-alleles 2 --max-alleles 2 \
--recode-INFO-all --recode --out $FRQ/pgx1000GP.VEP.biall

# Create a file with each sample in a given population and 
# Calculate allele frequencies of bi-allelic SNPS for individual populations
while read POPULATION; do
          grep $POPULATION $FRQ/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > $FRQ/$POPULATION.samples.list
          vcf-subset -c $FRQ/$POPULATION.samples.list $FRQ/pgx1000GP.VEP.biall.recode.vcf > $FRQ/$POPULATION.vcf
          vcftools --freq2 --vcf $FRQ/$POPULATION.vcf --out $FRQ/$POPULATION.Ind
done < $FRQ/phase3populations.txt

# Create a file with each sample in a given super-population
while read SUPERPOPULATION; do
          grep $SUPERPOPULATION $FRQ/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > $FRQ/$SUPERPOPULATION.samples.list
done < $FRQ/phase3superPopulations.txt

# Calculate the global allele frequency of bi-allelic SNPs
vcftools --freq2 --vcf $FRQ/pgx1000GP.VEP.biall.recode.vcf --out $FRQ/all

# Create file with allele count/number/freq INFO fields from vcf
vcftools --vcf $FRQ/pgx1000GP.VEP.biall.recode.vcf \
--get-INFO AC --get-INFO AN --get-INFO AF \
--out $FRQ/pgx1000GP.VEP.biall

# Create file with rsIDs for vcf (double hash saves header)
grep -v '##' $FRQ/pgx1000GP.VEP.biall.recode.vcf  | awk '{print $3}' \
> $FRQ/pgx1000GP.VEP.biall.rsID 

# Combine INFO and rsID files
paste $FRQ/pgx1000GP.VEP.biall.INFO $FRQ/pgx1000GP.VEP.biall.rsID \
> $FRQ/pgx1000GP.VEP.biall.INFO.rsID

rm $FRQ/*log

# List singleton markers
vcftools --vcf $FRQ/pgx1000GP.VEP.biall.recode.vcf --singletons --out $FRQ/pgx1000GP.VEP
# Number of singletons per sample
awk 'NR>1 {print $1}' $FRQ/integrated_call_samples_v3.20130502.ALL.panel > $FRQ/phase3samples.txt
echo "singletons" > $FRQ/singleton.count
while read SAMPLE; do grep $SAMPLE $FRQ/pgx1000GP.VEP.singletons | wc -l; done < $FRQ/phase3samples.txt \
>> $FRQ/singleton.count
paste $FRQ/integrated_call_samples_v3.20130502.ALL.panel $FRQ/singleton.count > $FRQ/phase3.singletons.per.sample 
rm $FRQ/singleton.count
