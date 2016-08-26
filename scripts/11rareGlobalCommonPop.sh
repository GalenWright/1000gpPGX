#!/bin/bash
ROOT_DIR=..
PGX=$ROOT_DIR/pgx_analyses
VCF=$PGX/vcf
FRQ=$PGX/Freq
POP=$PGX/popdiff

mkdir -p $PGX/popdiff

# Make a vcf with all rare variants (0.05%) in the global dataset
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --recode --max-maf 0.005 \
--recode-INFO-all --out $POP/pgx1000GP.VEP.rare

# Get variant IDs for all rare variants
grep -v ^# $POP/pgx1000GP.VEP.rare.recode.vcf | awk '{print $3}' > $POP/pgx1000GP.VEP.rare.IDs

# Make a vcf for each population for variants that occur at 5% in that population
# And are rare in global dataset
while read POPULATION; do
	vcf-subset -c $FRQ/$POPULATION.samples.list $POP/pgx1000GP.VEP.rare.recode.vcf |\
        vcftools --vcf - --recode --maf 0.05 --snps $POP/pgx1000GP.VEP.rare.IDs \
        --recode-INFO-all --out $POP/$POPULATION.pop.diff
done < $FRQ/phase3populations.txt

# Extract SNPs that are differentiated
grep -v ^# $POP/*pop.diff.recode.vcf > $POP/pop.diff.snps.genotypes
grep -v ^# $POP/*pop.diff.recode.vcf | awk '{print $3}'  > $POP/pop.diff.snps.rsIDs
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --recode --snps $POP/pop.diff.snps.rsIDs \
--recode-INFO-all --out $POP/pop.diff.snps	

while read POPULATION; do
	vcf-subset -c $FRQ/$POPULATION.samples.list $POP/pop.diff.snps.recode.vcf |\
	vcftools --freq2 --vcf - --out $POP/$POPULATION.pop.diff.snps.Ind
done < $FRQ/phase3populations.txt

awk '!/^($|#)/ {print $3"\t"$1":"$2":"$5}' $POP/pop.diff.snps.recode.vcf > $POP/rsIndexPopDiff.txt

# Extract populations
awk -F"/" '{print $4}' $POP/pop.diff.snps.genotypes | cut -c1-3 > $POP/temp.population.pop.diff.snps.genotypes
awk '{print $3}' $POP/pop.diff.snps.genotypes > $POP/temp.rs.pop.diff.snps.genotypes
paste $POP/temp.population.pop.diff.snps.genotypes  $POP/temp.rs.pop.diff.snps.genotypes > $POP/population.pop.diff.snps.genotypes

#Extract gene name and variant type
vcftools --vcf $POP/pop.diff.snps.recode.vcf --get-INFO CSQ --out $POP/pop.diff.snps
awk -F"|" 'NR>1 {print $4"\t"$2}' $POP/pop.diff.snps.INFO > $POP/temp.genes.pop.diff.snps.genotypes
paste $POP/temp.genes.pop.diff.snps.genotypes $POP/rsIndexPopDiff.txt > $POP/genes.pop.diff.snps.genotypes

rm $POP/*log $POP/temp*
