#!/bin/bash

ROOT_DIR=..
DATA=$ROOT_DIR/data
PGX=$ROOT_DIR/pgx_analyses
ENS=$PGX/gene_ens_details
VCF=$PGX/vcf
CSQ=$PGX/consequences
GENE=$CSQ/genes
MASK=$PGX/mask
VEP=~/Programs/ensembl-tools-release-83/scripts/variant_effect_predictor/

mkdir -p $PGX/mask
# Download 1000GP mask file
curl -O ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed
mv  *bed $MASK

# Make bed file with no "chr" from the 1000GP mask file
sed s/^chr//g $MASK/20141020.strict_mask.whole_genome.bed > $MASK/bed
mv $MASK/bed $MASK/20141020.strict_mask.whole_genome.bed

# Make vcf of masked/unmasked regions
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --bed $MASK/20141020.strict_mask.whole_genome.bed \
--out $MASK/pgx1000GP.VEP.accessible --recode --recode-INFO-all
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --exclude-bed $MASK/20141020.strict_mask.whole_genome.bed \
--out $MASK/pgx1000GP.VEP.inaccessible --recode --recode-INFO-all

# Look at segmental duplications
# Convert the following file to a bed
curl -O http://humanparalogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab
mv *tab $MASK
awk 'NR>1 {print $1"\t"$2"\t"$3}' $MASK/GRCh37GenomicSuperDup.tab | sed s/^chr//g > $MASK/GRCh37GenomicSuperDup.bed
rm $MASK/GRCh37GenomicSuperDup.tab 

# Make two vcf files, one within segmental dups, another outside of them
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --bed $MASK/GRCh37GenomicSuperDup.bed \
--out $MASK/pgx1000GP.VEP.segmental --recode --recode-INFO-all
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --exclude-bed $MASK/GRCh37GenomicSuperDup.bed \
--out $MASK/pgx1000GP.VEP.nonsegmental --recode --recode-INFO-all

# Make a vcf of SNPs found in both segmental duplications and inaccessible regions
vcftools --vcf $MASK/pgx1000GP.VEP.segmental.recode.vcf --exclude-bed $MASK/20141020.strict_mask.whole_genome.bed \
--out $MASK/pgx1000GP.VEP.segmental.inaccessible --recode --recode-INFO-all

# Calculate the number of variants per gene
echo total_variants > $GENE/totalVariantsPerGene.txt

while read gene;
do
grep -v ^# $GENE/$gene.vcf | wc -l >> $GENE/totalVariantsPerGene.txt
done < $ENS/final_pgx_ensID

# Look at number of variants in inaccessible regions
# Calculate the number of variants per gene
echo inaccessible_variants > $GENE/inaccessiblePerGene.txt

while read geneIn;
do
grep -v ^# $MASK/pgx1000GP.VEP.inaccessible.recode.vcf | grep -w $geneIn | wc -l >> $GENE/inaccessiblePerGene.txt
done < $ENS/final_pgx_ensID

# Look at number of variants in segmental duplications
# Calculate the number of variants per gene
echo segmental_variants > $GENE/segmentalPerGene.txt

while read geneSeg;
do
grep -v ^# $MASK/pgx1000GP.VEP.segmental.recode.vcf | grep -w $geneSeg | wc -l >> $GENE/segmentalPerGene.txt
done < $ENS/final_pgx_ensID

# Look at number of variants in inaccessible and segmental duplications
# Calculate the number of variants per gene
echo inaccess_segmental_variants > $GENE/inaccessSegmentalPerGene.txt

while read geneBoth;
do
grep -v ^# $MASK/pgx1000GP.VEP.segmental.inaccessible.recode.vcf | grep -w $geneBoth | wc -l \
>> $GENE/inaccessSegmentalPerGene.txt
done < $ENS/final_pgx_ensID

# Add heading to VIP genes gene symbol
echo gene > $GENE/symbol.txt
awk 'NR>1 {print $2}' $DATA/final_all_pharmacogenes_details >> $GENE/symbol.txt

# Calculate total length of each gene that was extracted
echo geneBedLengthTotal > $GENE/geneBedLengthTotal.txt
while read geneID
do
grep -w $geneID $ENS/pgx1000GPexomeSlop25.bed | \
sort -k1,1 -k2,2n | bedtools merge -i stdin | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' \
>> $GENE/geneBedLengthTotal.txt
done < $ENS/final_pgx_ensID

# Calculate transcript length of each gene that was extracted
echo geneBedLengthTranscript > $GENE/geneBedLengthTranscript.txt
while read geneID
do
grep -w $geneID $ENS/pgx1000GPexome.bed | \
sort -k1,1 -k2,2n | bedtools merge -i stdin | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' \
>> $GENE/geneBedLengthTranscript.txt
done < $ENS/final_pgx_ensID

# Combine all data together
paste $GENE/symbol.txt $GENE/totalVariantsPerGene.txt $GENE/inaccessiblePerGene.txt $GENE/segmentalPerGene.txt $GENE/inaccessSegmentalPerGene.txt $GENE/geneBedLengthTranscript.txt $GENE/geneBedLengthTotal.txt > $GENE/maskGeneSummary.txt

# To get summaries per gene, the following Rscript was run:
Rscript variantsPerGeneSummary.R
