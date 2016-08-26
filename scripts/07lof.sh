#!/bin/sh

ROOT_DIR=..
VEP=~/Programs/ensembl-tools-release-83/scripts/variant_effect_predictor/
PGX=$ROOT_DIR/pgx_analyses
VCF=$PGX/vcf
ENS=$PGX/gene_ens_details
LOF=$PGX/lof
in_VCF=$VCF/pgx1000GP.VEP.vcf

mkdir -p $PGX/lof

# Annotate with LOFTEE
# NB only use gencode
# ~/Resources contains the ancestral human sequence - change to path of your sequence
# Filter variants for only high quality LoF
# Remove CSQ flags that do not match
# Only select PGx genes
perl $VEP/variant_effect_predictor.pl \
--port 3337 --gencode_basic \
--input $in_VCF \
--output_file stdout \
--format vcf --vcf \
--species homo_sapiens \
--cache --fork 4 \
--compress "gunzip -c" \
--per_gene --check_existing \
-force --symbol --regulatory \
--sift b --polyphen b \
--plugin LoF,human_ancestor_fa:~/Resources/human_ancestor.fa.rz |\
perl $VEP/filter_vep.pl --force_overwrite \
--filter "LoF is HC" --force_overwrite \
--filter "Gene in $ENS/final_pgx_ensID" \
-o $LOF/pgx1000GP.VEP.hc.lof.gencode.vcf --only_matched
rm *html *log

# Remove multi-allelic variants 
# Remove variants with a Lof_flag
# Curate variants with >30% allele freq: remove rs11356919 (AF 69%) due to lack of evidence 
vcftools  --vcf $LOF/pgx1000GP.VEP.hc.lof.gencode.vcf \
--max-alleles 2 --recode-INFO-all --recode --stdout |\
grep -v 'SINGLE_EXON\|NAGNAG_SITE\|PHYLOCSF_WEAK\|PHYLOCSF_UNLIKELY_ORF\|PHYLOCSF_TOO_SHORT\|rs11356919' \
> $LOF/pgx1000GP.VEP.hc.lof.gencode.clean.vcf

# CALCULATE WHICH GENES HARBOUR LoF MUTATIONS and HOW MANY
# Extract information from vcf info file
vcftools --vcf $LOF/pgx1000GP.VEP.hc.lof.gencode.clean.vcf \
--get-INFO AC --get-INFO EAS_AF --get-INFO AMR_AF --get-INFO AFR_AF \
--get-INFO EUR_AF --get-INFO SAS_AF --get-INFO CSQ --out $LOF/LoFConsequences

# Create a unique ID for each of the LoF variants
awk 'NR>1 {print $1":"$2":"$4}' $LOF/LoFConsequences.INFO > $LOF/LoFConsequencesUniqueID
# Get allele counts for each of the genes
awk 'NR>1 {print $5,$6,$7,$8,$9,$10}' $LOF/LoFConsequences.INFO > $LOF/LoFConsequencesAC
# Get the relevant gene names for each of the LoF variants
awk -F\| 'NR>1 {print $4}' $LOF/LoFConsequences.INFO > $LOF/LoFConsequencesGene
# Combine the Unique ID, AC/AF and Gene files
paste $LOF/LoFConsequencesUniqueID $LOF/LoFConsequencesAC $LOF/LoFConsequencesGene \
> $LOF/LoFConsequencesIDacGENE
rm $LOF/LoFConsequencesUniqueID $LOF/LoFConsequencesAC $LOF/LoFConsequencesGene

# CALCULATE NUMBER OF LoF / SAMPLE
# Increase number of open files allowed (temp)
ulimit -n 4000

# Calculate the number of non-ref alleles per sample
vcftools --vcf $LOF/pgx1000GP.VEP.hc.lof.gencode.clean.vcf --012 --out $LOF/LofPerSample

# Return number of open files allowed to original
ulimit -n 256
rm $LOF/*log
