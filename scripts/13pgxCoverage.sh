#/bin/bash
ROOT_DIR=..
DATA=$ROOT_DIR/data
PGX=$ROOT_DIR/pgx_analyses
ENS=$PGX/gene_ens_details
CONS=$PGX/consequences
FRQ=$PGX/freq
EXOME=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/
COV=$PGX/coverage
# NB this script requires access to a multi-core server and takes several days to complete
mkdir -p $PGX/coverage

# Select sample IDs for checking exome coverage
awk 'NR>1 {print $1}' $FRQ/integrated_call_samples_v3.20130502.ALL.panel \
> $COV/integrated_call_samples_v3.20130502.ALL.IDs

# Download sample index for exome alignments
#curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.exome.alignment.index
#cut -f1 20130502.phase3.exome.alignment.index | grep mapped | grep -v unmapped \
#| grep -wf $COV/integrated_call_samples_v3.20130502.ALL.IDs \
#> $COV/integrated_call_samples_v3.20130502.ALL.IDs.index 
#rm 20130502.phase3.exome.alignment.index

# Generate a gene symbol file to loop through
# Remove header
awk 'NR>1 {print $2}' $DATA/final_all_pharmacogenes_details > $COV/gene_symbol

# Execute the coverage script for x samples in parallel
cat $COV/integrated_call_samples_v3.20130502.ALL.IDs | parallel -j 48 ./perSampleCoverage.sh {1}

