#/bin/bash
ROOT_DIR=..
DATA=$ROOT_DIR/data
PGX=$ROOT_DIR/pgx_analyses
ENS=$PGX/gene_ens_details
CONS=$PGX/consequences
FRQ=$PGX/freq
EXOME=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/
COV=$PGX/coverage

# $1 is the sample name being analysed

echo $1 > $COV/ind.$1.coverage.txt

while read gene
do
        grep -w ${gene} $ENS/pgx1000GPexomeSlop25.bed  | sort -k1,1 -k2,2n | bedtools merge -i stdin > $COV/${gene}.${1}.bed
        interval=$(cat $COV/${gene}.${1}.bed | awk -F'\t' 'BEGIN{SUM=0}{SUM+=$3-$2 }END{print SUM}')
        location=$(grep -w ${gene} $ENS/pgx1000GPexomeSlop25.bed | awk ' BEGIN {min=1E10} $2 < min \
        {min=$2; min_row=$2} max < $3 {max =$3; max_row=$3} END {print $1":"min_row"-"max_row;} ')
        idx=$(grep -w ${1} $COV/integrated_call_samples_v3.20130502.ALL.IDs.index)
        echo "Processing" ${1} ${pop} "and" ${gene} "(location:" $location "size" $interval" bp)"
        samtools view -b $EXOME/$idx ${location} > $COV/${gene}.${1}.bam
        # Test if file is bam file is valid with picard if it isn't, repeat download
        picard ValidateSamFile MODE=SUMMARY IGNORE=MATE_NOT_FOUND I=$COV/${gene}.${1}.bam | grep "No errors found"  > $COV/${gene}.${1}.validate
        [ -s $COV/${gene}.${1}.validate ] || (samtools view -b $EXOME/$idx ${location} > $COV/${gene}.${1}.bam)
        # Perform a second test to see if download worked (i.e. if samtools depth has no data) and, if not, repeat
        samtools depth -b $COV/${gene}.${1}.bed $COV/${gene}.${1}.bam > $COV/${gene}.${1}.depth
        [ -s $COV/${gene}.${1}.depth ] || (samtools view -b $EXOME/$idx ${location} > $COV/${gene}.${1}.bam)
        # Final test
        picard ValidateSamFile MODE=SUMMARY IGNORE=MATE_NOT_FOUND I=$COV/${gene}.${1}.bam | grep "No errors found"  > $COV/${gene}.${1}.validate
        [ -s $COV/${gene}.${1}.validate ] || (sleep 5s && samtools view -b $EXOME/$idx ${location} > $COV/${gene}.${1}.bam)
        samtools depth -b $COV/${gene}.${1}.bed $COV/${gene}.${1}.bam | awk -v awk_int=$interval '{sum+=$3} END {print sum/awk_int"X"}' >> $COV/ind.${1}.coverage.txt
        rm $COV/${gene}.${1}.bed $COV/${gene}.${1}.bam $COV/${gene}.${1}.validate $COV/${gene}.${1}.depth
done < $COV/gene_symbol
