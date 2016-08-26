#!/bin/sh

ROOT_DIR=..
PGX=$ROOT_DIR/pgx_analyses
FST=$PGX/fst
VCF=$PGX/vcf
FRQ=$PGX/freq
in_VCF=$VCF/pgx1000GP.VEP.vcf
out_VCF=$FST/pgx1000GP.VEP.fst.syn.vcf
VEP=~/Programs/ensembl-tools-release-83/scripts/variant_effect_predictor

mkdir -p $PGX/fst

# Filter for only synonymous variants
perl $VEP/variant_effect_predictor.pl \
--port 3337 \
--input $in_VCF \
--output_file stdout \
--format vcf --vcf \
--cache --fork 4 --no_stats \
--species homo_sapiens \
--compress "gunzip -c" \
--most_severe --check_existing \
-force | 
perl $VEP/filter_vep.pl --input $in_VCF \
--force -o $out_VCF \
--ontology --filter "Consequence is synonymous_variant"

# Calculate pair-wise Fst for synonymous variants in populations

for i in ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI
	do
		for j in ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI
			do
				vcftools \
				--vcf $FST/pgx1000GP.VEP.fst.syn.vcf \
				--weir-fst-pop $FRQ/$i.samples.list \
				--weir-fst-pop $FRQ/$j.samples.list \
				--out $FST/${i}_${j}.fst.txt
		done
done

# Get mean Fst estimates for each of the comparisons
grep "Weir and Cockerham weighted Fst estimate" $FST/*log | \
awk -F"/" '{print $NF}' | \
cut  -d'_' -f1 > $FST/weighted.fst.estimates.population1.txt
grep "Weir and Cockerham weighted Fst estimate" $FST/*log | \
awk -F"/" '{print $NF}' | \
cut  -d'_' -f2 | cut -c1-3 > $FST/weighted.fst.estimates.population2.txt

grep "Weir and Cockerham weighted Fst estimate" $FST/*log | \
awk '{print $NF}' > $FST/weighted.fst.estimates.values.txt 
paste $FST/weighted.fst.estimates.population1.txt  $FST/weighted.fst.estimates.population2.txt \
$FST/weighted.fst.estimates.values.txt > $FST/weighted.fst.estimates.txt
rm $FST/weighted.fst.estimates.population*.txt $FST/weighted.fst.estimates.values.txt

# For continental super-populations
for i in AFR AMR EAS EUR SAS
        do
                for j in AFR AMR EAS EUR SAS 
                        do
                                vcftools \
                                --vcf $VCF/pgx1000GP.VEP.vcf \
                                --weir-fst-pop $FRQ/$i.samples.list \
                                --weir-fst-pop $FRQ/$j.samples.list \
                                --out $FST/${i}_${j}.super.fst.txt
                done
done

rm $FST/*log
