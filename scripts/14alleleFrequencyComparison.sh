#!/bin/bash
ROOT_DIR=..
DATA=$ROOT_DIR/data
PGX=$ROOT_DIR/pgx_analyses
VCF=$PGX/vcf
CSQ=$PGX/consequences
FRQ=$PGX/freq
GENE=$CSQ/genes
COMP=$PGX/freq_comp

mkdir -p $PGX/freq_comp
# Extract the global minor allele frequencies from vcf for 1000GP and ExAc
# Position and 1000GP allele frequencies 
awk 'NR>1 {print $1":"$2}' $CSQ/pgx1000GP.VEP.CSQ.INFO > $COMP/freq.comp.pos.af
# ExAc
awk -F"|" 'NR>1 {print $36}' $CSQ/pgx1000GP.VEP.CSQ.INFO > $COMP/freq.comp.exac.af
# 1000GP
awk 'NR>1 {print $6}' $CSQ/pgx1000GP.VEP.CSQ.INFO > $COMP/freq.comp.1000gp.af 
# Make header
echo "variant	freq_1000GP	freq_exac" > $COMP/freq.comp.pos.1000gp.exac.af
# Add information for variants in both 1000GP and ExAC, remove multi-allelic
paste $COMP/freq.comp.pos.af $COMP/freq.comp.exac.af $COMP/freq.comp.1000gp.af | awk 'NF==3 {print}' | grep -v "," >> $COMP/freq.comp.pos.1000gp.exac.af
rm $COMP/freq.comp.pos.af $COMP/freq.comp.exac.af $COMP/freq.comp.1000gp.af

# Download HGDP QC'd data from the LASER: Locating Ancestry from SEquence Reads download page
curl -O http://csg.sph.umich.edu/chaolong/LASER/HGDP-938-632958.tar.gz
tar -xzvf HGDP-938-632958.tar.gz
mv HGDP $COMP/
mv $COMP/HGDP/* $COMP
rmdir $COMP/HGDP/
# Create top header for vcf
echo "##fileformat=VCFv4.1" > $COMP/HGDP_938.vcf

# Make column headers for vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT" > $COMP/HGDP_938.header
awk '{print $2}' $COMP/HGDP_938.geno | datamash transpose > $COMP/HGDP_938.header.samples
paste $COMP/HGDP_938.header $COMP/HGDP_938.header.samples >> $COMP/HGDP_938.vcf

# Create information field
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""100""\t"".""\t"".""\t""GT"}' $COMP/HGDP_938.site > $COMP/HGDP_938.information

# Create vcf coded genotypes
cat $COMP/HGDP_938.geno | datamash transpose |  awk 'NR>2' | sed 's/0/a/g' | sed 's/1/b/g' | sed 's/2/c/g' | sed 's/c/0\/0/g' | sed 's/b/0\/1/g' | sed 's/a/1\/1/g' | sed 's/-9/\.\/\./g' > $COMP/HGDP_938.genotypes

# Combine all files and make vcf
paste $COMP/HGDP_938.information $COMP/HGDP_938.genotypes >> $COMP/HGDP_938.vcf
rm $COMP/HGDP_938.information $COMP/HGDP_938.genotypes $COMP/HGDP_938.header.samples $COMP/HGDP_938.header

# Filter HGDP data for pharmacogenomic variants
grep -v ^# $VCF/pgx1000GP.VEP.vcf | awk '{print $3}' > $COMP/1000gpPgx.rs
vcftools --vcf $COMP/HGDP_938.vcf --snps $COMP/1000gpPgx.rs --recode --out $COMP/HGDP_938_pgx 
rm $COMP/HGDP_938.vcf
# Create list of rsIDs for later sub-setting with 1000GP
grep -v ^# $COMP/HGDP_938_pgx.recode.vcf | awk '{print $3}' > $COMP/HGDP_938_pgx.rsID
vcftools --vcf $VCF/pgx1000GP.VEP.vcf --snps $COMP/HGDP_938_pgx.rsID --recode --out $COMP/pgx1000GP.VEP.HGDP

# Create list of samples from each HGDP super population
for pop in africa america east-asia europe south-asia
do
grep -w $pop $DATA/hgdp.sample.info | awk '{print $2}' > $COMP/${pop}.hgdp.sample.id
vcf-subset -c $COMP/${pop}.hgdp.sample.id $COMP/HGDP_938_pgx.recode.vcf > $COMP/${pop}.hgpd.vcf
vcftools --vcf $COMP/${pop}.hgpd.vcf --freq2 --out $COMP/${pop}.hgdp
awk 'NR>1 {print $1":"$2,$6}' $COMP/${pop}.hgdp.frq > $COMP/${pop}.hgdp.frq.temp
echo "variant	freq_hgdp	freq_1000GP" > $COMP/${pop}.hgdp.1000GP.frq
rm $COMP/${pop}.hgpd.vcf
done

# Create list of samples from each 1000GP super population
for pop in AFR AMR EAS EUR SAS
do
vcf-subset -c $FRQ/${pop}.samples.list $COMP/pgx1000GP.VEP.HGDP.recode.vcf > $COMP/${pop}.1000GP.hgdp.vcf
vcftools --vcf $COMP/${pop}.1000GP.hgdp.vcf --freq2 --out $COMP/${pop}.1000GP.hgdp
awk 'NR>1 {print $6}' $COMP/${pop}.1000GP.hgdp.frq > $COMP/${pop}.1000GP.hgdp.temp
rm $COMP/${pop}.1000GP.hgdp.vcf
done

# Combine the relative super-populations from the HGDP and 1000GP
paste $COMP/africa.hgdp.frq.temp $COMP/AFR.1000GP.hgdp.temp >> $COMP/africa.hgdp.1000GP.frq
paste $COMP/america.hgdp.frq.temp $COMP/AMR.1000GP.hgdp.temp >> $COMP/america.hgdp.1000GP.frq
paste $COMP/east-asia.hgdp.frq.temp $COMP/EAS.1000GP.hgdp.temp >> $COMP/east-asia.hgdp.1000GP.frq
paste $COMP/europe.hgdp.frq.temp $COMP/EUR.1000GP.hgdp.temp >> $COMP/europe.hgdp.1000GP.frq
paste $COMP/south-asia.hgdp.frq.temp $COMP/SAS.1000GP.hgdp.temp >> $COMP/south-asia.hgdp.1000GP.frq

rm $COMP/*log $COMP/*temp
