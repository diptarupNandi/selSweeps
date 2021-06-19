plink --vcf \$vcfPath/driftOnly_500Ne_100ind_20Mb_biAll.vcf.gz --cm-map \$vcfPath/genetic_map_chr1_simulations_20Mb.txt 1 --recode --out \$vcfPath/driftOnly_500Ne_100ind_20Mb_biAll

sed '/^#/d' driftOnly_500kNe_100ind.vcf | awk '$2 < 10000000 {print $0}' >> driftOnly_chr1_500Ne_100ind.vcf

sed '/^#/d' driftOnly_5000kNe_100ind.vcf | awk '$2 < 10000000 {print $0}' >> driftOnly_chr1_5kNe_100ind.vcf

awk '/^#/ {print $0}' driftOnly_500kNe_100ind.vcf > driftOnly_chr2_500Ne_100ind.vcf

# Printing only the matched lines
sed '/^#/d' driftOnly_500kNe_100ind.vcf | awk -F'\t' -v OFS='\t' '$2 > 9999999 {$1=2; $2 =$2-10000000; print $0}' >> driftOnly_chr2_500Ne_100ind.vcf

bgzip -c driftOnly_chr1_500Ne_100ind.vcf > ./simul_VCFs/driftOnly_chr1_500Ne_100ind_20Mb.vcf.gz

tabix -p vcf ./simul_VCFs/driftOnly_chr1_500Ne_100ind_20Mb.vcf.gz

bcftools view -Oz -m 2 -M 2 -v snps simul_VCFs/driftOnly_chr1_5kNe_100ind_20Mb.vcf.gz -o simul_VCFs/VCFs_filtrd_20Mb/driftOnly_chr1_5kNe_100ind_20Mb_biAll.vcf.gz


plink --vcf \$vcfPath/driftOnly_500Ne_100ind_20Mb_biAll.vcf.gz --cm-map \$vcfPath/genetic_map_chr1_simulations_20Mb.txt 1 --recode --out \$vcfPath/driftOnly_500Ne_100ind_20Mb_biAll
