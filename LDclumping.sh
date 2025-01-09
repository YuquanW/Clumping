#!/bin/bash

# Script for LD Clumping
# Usage: ./LDclumping.sh <input_exp_path> <input_out_path> <output_exp_path> <output_out_path> <input_ld_path> <extract_snps_path>

# Ensure all arguments are passed
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_exp_path> <input_out_path> <output_exp_path> <output_out_path> <input_ld_path> <extract_snps_path>"
    exit 1
fi

# Input parameters
input_exp_path=$1  # Input exposure VCF file
input_out_path=$2  # Input outcome VCF file
output_exp_path=$3  # Output exposure results file
output_out_path=$4  # Output outcome results file
input_ld_path=$5  # Input LD genotype file (plink binary files)
extract_snps_path=$6  # Path to extract_snps.py script

# Temporary directory for intermediate files
temp_dir="./common"
mkdir -p $temp_dir

echo "Processing exposure and outcome data..."

# Step 1: Find common variants between input exposure and outcome files
bcftools isec -n=2 -c all -p $temp_dir $input_exp_path $input_out_path

# Step 2: Normalize and remove multiallelic SNPs
bcftools norm -m+ $temp_dir/0000.vcf -o $temp_dir/tmp.vcf
awk '/#/{print;next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' $temp_dir/tmp.vcf > $temp_dir/tmp1.vcf

# Step 3: Remove ambiguous SNPs
awk '/#/{print;next}{if(!($4 ~ /A/ && $5 ~ /T/) && !($4 ~ /T/ && $5 ~ /A/) && !($4 ~ /G/ && $5 ~ /C/) && !($4 ~ /C/ && $5 ~ /G/)){print}}' $temp_dir/tmp1.vcf > $temp_dir/tmp.vcf

# Step 4: Filter by allele frequency
bcftools view -i 'FMT/AF > 0.01 & FMT/AF < 0.99' $temp_dir/tmp.vcf -o $temp_dir/Exposure.vcf

# Step 5: Extract SNPs using the Python script
python $extract_snps_path -i $temp_dir/0001.vcf -r $temp_dir/Exposure.vcf -o $temp_dir/Outcome.vcf

# Step 6: Query clumping data
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%LP\t%SE]\n' $temp_dir/Exposure.vcf -o $temp_dir/clump-1.txt
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%LP\t%SE]\n' $temp_dir/Outcome.vcf -o $temp_dir/clump-2.txt

# Step 7: Calculate P values and merge clumping data
awk 'NR==FNR{C2[NR]=$NF; next} {PS = 2/(1/$NF^2 + 1/C2[FNR]^2); NF--; OFS="\t"; print $1, $2, $3, $4, $5, $6, PS}' \
    $temp_dir/clump-2.txt $temp_dir/clump-1.txt > $temp_dir/clump.txt
sed -i '1s/^/CHR\tPOS\tSNP\tother_allele.exposure\teffect_allele.exposure\tpval.exposure\tP\n/' $temp_dir/clump.txt

# Step 8: Perform LD clumping with plink2
plink2 --bfile $input_ld_path --clump $temp_dir/clump.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.001 --clump-kb 10000 --out $temp_dir/filtered

# Step 9: Filter and save clumped exposure and outcome data
bcftools view -R $temp_dir/filtered.clumped $input_exp_path -o $temp_dir/clumped_exposure.vcf
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP\t%AF\t%SS]\n' $temp_dir/clumped_exposure.vcf -o $output_exp_path
sed -i '1s/^/CHR\tPOS\tSNP\tother_allele.exposure\teffect_allele.exposure\tbeta.exposure\tse.exposure\tpval.exposure\teaf.exposure\tsamplesize.exposure\n/' $output_exp_path

bcftools view -R $temp_dir/filtered.clumped $input_out_path -o $temp_dir/clumped_outcome.vcf
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP\t%AF\t%SS]\n' $temp_dir/clumped_outcome.vcf -o $output_out_path
sed -i '1s/^/CHR\tPOS\tSNP\tother_allele.outcome\teffect_allele.outcome\tbeta.outcome\tse.outcome\tpval.outcome\teaf.outcome\tsamplesize.outcome\n/' $output_out_path

# Clean up
rm -rf $temp_dir

# Output results
echo "Exposure results saved to $output_exp_path."
echo "Outcome results saved to $output_out_path."
