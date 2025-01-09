import argparse

def extract_subset_vcf(input_vcf, reference_vcf, output_vcf):
    # Step 1: Read the reference VCF and collect CHROM and POS pairs
    snps_set = set()
    with open(reference_vcf, 'r') as ref_file:
        for line in ref_file:
            if not line.startswith('#'):  # Skip header lines
                cols = line.strip().split('\t')
                chrom = cols[0]  # CHROM is the first column
                pos = cols[1]    # POS is the second column
                snps_set.add((chrom, pos))
    
    # Step 2: Read the input VCF file and filter by CHROM and POS
    with open(input_vcf, 'r') as in_file, open(output_vcf, 'w') as out_file:
        for line in in_file:
            if line.startswith('#'):  # Copy header lines as is
                out_file.write(line)
            else:
                cols = line.strip().split('\t')
                chrom = cols[0]  # CHROM is the first column
                pos = cols[1]    # POS is the second column
                if (chrom, pos) in snps_set:
                    out_file.write(line)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Extract a subset of VCF based on another VCF.")
    parser.add_argument('-i', '--input', required=True, help="Input VCF file to filter.")
    parser.add_argument('-r', '--reference', required=True, help="Reference VCF file with CHROM and POS to filter by.")
    parser.add_argument('-o', '--output', required=True, help="Output VCF file with the filtered data.")
    
    args = parser.parse_args()

    # Call the function to perform the extraction
    extract_subset_vcf(args.input, args.reference, args.output)

if __name__ == "__main__":
    main()
