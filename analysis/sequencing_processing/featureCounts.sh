#!/bin/bash

# Extract BAM file names from database.csv
bam_files=$(awk -F"," 'NR>1{print $1".bam"}' database.csv)

# Define species (you may need to set this manually or dynamically)
species=$1  # Replace with your species

# Set the reference GTF file
gtf_file="/reference/$species/ref_nm.gtf"

# Initialize featureCounts options
options="-T 56 -O -s 2 -a $gtf_file -t transcript -g gene_id"

# Check the first BAM file for type
first_bam=$(echo $bam_files | awk '{print $1}')

if [[ -f $first_bam ]]; then
    echo "Checking BAM file: $first_bam"

    # Use samtools to determine if the first BAM file is paired-end
    if samtools flagstat $first_bam | grep -q "0 + 0 paired in sequencing"; then
        echo "$first_bam is single-end."
        paired=""
    else

        echo "$first_bam is paired-end."
        paired="-p"
    fi
else
    echo "First BAM file $first_bam does not exist. Exiting."
    exit 1
fi

# Run featureCounts for all BAM files
echo "Running featureCounts for all BAM files..."
featureCounts $options $paired -o featureCounts.txt $bam_files
