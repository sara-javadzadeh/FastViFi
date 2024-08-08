#!/bin/bash

if [ ! $# -eq 2 ]; then
    echo "Usage: bash $0 <input_bam> <output_prefix>"
    exit 1
fi
input=$1
output=$2
output_dir=$(dirname $output)
echo "Running filtering with input $input and output ${output}.bam"

echo "Finding filtered read ids"
samtools view -h $input | awk -e '{ if(substr($1,1,1) != "@") {if (and($2,12) > 0  || $3 !~ "chr*") print $1}}' | sort --uniq > ${output_dir}/filtered_read_ids.txt
echo "Copying the header"
samtools view -H $input > ${output}.sam
echo "Copying the target reads (and mates)"
samtools view $input | grep -f ${output_dir}/filtered_read_ids.txt >> ${output}.sam
echo "Generating bam file"
samtools view -h -Sb ${output}.sam > ${output}.bam
echo "Computing the number of reads"
num_reads=$(samtools view -c ${output}.bam)
echo "Result has $num_reads reads"
