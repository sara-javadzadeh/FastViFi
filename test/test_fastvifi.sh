#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: bash $0 <path to kraken directory> <path to where kraken databases are stored> <path to vifi directory>"
    exit 1
fi

kraken_path=$1
kraken_bin_path="${kraken_path}/kraken2"
kraken_db_path=$2
vifi_path=$3
vifi_script_path="${vifi_path}/scripts/run_vifi.py"

input_file_1="$(pwd)/test_reads_1.fq"
input_file_2="$(pwd)/test_reads_2.fq"
output_file="$(pwd)/test_output.txt"

echo "Running toy sample with input files $input_file_1 and $input_file_2" >> $output_file
output_dir="$(pwd)/fastvifi_output_files"

if [ ! -d $output_dir ]; then
	echo "Creating the directory for FastViFi outputs: $output_dir" > $output_file
	mkdir $output_dir
else
	echo "Output directory exists: $output_dir" > $output_file
fi

human_chr_list="$(pwd)/human_chr_list.txt"

/usr/bin/time -v python ../run_kraken_vifi_pipeline.py --output-dir $output_dir --input-file $input_file_1 --input-file-2 $input_file_2 --level sample-level-validation-intermediate --kraken-path $kraken_bin_path --kraken-db-path $kraken_db_path --vifi-path $vifi_script_path --virus hpv --human-chr-list $human_chr_list --skip-bwa-filter --keep-intermediate-files &>> $output_file

