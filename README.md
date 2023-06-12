# FastViFi
A software that detects viral infection and integration sites for different viruses. FastViFi relies on ViFi and Kraken2 tools. Please use the forked repositories referenced in the installation section which are modified to integrate into FastViFi pipeline.

Manuscript for this software is in preparation.

# Download Data Files
In order to run FastViFi, the followind data files should be downloaded.

### Pre-built Kraken2 Datasets
Kraken datasets for **sample-level** FastViFi for HPV virus are available on https://drive.google.com/file/d/1QYn5lDWjvhtIWCrwmzDc_1fy8ANrXWz1/view?usp=sharing. After downloading the compressed file, run the following command to unpack:
`tar -xzvf kraken_datasets.tar.gz`

We recommend running FastViFi with **sample-level** first. It is possible to customize the level of sensitivity (with the cost of runtime) with a different set of parameters to FastViFi which require different Kraken2 datasets. More information on customizing FastViFi is provided in the FastViFi manuscript (in preparation). For more information, contact: saraj@eng.ucsd.edu.

If Kraken2 is installed locally, move both datasets to the path where Kraken2 is installed. Otherwise, to run FastViFi based on the Docker image, the directory to the Kraken2 datasets should be provided to `run_kraken_vifi_docker.py` with the `--kraken-db-path` argument.
To create custom datasets, see the "Custom Kraken datasets" section below.

### ViFi Datasets
If ViFi is installed locally, the viral references is downloaded at the time of cloning ViFi. Otherwise, to run FastViFi based on the Docker image, download the viral references from the ViFi repository (https://github.com/sara-javadzadeh/ViFi). After cloning the repository and un-compressing the file `viral_data.tar.gz`, the viral reference for each of the desired viruses can be found on `viral_data/<VIRUS>/<VIRUS>.unaligned.fas*`  where `<VIRUS>` can be either of the following: hpv, hbv, hcv or ebv. Note that to run FastViFi Docker image, the viral reference directory should be provided to `run_kraken_vifi_docker.py` with the `--vifi-viral-ref-dir` argument.

Whether or not ViFi is installed locally, the human reference data file should be downloaded from https://drive.google.com/file/d/1XBZbwgcV1n2AWWAyt2RWfSKKxzssRFBo/view?usp=share_link. If ViFi is locally installed, the downloaded file should be placed inside the ViFi directory. Otherwise, to run FastViFi based on the Docker image, the human reference directory should be passed to `run_kraken_vifi_docker.py` with the `--vifi-human-ref-dir` argument.

# Installation
There are two ways to run FastViFi: using the Docker image or locally installing dependencies. The Docker image includes pre-installed dependencies (Kraken and ViFi). However, to create a smaller sized Docker image, the data files are not included in the image. To download the necessary data files, refer to the "Download Data Files" section. Follow the instructions on either section "Running FastViFi Docker Image" or "Running FastViFi by Locally Installing Dependencies". Correct installation can be verified by running the toy example (see section "Toy Example").

## Running FastViFi Docker Image
To run the Docker image for FastViFi, first download the image from DockerHub:
`docker pull sarajava/fastvifi`
Next, make sure the necessary datasets are installed (see section "Download Data Files").

To run FastViFi with the Docker image, run:
`python run_kraken_vifi_docker.py --input-file <input file (BAM or FASTQ)> [--input-file-2 <input file 2 (for the second FASTQ file)>] --output-dir <output directory> --virus <virus> --kraken-db-path <path to where the downloaded kraken datasets> --vifi-viral-ref-dir <viral referece directory> --vifi-human-ref-dir <path to the human reference directory>`

If the input is a bam file, provide the path to the file with `--input-file` argument. Ignore the --input-file-2 argument. However, if the input is in form of a pair of FASTA/FASTQ files (for paired end reads), use the `--input-file-2` argument to provide the path to the second file and use `--skip-bwa-filter`. The `--skip-bwa-filter` refers to filtering reads based on alignment into the human genome in an already aligned set of reads, which is only usable when the input is a BAM file.

Output directory should be writable by other users (i.e., the user in the Docker container). `run_kraken_vifi_docker.py` will change the access permissions if that is not the case.

For the `--virus` argument, we currently provide datasets for all of the following: hpv, hbv, hcv or ebv. Note that multiple viruses can be provided to FastViFi by adding additional `--virus` arguments.

Arguments `--kraken-db-path`, `--vifi-viral-ref-dir` and `--vifi-human-ref-dir` are dedicated to providing the path to the corresponding dataset. See section "Download Data Files" for instructions on how to download the datasets.

## Running FastViFi by Locally Installing Dependencies
Follow the installation guides on https://github.com/sara-javadzadeh/ViFi and https://github.com/sara-javadzadeh/kraken2. Install the tools in any directory. To run FastViFi, the path to ViFi and Kraken should be passed through flags.


FastViFi can be used on WGS or RNA-seq files. The recommended process for running FastViFi is to run **sample-level** FastViFi on all samples (lower runtime but less sensitive) followed by running **read-level** FastViFi on samples where viral reads were detected. Use the following command to run FastViFi:

```
python run_kraken_vifi_pipeline.py --output-dir $OUTPUT_DIR --input-file $INPUT_BAM --level $FASTVIFI_LEVEL --kraken-path $KRAKEN_PATH --vifi-path $VIFI_PATH --human-chr-list $human_chr_list.txt --virus $VIRUSES
```
- `OUTPUT_DIR` is the directory to store all the output files.
- `$INPUT_BAM` is the path to the input BAM file. FastViFi also works on a pair of FASTQ files for unaligned paired-end reads. See the manual or contact the author for more information.
- `$FASTVIFI_LEVEL` sets the configurations based on the level of sensitivity in FastViFi. Set the `--level` to `sample-level-validation-intermediate` for **sample-level** FastViFi and `sensitive-level-validation-intermediate` for **read-level** FastViFi for each FastViFi call. It is possible to add customized configuration values. For more information, please contact the author.
- `$KRAKEN_PATH` is the path to the Kraken2 executable file named kraken2 within the Kraken2 directory. See `Installation` section for installing Kraken2 and setting up customized datasets prior to running FastViFi.
- `$VIFI_PATH` is the path to ViFi main python script named `run_vifi.py`. See `Installation` setion for installing ViFi prior to running FastViFi.

When input is in form of reads aligned to the human reference, use `--human-chr-list human_chr_list.txt` to indicate the list of reference names in the `human_chr_file.txt` file (one reference per line). The list of reference names can be extracted by calling `samtools view -H $INPUT_BAM` where `$INPUT_BAM` is the input BAM file.

Indicate the target viruses listed after the flag `--virus`. E.g.: `--virus hbv hcv`. Currently 4 viruses are supported: hpv, hbv, hcv and ebv.To create references for other viruses, use the manual on ViFi https://github.com/sara-javadzadeh/ViFi.


To speedup the runtime, use the flag `--keep-intermediate-files` when running sample-level FastViFi on a sample. Once the sample-level FastViFi is completed, run FastViFi on read-level with the flag `--skip-bwa-filter` on the same sample. Note that the bwa-filter (when `--skip-bwa-filter` is not set) is in fact alignment based filter. Rather than running BWA, it depends on alignment information on the input BAM file.

FastViFi will remove all the intermediate and output files if no viral reads were detected. To avoid removing all the files, use the flag `--keep-all-virus-files` to manually inspect the files and troubleshoot in case of an error.

## Toy Example
To test FastViFi against a toy example, change the directory to `test` and run `bash ./test_fastvifi.sh <path to the kraken2 directory> <path to ViFi directory>`.

Input files for the toy sample inlude:
- 6 paired end reads originated from HPV16,
- 4 Paired end reads originated from the human genome, and
- 4 Paired end reads where one mate is originated from the human genome and the other is originated from HPV185 (not in the viral references).
The origin of each read in the toy example is indicated in the read name.

The expected behaviour for FastViFi is to filter out the human originated reads, report reads originated from HPV16 as viral and report a hybrid human-viral junction with the HPV185 and human genome reads. To check the output, change the directory to `fastvifi_output_files`.
- The file `output_hpv.viral.bam` stores all the reads where both mates were mapped to the reference viral genomes. If Samtools is installed, Run `samtools view output_hpv.viral.bam` to see 12 entries (6 paired end reads) mapped to the viral reference.
- The file `output_hpv.fixed.trans.bam` stores the reads where one mate is mapped to the human genome and the other is mapped to viral genomes. Run `samtools view output_hpv.fixed.trans.bam` and look for 8 entries (4 paired end reads) that represent the presence of hybrid human-viral DNA.
- The file `output_hpv.clusters.txt` stores all the hybrid human-viral DNA/RNA junctions clusters. There should be a single cluster with the 4 paired end reads.

## Custom Kraken datasets

Touild Kraken datasets based on the **sample-level** or **read-level** configurations on FastViFi, the following datasets are necessary:
- One Kraken dataset including viral and human references with k=25 and
Plus, for **sample-level**:
- One Kraken dataset including only viral references with k=22
or for **read-level**:
- One Kraken dataset including only viral references with k=18
Each Kraken dataset typically takes 30GB-40GB space.

To build the custom dataset, use the script `build_custom_kraken_database.sh` on the GitHub repository https://github.com/sara-javadzadeh/kraken2. The script automatically creates all the three datasets for a given virus name and corresponding reference file. It takes time to download the Kraken taxonomy and genomic information. We recommend using the same set of viruses for Kraken2 datasets and ViFi to simplify the analysis.

Before proceeding with the custom Kraken dataset, check the following in the Kraken datasets:

- The name of the Kraken datasets should follow the format:
  - `Kraken2StandardDB_k_<KMER_LEN>_<VIRUS>_hg` for the first kraken filter including human reference and
  - `Kraken2StandardDB_k_<KMER_LEN>_<VIRUS>` for the next two kraken filters with the viral references only.
  -  By replacing `<KMER_LEN>` by the respective k value mentioned above and `<VIRUS>` with the target virus name among the list: (hpv, hbv, hcv, ebv).
- In each dataset directory, check for presence of the following three index files: hash.k2d, opts.k2d and taxo.k2d
You can further inspect each dataset with `kraken2-inspect` script in Kraken2 directory by running `./kraken2-inspect --db <dataset name>`

## Publication

Javadzadeh, Sara, et al. "FastViFi: Fast and accurate detection of (Hybrid) Viral DNA and RNA." NAR genomics and bioinformatics 4.2 (2022): lqac032.

https://academic.oup.com/nargab/article/4/2/lqac032/6574445
