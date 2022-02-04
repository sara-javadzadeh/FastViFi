# FastViFi
A software to detect viral infection and integration sites for different viruses. FastViFi relies on ViFi and Kraken2 tools. Please use the forked repositories referenced in the installation section which are modified to integrate into FastViFi pipeline.

Manuscript for this software is in preparation.

# Installation
Follow the installation guides on https://github.com/sara-javadzadeh/ViFi and https://github.com/sara-javadzadeh/kraken2. Install the tools in any directory. To run FastViFi, the path to ViFi and Kraken should be passed through flags.

### Pre-built Kraken databases

Once Both tools are installed, build Kraken databases based on the **sample-level** and **read-level** configurations on FastViFi, namely:
- One Kraken database including viral and human references with k=25 and
Plus, for **sample-level**:
- One Kraken datasets including only viral references with k=22
and for **read-level**:
- One Kraken datasets including only viral references with k=18

Each Kraken dataset typically takes 30GB-40GB space.

Kraken databses for **sample-level** FastViFi for HPV are available on https://drive.google.com/file/d/1QYn5lDWjvhtIWCrwmzDc_1fy8ANrXWz1/view?usp=sharing. After downloading the compressed file, run the following command to unpack:
`tar -xzvf kraken_datasets.tar.gz`

Then, move both datasets to the path where Kraken is installed.


### Custom Kraken databases

To build the custom databases, use the script `build_custom_kraken_database.sh` on the GitHub repository https://github.com/sara-javadzadeh/kraken2. The script automatically creates all the three databases for a given virus name and corresponding reference file. It takes time to download the Kraken taxonomy and genomic information. We recommend using the same set of viruses for Kraken2 databases and ViFi to simplify the analysis. The viral reference for each virus can be obtained from the ViFi repository (https://github.com/sara-javadzadeh/ViFi). After cloning the repository and un-compressing the file `viral_data.tar.gz`, the viral reference for each of the desired viruses can be found on `viral_data/<VIRUS>/<VIRUS>.unaligned.fas*`  where `<VIRUS>` can be either of the following: hpv, hbv, hcv or ebv.

Before proceeding with the custom Kraken dataset, check the following in the Kraken databases:

- The name of the Kraken databases should follow the format:
  - `Kraken2StandardDB_k_<KMER_LEN>_<VIRUS>_hg` for the first kraken filter including human reference and
  - `Kraken2StandardDB_k_<KMER_LEN>_<VIRUS>` for the next two kraken filters with the viral references only.
  -  By replacing `<KMER_LEN>` by the respective k value mentioned above and `<VIRUS>` with the target virus name among the list: (hpv, hbv, hcv, ebv).
- In each database file, check for presence of the following three index files: hash.k2d, opts.k2d and taxo.k2d
You can further inspect each database with `kraken2-inspect` script in Kraken2 directory by running `./kraken2-inspect --db <database name>`

# Toy Example
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

# Running
FastViFi can be used on WGS or RNA-seq files. The recommended process for running FastViFi is to run **sample-level** FastViFi on all samples (lower runtime but less sensitive) followed by running **read-level** FastViFi on samples where viral reads were detected. Use the following command to run FastViFi:

```
python run_kraken_vifi_pipeline.py --output-dir $OUTPUT_DIR --input-file $INPUT_BAM --level $FASTVIFI_LEVEL --kraken-path $KRAKEN_PATH --vifi-path $VIFI_PATH --human-chr-list $human_chr_list.txt --virus $VIRUSES
```
- `OUTPUT_DIR` is the directory to store all the output files.
- `$INPUT_BAM` is the path to the input BAM file. FastViFi also works on a pair of FASTQ files for unaligned paired-end reads. See the manual or contact the author for more information.
- `$FASTVIFI_LEVEL` sets the configurations based on the level of sensitivity in FastViFi. Set the `--level` to `sample-level-validation-intermediate` for **sample-level** FastViFi and `sensitive-level-validation-intermediate` for **read-level** FastViFi for each FastViFi call. It is possible to add customized configuration values. For more information, please contact the author.
- `$KRAKEN_PATH` is the path to the Kraken2 executable file named kraken2 within the Kraken2 directory. See `Installation` section for installing Kraken2 and setting up customized databases prior to running FastViFi.
- `$VIFI_PATH` is the path to ViFi main python script named `run_vifi.py`. See `Installation` setion for installing ViFi prior to running FastViFi.

When input is in form of reads aligned to the human reference, use `--human-chr-list human_chr_list.txt` to indicate the list of reference names in the `human_chr_file.txt` file (one reference per line). The list of reference names can be extracted by calling `samtools view -H $INPUT_BAM` where `$INPUT_BAM` is the input BAM file.

Indicate the target viruses listed after the flag `--virus`. E.g.: `--virus hbv hcv`. Currently 4 viruses are supported: hpv, hbv, hcv and ebv.To create references for other viruses, use the manual on ViFi https://github.com/sara-javadzadeh/ViFi.


To speedup the runtime, use the flag `--keep-intermediate-files` when running sample-level FastViFi on a sample. Once the sample-level FastViFi is completed, run FastViFi on read-level with the flag `--skip-bwa-filter` on the same sample. Note that the bwa-filter (when `--skip-bwa-filter` is not set) is in fact alignment based filter. Rather than running BWA, it depends on alignment information on the input BAM file.

FastViFi will remove all the intermediate and output files if no viral reads were detected. To avoid removing all the files, use the flag `--keep-all-virus-files` to manually inspect the files and troubleshoot in case of an error.
