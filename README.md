# FastViFi
A software to detect viral infection and integration sites for different viruses. FastViFi relies on ViFi and Kraken2 softwares originally available on https://github.com/namphuon/ViFi and https://github.com/DerrickWood/kraken2. Please use the forked repositories referenced in the installation section which are modified to integrate into FastViFi pipeline.

Manuscript for this software is in preparation.

# Installation
Follow the installation guides on https://github.com/sara-javadzadeh/ViFi and https://github.com/sara-javadzadeh/kraken2.
To get the reference viral dataset for ViFi, please contact the author by sending an email to (saraj [at] eng [ dot ] ucsd [ dot ] edu).
Once Both tools are installed, build Kraken databases based on the **sample-level** and **read-level** configurations on FastViFi, i.e. one Kraken database including viral and human references with k=25 and two kraken datasets including only viral references with k=18 and k=22. The name of the Kraken databases should follow the format `Kraken2StandardDB_k_<KMER_LEN>_<VIRUS>_hg`, replacing `<KMER_LEN>` by the respective k value mentioned above and `<VIRUS>` with the target virus name among the list: (hpv, hbv, hcv, ebv).

# Running
FastViFi can be used on WGS or RNA-seq files. The recommended process for running FastViFi is to run **sample-level** FastViFi on all samples (lower runtime but less sensitive) followed by running **read-level** FastViFi on samples where viral reads were detected.

```
python run_kraken_vifi_pipeline.py --output-dir $OUTPUT_DIR --input-file $INPUT_BAM --level $FASTVIFI_LEVEL --kraken-path $KRAKEN_PATH --vifi-path $VIFI_PATH --human-chr-list $human_chr_list.txt --virus $VIRUSES
```
`OUTPUT_DIR` is the directory to store all the output files.
`$INPUT_BAM` is the path to the input BAM file. FastViFi also works on a pair of FASTQ files for unaligned paired-end reads. See the manual or contact the author for more information.
`$FASTVIFI_LEVEL` sets the configurations based on the level of sensitivity in FastViFi. Set the `--level` to `sample-level-validation-intermediate` for **sample-level** FastViFi and `sensitive-level-validation-intermediate` for **read-level** FastViFi for each FastViFi call. It is possible to add customized configuration values. For more information, please contact the author.
`$KRAKEN_PATH` and `$VIFI_PATH` are the path to the Kraken and ViFi tool (should be installed before running FastViFi).

When input is in form of reads aligned to the human reference, use `--human-chr-list human_chr_list.txt` to indicate the list of reference names in the `human_chr_file.txt` file (one reference per line). The list of reference names can be extracted by calling `samtools view -H $INPUT_BAM` where `$INPUT_BAM` is the input BAM file.

Indicate the target viruses listed after the flag `--virus`. E.g.: `--virus hbv hcv`. Currently 4 viruses are supported: hpv, hbv, hcv and ebv.To create references for other viruses, use the manual on ViFi https://github.com/sara-javadzadeh/ViFi.


To speedup the runtime, use the flag `--keep-intermediate-files` when running sample-level FastViFi on a sample. Once the sample-level FastViFi is completed, run FastViFi on read-level with the flag `--skip-bwa-filter` on the same sample.

FastViFi will remove all the intermediate and output files if no viral reads were detected. To avoid removing all the files, use the flag `--keep-all-virus-files` to manually inspect the files and troubleshoot in case of an error.
