import argparse
import sys
import os
import stat
from datetime import timedelta
from shutil import copyfile
import subprocess
import time
import pysam

virus_detection_mode = ('sample-level',
                        'sensitive-level')
virus_types = ("all", "hpv", "hbv", "hcv", "ebv", "hpv_655", "hbv_2012")
default_config = {
    virus_detection_mode[0]: {'k1': 25, 'k2': 22, 'u1': 0.8, 'u2': 0.9, 't1': 0.4},
    virus_detection_mode[1]: {'k1': 25, 'k2': 18, 'u1': 0.6, 'u2': 0.8, 't1': 0.8}}

def log_time(log_file):
    formatted_time = subprocess.check_output('date -u "+DATE: %Y-%m-%d TIME: %H:%M:%S"', shell=True)
    log_file.write(formatted_time.decode('utf-8') + os.linesep)

def get_formatted_time(start_time):
    return str(timedelta(seconds=time.time() - start_time))

def is_aligned_reads(input_file):
    return input_file.endswith(".bam") or input_file.endswith(".sam")

def log_error_and_exit(error_message, log_file_pipeline):
    log_file_pipeline.write(error_message + os.linesep)
    print(error_message)
    exit(1)

def check_arguments_validity(args):
    # Set default configuration of parameters if not set by the user.
    if args.k1 is None:
        args.k1 = default_config[args.level]['k1']
    if args.k2 is None:
        args.k2 = default_config[args.level]['k2']
    if args.u1 is None:
        args.u1 = default_config[args.level]['u1']
    if args.u2 is None:
        args.u2 = default_config[args.level]['u2']
    if args.t1 is None:
        args.t1 = default_config[args.level]['t1']
    # Check if the arguments set by the user is in the valid range.
    message = None
    if args.k1 < 0:
        message = "Error: Argument k1 should be a positive number"
    elif args.k2 < 0:
        message = "Error: Argument k2 should be a positive number"
    elif args.t1 < 0 or args.t1 > 1:
        message = "Error: Argument t1 should be a float in [0, 1]"
    elif args.u1 < 0 or args.u1 > 1:
        message = "Error: Argument u1 should be a float in [0, 1]"
    elif args.u2 < 0 or args.u2 > 1:
        message = "Error: Argument u2 should be a float in [0, 1]"
    elif (not is_aligned_reads(args.input_file)) and args.input_file_2 is None:
        message = ("Error: If the input file is not a BAM/SAM file, "
            "two FASTQ files should be provided using --input-file and input-file-2")
    if message is not None:
        print(message)
        exit(1)

def check_kraken_db(db_path, log_file_pipeline):
    if not os.path.exists(db_path):
        error_message ="Error: Kraken database {} not created. ".format(db_path) +\
              "Create the database with corresponding k-mer length " +\
              "and indices."
        log_error_and_exit(error_message, log_file_pipeline)

def remove_if_exists(filenames):
    for filename in filenames:
        if os.path.exists(filename):
            os.remove(filename)

def remove_vifi_output_files(output_dir, prefix):
    vifi_file_names_suffix = [".bam", ".clusters.txt.range", ".fixed.trans.cs.bam",
                              ".misc.bam", ".unknown.bam", ".viral.cs.bam",
                              ".fixed.trans.bam", ".fixed.trans.cs.bam.bai",
                              ".viral.cs.bam.bai", ".clusters.txt",
                              ".trans.bam", ".viral.bam"]
    for file_name_suffix in vifi_file_names_suffix:
        file_name = prefix + file_name_suffix
        remove_if_exists([os.path.join(args.output_dir, file_name)])
    remove_if_exists([os.path.join(args.output_dir, "hmms.txt")])

def remove_kraken_intermediate_files_if_exist(keep_intermediate_files_flag,
                                              first_level_kraken_output, final_level_kraken_output,
                                              first_level_kraken_report, final_level_kraken_report):
    # If args.keep_intermediate_files is false, then the next files are set as os.devnull
    if keep_intermediate_files_flag:
        remove_if_exists([first_level_kraken_output, final_level_kraken_output,
                         first_level_kraken_report, final_level_kraken_report])

def parse_input_args(docker_run=False):
    # Some arguments are required for a regular run, but not required when running on Docker.
    # Note that the same function is used to parse the arguments in both cases.
    # required=docker_run means that if running on Docker mode, the argument is required, otherwise not required.
    # required= not docker_run means that if running on Docker mode, the argument is not required, otherwise is required.
    # For example, kraken_path is required on a regular to run kraken. But is not required on the Docker run because
    # kraken is installed in the prebuilt image.
    parser = argparse.ArgumentParser(
        description="Rapidly detect viral reads from genomic data.")

    dependent_tools_args = parser.add_argument_group('Arguments related to the reference files and dependent tools')
    dependent_tools_args.add_argument('--kraken-path', required= not docker_run,
        help='Path to the Kraken tool.')
    dependent_tools_args.add_argument('--kraken-db-path', required=docker_run,
        help='Path to the directory with the Kraken databases.')
    dependent_tools_args.add_argument('--vifi-path', required= not docker_run,
        help='Path to the ViFi python run_vifi.py script.')
    dependent_tools_args.add_argument('--vifi-human-ref-dir', required=docker_run,
        help='Path to the directory including ViFi human reference. It is the <data_repo> directory when downloaded based on instructions on the README.md file.')
    dependent_tools_args.add_argument('--human-chr-list', required= not docker_run,
        help='Path to the file containing all the human chromosome names used to align the input bam file.\n' +
             'Should be exactly the same name as headers in the input bam file, one reference name per line.\n' +
             'For an example, see test/human_chr_list.txt')
    dependent_tools_args.add_argument('--vifi-viral-ref-dir', required=docker_run,
        help='Path to the directory including ViFi viral reference. \n' +
             'The pre-built viral reference files for HPV, HBV, HCV and EBV are available at the GitHub repository https://github.com/sara-javadzadeh/ViFi.')


    input_output_args = parser.add_argument_group('Arguments related to input and output')
    input_output_args.add_argument('--output-dir', required=True,
        help='The path for the intermediate and output files')
    input_output_args.add_argument("--input-file", required=True,
        help='Path to the input file.\n' +
             'Input file can be either of the following:\n' +
             '- A BAM/SAM file storing the aligned reads to human reference (and possibly the target viral strains)\n' +
             '- A pair of FASTQ files storing the paired end reads with the corresponding reads appearing in the same order.' +
             ' to use this option, provide the first FASTQ paired end read using --input-file and the secone one using --input-file-2')
    input_output_args.add_argument('--input-file-2', default=None,
        help='Path to one of the FASTQ files storing the paired end reads. ' +
             'The FASTQ file containing corresponding read mates should be provided using --input-file argument. ' +
             'Both FASTQ files are required in case FASTQ is the chosen input file.')


    container_args = parser.add_argument_group('Arguments related to running with a container')
    container_args.add_argument('--docker', action="store_true", default=False,
        help='Run FastViFi in Docker mode [False]. This argument is intended for internal use. \n' +
             'To run containerized FastViFi, run run_kraken_vifi_container.py with --docker.')
    container_args.add_argument('--singularity', action="store_true", default=False,
        help='Run FastViFi in container mode with Singularity [False]. This argument is intended for internal use. \n' +
             'To run containerized FastViFi, run run_kraken_vifi_container.py with --singularity.')

    virus_args = parser.add_argument_group('Arguments related to viruses')
    virus_args.add_argument('--virus', default=None, nargs="+",
        help='The virus name used for this experiment e.g. hpv.\n' +
        "Alternatively you can add a custom virus with corresponding Kraken databases and ViFi HMMs.")
    virus_args.add_argument('--keep-all-virus-files', default=False, action="store_true",
        help='Keep all files for viruses where FastViFi could not find any viral reads.\n' +
             'Enabling this flag increases the number of output files signifincantly ' +
             'when searching among all viral databases [False].')


    steps_args = parser.add_argument_group('Arguments related to steps of the pipeline and intermediate files')
    steps_args.add_argument('--keep-intermediate-files', action="store_true", default=False,
        help='Keep the FASTQ files filtered for each filtering step as well as kraken output and report [False].')
    steps_args.add_argument('--skip-bwa-filter', action="store_true", default=False,
        help='Skip the BWA filtering on aligned reads [False]')
    steps_args.add_argument('--skip-kraken-filters', action="store_true", default=False,
        help='Skip the Kraken filtering steps and direct all input to ViFi.\n' +
             'Using this flag drastically increases the runtime [False].')
    steps_args.add_argument('--one-level-kraken', action="store_true", default=False,
        help='Skip the second Kraken filtering step [False].')
    steps_args.add_argument('--skip-vifi', action="store_true", default=False,
        help='Skip the ViFi call. Report the results of Kraken filtering.\n' +
             'This flag results in a high number of false positive and non-viral reads[False].')
    steps_args.add_argument('--skip-vifi-hmms', action="store_true", default=False,
        help='Do not use HMMs when calling ViFi. This does not disable ViFi, only disables the call to HMM alignment.\n' +
             'This flag is useful when there are not many viral strains available,\n' +
             'therefore, the HMMs are not likely to capture conserved or highly variant regions in the viral genome[False].')

    filtering_args = parser.add_argument_group('Arguments related to filtering parameters')
    filtering_args.add_argument('--level', '-l', default=virus_detection_mode[0], choices=virus_detection_mode,
        help='The level of virus detection [{}]'.format(virus_detection_mode[0]))
    filtering_args.add_argument('--k1', type=int, default=None,
        help='The length of k-mers used in the first Kraken filter including human and viral references [{} for {} and {} for {}]'
        .format(default_config[virus_detection_mode[0]]['k1'], virus_detection_mode[0],
                default_config[virus_detection_mode[1]]['k1'], virus_detection_mode[1]))
    filtering_args.add_argument('--k2', type=int, default=None,
        help='The length of k-mers used in the second Kraken filter including only viral references [{} for {} and {} for {}]'
        .format(default_config[virus_detection_mode[0]]['k2'], virus_detection_mode[0],
                default_config[virus_detection_mode[1]]['k2'], virus_detection_mode[1]))
    filtering_args.add_argument('--u1', type=float, default=None,
        help='The parameter u1 reflecting the ratio of unmapped k-mers to consider the read unmapped (and pass) ' +
             'in the first filter. The value should be between 0 and 1 [{} for {} and {} for {}]'
        .format(default_config[virus_detection_mode[0]]['u1'], virus_detection_mode[0],
                default_config[virus_detection_mode[1]]['u1'], virus_detection_mode[1]))
    filtering_args.add_argument('--u2', type=float, default=None,
        help='The parameter u2 reflecting the ratio of unmapped k-mers to consider the read unmapped (and discard) ' +
             'in the first filter. The value should be between 0 and 1 [{} for {} and {} for {}]'
        .format(default_config[virus_detection_mode[0]]['u2'], virus_detection_mode[0],
                default_config[virus_detection_mode[1]]['u2'], virus_detection_mode[1]))
    filtering_args.add_argument('--t1', type=float, default=None,
        help='The threshold on the score of the read to be considered viral. ' +
             'The value should be between 0 and 1 [{} for {} and {} for {}]'
        .format(default_config[virus_detection_mode[0]]['t1'], virus_detection_mode[0],
                default_config[virus_detection_mode[1]]['t1'], virus_detection_mode[1]))
    filtering_args.add_argument('--min-hybrid-support', type=int, default=4,
        help="Minimum number of supporting paired end reads to report a a hybrid junction.")

    dev_args = parser.add_argument_group('Arguments related to development of FastViFi')
    dev_args.add_argument('--gt-viral-path', default=None,
        help='The path to the file containing the read ids for the ground truth viral reads. ' +
             'This is for the experiments on the simulated reads where the ground truth ' +
             'viral reads are known. By providing this option, the true positive, false positive, ' +
             'and false negative counts will be reported.')
    dev_args.add_argument('--kraken-grid-search', action='store_true',
        help='Compute precision and recall based on kraken output Fastq files instead of ViFi output files. ' +
             'Use with --gt-viral-path to get the precision and recall.')

    running_args = parser.add_argument_group('Running arguments')
    running_args.add_argument("--threads", type=int, default=1,
        help='Number of threads to use when running Kraken and ViFi [1]')
    running_args.add_argument('--verbose', '-v', action='store_true',
        help='Print verbose debugging information.')

    args = parser.parse_args()

    check_arguments_validity(args)
    return args

def run_kraken_vifi(virus, args, log_file_pipeline, log_file_pipeline_shell,
                    bwa_filtered_fq_filename_1, bwa_filtered_fq_filename_2):
    # Defining file paths
    if args.kraken_db_path is not None:
        kraken_db_path = args.kraken_db_path
    else:
        kraken_db_path = os.path.dirname(args.kraken_path)
    kraken_db_1 = os.path.join(kraken_db_path, "Kraken2StandardDB_k_{}_{}_hg".format(args.k1, virus))
    kraken_db_2 = os.path.join(kraken_db_path, "Kraken2StandardDB_k_{}_{}".format(args.k2, virus))
    first_level_kraken_filtered_code = os.path.join(args.output_dir, "reads_passing_kraken_first_level_for_virus_{}#.fq".format(virus))
    first_level_kraken_filtered_fq_1 = os.path.join(args.output_dir, "reads_passing_kraken_first_level_for_virus_{}_1.fq".format(virus))
    first_level_kraken_filtered_fq_2 = os.path.join(args.output_dir, "reads_passing_kraken_first_level_for_virus_{}_2.fq".format(virus))
    final_level_kraken_filtered_code = os.path.join(args.output_dir, "reads_passing_kraken_filter_for_virus_{}#.fq".format(virus))
    vifi_input_fq_1 = os.path.join(args.output_dir, "reads_passing_kraken_filter_for_virus_{}_1.fq".format(virus))
    vifi_input_fq_2 = os.path.join(args.output_dir, "reads_passing_kraken_filter_for_virus_{}_2.fq".format(virus))

    log_file_pipeline.write("Running kraken vifi pipline on virus {} with output dir {} with DB {} {}".format(
        virus, args.output_dir, kraken_db_1, kraken_db_2) + os.linesep)

    if args.keep_intermediate_files:
        first_level_kraken_output = os.path.join(args.output_dir, "kraken_output_classify_first_level_for_virus_{}".format(virus))
        final_level_kraken_output = os.path.join(args.output_dir, "kraken_output_classify_final_level_for_virus_{}".format(virus))
        first_level_kraken_report = os.path.join(args.output_dir, "kraken_report_first_level_for_virus_{}".format(virus))
        final_level_kraken_report = os.path.join(args.output_dir, "kraken_report_final_level_for_virus_{}".format(virus))
    else:
        first_level_kraken_output = final_level_kraken_output = first_level_kraken_report = final_level_kraken_report = os.devnull


    start_timer = time.time()
    # Kraken step
    if not args.skip_kraken_filters:
        check_kraken_db(kraken_db_1, log_file_pipeline)
        # T_2 does not affect the results of the kraken filter as the reference is viral only.
        t_2 = 0.5
        log_file_pipeline.write("Running Kraken2 to classify the reads in {} and {}".format(
            bwa_filtered_fq_filename_1, bwa_filtered_fq_filename_2) + os.linesep)
        shell_output = subprocess.check_output(
            "/usr/bin/time -v {} --use-names ".format(args.kraken_path) +
            "--report {} ".format(first_level_kraken_report) +
            "--db {} --threads {} --paired ".format(kraken_db_1, args.threads) +
            "--f-threshold {} --keep-unmapped-reads ".format(args.t1) +
            "--unmapped-threshold {} ".format(args.u1) +
            "--classified-out {} {} {} ".format(first_level_kraken_filtered_code,
                bwa_filtered_fq_filename_1, bwa_filtered_fq_filename_2) +
            " --output {} ".format(first_level_kraken_output), shell=True)
        log_file_pipeline_shell.write(shell_output)

        # Skip the second level filtering if indicated in the arguments.
        # One level kraken is not the most optimal filter.
        # This flag is for recording, comparing and running experiments.
        if args.one_level_kraken:
            copyfile(first_level_kraken_filtered_fq_1, vifi_input_fq_1)
            copyfile(first_level_kraken_filtered_fq_2, vifi_input_fq_2)
        else:
            check_kraken_db(kraken_db_2, log_file_pipeline)
            shell_output = subprocess.check_output(
                "/usr/bin/time -v {} --use-names ".format(args.kraken_path) +
                "--report {} ".format(final_level_kraken_report) +
                "--db {} --threads {} --paired ".format(kraken_db_2, args.threads) +
                "--f-threshold {} ".format(t_2) +
                "--unmapped-threshold {} ".format(args.u2) +
                "--classified-out {} ".format(final_level_kraken_filtered_code) +
                "{} {} ".format(first_level_kraken_filtered_fq_1, first_level_kraken_filtered_fq_2) +
                "--output {} ".format(final_level_kraken_output), shell=True)
            log_file_pipeline_shell.write(shell_output)
        # if there are no viral reads passing kraken filter, skip the next steps.
        if not os.path.exists(vifi_input_fq_1) and not args.keep_all_virus_files:
            log_file_pipeline.write("Warning: No reads passed kraken filter for virus {}. ".format(virus) +
                  "Skipping ViFi step. Consider changing the virus type or kraken parameters." + os.linesep)
            remove_kraken_intermediate_files_if_exist(keep_intermediate_files_flag=args.keep_intermediate_files,
                      first_level_kraken_output=first_level_kraken_output,
                      final_level_kraken_output=final_level_kraken_output,
                      first_level_kraken_report=first_level_kraken_report,
                      final_level_kraken_report=final_level_kraken_report)
            return
        num_lines = sum([1 for line in open(vifi_input_fq_1)])
        log_file_pipeline.write('In total {} paired end reads passed kraken filtering in {} time'.format(
            num_lines/4, get_formatted_time(start_timer)) + os.linesep)
        log_time(log_file_pipeline)

    else:
        if not os.path.exists(vifi_input_fq_1) or not os.path.exists(vifi_input_fq_2):
            vifi_input_fq_1 = bwa_filtered_fq_filename_1
            vifi_input_fq_2 = bwa_filtered_fq_filename_2
            #copyfile(bwa_filtered_fq_filename_1, vifi_input_fq_1)
            #copyfile(bwa_filtered_fq_filename_2, vifi_input_fq_2)

    start_timer = time.time()

    # ViFi step
    vifi_output_prefix = "output_" + str(virus)
    vifi_docker_arg = " --docker "
    if args.docker:
        vifi_docker_arg = ""

    command = "/usr/bin/time -v python {} ".format(args.vifi_path) +\
              "{}".format(vifi_docker_arg) +\
              "-f {} -r {} ".format(vifi_input_fq_1, vifi_input_fq_2) +\
              "-o {} --virus {} ".format(args.output_dir, virus) +\
              "-c {} ".format(args.threads) +\
              "--prefix {} ".format(vifi_output_prefix) + \
              "--min-support {} ".format(args.min_hybrid_support)
    if not args.skip_vifi:
        log_file_pipeline.write("Starting ViFi classification" + os.linesep)
        if args.vifi_human_ref_dir is not None:
            command += ' --hg_data_dir {} '.format(args.vifi_human_ref_dir)
        if args.level == 'sensitive-level-validation-intermediate':
            command += " --sensitive"
        if args.skip_vifi_hmms:
            command += "  --disable_hmms"
        if args.verbose:
            log_file_pipeline.write(command + os.linesep)
        shell_output = subprocess.check_output(
            command, shell=True)
        log_file_pipeline_shell.write(shell_output)

        log_file_pipeline.write('ViFi completed in {} time'.format(get_formatted_time(start_timer)) + os.linesep)
        log_time(log_file_pipeline)
    else:
        message = "Warning! Running ViFi is skipped. The final report of viral reads may be incomplete unless kraken-grid-search is set and a ground truth file is provided."
        print(message)
        log_file_pipeline.write(message + os.linesep)

    # Report detected viral reads
    num_viral_reads = write_viral_reads_report(args=args, virus=virus,
                                               prefix=vifi_output_prefix,
                                               log_file_pipeline=log_file_pipeline,
                                               vifi_input_fq_1=vifi_input_fq_1)
    if not args.skip_vifi:
        # Remove intermediate files
        if not args.keep_intermediate_files:
            remove_if_exists([vifi_input_fq_1, vifi_input_fq_2,
                             first_level_kraken_filtered_fq_1,
                             first_level_kraken_filtered_fq_2])

        # Remove all ViFi files if no viral read was found.
        # This is a helpful feature to reduce output files when searching among all viral databases.
        if num_viral_reads == 0 and not args.keep_all_virus_files:
            log_file_pipeline.write("Removing Vifi and kraken files for virus {} as no viral reads were detected".format(virus) + os.linesep)
            remove_vifi_output_files(output_dir=args.output_dir, prefix=vifi_output_prefix)
            remove_if_exists([vifi_input_fq_1, vifi_input_fq_2,
                             first_level_kraken_filtered_fq_1,
                             first_level_kraken_filtered_fq_2])
            remove_kraken_intermediate_files_if_exist(keep_intermediate_files_flag=args.keep_intermediate_files,
                      first_level_kraken_output=first_level_kraken_output,
                      final_level_kraken_output=final_level_kraken_output,
                      first_level_kraken_report=first_level_kraken_report,
                      final_level_kraken_report=final_level_kraken_report)


def write_viral_reads_report(args, virus, prefix, log_file_pipeline, vifi_input_fq_1):
    if args.kraken_grid_search:
        all_viral_read_ids = []
        kraken_output_file = pysam.FastxFile(vifi_input_fq_1)
        for read in kraken_output_file:
            all_viral_read_ids.append(read.name.replace("/1", ""))
        #print(all_viral_read_ids)
    else:
        # Extract the viral reads
        viral_aligned_file = pysam.AlignmentFile(os.path.join(args.output_dir, "{}.viral.bam".format(prefix)), 'rb')
        trans_viral_file = pysam.AlignmentFile(os.path.join(args.output_dir, "{}.trans.bam".format(prefix)), 'rb')
        fixed_trans_viral_read_ids = []
        viral_aligned_read_ids = []
        trans_viral_read_ids = []
        for read in viral_aligned_file.fetch(until_eof=True):
            viral_aligned_read_ids.append(read.query_name)
        for read in trans_viral_file.fetch(until_eof=True):
            trans_viral_read_ids.append(read.query_name)
        fixed_trans_filename = os.path.join(args.output_dir, "{}.fixed.trans.bam".format(prefix))
        if os.path.exists(fixed_trans_filename):
            fixed_trans_viral_file = pysam.AlignmentFile(fixed_trans_filename, 'rb')
            for read in fixed_trans_viral_file.fetch(until_eof=True):
                fixed_trans_viral_read_ids.append(read.query_name)

        all_viral_read_ids = list(set(trans_viral_read_ids + viral_aligned_read_ids + fixed_trans_viral_read_ids))

        log_file_pipeline.write("In total {} paired end viral reads are detected for {}".format(len(all_viral_read_ids), virus) + os.linesep)
        log_file_pipeline.write("Of which {} single reads are fixed trans read_ids, {} viral aligned read_ids, and {} trans read_ids".format(
            len(fixed_trans_viral_read_ids), len(viral_aligned_read_ids), len(trans_viral_read_ids)) + os.linesep)
        log_file_pipeline.write("Of which {} are unique to trans, {} are unique to fixed_trans".format(
            len(set(trans_viral_read_ids).difference(set(viral_aligned_read_ids).union(fixed_trans_viral_read_ids))),
            len(set(fixed_trans_viral_read_ids).difference(set(viral_aligned_read_ids).union(trans_viral_read_ids)))) + os.linesep + os.linesep)

    if args.gt_viral_path is not None:

        gt_viral_file = open(args.gt_viral_path, 'r')
        gt_viral_read_ids = list(set([line.strip() for line in gt_viral_file.readlines()]))
        tp = len([read_id for read_id in gt_viral_read_ids if read_id in all_viral_read_ids])
        fp = len([read_id for read_id in all_viral_read_ids if read_id not in gt_viral_read_ids])
        fn = len([read_id for read_id in gt_viral_read_ids if read_id not in all_viral_read_ids])
        precision = (tp) / (tp + fp +sys.float_info.epsilon)
        recall = (tp) / (tp + fn +sys.float_info.epsilon)
        log_file_pipeline.write("TP: {}, FP: {}, FN: {}, precision: {}, recall: {}".format(
        tp, fp, fn, precision, recall) + os.linesep)
    return len(all_viral_read_ids)


def run_pipeline(args):
    bwa_filtered_filename_prefix = "reads_passing_bwa_filter"
    bwa_filtered_fq_filename_1 = os.path.join(args.output_dir, "reads_passing_bwa_filter_1.fq")
    bwa_filtered_fq_filename_2 = os.path.join(args.output_dir, "reads_passing_bwa_filter_2.fq")

    if os.path.isdir(args.output_dir):
        message = "Warning! output directory already exists. Rewriting previous outputs."
        print(message)
    else:
        os.makedirs(args.output_dir)
    log_file_path = os.path.join(args.output_dir, "log_fastvifi.txt")
    verbose_log_file_path = os.path.join(args.output_dir, "log_fastvifi_verbose.txt")
    log_file_pipeline = open(log_file_path, 'w+')
    log_file_pipeline_shell = open(verbose_log_file_path, 'wb+')
    log_file_pipeline.write("Processing file {} to output {}".format(
        args.input_file, args.output_dir) + os.linesep)
    log_file_pipeline.write("Arguments and values:" + os.linesep)
    for arg, value in sorted(vars(args).items()):
        log_file_pipeline.write("{} : {}".format(arg, value) + os.linesep)

    log_time(log_file_pipeline)
    start_timer = time.time()
    human_chr_list = args.human_chr_list

    # Filter aligned reads and discard reads mapping to human reference early on.
    if not args.skip_bwa_filter:
        if is_aligned_reads(args.input_file):
            log_file_pipeline.write("Filtering reads from input file {}".format(args.input_file) + os.linesep)
            bwa_filter_file = "filter_reads_bwa_efficient.py"
            if args.docker:
                bwa_filter_file = "/home/fastvifi/" + bwa_filter_file

            shell_output = subprocess.check_output("/usr/bin/time -v python {} {} {} {} {} &>> output".format(
                bwa_filter_file,
                args.input_file, args.output_dir, human_chr_list, bwa_filtered_filename_prefix), shell=True)
            log_file_pipeline_shell.write(shell_output)

            num_lines = sum([1 for line in open(bwa_filtered_fq_filename_1)])
            log_file_pipeline.write('In total {} paired end reads passed bwa alignment filtering in {} time'.format(
                num_lines/4, get_formatted_time(start_timer)) + os.linesep)
            log_time(log_file_pipeline)

        else:
            log_file_pipeline.write("Skipping BWA filtering step with FASTQ input reads")
            pass
    else:
        # This is if we want to skip BWA and feed the two input files to kraken (e.g. test datasets)
        # Otherwise we don't copy the files so that same files in previous runs are kept which leads to
        # faster debugging process.
        if (args.input_file_2 is not None and
                bwa_filtered_fq_filename_1 != args.input_file and
                bwa_filtered_fq_filename_2 != args.input_file_2):
            bwa_filtered_fq_filename_1 = os.path.join(args.output_dir, "input_file_1.fq")
            bwa_filtered_fq_filename_2 = os.path.join(args.output_dir, "input_file_2.fq")
            copyfile(args.input_file, bwa_filtered_fq_filename_1)
            copyfile(args.input_file_2, bwa_filtered_fq_filename_2)

    selected_viruses = []
    # If virus type is not chosen, run Kraken and ViFi on all viruses.
    # Otherwise, run only on the specified viruses.
    selected_viruses = args.virus
    if args.virus == None:
        selected_viruses = virus_types[1:]

    for virus in selected_viruses:
        run_kraken_vifi(virus=virus, args=args, log_file_pipeline=log_file_pipeline,
                        log_file_pipeline_shell=log_file_pipeline_shell,
                        bwa_filtered_fq_filename_1=bwa_filtered_fq_filename_1,
                        bwa_filtered_fq_filename_2=bwa_filtered_fq_filename_2)
    if not args.keep_intermediate_files:
        remove_if_exists([bwa_filtered_fq_filename_1, bwa_filtered_fq_filename_2])
    if args.docker:
        os.system("chmod 777 -R {}".format(os.path.join(args.output_dir, "tmp")))

if __name__ == "__main__":
    args = parse_input_args()

    run_pipeline(args)
