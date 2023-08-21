import os
import subprocess
from run_kraken_vifi_pipeline import parse_input_args

def check_arguments(args):
    pass

def call_fastvifi_pipeline(args):
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    print("Changing the mode of the output directory to be writable by other users (i.e., the user in the docker container)")
    os.system("chmod 777 {}".format(args.output_dir))
    if args.kraken_db_path is None:
        print("Error: Kraken database path is not set." + \
        "See https://github.com/sara-javadzadeh/FastViFi for " + \
        "instructions on how to download the databases and provide " + \
        "the path to the program using --kraken-db-path argument.")
        exit(1)
    kraken_db_path = args.kraken_db_path
    vifi_hg_data_path = os.path.join(args.vifi_human_ref_dir)
    vifi_viral_data_path = os.path.join(args.vifi_viral_ref_dir)
    docker_image_tag = "sarajava/fastvifi"
    # Command for Singularity
    if args.singularity:
        # Running as a user
        singularity_image_tag = "fastvifi_latest.sif"
        docker_image_tag = " docker://{}".format(docker_image_tag)
        if not os.path.exists(singularity_image_tag):
            # Build a singularity .sif file from the docker image, if doesn't exist.
            command = "singularity pull {};".format(docker_image_tag) + \
            # This is for building a .sif file from a local docker image.
            #command = "singularity build {} {}; ".format(singularity_image_tag, docker_image_tag)
        else:
            # Once the .sif file is created, it can be re-used for future runs.
            # No need to run the build command again.
            command = ""
        command += " singularity run " + \
        "--bind {}:/home/input/{} ".format(args.input_file, os.path.basename(args.input_file))

        if args.input_file_2 is not None:
            command += "--bind {}:/home/input/{} ".format(args.input_file_2, os.path.basename(args.input_file_2))
        command += "--bind {}:/home/kraken2-db ".format(kraken_db_path) + \
            "--bind {}:/home/data_repo/ ".format(vifi_hg_data_path) + \
            "--bind {}:/home/repo/data/ ".format(vifi_viral_data_path) + \
            "--bind {}:/home/output/ ".format(args.output_dir) + \
            "--bind ./run_kraken_vifi_pipeline.py:/home/fastvifi/run_kraken_vifi_pipeline.py " + \
            "{} ".format(singularity_image_tag) + \
            "python /home/fastvifi/run_kraken_vifi_pipeline.py " + \
            "--kraken-path /home/kraken2/kraken2 " + \
            "--vifi-path /home/ViFi/scripts/run_vifi.py " + \
            "--output /home/output " + \
            "--human-chr-list /home/data_repo/GRCh38/chrom_list.txt " + \
            "--kraken-db-path /home/kraken2-db " + \
            "--docker "
    # Command for Docker
    if args.docker:
        docker_image_tag = docker_image_tag + ":latest"
        command = "docker run --rm " + \
        "--read-only -v {}:/home/input/{} ".format(args.input_file, os.path.basename(args.input_file))
        if args.input_file_2 is not None:
            command += "--read-only -v {}:/home/input/{} ".format(args.input_file_2, os.path.basename(args.input_file_2))

        command += "--read-only -v {}:/home/kraken2-db ".format(kraken_db_path) + \
            "--read-only -v {}:/home/data_repo/ ".format(vifi_hg_data_path) + \
            "--read-only -v {}:/home/repo/data/ ".format(vifi_viral_data_path) + \
            "-v {}:/home/output ".format(args.output_dir) + \
            "{} ".format(docker_image_tag) + \
            "python /home/fastvifi/run_kraken_vifi_pipeline.py " + \
            "--kraken-path /home/kraken2/kraken2 " + \
            "--vifi-path /home/ViFi/scripts/run_vifi.py " + \
            "--output /home/output " + \
            "--human-chr-list /home/data_repo/GRCh38/chrom_list.txt " + \
            "--kraken-db-path /home/kraken2-db " + \
            "--docker "
    if args.virus is None:
        print("Error: At least one virus should be provided with the --virus argument.")
        exit(1)
    for virus in args.virus:
        command += "--virus {} ".format(virus)
    command += "--input-file /home/input/{} ".format(os.path.basename(args.input_file))
    if args.input_file_2 is not None:
        command += "--input-file-2 /home/input/{} ".format(os.path.basename(args.input_file_2))
    if args.skip_bwa_filter:
        command += " --skip-bwa-filter "
    print(command)
    shell_output = subprocess.check_output(
        command, shell=True)
    print(shell_output)

if __name__ == "__main__":
    args = parse_input_args(docker_run=True)
    check_arguments(args)
    call_fastvifi_pipeline(args)
