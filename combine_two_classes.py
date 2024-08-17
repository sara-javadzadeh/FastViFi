import os
import argparse
import subprocess
import pysam

class Cluster:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.reads_class_1 = []
        self.reads_class_2 = []

    def add_read_class_1(self, name, chrom, position, umi):
        self.reads_class_1.append((name, chrom, int(position), umi))

    def add_read_class_2(self, name, chrom, position, umi):
        self.reads_class_2.append((name, chrom, int(position), umi))

    def get_supporting_reads(self):
        return len(self.reads_class_1) + len(self.reads_class_2)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Combine class 1 and class 2 hybrid reads. " + \
        "Class 1 are single reads where the read is both aligned to the human and viral genomes " + \
        "which is captured based on primary and supplementary alignments. " + \
        "Class 2 is if one mate of the paired end read is mapped to human genome and the other is mapped to the viral genome.")
    parser.add_argument("--class-1-read-ids-1",
                        type=str, required=True,
                        help="The output file generated after running find_class_1_reads.py. The filename is class_1_read_ids_1.txt. " + \
                        "These are single hybrid reads that are first in pair.")
    parser.add_argument("--class-1-read-ids-2",
                        type=str, required=True,
                        help="The output file generated after running find_class_1_reads.py. The filename is class_1_read_ids_2.txt. " + \
                        "These are single hybrid reads that are second in pair.")
    parser.add_argument("--fastvifi-cluster",
                        type=str, required=True,
                        help="The output file generated after running the FastViFi pipeline. The filename is output_<virus>.clusters.txt"
                        )
    parser.add_argument("--fastvifi-bam",
                        type=str, required=True,
                        help="The output bam file generated after running the FastViFi pipeline. The filename is output_<virus>.fixed.trans.bam"
                        )
    parser.add_argument("--input-bam",
                        type=str, required=True,
                        help="The input bam file where find_class_1_reads.py was called on. " + \
                        "This is used to get the full read information from the output of find_class_1_reads.py."
                        )
    parser.add_argument("--output-dir",
                         type=str, required=True,
                         help="directory where all the output files are written into.")
    args = parser.parse_args()
    return args

def combine_bam_files(class_1_filename, class_2_filename, output_dir):
    combined_filename = os.path.join(output_dir, "reads_supporting_junction.sam")
    command = "samtools merge -f -O BAM {combined} {input_1} {input_2}".format(
                    combined=combined_filename,
                    input_1=class_1_filename,
                    input_2=class_2_filename)
    shell_output = subprocess.check_output(command, shell=True)
    logfile.write(b"Combining reads in two classes into a single bam file\n")
    logfile.write(shell_output)
    return combined_filename

def extract_class_1_reads(input_bam, class_1_read_ids_1, class_1_read_ids_2, logfile, output_dir):
    #samtools view ../antoine_hiv/output_hiv.trans.bam | gawk '{if (and($2, 64) > 0) print $0}'| grep -f hybrid_single_read_ids_1.txt
    class_1_reads_1 = os.path.join(output_dir, "class_1_reads_first_in_pair.sam")
    class_1_reads_2 = os.path.join(output_dir,"class_1_reads_second_in_pair.sam")
    class_1_reads_all = os.path.join(output_dir,"class_1_reads.bam")

    # Extract the reads from read ids
    command = "samtools view -H {input_file} > {reads}; samtools view {input_file} | ".format(
                    input_file=input_bam,
                    reads=class_1_reads_1,
                        ) + \
                "awk -e '{if (and($2, 64) > 0) print $0}' | " + \
                "grep -f {ids} >> {reads}".format(
                    ids=class_1_read_ids_1,
                    reads=class_1_reads_1)
    shell_output = subprocess.check_output(command, shell=True)
    logfile.write(b"Extracting class 1 reads (first in pair) from the input bam\n")
    logfile.write(shell_output)

    # Repeat for class 1 reads that are second in pair
    command = "samtools view -H {input_file} > {reads}; samtools view {input_file} | ".format(
                    input_file=input_bam,
                    reads=class_1_reads_2,
                        ) + \
                "awk -e '{if (and($2, 64) > 0) print $0}' | " + \
                "grep -f {ids} >> {reads}".format(
                    input_file=input_bam,
                    ids=class_1_read_ids_2,
                    reads=class_1_reads_2)
    shell_output = subprocess.check_output(command, shell=True)
    logfile.write(b"Extracting class 1 reads (second in pair) from the input bam\n")
    logfile.write(shell_output)

    # Combine the two files into one
    command = "samtools merge -f -O BAM {combined} {input_1} {input_2}".format(
                    combined=class_1_reads_all,
                    input_1=class_1_reads_1,
                    input_2=class_1_reads_2)
    shell_output = subprocess.check_output(command, shell=True)
    logfile.write(b"Combining two class 1 sam files into one\n")
    logfile.write(shell_output)

    return class_1_reads_all

def get_reads(filename):
    aln_file = pysam.AlignmentFile(filename)
    reads = []
    id_to_read = {}
    for read in aln_file.fetch(until_eof=True):
        reads.append(read)
        if read.query_name in id_to_read:
            id_to_read[read.query_name].append(read)
        else:
            id_to_read[read.query_name] = []
    return reads, id_to_read

def get_umi(read_id, id_to_reads, verbose=False):
    #TODO: fix the issue with UMI
    reads = id_to_reads[read_id]
    if verbose:
        print("get tags")
        print(reads[0].get_tags(with_value_type=True))
    #print("get tag")
    #print(reads[0].get_tag("UB", with_value_type=True))
    # CB:Z: + UB:Z:
    return "-"

def extract_class_2_clusters(cluster_filename, id_to_reads):
    cluster_file = open(cluster_filename, "r")
    clusters = []
    for line in cluster_file.readlines():
        # Skip header lines
        if line.startswith("#chr") or line.startswith("##=="):
            continue
        # Cluster info
        if not line.startswith("##"):
            chrom, start, end, total, forward, reverse = line.strip().split()
            # Only store the chromosome, start and end, as a tuple
            current_cluster = Cluster(chrom=chrom, start=start, end=end)
            clusters.append(current_cluster)
        # Read info within each cluster
        if line.startswith("##"):
            read_id, chrom, position, reverse, read_1 = line.strip().split()
            # First two characters are ##
            read_id = read_id[2:]
            umi = get_umi(read_id, id_to_reads)
            current_cluster.add_read_class_2(name=read_id,
                                            chrom=chrom,
                                            position=position,
                                            umi=umi)
            print("read class 2 added")
    print("len read class 2 ", len(clusters[0].reads_class_2))

    return clusters

def find_overlaps(class_1_reads, id_to_reads, clusters):
    unmerged_class_1_reads = []
    for read in class_1_reads:
        added_to_cluster = False
        for cluster in clusters:
            if read.reference_name == cluster.chrom and \
                read.query_alignment_start >= cluster.start and \
                read.query_alignment_start <= cluster.end:
                umi = get_umi(read.query_name, id_to_reads)
                cluster.add_read_class_1(name=read.query_name,
                                         chrom=read.reference_name,
                                         position=read.query_alignment_start,
                                         umi=umi)
                added_to_cluster = True
                break
        if not added_to_cluster:
            unmerged_class_1_reads.append(read)

    return clusters, unmerged_class_1_reads

def get_read_string_in_cluster(read, read_class):
    name, chrom, position, umi = read
    return "{rclass}\t{rid}\t{chrom}\t{position}\t{umi}\n".format(
            rclass=read_class,
            rid=name,
            chrom=chrom,
            position=position,
            umi=umi)

def get_read_string(read, id_to_reads, read_class):
    return "{rclass}\t{rid}\t{chrom}\t{position}\t{umi}\n".format(
            rclass=read_class,
            rid=read.query_name,
            chrom=read.reference_name,
            position=read.query_alignment_start,
            umi=get_umi(read.query_name, id_to_reads))

def write_results(clusters, unmerged_class_1_reads, id_to_reads, output_dir):
    clusters = sorted(clusters, key=lambda cluster: cluster.get_supporting_reads(), reverse=True)
    out_file = open(os.path.join(output_dir, "clusters_summary.txt"), "w+")
    out_file.write("### Clusters where both/either class 1 and class 2 reads are present\n")
    out_file.write("#chrom\tstart\tend\tnum_class_1_reads\tnum_class_2_reads\n")
    # First write clusters where both class 1 and class 2 are present
    for cluster in clusters:
        out_file.write("{chrom}\t{start}\t{end}\t{num_c_1}\t{num_c_2}\n".format(
                        chrom=cluster.chrom,
                        start=cluster.start,
                        end=cluster.end,
                        num_c_1=len(cluster.reads_class_1),
                        num_c_2=len(cluster.reads_class_2),
                        ))
        out_file.write("#read_class\tread_id\tchr\tposition\tUMI\n")
        for read in cluster.reads_class_1:
            out_file.write(get_read_string(read, "class_1"))
        for read in cluster.reads_class_2:
            out_file.write(get_read_string_in_cluster(read, "class_2"))

    out_file.write("### Remaining class 1 reads\n")
    for read in unmerged_class_1_reads:
        out_file.write(get_read_string(read, id_to_reads, "class_1"))

    out_file.close()

if __name__ == "__main__":
    args = parse_arguments()

    # Create log file
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    log_filename = os.path.join(args.output_dir, "log_combine_results.log")
    logfile = open(log_filename, "wb+")

    # Extract a file with all read info from class 1
    class_1_filename = extract_class_1_reads(input_bam=args.input_bam,
                                          class_1_read_ids_1=args.class_1_read_ids_1,
                                          class_1_read_ids_2=args.class_1_read_ids_2,
                                          logfile=logfile,
                                          output_dir=args.output_dir)
    # Extract actual reads
    class_1_reads, id_to_read_1 = get_reads(class_1_filename)
    class_2_reads, id_to_read_2 = get_reads(args.fastvifi_bam)
    all_reads = class_1_reads + class_2_reads
    id_to_reads = id_to_read_1.copy()
    id_to_reads.update(id_to_read_2)

    # Read clusters from file
    clusters = extract_class_2_clusters(args.fastvifi_cluster, id_to_read_2)

    # Create a single bam file with all class_1 or class_2 reads.
    # This way, we have all supporting reads in a single bam file and we don't need to
    # go back to the input bam file for more information.
    combined_bam = combine_bam_files(
            class_1_filename=class_1_filename,
            class_2_filename=args.fastvifi_bam,
            output_dir=args.output_dir)

    clusters, unmerged_class_1_reads = find_overlaps(
            id_to_reads=id_to_read_1,
            class_1_reads=class_1_reads,
            clusters=clusters)
    #TODO: add second alignment and cigar for class 1 reads
    write_results(clusters=clusters,
                unmerged_class_1_reads=unmerged_class_1_reads,
                id_to_reads=id_to_reads,
                output_dir=args.output_dir)
