import argparse
import pysam
import os
import sys
import time
from collections import defaultdict
import argparse

class HybridRead:
    def __init__(self, query_name, human_chr_names):
        self.first_mate_ref = set()
        self.second_mate_ref = set()
        self.query_name = query_name
        self.human_chr_names = human_chr_names

    def add_ref(self, ref, flag):
        if ref in self.human_chr_names:
            ref_code = "human"
        else:
            ref_code = "viral"
        if flag & 64 == 64:  # First in pair
            self.first_mate_ref.add(ref_code)
        else: # Otherwise, it's 0 and it's second in pair
            self.second_mate_ref.add(ref_code)

    def is_hybrid(self):
        if "human" in self.first_mate_ref and \
           "viral" in self.first_mate_ref:
            return 1
        if "human" in self.second_mate_ref and \
           "viral" in self.second_mate_ref:
            return 2
        return 0

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="This script extracts class 1 reads. " + \
        "Class 1 are single hybrid reads where the read is both aligned to the human and viral genomes " + \
        "which is captured based on primary and supplementary alignments. " + \
          "Two text files are generated with class 1 single reads, depending on the read being the first or second in the pair.\n")
    parser.add_argument("--input",
                        type=str, required=True,
                        help="input bam/sam file where class 1 is to be extracted")
    parser.add_argument("--output-dir",
                        type=str, required=True,
                        help="output directory where output files will be written into.")
    parser.add_argument("--human-chr-list",
                        type=str, required=True,
                        help="a file with each line representing human chromosome list. " + \
                        "The list is used to identify human and non-human (i.e. viral) alignment.")
    args = parser.parse_args()
    return args

def get_human_chr_list(filename):
    human_chr_names = {}
    chrnames_file = open(filename, "r")
    print("reading chr names")
    for line in chrnames_file.readlines():
        human_chr_names[line.strip().replace(">", "")] = True
    chrnames_file.close()
    # deprecated chr names
    #    human_chr_list = ["chr" + str(i) for i in range(1, 23)]
    #    human_chr_list.extend(["chrX", "chrY", "chrM"])
    #    for chrom in human_chr_list:
    #        human_chr_names[chrom] = True
    return human_chr_names

def read_input_file(input_bamfile, human_chr_names_filename):
    human_chr_names = get_human_chr_list(human_chr_names_filename)
    bamfile_references = input_bamfile.references
    # Separating read_1 from read_2 because we want to treat each single read separately.
    hybrid_single_reads_1 = []
    hybrid_single_reads_2 = []
    candidate = {}
    reading_start_time = time.time()
    for i, read in enumerate(input_bamfile.fetch(until_eof=True)):
        if i % 100000000 == 0:
            print("Processing read # {}, time spent on reading {} seconds".format(
                    i+1, time.time() - reading_start_time))
        if read.query_name not in candidate:
            candidate[read.query_name] = HybridRead(read.query_name, human_chr_names)
        candidate[read.query_name].add_ref(
                ref=read.reference_name,
                flag=read.flag)
        if candidate[read.query_name].is_hybrid() == 1:
            hybrid_single_reads_1.append(read.query_name)
        elif candidate[read.query_name].is_hybrid() == 2:
            hybrid_single_reads_2.append(read.query_name)
        else:
            # Not a hybrid
            pass

    # Remove duplicates
    hybrid_single_reads_1 = list(set(hybrid_single_reads_1))
    hybrid_single_reads_2 = list(set(hybrid_single_reads_2))

    print("num hybrid single reads (first mate): ", len(hybrid_single_reads_1))
    print("num hybrid single reads (second mate): ", len(hybrid_single_reads_2))
    print("time of reading: {}".format(time.time() - reading_start_time))
    return hybrid_single_reads_1, hybrid_single_reads_2


if __name__ == "__main__":
    args = parse_arguments()
    input_filename = args.input
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    human_chr_names_filename = args.human_chr_list

    print("Finding hybrid single reads reads from {}".format(input_filename))
    if input_filename.endswith(".bam"):
        input_bamfile = pysam.AlignmentFile(input_filename, "rb")
    elif input_filename.endswith(".cram"):
        input_bamfile = pysam.AlignmentFile(input_filename, "rc")
    else:
        print("ERROR: Cannot open file {}\nPlease provide a BAM or CRAM file as input.".format(input_filename))
        exit(1)
    # To save space and time, output bam file not created for now.
    #output_fq_file_1 = open(output_fq_filename_1, "w")
    #output_fq_file_2 = open(output_fq_filename_2, "w")

    hybrid_single_reads_1, hybrid_single_reads_2 = read_input_file(input_bamfile, human_chr_names_filename)

    output_filename_1 = os.path.join(output_dir, "class_1_read_ids_1.txt")
    output_filename_2 = os.path.join(output_dir, "class_1_read_ids_2.txt")

    # Write single hybrid reads where the hybrid read was first mate
    with open(output_filename_1, "w+") as output_file:
        for read_id in hybrid_single_reads_1:
            output_file.write(read_id + "\n")

    # Write single hybrid reads where the hybrid read was second mate
    with open(output_filename_2, "w+") as output_file:
        for read_id in hybrid_single_reads_2:
            output_file.write(read_id + "\n")
