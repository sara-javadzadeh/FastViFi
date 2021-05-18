import pysam
import os
import sys
import time
from collections import defaultdict

def get_human_chr_list():
    human_chr_list = ["chr" + str(i) for i in range(1, 23)]
    human_chr_list.extend(["chrX", "chrY", "chrM"])
    human_chr_names = {}
    for chrom in human_chr_list:
        human_chr_names[chrom] = True
    return human_chr_names

def is_read_or_pair_unmapped(read):
    if read.is_paired and (read.is_unmapped or read.mate_is_unmapped):
        return True
    return False

def is_potentially_viral(read, human_chr_names):
    # Assume the read is paired. (should be checked before running this function)
    # If one end is unmapped, the read is potentially viral
    if read.is_unmapped or read.mate_is_unmapped:
        return True
    # If both the read and the mate are mapped to human genome, they are both non-viral
    if read.reference_name in human_chr_names and read.next_reference_name in human_chr_names:
            return False
    return True

def write_to_fastq(fq_file, read_name, sequence, quality):
    # quality in unsigned char (not ascii)
    fq_file.write('@{}\n{}\n+\n{}\n'.format(read_name, ''.join(sequence), ''.join("%c" % (char + 33) for char in quality)))

def write_sequence(read, output_fq_file_1, output_fq_file_2):
    # To save space and time, output bam file not created for now.
    #output_bamfile.write(read)
    if read.is_read1:
        write_to_fastq(output_fq_file_1, read.query_name + "/1", read.query_sequence, read.query_qualities)
    elif read.is_read2:
        write_to_fastq(output_fq_file_2, read.query_name + "/2", read.query_sequence, read.query_qualities)
    else:
        print("WARNING: read {} is neither read1 nor read2".format(read.query_name))

def duplicate_single_read(read_1, read_list):
    if len(read_list) == 0:
        return False
    read_2 = read_list[-1]
    return read_1.query_name == read_2.query_name and\
           read_1.query_sequence == read_2.query_sequence

def read_input_file_flag_names(input_bamfile):
    human_chr_names = get_human_chr_list()
    print("human_chr_names: ", human_chr_names)
    bamfile_references = input_bamfile.references
    for reference in bamfile_references:
        if reference not in human_chr_names:
            has_only_human_references = False
            break
    reads_passing_filter = defaultdict(list)
    num_unmapped = 0
    num_unmapped_after_checks = 0
    reading_start_time = time.time()
    secondary_counter = 0
    for read in input_bamfile.fetch(until_eof=True):
        if has_only_human_references and (not read.is_unmapped) and (not read.mate_is_unmapped):
            continue
        elif (not read.is_unmapped) and read.reference_name in human_chr_names and \
            (not read.mate_is_unmapped) and read.next_reference_name in human_chr_names:
            continue
        elif (not read.is_paired) or read.is_secondary or read.is_supplementary or read.is_duplicate:
            secondary_counter += 1
            continue
        #if read.is_unmapped or read.mate_is_unmapped:
            # unmapped reads don't appear next to each other in the bam file
            # keeping them in a separate dict to write in fq files at the end
            # in order to keep both fq files having reads in order.
        reads_passing_filter[read.query_name].append(read)

    print("num secondary, supplementary or duplicate : ", secondary_counter)
    print("num_unmapped: ", num_unmapped)
    print("num_unmapped_after_checks: ", num_unmapped_after_checks)
    print("num reads passing filter: ", len(reads_passing_filter))
    print("time of reading: {}".format(time.time() - reading_start_time))
    return reads_passing_filter

def read_input_file(input_bamfile):
    human_chr_names = get_human_chr_list()
    #print("human_chr_names: ", human_chr_names)
    bamfile_references = input_bamfile.references
    #print("bamfile_references", bamfile_references)
    reads_passing_filter = defaultdict(list)
    reading_start_time = time.time()
    secondary_counter = 0
    for i, read in enumerate(input_bamfile.fetch(until_eof=True)):
        if i % 100000000 == 0:
            print("Processing read # {}, time spent on reading so far: {} seconds".format(i+1, time.time() - reading_start_time))
        if not (read.flag & 15 == 3 and read.reference_name in human_chr_names and
                read.next_reference_name in human_chr_names):
            # Skip this read if the read is paired and mapped in proper pair
            # and both the read and the mate are mapped
            # and both the read and mate are mapped to human reference,
            if (read.flag ^ 1) & 3841 > 0:
                # flip the paired end read flag. Skip this iteration if:
                # not paired, not primary, low quality, duplicate or supplementary
                secondary_counter += 1
            else:
                reads_passing_filter[read.query_name].append(read)

    print("num secondary, supplementary or duplicate : ", secondary_counter)
    print("num reads passing filter: ", len(reads_passing_filter))
    print("time of reading: {}".format(time.time() - reading_start_time))
    return reads_passing_filter

def write_selected_reads(reads_passing_filter):
    writing_start_time = time.time()
    for query_name, read_mates in reads_passing_filter.items():
        #if i % 1000000 == 0:
        #    print("Processing potentially viral read # {}, time spent on writing ".format(i, time.time() - writing_start_time))
        write_sequence(read_mates[0], output_fq_file_1, output_fq_file_2)
        write_sequence(read_mates[1], output_fq_file_1, output_fq_file_2)
        if not ((read_mates[0].is_read1 and read_mates[1].is_read2) or \
            (read_mates[1].is_read1 and read_mates[0].is_read2)):
            print("WARNING: Read mates read_1 and read_2 are not as expected for read {} with flag {}".format(read_mates[0].query_name, read_mates[0].flag))
    print("time of writing: {}".format(time.time() - writing_start_time))


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: {} <path to input bam file> <path to output directory> <prefix>\n".format(sys.argv[0]) +
              "Two Fastq files for paired end reads are created.\n" +
              "Note: human chromosomes reference names are assumed to be in the format "+
              "of chr1-chr22, chrX, chrY, chrM. If that is not the case for the BAM file, " +
              "update the function get_human_chr_list() within the code. Easier way to " +
              "feed in the human chromosome reference names to be provided in the future.")
        exit(1)
    input_filename = sys.argv[1]
    output_dir = sys.argv[2]
    prefix = sys.argv[3]
    output_fq_filename_1 = os.path.join(output_dir, prefix + "_1.fq")
    output_fq_filename_2 = os.path.join(output_dir, prefix + "_2.fq")
    print("Filtering reads from {}".format(input_filename))
    if input_filename.endswith(".bam"):
        input_bamfile = pysam.AlignmentFile(input_filename, "rb")
    elif input_filename.endswith(".cram"):
        input_bamfile = pysam.AlignmentFile(input_filename, "rc")
    else:
        print("ERROR: Cannot open file {}\nPlease provide a BAM or CRAM file as input.".format(input_filename))
        exit(1)
    # To save space and time, output bam file not created for now.
    output_fq_file_1 = open(output_fq_filename_1, "w")
    output_fq_file_2 = open(output_fq_filename_2, "w")

    reads_passing_filter = read_input_file(input_bamfile)
    write_selected_reads(reads_passing_filter)
