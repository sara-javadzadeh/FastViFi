import sys

max_num_reads = 100*1000

def represent_nums(num):
    if num < max_num_reads:
        return "%0.2f" % num
    return "infinity"

if __name__ == "__main__":
    #mode, k_1, k_2, f_thr_1, unmapped_thr_1, unmapped_thr_2, precision, recall, \
    #num reads passing first kraken filter, num viral reads passing first kraken filter, and human, and contaminants, \
    #num reads passing second kraken filter, num viral reads passing second kraken filter, and human, and contaminants, \
    #num reads pre filtering, runtime, experiment name
    if len(sys.argv) < 6:
        print("Usage: {} summary_file_name min_recall min_efficiency min_human_efficiency min_contaminant_efficiency".format(sys.argv[0]))
        exit(1)
    filename = sys.argv[1]
    num_human_reads = 100*1000
    num_contaminant_reads = 5250
    min_recall = float(sys.argv[2])
    min_efficiency = float(sys.argv[3])
    min_human_efficiency = int(sys.argv[4])
    min_contaminant_efficiency = int(sys.argv[5])

    for line in open(filename, "r"):
        if line[0] == "#":
            continue
        line = line.strip()
        words = line.split(',')
        recall = float(words[7])
        precision = float(words[6])
        num_reads_past_filter = int(words[12])
        num_reads_pre_filter = int(words[16])
        efficiency = float(num_reads_pre_filter) / (num_reads_past_filter + sys.float_info.epsilon)
        efficiency_human = float(num_human_reads) / (int(words[14]) + sys.float_info.epsilon)
        efficiency_contaminants = float(num_contaminant_reads) / (int(words[15]) + sys.float_info.epsilon)
        k_1 = int(words[1])
        k_2 = int(words[2])

        if recall >= min_recall and efficiency >= min_efficiency and efficiency_human >= min_human_efficiency and efficiency_contaminants >= min_contaminant_efficiency:
            efficiency_human = represent_nums(efficiency_human)
            efficiency_contaminants = represent_nums(efficiency_contaminants)
            print("Recall: {},\tefficiency: {},\thuman_eff: {},\tcontaminant_eff: {},\tprecision {},\texp: {}"\
                    .format("%.4f" % recall, "%.2f" % efficiency, efficiency_human,
                        efficiency_contaminants, "%.2f" % precision, words[-1].split("/")[-1]))
