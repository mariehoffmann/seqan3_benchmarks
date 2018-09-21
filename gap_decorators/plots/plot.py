import csv
import matplotlib.pyplot as plt
import sys

# format [benchmark_idx: %u, datastructure: %s, alphabet:%s, gapflag:%u, seqlen: %u, runtime: %f, runtime_stddev: %f, heap_size: %u]
BM = 0  # position of benchmark index
DS = 1  # position of datastructure
AL = 2  # position of alphabet type
GF = 3  # position of gap flag
SL = 4  # position of sequence length
RTavg = 5   # position of average runtime
RTstddev = 6    # position of runtime deviation
HS = 7  # position of heap size

result_files = ['../results/benchmark_1_1_2_4.csv']

def plot(benchmark_idx):
    R = []  # Result matrix
    with open(result_files[benchmark_idx-1], 'r') as f:
        csvreader = csv.reader(f, delimiter=',')
        next(csvreader)
        for row in csvreader:
            R.append([int(row[0]), row[1], row[2], int(row[3]), int(row[4]), float(row[5]), float(row[6]), int(row[7])])
    for result in R:
        print(result)

    # number of different data structures
    data_structures = sorted(list(set([row[DS] for row in R])))
    num_DS = len(data_structures)
    # number of unique alphabet types
    alphabet_types = sorted(list(set([row[AL] for row in R])))
    num_AL = len(alphabet_types)
    # gap_flags tested
    gap_flags = set([row[GF] for row in R])
    # sequence lengths
    seq_lens = sorted(list(set([row[SL] for row in R])))

    for alphabet_type in alphabet_types[:1]:
        plot_idx = 1
        for gap_flag in gap_flags:
            # row1: gap_flag=0, row2: gap_flag=1, col1: runtimes, col2: heap size
            for data_structure in data_structures:
                print("alphabet_type = {}, gap_flag = {}, data_structure = {}".format(alphabet_type, gap_flag, data_structure))
                R_sub =  [row for row in R if row[GF] == gap_flag and row[AL] == alphabet_type and row[DS] == data_structure]
                data = sorted([(row[SL], row[RTavg], row[HS]) for row in R_sub], key=lambda t: t[0])

                if (len(seq_lens) != len(data)):
                    print("ERROR: missing data for alphabet_type = {}, gap_flag = {}, data_structure = {}".\
                    format(alphabet_type, gap_flag, data_structure))
                    sys.exit(-1)

                plt.subplot(len(gap_flags), 2, plot_idx)
                print(seq_lens)
                print([t[1] for t in data])
                plt.plot(seq_lens, [t[1] for t in data])
                plt.title('Runtimes for ALPHABET = {}, GAP_FLAG = {}'.format(alphabet_type, gap_flag))
                plt.xlabel('sequence lengths')
                plt.ylabel('Runtimes in nanoseconds')

                plt.subplot(len(gap_flags), 2, plot_idx+1)
                plt.plot(seq_lens, [t[2] for t in data])
                plt.title('Resident set size for ALPHABET = {}, GAP_FLAG = {}'.format(alphabet_type, gap_flag))
                plt.xlabel('sequence lengths')
                plt.ylabel('Maximum Resident Set Sizes')

            plot_idx += 2

if __name__ == "__main__":
    if len(sys.argv) == 1:
        plot(1)
    else:
        plot(int(sys.argv[1]))
