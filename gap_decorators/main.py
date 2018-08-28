#!/usr/local/bin/python3
'''
    Script for generating 2*2*|methods| c++ source files for each benchmark
    (read, write random or deterministic), method (gap presentation), datatype
    (dna15 or compressed dna4), and sequence state (gapped or ungapped).
    It is recommended to create environmental variables for the home and
    installation directories.
    Exemplary call:

            python main.py 1 $HOME_DIR $SEQAN3_DIR $RANGEV3_DIR $SDSL_DIR 32 2 8

    to produce, run, and evaluate files for benchmark1 with sequence length
    ranges in [2^2, 2^3, ..., 2^8] and 32 experiment repetitions.
'''
from collections import namedtuple
import copy
from datetime import datetime
import multiprocessing as mp
import os
import random
import re
from shutil import copyfile
import sys

# directories for compilation
dirs = namedtuple("dirs", "home seqan3 rangev3 sdsl")
# parameter settings for benchmarks
# REPEAT: number of experiments repetitions
# POW1, POW2: sequence length ranges from [2**POW1, 2**(POW1+1), ..., 2**(POW2)]
param = namedtuple("param", "REPEAT POW1 POW2")

# result collector for single runs
result_collector = namedtuple("result_collector", "binary seq_len runtime_avg runtime_stddev max_resident_set_size code")

# seqan3 data types to test, dna4 = seqan3::dna4, dna4_compressed = seqan3::bitcompressed_vector<dna4>
base_types = ['dna4', 'dna4_compressed']  # 'dna15', 'dna4',

# local directory to place generated cpp files
src_dir = "./src"

# local directory to collect time and heap measurements
result_dir = "./results"

# local directory to collect terminal outputs of single runs
log_dir = "./log"

# compilation string for generated cpp source files
compile_str = "g++ -std=c++17 -Wall -fconcepts -I<seqan3_dir>/include -I<rangev3_dir>/include -I<sdsl_dir>/include {0}/<cpp_file> -o {0}/<benchmark>".format(src_dir)

#: heap profiler with positional arguments 0 for output file name, and 1 for binary name
valgrind_cmd = "valgrind --tool=massif --massif-out-file={0} {1}"

# measure heap consumption (MacOS, FreeBSD), 0: binary, 1: terminal output, 2: time output for heap profiling
space_cmd = "(/usr/bin/time -l {0} 0) > {1} 2> {2}"

# measure runtime by setting csv_flag for binary call, 0: binary, 1: terminal output for logging
time_cmd = "{0} 1"

# benchmark names used for binary and result generation
benchmarks = [["benchmark" + str(i+1) + meth for meth in ["_GS", "_GVsd", "_GVbit", "_AL", "_AS"]] for i in range(3)]

# the different approaches to represent gaps, see gapped_sequence.cpp, etc.
gap_decorators = ['gapped_sequence', 'gap_vector_sd', 'gap_vector_bit', 'anchor_list', 'anchor_set']

# info strings for the different approaches
info = ["Gapped Sequence container<gapped>", "Gap Vector sdsl::sd_vector", "Gap Vector sdsl::bit_vector", "Anchor List", "Anchor Set"]

# sed command for in-place text substitutions
sed = "sed -i.bak -e 's|<\[{}\]>|{}|g' -- ./src/{}"

# compile single binary
def compile(cmd):
    os.system(cmd)

# compile in parallel
def compile_parallel(compile_strings):
    print("compile_strings = " + str(compile_strings))
    binaries = []
    procs = []
    for cmd in compile_strings:
        print("next cmd: ")
        print(cmd)
        p = mp.Process(target=compile, args=(str(cmd),))
        p.start()
        procs.append(p)
        print("proc " + str(p) + " with cmd = " + cmd)
    print(str(procs))
    for id, p in enumerate(procs):
        p.join()
        # extract binary name from compilation command
        binary_rx = re.compile("-o {}/(.+)$".format(src_dir))
        search_obj = binary_rx.search(compile_strings[id])
        if search_obj is None:
            print("STATUS: Error - could not extract binary name")
            sys.exit(-1)
        print("search_obj = " + search_obj.group(1))
        binaries.append(search_obj.group(1))
    return binaries

# run single binary twice, once for runtime performance and once for heap profiling
def run(binary, id, return_dict):
    #cmd = "{}".format(str(os.path.join(src_dir, binary)))  #, param.REPEAT, param.POW1, param.POW2)
    #cmd = valgrind_cmd.format(str(os.path.join(result_dir, binary + ".massif.out")), str(os.path.join(src_dir, binary)))
    path_to_binary = os.path.join(src_dir, binary)
    path_to_space_log = os.path.join(log_dir, binary + ".space.log")
    path_to_space_out = os.path.join(result_dir, binary + ".space.out")
    #path_to_time_log = os.path.join(log_dir, binary + ".time.log")
    path_to_time_out = os.path.join(result_dir, binary + ".csv")
    if os.path.isfile(path_to_binary) is False:
        print("STATUS: binary '" + binary + "' not found.")
        return
    cmd = space_cmd.format(str(path_to_binary), str(path_to_space_log), str(path_to_time_out))
    # A. Profile heap
    code = os.system(cmd)
    # collect heap profiling output
    rc = copy.copy(result_collector)
    rc.binary = binary
    # extract 'maximum resident set size'
    with open(path_to_time_out, 'r') as f:
        line = f.readlines()[1]
        mobj = re.compile("^\s*(\d+)\s+maximum resident set size\s*$").match(line)
        if mobj is None:
            print("STATUS: Error - could not extract resident size from '" + str(path_to_time_out) + "'")
            sys.exit(0)
    rc.max_resident_set_size = int(mobj.group(1))
    print("max set size: " + str(rc.max_resident_set_size))
    # B. Measure runtimes
    cmd = time_cmd.format(str(path_to_binary))
    # high byte: binary exit code, low byte: signal num
    code = os.system(cmd) >> 8
    rc.code = code
    # TODO: extract avg runtime and stddev
    with open(path_to_time_out, 'r') as f:
        line = f.readlines()[1].split(',')
        rc.seq_len = int(line[0])
        rc.runtime_avg = float(line[1])
        rc.runtime_stddev = float(line[2])
    print(rc.binary + ": \t" + "\t".join([str(rc.seq_len), str(rc.runtime_avg) + "\u00B1" + str(rc.runtime_stddev), str(rc.max_resident_set_size), str(rc.code)]))
    sys.exit(0)
    return_dict[id] = rc

# run binaries in parallel
def run_parallel(binaries):
    procs = []
    manager = mp.Manager()
    return_dict = manager.dict()
    print(binaries)
    for id, cmd in enumerate(binaries):
        p = mp.Process(target=run, args=(cmd, id, return_dict))
        print("STATUS: start " + cmd)
        p.start()
        procs.append(p)
    for p in procs:
        p.join()
    return return_dict

# TODO: create one source files per sequence length!
def generate_src_files(idx):
    # create seed for random number generator in benchmarks
    random.seed(datetime.now())
    seed = random.random()
    src_files = []
    for i, benchmark in enumerate(benchmarks[idx-1]):
        for base_type in base_types:
            for GAP_FLAG in range(2):
                benchmark_name = benchmark + "_".join(["", base_type, "GAPFLAG", str(GAP_FLAG), "REPEAT", str(param.REPEAT), "POW1", str(param.POW1), "POW2", str(param.POW2)])
                print(benchmark_name)
                b = compile_str.replace("<home_dir>", dirs.home)
                b = b.replace("<seqan3_dir>", dirs.seqan3)
                b = b.replace("<rangev3_dir>", dirs.rangev3)
                b = b.replace("<sdsl_dir>", dirs.sdsl)
                b = b.replace("<benchmark_file>", benchmark)
                b = b.replace("<benchmark>", benchmark_name)

                template_file = "benchmark{}.template.cpp".format(idx)
                # substitute datatype in template cpp file
                if (os.path.isfile(template_file) == False):
                    print("Error: file '" + template_file + "' not in current directory.")
                    sys.exit(1)
                print("GAP flag = " + str(GAP_FLAG))
                print("benchmark = " + benchmark + "\nbase_type = " + base_type)
                #benchmark_name = benchmark + "_" + alphabet_type + "_" + str(GAP_FLAG)
                # copy template
                cpp_file = benchmark_name + ".cpp"
                copyfile(template_file, os.path.join(src_dir, cpp_file))
                print("cpp_file = " + cpp_file)
                b = b.replace("<cpp_file>", cpp_file)
                print(b)
                # print info string
                sed_cmd = sed.format("info", "auto-generated benchmark", cpp_file)
                print(sed_cmd)
                os.system(sed_cmd)
                # set seed
                sed_cmd = sed.format("seed", seed, cpp_file)
                os.system(sed_cmd)

                # substitute datatype in file
                alphabet_type_rx = re.compile("(\w+\d+)_?.*?")
                alphabet_type = alphabet_type_rx.search(base_type).group(1)
                print(alphabet_type)
                letter_A = "alphabet_type::A"
                if i is 0:  #  ="Gapped Sequence":
                    letter_A = "alphabet_type{{{}::A}}".format(alphabet_type)
                sed_cmd = sed.format("letter_A", letter_A, cpp_file)
                print(sed_cmd)
                os.system(sed_cmd)
                # substitute underlying sequence data type
                sed_cmd = sed.format("alphabet_type", alphabet_type, cpp_file)
                if i is 0:  # = "Gapped Sequence"
                    sed_cmd = sed.format("alphabet_type", "gapped<{}>".format(alphabet_type), cpp_file)
                print(sed_cmd)
                os.system(sed_cmd)
                # substitute value_type
                #case gapped sequence: vector<alphabet_type>, else: vector<gapped<alphabet_type>>
                value_type = alphabet_type
                if i == 0:
                    value_type = "gapped<{}>".format(value_type)
                print("value_type = " + value_type)
                sed_cmd = sed.format("value_type", value_type, cpp_file)
                print(sed_cmd)
                os.system(sed_cmd)

                # substitute container type of underlying sequence
                container_type = "std::vector<{}>"
                if base_type.endswith("_compressed") == True:
                    container_type = "seqan3::bitcompressed_vector<{}>"
                print("container_type = " + container_type.format(value_type))
                sed_cmd = sed.format("container_type", container_type.format(value_type), cpp_file)
                os.system(sed_cmd)

                # substitute gap_decorator
                gap_decorator = gap_decorators[i]
                sed_cmd = sed.format("gap_decorator", gap_decorator, cpp_file)
                os.system(sed_cmd)

                # substitute output filename
                result_file = os.path.join(result_dir, benchmark_name + ".csv")
                print("result file name: " + result_file)
                sed_cmd = sed.format("benchmark.csv", result_file, cpp_file)
                os.system(sed_cmd)
                print(sed_cmd)

                # substitute LOG_LEVEL macro
                sed_cmd = sed.format("LOG_LEVEL", abs(hash(benchmark_name)) % 10**8, cpp_file)
                os.system(sed_cmd)
                print(sed_cmd)

                # substitute benchmark name
                sed_cmd = sed.format("benchmark", benchmark_name, cpp_file)
                os.system(sed_cmd)
                print(sed_cmd)

                # substitute sequence length range as power of 2 and number of experiment repetitions
                sed_cmd = sed.format("POW1", param.POW1, cpp_file)
                os.system(sed_cmd)
                sed_cmd = sed.format("POW2", param.POW2, cpp_file)
                os.system(sed_cmd)
                sed_cmd = sed.format("REPEAT", param.REPEAT, cpp_file)
                os.system(sed_cmd)

                # set GAP_FLAG
                sed_cmd = sed.format("GAP_FLAG", GAP_FLAG, cpp_file)
                os.system(sed_cmd)
                src_files.append(b)
                break
            break
        break
    return src_files

if __name__ == "__main__":
    if len(sys.argv) != 9:
        print("Usage: python main.py benchmark_idx[1:3] home_dir seqan3_dir rangev3_dir sdsl_dir repeat:int pow1:int pow2:int")
    else:
        dirs.home = sys.argv[2]
        dirs.seqan3 = sys.argv[3]
        dirs.rangev3 = sys.argv[4]
        dirs.sdsl = sys.argv[5]
        param.REPEAT = int(sys.argv[6])
        param.POW1 = int(sys.argv[7])
        param.POW2 = int(sys.argv[8])
        print("STATUS: Generate source files ...")
        src_files = generate_src_files(int(sys.argv[1]))
        for src_file in src_files:
            print("\t" + src_file)
        print("STATUS: Done\nSTATUS: Compile source files ...")
        binaries = compile_parallel(src_files)
        print(binaries)
        #sys.exit(0)
        print("STATUS: Done\nSTATUS: Execute binaries in parallel ...")
        return_codes = run_parallel(binaries)
        for id, code in return_codes.items():
            if code != 1:
                print("STATUS: " + binaries[id] + " aborted with return code " + str(code) + ".")
            else:
                print("STATUS: Results for " + binaries[id] + " are written to ./results/" + binaries[id] + ".csv")
        # delete backup files
        for fn in [fn for fn in os.listdir(src_dir) if fn.endswith(".bak")]:
            os.remove(os.path.join(src_dir, fn))
