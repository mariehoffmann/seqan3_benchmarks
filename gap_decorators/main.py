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

# regex for extracting maximum resident heap size from 'time' command above
size_rx = re.compile("^\s*(\d+)\s+maximum resident set size\s*$")

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

# input: dictionary of type <class 'multiprocessing.managers.DictProxy'>
# 0: binaries
# 1: seq_len_list
# 2: runtime_avg_list
# 3: runtime_stddev_list
# 4: max_resident_set_sizes
# 5: return code
def print_results(result_dict):
    print(type(result_dict))
    d = dict(result_dict)
    print(type(d))
    print(d.items())
    print(iter(d.keys()))
    for id, item in d.items():
        print(id)
        print(item)
        print("RESULT: " + item[0][0])
        header = "\t".join(["seq_len", "runtime avg [ns]", "runtime stddev [ns]", "maximum_resident_set_size"])
        print(header + "\n" + "-"*len(header))
        for seq_len, runtime_avg, runtime_stddev, max_resident_set_size in zip(item[1], item[2], item[3], item[4]):
            print("\t\t".join([str(seq_len), str(runtime_avg), str(runtime_stddev), str(max_resident_set_size)]))

# compile single binary
def compile(cmd):
    os.system(cmd)

# compile in parallel, input [[str_run, strs_heap]], output [[binary_run, binaries_heap]]
def compile_parallel(compile_string_llist):
    print("compile_strings = " + str(compile_string_llist))
    binaries = []
    procs = []
      # flatten
    for i, compile_strings in enumerate(compile_string_llist):
        compile_strings = [compile_strings[0]] + compile_strings[1:]
        print("next cmd: ")
        procs_i = []
        for cmd in compile_strings:
            print(cmd)
            p = mp.Process(target=compile, args=(str(cmd),))
            p.start()
            procs_i.append(p)
            print("proc " + str(p) + " with cmd = " + cmd)
        procs.append(procs_i)
    print(str(procs))
    for i, procs_i in enumerate(procs):
        binaries_i = []
        for id, p in enumerate(procs_i):
            p.join()
            # extract binary name from compilation command
            binary_rx = re.compile("-o {}/(.+)$".format(src_dir))
            search_obj = binary_rx.search(compile_string_llist[i][id])
            if search_obj is None:
                print("STATUS: Error - could not extract binary name")
                sys.exit(-1)
            print("search_obj = " + search_obj.group(1))
            binaries_i.append(search_obj.group(1))
        binaries.append(binaries_i)
    for binary in binaries:
        print(binary)
    return binaries

# runtime and heap profiling input [binary_run, binaries_heap]
def run(binaries, id, return_dict):
    # collect heap profiling output
    # 0: binaries
    # 1: seq_len_list
    # 2: runtime_avg_list
    # 3: runtime_stddev_list
    # 4: max_resident_set_sizes
    # 5: return code
    rc = [[] for _ in range(6)]
    for i, binary in enumerate(binaries):
        rc[0].append(binary)
        path_to_binary = os.path.join(src_dir, binary)
        if os.path.isfile(path_to_binary) is False:
            print("ERROR: binary '" + binary + "' not found.")
            return
        if i == 0:  # case: profile runtimes
            path_to_time_out = os.path.join(result_dir, binary + ".csv")
            cmd = time_cmd.format(str(path_to_binary))
            # high byte: binary exit code, low byte: signal num
            code = os.system(cmd) >> 8
            rc[5] = code
            # TODO: extract avg runtime and stddev
            with open(path_to_time_out, 'r') as f:
                for line in f.readlines()[1:]:
                    line = line.strip().split(',')
                    rc[1].append(int(line[0]))
                    rc[2].append(float(line[1]))
                    rc[3].append(float(line[2]))
            print(rc[0][0] + ": \t" + "\t".join([str(rc[1]), str(rc[2]) + "\u00B1" + str(rc[3]), str(rc[4]), str(rc[5])]))

        else:   # case: profile heap
            # space_cmd = "(/usr/bin/time -l {0} 0) > {1} 2> {2}"
            path_to_space_log = os.path.join(log_dir, binary + ".space.log")
            path_to_space_out = os.path.join(result_dir, binary + ".space.out")
            cmd = space_cmd.format(str(path_to_binary), str(path_to_space_log), str(path_to_space_out))
            print(cmd)
            code = os.system(cmd)
            print(code >> 8)
            # extract 'maximum resident set size'
            with open(path_to_space_out, 'r') as f:
                line = f.readlines()[1]
                mobj = size_rx.match(line)
                if mobj is None:
                    print("ERROR: Could not extract resident size from '" + str(path_to_time_out) + "'")
                    sys.exit(0)
            rc[4].append(int(mobj.group(1)))
            print("max set size: " + str(rc[4][-1]))

    print(rc[0])
    print(rc[1])
    print(rc[2])
    print(rc[3])
    print(rc[4])
    return_dict[id] = rc

# run binaries in parallel
# input: [[binary_run, binaries_heap]]
def run_parallel(binary_llist):
    procs = []
    manager = mp.Manager()
    return_dict = manager.dict()

    print(binary_llist)
    for binary_list in binary_llist:
        p = mp.Process(target=run, args=(binary_list, id, return_dict))
        for binary in binary_list:
            print("STATUS: start " + binary)
        p.start()
        procs.append(p)
        for p in procs:
            p.join()
    return return_dict

# create cpp files from template and compilation strings,
# return format: [[compile_string_run, compile_strings_heap]]
def generate_src_files(idx):
    # create seed for random number generator in benchmarks
    random.seed(datetime.now())
    seed = random.random()
    compile_str_list = []
    for i, benchmark in enumerate(benchmarks[idx-1]):
        for base_type in base_types:
            for GAP_FLAG in range(2):
                benchmark_name_run = benchmark + "_".join(["", base_type, "GAPFLAG", str(GAP_FLAG), "RUN", "REPEAT", str(param.REPEAT), "POW1", str(param.POW1), "POW2", str(param.POW2)])
                # one for each power of two with REPEAT=1
                benchmark_names_heap = [benchmark + "_".join(["", base_type, "GAPFLAG", str(GAP_FLAG), "HEAP", "POW", str(pow)]) for pow in range(param.POW1, param.POW2+1)]
                # TODO: use better run file for copying
                print(benchmark_name_run)
                print(benchmark_names_heap)

                # cpp file names
                cpp_file_run = benchmark_name_run + ".cpp"
                cpp_files_heap = [benchmark_name + ".cpp" for benchmark_name in benchmark_names_heap]

                # build individual compile strings
                compile_str_new = compile_str.replace("<home_dir>", dirs.home)
                compile_str_new = compile_str_new.replace("<seqan3_dir>", dirs.seqan3)
                compile_str_new = compile_str_new.replace("<rangev3_dir>", dirs.rangev3)
                compile_str_new = compile_str_new.replace("<sdsl_dir>", dirs.sdsl)
                compile_str_new = compile_str_new.replace("<benchmark_file>", benchmark)
                compile_str_run = compile_str_new.replace("<benchmark>", benchmark_name_run)
                compile_str_run = compile_str_run.replace("<cpp_file>", cpp_file_run)
                compile_strs_heap = [compile_str_new.replace("<benchmark>", benchmark_name) for benchmark_name in benchmark_names_heap]
                compile_strs_heap = [compile_str_heap.replace("<cpp_file>", cpp_file) for compile_str_heap, cpp_file in zip(compile_strs_heap, cpp_files_heap)]

                print(compile_str_run)
                for c in compile_strs_heap:
                    print(c)
                compile_str_list.append([compile_str_run] + compile_strs_heap)

                template_file = "benchmark{}.template.cpp".format(idx)
                # substitute datatype in template cpp file
                if (os.path.isfile(template_file) == False):
                    print("Error: file '" + template_file + "' not in current directory.")
                    sys.exit(1)

                # copy template
                template_file_new = cpp_file_run
                copyfile(template_file, os.path.join(src_dir, template_file_new))

                # print info string
                sed_cmd = sed.format("info", "auto-generated benchmark", template_file_new)
                print(sed_cmd)
                os.system(sed_cmd)
                # set seed
                sed_cmd = sed.format("seed", seed, template_file_new)
                os.system(sed_cmd)

                # substitute datatype in file
                alphabet_type_rx = re.compile("(\w+\d+)_?.*?")
                alphabet_type = alphabet_type_rx.search(base_type).group(1)
                print(alphabet_type)
                letter_A = "alphabet_type::A"
                if i is 0:  #  ="Gapped Sequence":
                    letter_A = "alphabet_type{{{}::A}}".format(alphabet_type)
                sed_cmd = sed.format("letter_A", letter_A, template_file_new)
                print(sed_cmd)
                os.system(sed_cmd)
                # substitute underlying sequence data type
                sed_cmd = sed.format("alphabet_type", alphabet_type, template_file_new)
                if i is 0:  # = "Gapped Sequence"
                    sed_cmd = sed.format("alphabet_type", "gapped<{}>".format(alphabet_type), template_file_new)
                print(sed_cmd)
                os.system(sed_cmd)
                # substitute value_type
                #case gapped sequence: vector<alphabet_type>, else: vector<gapped<alphabet_type>>
                value_type = alphabet_type
                if i == 0:
                    value_type = "gapped<{}>".format(value_type)
                print("value_type = " + value_type)
                sed_cmd = sed.format("value_type", value_type, template_file_new)
                print(sed_cmd)
                os.system(sed_cmd)

                # substitute container type of underlying sequence
                container_type = "std::vector<{}>"
                if base_type.endswith("_compressed") == True:
                    container_type = "seqan3::bitcompressed_vector<{}>"
                print("container_type = " + container_type.format(value_type))
                sed_cmd = sed.format("container_type", container_type.format(value_type), template_file_new)
                os.system(sed_cmd)

                # substitute gap_decorator
                gap_decorator = gap_decorators[i]
                sed_cmd = sed.format("gap_decorator", gap_decorator, template_file_new)
                os.system(sed_cmd)

                # substitute LOG_LEVEL macro
                sed_cmd = sed.format("LOG_LEVEL", abs(hash(benchmark_name_run)) % 10**8, template_file_new)
                os.system(sed_cmd)
                print(sed_cmd)

                # set GAP_FLAG
                sed_cmd = sed.format("GAP_FLAG", GAP_FLAG, template_file_new)
                os.system(sed_cmd)

                # distribute current copy (note: cpp_file_run is identical with template_file_new up to here)
                for cpp_file in cpp_files_heap:
                    copyfile(os.path.join(src_dir, template_file_new), os.path.join(src_dir, cpp_file))
                    print("copied " + str(cpp_file))

                # substitute output filename
                result_file_run = os.path.join(result_dir, benchmark_name_run + ".csv")
                print("result file name: " + result_file_run)
                sed_cmd = sed.format("benchmark.csv", result_file_run, template_file_new)
                os.system(sed_cmd)
                print(sed_cmd)

                # substitute benchmark name
                sed_cmd = sed.format("benchmark", benchmark_name_run, cpp_file_run)
                os.system(sed_cmd)
                print(sed_cmd)
                for cpp_file, benchmark_name in zip(cpp_files_heap, benchmark_names_heap):
                    sed_cmd = sed.format("benchmark", benchmark_name, cpp_file)
                    os.system(sed_cmd)
                    print(sed_cmd)

                # substitute sequence length range as power of 2 and number of experiment repetitions
                sed_cmd = sed.format("POW1", param.POW1, cpp_file_run)
                os.system(sed_cmd)
                sed_cmd = sed.format("POW2", param.POW2, cpp_file_run)
                os.system(sed_cmd)

                # note: for heap runs pow1 = pow2
                sed_cmds = [sed.format("POW1", p, cpp_file) for p, cpp_file in zip(range(param.POW1, param.POW2 + 1), cpp_files_heap)]
                [os.system(sed_cmd) for sed_cmd in sed_cmds]
                sed_cmds = [sed.format("POW2", p, cpp_file) for p, cpp_file in zip(range(param.POW1, param.POW2 + 1), cpp_files_heap)]
                [os.system(sed_cmd) for sed_cmd in sed_cmds]

                # set REPEAT=1 for heap runs, REPEAT range for runtime version
                sed_cmd = sed.format("REPEAT", param.REPEAT, cpp_file_run)
                os.system(sed_cmd)

                sed_cmds = [sed.format("REPEAT", 1, cpp_file) for cpp_file in cpp_files_heap]
                [os.system(sed_cmd) for sed_cmd in sed_cmds]
            break
        break
    os.system("rm ./src/*.bak")
    return compile_str_list

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
        compile_string_llist = generate_src_files(int(sys.argv[1]))
        for compile_string_list in compile_string_llist:
            print(compile_string_list)
            for f in compile_string_list:
                print(f)
            print("\n")
        print("STATUS: Done\nSTATUS: Compile source files ...")
        binaries = compile_parallel(compile_string_llist)
        print(binaries)
        #sys.exit(0)
        print("STATUS: Done\nSTATUS: Execute binaries in parallel ...")
        result_dict = run_parallel(binaries)
        print_results(result_dict)
        # delete backup files
        for fn in [fn for fn in os.listdir(src_dir) if fn.endswith(".bak")]:
            os.remove(os.path.join(src_dir, fn))
