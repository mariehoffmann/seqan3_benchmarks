#!/usr/local/bin/python3
'''
    Script for generating 2*2*|methods| c++ source files for each benchmark
    (read, write random or deterministic), method (gap presentation), datatype
    (dna15 or compressed dna4), and sequence state (gapped or ungapped).
    It is recommended to create environmental variables for the home and
    installation directories.
    Exemplary call:

            python3 main.py 1 $HOME_DIR $SEQAN3_DIR $RANGEV3_DIR $SDSL_DIR 10 3 12

    to produce, run, and evaluate files for benchmark1 with sequence length
    ranges in [2^2, 2^3, ..., 2^8] and 32 experiment repetitions.
'''

from collections import namedtuple
import copy
from datetime import datetime
import math
import multiprocessing as mp
import os
import random
import re
from shutil import copyfile
import sys
import timeit

# binaries produce their own csv files with runtime results
CSV_FLAG = 0

# directories for compilation
dirs = namedtuple("dirs", "home seqan3 rangev3 sdsl")
# parameter settings for benchmarks
# REPEAT: number of experiments repetitions
# POW1, POW2: sequence length ranges from [2**POW1, 2**(POW1+1), ..., 2**(POW2)]
param = namedtuple("param", "REPEAT POW1 POW2")

# seqan3 data types to test, dna4 = seqan3::dna4, dna4_compressed = seqan3::bitcompressed_vector<dna4>
base_types = ['dna4']  #, 'dna4compressed', 'dna15']  # 'dna15', 'dna4',

# local directory to place generated cpp files
src_dir = "./src"
if os.path.isdir(src_dir) is False:
    os.mkdir(src_dir)

# local directory to collect time and heap measurements
result_dir = "./results"
#result_dir = "./results_AB"
if os.path.isdir(result_dir) is False:
    os.mkdir(result_dir)

# local directory to collect terminal outputs of single runs
log_dir = "./log"
if os.path.isdir(log_dir) is False:
    os.mkdir(log_dir)

# compilation string for generated cpp source files
compile_str = "g++-8 -std=c++17 -Wall -fconcepts -I<seqan3_dir>/include -I<rangev3_dir>/include -I<sdsl_dir>/include {0}/<cpp_file> -o {0}/<benchmark>".format(src_dir)
if os.uname()[0] == 'Darwin':
    compile_str = "g++ -std=c++17 -Wall -fconcepts -I<seqan3_dir>/include -I<rangev3_dir>/include -I<sdsl_dir>/include {0}/<cpp_file> -o {0}/<benchmark>".format(src_dir)

#: heap profiler with positional arguments 0 for output file name, and 1 for binary name
valgrind_cmd = "valgrind --tool=massif --massif-out-file={0} {1}"

# measure heap consumption (MacOS, FreeBSD), 0: binary, 1: terminal output, 2: time output for heap profiling
space_cmd = "(time {0} 0) > {1} 2> {2}"
if os.uname()[0] == 'Darwin':
    space_cmd = "(/usr/bin/time -l {0} 0) > {1} 2> {2}"

# regex for extracting maximum resident heap size from 'time' command above
size_rx = re.compile(".*?\s*(\d+)\s+maximum resident set size.*?")

# measure runtime by setting csv_flag for binary call, 0: binary, 1: terminal output for logging
time_cmd = "{0} 1"

# benchmark names used for binary and result generation
benchmarks = [["benchmark" + str(i+1) + meth for meth in ["_GS", "_GVsd", "_GVbit", "_AL", "_AS2", "_AB"]] for i in range(3)]
#benchmarks = [["benchmark" + str(i+1) + meth for meth in ["_GS", "_AL", "_AS2", "_AB"]] for i in range(3)]

# the different approaches to represent gaps, see gapped_sequence.cpp, etc.
gap_decorators = ['gapped_sequence', 'gap_vector_sd', 'gap_vector_bit', 'anchor_list', 'anchor_set2', 'anchor_blocks']
#gap_decorators = ['gapped_sequence', 'anchor_list', 'anchor_set2', 'anchor_blocks']
#gap_decorators = ['anchor_blocks_8', 'anchor_blocks_16', 'anchor_blocks_32', 'anchor_blocks_64', 'anchor_blocks_128', 'anchor_blocks_256']

# black list
ignore_binary_rx = re.compile('XX')  #'.+?GV(sd|bit)_dna.+?_GAPFLAG_\d_(RUN|HEAP).+?')

# info strings for the different approaches
info = ["Gapped Sequence container<gapped>", "Gap Vector sdsl::sd_vector", "Gap Vector sdsl::bit_vector", "Anchor List", "Anchor Set2", "Anchor Block"]
#info = ["Gapped Sequence container<gapped>", "Anchor List", "Anchor Set2", "Anchor Block"]

# sed command for in-place text substitutions
sed = "sed -i.bak -e 's|<\[{}\]>|{}|g' -- ./src/{}"

# input: dictionary of type <class 'multiprocessing.managers.DictProxy'>
# 0: binaries
# 1: seq_len_list
# 2: runtime_avg_list
# 3: runtime_stddev_list
# 4: max_resident_set_sizes
# 5: return code
def print_results(result_dict, outfile="out.csv"):
    d = dict(result_dict)
    print(d.items())
    benchmark_rx = re.compile("^benchmark(\d+)_(\w+)_([\w\d]+)_GAPFLAG_([0,1]).+?$")
    with open(os.path.join(result_dir, outfile), 'w') as f:
        f.write(",".join(["benchmark_id", "gap_decorator", "alphabet_type", "GAPFLAG", "seq_len", "runtime avg [ns]", "runtime stddev [ns]", "maximum_resident_set_size"]) + "\n")
        for id, item in d.items():
            print("RESULT: " + item[0][0])
            max_rows = max([len(item[i]) for i in range(1,5)])
            if max_rows != min ([len(item[i]) for i in range(1,5)]):
                print("WARNING: Results are missing ")
            header = "\t".join(["seq_len", "runtime avg [ns]", "runtime stddev [ns]", "maximum_resident_set_size"])
            print(header + "\n" + "-"*len(header))
            mobj = benchmark_rx.match(item[0][0])
            benchmark_id = mobj.group(1)
            gap_decorator = mobj.group(2)
            alphabet_type = mobj.group(3)
            GAPFLAG = mobj.group(4)

            for seq_len, runtime_avg, runtime_stddev, max_resident_set_size in \
            zip(item[1] + [-1 for _ in range(max_rows-len(item[1]))], \
            item[2] + [-1 for _ in range(max_rows-len(item[2]))], \
            item[3] + [-1 for _ in range(max_rows-len(item[3]))], \
            item[4] + [-1 for _ in range(max_rows-len(item[4]))]):
                print("\t\t".join([str(seq_len), str(runtime_avg), str(runtime_stddev), str(max_resident_set_size)]))
                f.write(",".join([benchmark_id, gap_decorator, alphabet_type, GAPFLAG, str(seq_len), str(runtime_avg), str(runtime_stddev), str(max_resident_set_size)]) + "\n")
            print("")
            #f.write("\n")
        print("STATUS: output written to " + str(os.path.join(result_dir, outfile)))

# delete intermediate result and backup files
def cleanup():
    for fn in [fn for fn in os.listdir(src_dir) if fn.endswith(".bak")]:
        os.remove(os.path.join(src_dir, fn))
    for fn in [fn for fn in os.listdir(result_dir) if fn.endswith(".space.out")]:
        os.remove(os.path.join(result_dir, fn))

# compile single binary
def compile(cmd):
    code = os.system(cmd) >> 8
    if code != 0:
        print("ERROR: compilation failed with code " + str(code) + " " + cmd)
        sys.exit()

# compile in parallel, input [[str_run, strs_heap]], output [[binary_run, binaries_heap]]
def compile_parallel(compile_string_llist, continue_flag=False):
    #print("compile_strings = " + str(compile_string_llist))
    binaries = []
    binary_rx = re.compile("-o {}/(.+)$".format(src_dir))
    num_workers = mp.cpu_count() - 2

    for llist in compile_string_llist:
        for l in llist:
            print(l[-50:])
        print("---------------------------")
    for i, compile_strings in enumerate(compile_string_llist):
        compile_strings = [compile_strings[0]] + compile_strings[1:]
        for cs in compile_strings:
            print(cs)
        procs_i = set()
        # launch warp-wise, otherwise OSError in popen_fork.py: Too many open files
        k = 0
        # extract binary name from compilation command
        binaries_i = []
        for compile_string in compile_strings:
            search_obj = binary_rx.search(compile_string)
            if search_obj is None:
                print("STATUS: Error - could not extract binary name")
                sys.exit(-1)
            print("DEBUG: append new binary to sublist: {}".format(search_obj.group(1)))
            binaries_i.append(search_obj.group(1))
            if continue_flag == True and os.path.isfile(os.path.join(src_dir, binaries_i[-1])) is False:
                print("DEBUG: binary not found, compile: '{}'".format(binaries_i[-1]))
                os.system(compile_string)
        for binary in binaries_i:
            print(binary)

        if continue_flag == False:
            for j in range(max(1, math.ceil(len(compile_strings)/num_workers))):
                for cmd in compile_strings[j*num_workers : (j+1)*num_workers]:
                    print("DEBUG: launch cmd = " + cmd[-50:])
                    p = mp.Process(target=compile, args=(str(cmd),))
                    p.start()
                    procs_i.add(p)
                    k += 1
                for id, p in enumerate(procs_i):
                    p.join()
                procs_i.clear()

        print("DEBUG: compiled " + str(k) + " binaries")
        print("DEBUG: binaries returned after compilation: ")
        binaries.append(binaries_i)
        for b in binaries_i:
            print(b)
        #sys.exit()
        for binary_list in binaries:
            for binary in binary_list:
                print(binary)
            print("---------------------------")
        #sys.exit()
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
    #print("DEBUG: run with {} binaries".format(len(binaries)))
    for i, binary in enumerate(binaries):
        rc[0].append(binary)
        path_to_binary = os.path.join(src_dir, binary)
        if os.path.isfile(path_to_binary) is False:
            print("ERROR: binary '" + binary + "' not found.")
            return

        if ignore_binary_rx.match(binary) is not None:
            print("STATUS: ignore {}".format(binary))
            rc[1].append(-1)
            rc[2].append(-1)
            rc[3].append(-1)
            rc[4].append(-1)
            rc[5] = -1

        else:
            if i == 0:  # case: profile runtimes
                path_to_time_out = os.path.join(result_dir, binary + ".csv")
                # launch in case no runtime result file exists
                if continue_flag == False or os.path.isfile(path_to_time_out) == False and continue_flag == True:
                    print("DEBUG: run2 '{}'".format(binary))
                    cmd = time_cmd.format(str(path_to_binary))
                    # high byte: binary exit code, low byte: signal num
                    code = os.system(cmd) >> 8
                    rc[5] = code
                with open(path_to_time_out, 'r') as f:
                    for line in f.readlines()[1:]:
                        line = line.strip().split(',')
                        rc[1].append(int(line[0]))
                        rc[2].append(float(line[1]))
                        rc[3].append(float(line[2]))
                print(rc[0][0] + ": \t" + "\t".join([str(rc[1]), str(rc[2]) + "\u00B1" + str(rc[3]), str(rc[4]), str(rc[5])]))

            else:   # case: profile heap
                # average over x runs
                runs = 2
                path_to_space_log = os.path.join(log_dir, binary + ".space.log")
                path_to_space_out = os.path.join(result_dir, binary + ".space.out")
                cmd = space_cmd.format(str(path_to_binary), str(path_to_space_log), str(path_to_space_out))
                sum_max_rss = 0
                path_to_space_out_isvalid = True if os.path.isfile(path_to_space_out) == True else False
                if path_to_space_out_isvalid == True:
                    with open(path_to_space_out, 'r') as f:
                        line = ''.join(f.readlines()[1:3])
                        mobj = size_rx.match(line)
                        if mobj is None:
                            path_to_space_out_isvalid = False

                for _ in range(runs):
                    print("{}: {}".format(path_to_space_out, os.path.isfile(path_to_space_out)))
                    if continue_flag == False or (continue_flag == True and path_to_space_out_isvalid == False):
                        print("STATUS: {} ...".format(cmd))
                        code = os.system(cmd)

                    # extract 'maximum resident set size'
                    with open(path_to_space_out, 'r') as f:
                        line = ''.join(f.readlines()[1:3])
                        mobj = size_rx.match(line)
                        if mobj is None:
                            print("ERROR: Could not extract resident size from '" + str(path_to_space_out) + "', re-run")
                            os.system("cat {}".format(path_to_space_out))
                            #sys.exit(0)
                            sum_max_rss = -1
                        else:
                            sum_max_rss += int(mobj.group(1))
                rc[4].append(int(sum_max_rss/runs))
    print("DEBUG: add result_collector with id = " + str(id))
    return_dict[id] = rc

# run binaries in parallel
# input: [[binary_run, binaries_heap]]
def run_parallel(binary_llist, continue_flag=False):
    manager = mp.Manager()
    return_dict = manager.dict()
    num_workers = mp.cpu_count() - 2

    for i, binary_list in enumerate(binary_llist):
        print("STATUS: launch swap " + str(i))
        p = mp.Process(target=run, args=(binary_list, i, return_dict))
        p.start()
        p.join()
        print("STATUS: swap " + str(i) + " done")
    return return_dict

# create cpp files from template and compilation strings,
# return format: [[compile_string_run, compile_strings_heap]]
def generate_src_files(idx, continue_flag=False):
    # create seed for random number generator in benchmarks
    random.seed(datetime.now())
    seed = random.random()
    compile_str_list = []
    print(benchmarks[idx-1])
    #i = 0
    #benchmark = benchmarks[idx-1][-1]
    for i, benchmark in enumerate(benchmarks[idx-1]):
        for base_type in base_types:
            for GAP_FLAG in range(2):
                benchmark_name_run = benchmark + "_".join(["", base_type, "GAPFLAG", str(GAP_FLAG), "RUN", "REPEAT", str(param.REPEAT), "POW1", str(param.POW1), "POW2", str(param.POW2)])
                # one for each power of two with REPEAT=1
                benchmark_names_heap = [benchmark + "_".join(["", base_type, "GAPFLAG", str(GAP_FLAG), "HEAP", "POW", str(pow)]) for pow in range(param.POW1, param.POW2+1)]

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

                compile_str_list.append([compile_str_run] + compile_strs_heap)
                # case: just create file list, source files already exist
                if continue_flag == True:
                    continue
                template_file = "benchmark{}.template.cpp".format(idx)
                # substitute datatype in template cpp file
                if (os.path.isfile(template_file) == False):
                    print("ERROR: file '" + template_file + "' not in current directory.")
                    sys.exit(1)

                # copy template
                template_file_new = cpp_file_run
                print("src_dir = {}".format(src_dir))
                print("template_file_new = {}".format(template_file_new))

                copyfile(template_file, os.path.join(src_dir, template_file_new))

                # print info string
                sed_cmd = sed.format("info", "auto-generated benchmark", template_file_new)
                os.system(sed_cmd)
                # set seed
                sed_cmd = sed.format("seed", seed, template_file_new)
                os.system(sed_cmd)

                # substitute datatype in file
                alphabet_type_rx = re.compile("(\w+\d+)_?.*?")
                alphabet_type = alphabet_type_rx.search(base_type).group(1)
                letter_A = "alphabet_type::A"
                if i is 0:  #  ="Gapped Sequence":
                    letter_A = "alphabet_type{{{}::A}}".format(alphabet_type)
                sed_cmd = sed.format("letter_A", letter_A, template_file_new)
                os.system(sed_cmd)
                # substitute underlying sequence data type
                sed_cmd = sed.format("alphabet_type", alphabet_type, template_file_new)
                if i is 0:  # = "Gapped Sequence"
                    sed_cmd = sed.format("alphabet_type", "gapped<{}>".format(alphabet_type), template_file_new)
                os.system(sed_cmd)
                # substitute gapped alphabet type
                sed_cmd = sed.format("gapped_alphabet_type", "union_composition<{}, gap>".format(alphabet_type), template_file_new)
                os.system(sed_cmd)
                # substitute value_type
                #case gapped sequence: vector<alphabet_type>, else: vector<gapped<alphabet_type>>
                value_type = alphabet_type
                if i == 0:
                    value_type = "gapped<{}>".format(value_type)
                sed_cmd = sed.format("value_type", value_type, template_file_new)
                os.system(sed_cmd)

                # substitute container type of underlying sequence
                container_type = "std::vector<{}>"
                if base_type.endswith("_compressed") == True:
                    container_type = "seqan3::bitcompressed_vector<{}>"
                sed_cmd = sed.format("container_type", container_type.format(value_type), template_file_new)
                os.system(sed_cmd)

                # substitute gap_decorator
                print("i = {}".format(i))
                print(gap_decorators)
                gap_decorator = gap_decorators[i]
                sed_cmd = sed.format("gap_decorator", gap_decorator, template_file_new)
                os.system(sed_cmd)

                # substitute LOG_LEVEL macro
                sed_cmd = sed.format("LOG_LEVEL", abs(hash(benchmark_name_run)) % 10**8, template_file_new)
                os.system(sed_cmd)

                # set GAP_FLAG
                sed_cmd = sed.format("GAP_FLAG", GAP_FLAG, template_file_new)
                os.system(sed_cmd)

                # distribute current copy (note: cpp_file_run is identical with template_file_new up to here)
                for cpp_file in cpp_files_heap:
                    copyfile(os.path.join(src_dir, template_file_new), os.path.join(src_dir, cpp_file))

                # substitute output filename
                result_file_run = os.path.join(result_dir, benchmark_name_run + ".csv")
                sed_cmd = sed.format("benchmark.csv", result_file_run, template_file_new)
                os.system(sed_cmd)

                # substitute benchmark name
                sed_cmd = sed.format("benchmark", benchmark_name_run, cpp_file_run)
                os.system(sed_cmd)
                for cpp_file, benchmark_name in zip(cpp_files_heap, benchmark_names_heap):
                    sed_cmd = sed.format("benchmark", benchmark_name, cpp_file)
                    os.system(sed_cmd)

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

    os.system("rm ./src/*.bak")
    print("STATUS: binary list")
    for s in compile_str_list:
        for item in s:
            print(item.split(" ")[-1])
    return compile_str_list

'''
 [continue_flag] assumes running has been interrupted, i.e. all binaries exist,
 execution will be continued for the first file in the list for which no valid
 result files are found
'''
if __name__ == "__main__":
    if len(sys.argv) not in [9, 10]:
        print("Usage: python main.py benchmark_idx[1:3] home_dir seqan3_dir rangev3_dir sdsl_dir repeat:int pow1:int pow2:int [continue]")
        for arg in sys.argv:
            print(arg)
    else:
        start = timeit.default_timer()
        dirs.home = sys.argv[2]
        dirs.seqan3 = sys.argv[3]
        dirs.rangev3 = sys.argv[4]
        dirs.sdsl = sys.argv[5]
        param.REPEAT = int(sys.argv[6])
        param.POW1 = int(sys.argv[7])
        param.POW2 = int(sys.argv[8])
        continue_flag = True if len(sys.argv) == 10 else False
        print("STATUS: Generate source files ...")
        compile_string_llist = generate_src_files(int(sys.argv[1]), continue_flag)
        print("STATUS: Source file generation done\nSTATUS: Compile source files ...")
        binaries = compile_parallel(compile_string_llist, continue_flag)
        print(binaries)
        print("STATUS: Source file compilation done\nSTATUS: Execute binaries in parallel ...")
        result_dict = run_parallel(binaries, continue_flag)
        print("STATUS: Binary execution done")
        stop = timeit.default_timer()
        print_results(result_dict, "out_{}_{}_{}_{}.csv".format(sys.argv[1], param.REPEAT, param.POW1, param.POW2))
        print('Time: ' + str(stop - start))
        #cleanup()
