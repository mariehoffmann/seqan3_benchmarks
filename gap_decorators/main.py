#!/usr/local/bin/python3
# script for generating 3*2*|methods| cpp files for each benchmark, method and datatype

from collections import namedtuple
import os
import random
from datetime import datetime
from multiprocessing import Process
import re
from shutil import copyfile
import sys

# directories for compilation
dirs = namedtuple("dirs", "home seqan3 rangev3")
# parameter settings for benchmarks
# REPEAT: number of experiments repetitions
# POW1, POW2: sequence length ranges from [2**POW1, 2**(POW1+1), ..., 2**(POW2)]
param = namedtuple("param", "REPEAT POW1 POW2")

alphabet_types = ['dna15', 'dna4']

src_dir = "./src"

compile_str = "g++ -std=c++17 -DNDEBUG -O3 -msse4.2 -I<home_dir>/include -L<home_dir>/lib -lsdsl -ldivsufsort -ldivsufsort64 -I<seqan3_dir>/seqan3/include/ -I. -I<rangev3_dir>/range-v3/include/ -fconcepts -Wall -Wextra {0}/<cpp_file> -o {0}/<benchmark>".format(src_dir)
info = ["Gapped Sequence", "Gap Vector sdsl::sd_vector", "Gap Vector sdsl::bit_vector", "Anchor List", "Anchor Set"]
benchmarks = [["benchmark" + str(i+1) + meth for meth in ["_GS", "_GVsd", "_GVbit", "_AL", "_AS"]] for i in range(3)]
print(benchmarks)
gap_decorators = ['std::vector<gapped<alphabet_type>>', 'gap_vector_sd', 'gap_vector_bit', 'anchor_list', 'anchor_set']

# sed command for in-place text substitutions
sed = "sed -i.bak -e 's|<\[{}\]>|{}|g' -- ./src/{}"

# compile single binary
def compile(cmd):
    os.system(cmd)

# compile in parallel
def compile_parallel(compile_strings):
    procs = []
    
    for cmd in compile_strings:
        print("next cmd: ")
        print(cmd)
        print(len(cmd))
        p = Process(target=compile, args=(str(cmd),))
        p.start()
        procs.append(p)
    for p in procs:
        p.join()

# run single binary
def run(binary):
    cmd = "{} {} {} {}".format(binary, param.REPEAT, param.POW1, param.POW2)
    os.system(cmd)

# run binaries in parallel
def run_parallel(binaries):
    procs = []
    for cmd in binaries:
        p = Process(target=run, args=(cmd))
        p.start()
        procs.append(p)
    for p in procs:
        p.join()

def generate_src_files(idx):
    # create seed for random number generator in benchmarks
    random.seed(datetime.now())
    seed = random.random
    src_files = []
    for i, benchmark in enumerate(benchmarks[idx-1]):
        for alphabet_type in alphabet_types:
            b = compile_str.replace("<home_dir>", dirs.home)
            b = b.replace("<seqan3_dir>", dirs.seqan3)
            b = b.replace("<rangev3_dir>", dirs.rangev3)
            b = b.replace("<benchmark_file>", benchmark)
            b = b.replace("<benchmark>", benchmark + "_".join(["", alphabet_type, "REPEAT", str(param.REPEAT), "POW1", str(param.POW1), "POW2", str(param.POW2)]))

            template_file = "benchmark{}.template.cpp".format(idx)
            # substitute datatype in template cpp file
            if (os.path.isfile(template_file) == False):
                print("Error: file '" + template_file + "' not in current directory.")
                sys.exit(1)
            benchmark_name = benchmark + "_" + alphabet_type
            # copy template
            cpp_file = benchmark_name + ".cpp"
            copyfile(template_file, os.path.join(src_dir, cpp_file))
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
            sed_cmd = sed.format("alphabet_type", alphabet_type, cpp_file)
            print(sed_cmd)
            os.system(sed_cmd)
            # substitute underlying sequence data type
            letter_A = "alphabet_type::A"
            if info[i] == "Gapped Sequence":
                letter_A = "gapped<alphabet_type>{alphabet_type::A}"
            sed_cmd = sed.format("letter_A", letter_A, cpp_file)
            print(sed_cmd)
            os.system(sed_cmd)
            # substitute gap_decorator
            gap_decorator = gap_decorators[i]
            sed_cmd = sed.format("gap_decorator", gap_decorator, cpp_file)
            os.system(sed_cmd)
            # substitute output filename
            sed_cmd = sed.format("benchmark.csv", benchmark_name, cpp_file)
            os.system(sed_cmd)
            print(sed_cmd)
            # substitute LOG_LEVEL macro
            sed_cmd = sed.format("benchmark", benchmark_name, cpp_file)
            os.system(sed_cmd)
            print(sed_cmd)
            src_files.append(b)
    
        break
    return src_files



#python main.py 1 $HOME_DIR $SEQAN3_DIR $RANGEV3_DIR 8 2 8

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python main.py benchmark_idx[1:3] home_dir seqan3_dir rangev3_dir repeat:int pow1:int pow2:int")
    else:
        dirs.home = sys.argv[2]
        dirs.seqan3 = sys.argv[3]
        dirs.rangev3 = sys.argv[4]
        param.REPEAT = int(sys.argv[5])
        param.POW1 = int(sys.argv[6])
        param.POW2 = int(sys.argv[7])
        src_files = generate_src_files(int(sys.argv[1]))
        for src_file in src_files:
            print(src_file)
        binaries = compile_parallel(src_files)
# run_parallel(binaries)

