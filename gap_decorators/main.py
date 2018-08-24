#!/usr/local/bin/python3
# script for generating 3*2*|methods| cpp files for each benchmark, method and datatype

from collections import namedtuple
import os
import random
from datetime import datetime
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

compile_str = "g++ -std=c++17 -DNDEBUG -O3 -msse4.2 -I<home_dir>/include -L<home_dir>/lib -lsdsl -ldivsufsort -ldivsufsort64 -I<seqan3_dir>/seqan3/include/ -I. -I<rangev3_dir>/range-v3/include/ -fconcepts -Wall -Wextra <cpp_file> -o <benchmark>"
info = ["Gapped Sequence", "Gap Vector sdsl::sd_vector", "Gap Vector sdsl::bit_vector", "Anchor List", "Anchor Set"]
benchmarks = [["benchmark" + str(i+1) + meth for meth in ["_GS", "_GVsd", "_GVbit", "_AL", "_AS"]] for i in range(3)]
print(benchmarks)
structures = ['std::vector<gapped<<[alphabet_type]>>>', 'gap_vector_sd', 'gap_vector_bit', 'anchor_list', 'anchor_set']
src_dir = './src'

def launch(idx):
    # create seed for random number generator in benchmarks
    random.seed(datetime.now())
    seed = random.random
    compile_list = []
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
            cpp_file = os.path.join(src_dir, benchmark_name + ".cpp")
            copyfile(template_file, cpp_file)
            b = b.replace("<cpp_file>", cpp_file)
            print(b)
            # print info string
            sed = "sed -i '.bak' 's/<[info]>/auto-generated benchmark/g' {}".format(cpp_file)
            print(sed)
            os.system(sed)
            sys.exit(0)
            # set seed
            sed = "sed -i '.bak' 's/<[seed]>/{}/g' {}".format(seed, cpp_file)
    
            # substitute datatype in file
            sed = "sed -i '.bak' 's/<[alphabet_type]>/{}/g' {}".format(alphabet_type, cpp_file)
            print(sed)
            # substitute underlying sequence data type
            letter_A = "alphabet_type::A"
            if info[i] == "Gapped Sequence":
                letter_A = "gapped<alphabet_type>{alphabet_type::A}"
            print(sed)
            sed = "sed -i '.bak' 's/<[letter_A]>/{}/g' {}".format(letter_A, cpp_file)
            # substitute structure
            structure = structures[i]
            sed = "sed -i '.bak' 's/<[gap_decorator]>/{}/g' {}".format(structure, cpp_file)
            # substitute output filename
            sed = "sed -i '.bak' 's/<[benchmark.csv]>/{}.csv/g' {}".format(benchmark_name, cpp_file)
            print(sed)
            # substitute LOG_LEVEL macro
            sed = "sed -i '.bak' 's/<[benchmark]>/{}/g' {}".format(benchmark_name, cpp_file)
            print(sed)
            compile_list.append(b)

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
        launch(int(sys.argv[1]))


