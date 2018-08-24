// <[info]>

#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include <iostream>
#include <random>
#include <unordered_map>
#include <utility>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

#include "anchor_list.hpp"
#include "anchor_set.hpp"
#include "gap_vector.hpp"
#include "utilities.hpp"

// g++ -std=c++17 -DNDEBUG -O3 -msse4.2 -I/<path_to>/include -L/<path_to>/lib -lsdsl -ldivsufsort -ldivsufsort64 -I/<path_to>/seqan3/include/ -I. -I/<path_to>/range-v3/include/ -fconcepts -Wall -Wextra benchmark1.cpp -o <[benchmark]>

#define LOG_LEVEL_<[benchmark]> 0
#define SEED <[seed]>

using namespace seqan3;

void benchmark1(String binary_name, long long int REPEAT, short unsigned POW1, short unsigned POW2)
{
    
    std::cout << "Starting " << binary_name << " with REPEAT=" << REPEAT << " ..." << std::endl;
    
    using alphabet_type = <[alphabet_type]>; // to be replace by parser
    using inner_type = typename std::vector<alphabet_type>;
    using size_type = typename ranges::v3::size_type_t<inner_type>;
    using time_t = long long int;
    
    // setup
    std::cout << "set up seq_length ...\n";
    std::vector<size_type> seq_lengths;
    for (short unsigned int p = POW1-1; p < POW2; ++p) seq_lengths.push_back(1 << (p+1));
    std::vector<alphabet_type> seq;    // gap-free sequences for all 3 implementations
    std::vector<gapped<alphabet_type>> gs;          // gapped sequence
    
    std::vector<size_type> gaps(0);
    // store accumulated statistics
    std::vector<stats<time_t, size_type>> results;
    
    // seed random number generator
    std::default_random_engine generator;
    generator.seed(seed);
    
    std::array<long long int, N*R> durations;
    time_t avg_duration, stddev_duration;
    for (auto seq_len : seq_lengths)
    {
        // reset durations
        if (LOG_LEVEL_<[benchmark]>) std::cout << "reset durations ...\n";
        std::fill(durations[i].begin(),durations[i].end(), 0);
        if (LOG_LEVEL_<[benchmark]>) std::cout << "done, continue benchmark with seq_len = " << seq_len << std::endl;
        for (long unsigned int round = 0; round < N; ++round)
        {
            if (LOG_LEVEL_<[benchmark]>) std::cout << "start resizing ...\n";
            seq.resize(seq_len);
            if (LOG_LEVEL_<[benchmark]>) std::cout << "done resizing.\n";
            
            //fill vector with A
            if (LOG_LEVEL_<[benchmark]>) std::cout << "start filling sequence ...\n";
            std::fill(seq.begin(), seq.end(), <[letter_A]>);
            if (LOG_LEVEL_<[benchmark]>) std::cout << "done filling.\n";
            
            if (LOG_LEVEL_<[benchmark]>) std::cout << "init seq length is : " << seq.size() << std::endl;
            // sample gap lengths
            gaps.resize(seq_len);
            sample<size_type>(&gaps, seq_len);
            if (LOG_LEVEL_<[benchmark]>) std::cout << "sampling done" << std::endl;
            
            // i) anchor list implementation
            if (LOG_LEVEL_<[benchmark]>) std::cout << "init gap decorator al ...\n";
            
            <[gap_decorator]> gap_decorator(&seq);
            
            // insert gaps
            size_type gap_acc = 0;
            if (GAP_FLAG)
            {
                for (size_type i = gaps.size()-1; i != 0; --i)
                {
                    if (LOG_LEVEL_<[benchmark]>) std::cout << "gap_len = " << gaps[i] << std::endl;
                    if (gaps[i] > 0)
                    {
                        if (LOG_LEVEL_<[benchmark]>) std::cout << "insert gap (" << i << ", " << gaps[i] << ") into structure ...\n";
                        gap_decorator.insert_gap(i, gaps[i]);
                        if (LOG_LEVEL_<[benchmark]>) print_sequence<anchor_set<std::vector<alphabet_type>>>(as);
                        
                        gap_acc += gaps[i];
                    }
                }
            }
            if (LOG_LEVEL_<[benchmark]>) std::cout << "num gaps inserted: " << gap_acc << std::endl;
            // shorten container for not exceeding target sequence length
            if (LOG_LEVEL_<[benchmark]>) std::cout << "resize with final seq_len = " << seq_len << std::endl;
            gap_decorator.resize(seq_len);
            if (LOG_LEVEL_<[benchmark]>) std::cout << "... done\n";
            
            int gap_ctr = 0;
            for (auto s : gs){ if (s == gap::GAP) ++gap_ctr;}
            std::cout << "gap proportion: " << (float)gap_ctr/(float)seq_len << std::endl;
            
            if (LOG_LEVEL_<[benchmark]>) {std::cout <<"aligned sequence = "; print_sequence<<[gap_decorator>]>>(gap_decorator);}
            
            // perform reading at random positions R times
            std::mt19937 generator(seed); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<> uni_dis(0.0, static_cast<double>(seq_len));
            
            size_type pos;
            std::chrono::high_resolution_clock::time_point t1, t2;
            
            std::vector<gapped<alphabet_type>> aux(10);
            for (size_type j = 0; j < R>>1; ++j)
            {
                // sample read position
                pos = uni_dis(generator);
                
                // read on anchor list (AS)
                t1 = std::chrono::high_resolution_clock::now();
                [[maybe_unused]] auto val = gap_decorator[pos];
                aux[j % aux.size()] = val;
                t2 = std::chrono::high_resolution_clock::now();
                durations[round*R + j] = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
                
                //std::cout << "read at " << pos << ": " << durations[round*R + j] << std::endl;
                
            } // R read ops
            for (auto a : aux)  std::cout << a << std::endl;
        } // N experiment repetitions
        
        std::cout << seq_len << " done\n";
        time_t quantile = 99;
        avg_duration = avg<R*N>(&durations, quantile);
        stddev_duration = stddev<R*N>(&durations, quantile);
        if (LOG_LEVEL_<[benchmark]>) std::cout << "avg duration time of reading of approach " << i << ": " << avg_duration << "\u00B1" << stddev_duration << " ns\n";
        results.push_back(stats<time_t, size_type>{seq_len, N*R, avg_duration, stddev_duration});
        
    } // seq_len
    std::cout << "Benchmark 1: read-only at random positions on GAP_FLAG= " << GAP_FLAG << " sequence\nwith " << N << " repetitions and " << R << " read operations per experiment\n";
    std::cout << "seq_len\t\tavg\t\tstddev\n---------------------------------------\n";
    
    for (auto it = results.begin(); it < results.end(); ++it)
    {
        std::cout << (*it).seq_len << "\t\t";
        std::cout << (*it).avg << "\u00B1" << (*it).stddev << "ns\t";
        std::cout << std::endl;
    }
    // write to csv file
    std::ofstream outfile;
    outfile.open ("<[benchmark.csv]>");
    outfile << "<[benchmark]>: " << "Read-only at random positions on GAP_FLAG= " << GAP_FLAG << " sequence\nwith " << N << " repetitions and " << R << " read operations per experiment\n";
    
    outfile << "#seq_len,avg[ns],stddev[ns]\n";
    for (auto it = results.begin(); it < results.end(); ++it)
        outfile << (*it).seq_len << "," << (*it).avg << "," << (*it).stddev << "\n";
    
    outfile.close();
}

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cout << "Usage: <[benchmark]> <REPEAT> <POW1> <POW2>\n";
        sys.exit(-1);
    }
    benchmark1(argv[0], atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
}
