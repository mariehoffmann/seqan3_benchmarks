// <[info]>
// Insert/Erase gap from left to right. Gaps of length 2 will be inserted and alternating head or tail erased.

#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <fstream>
#include <numeric>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

#include "../anchor_blocks.hpp"
#include "../anchor_blocks2.hpp"
#include "../anchor_list.hpp"
#include "../anchor_set.hpp"
#include "../gap_vector_bit.hpp"
#include "../gap_vector_sd.hpp"
#include "../gapped_sequence.hpp"
#include "../utilities.hpp"

// g++ -std=c++17 -fconcepts -Wall -Wextra -I $SEQAN3_DIR/include -I $RANGEV3_DIR/include -I $SDSL_DIR/include ./src/<[benchmark]>.cpp -o ./src/<[benchmark]>

#define LOG_LEVEL_<[LOG_LEVEL]> 0
#define NUM_OP 1024     // number of operations performed per experiment
#define REPEAT <[REPEAT]>
#define POW1 <[POW1]>
#define POW2 <[POW2]>
#define GAP_FLAG <[GAP_FLAG]>   // 0: operate on ungapped sequence, 1: operate on already gapped sequence

using namespace seqan3;

void benchmark3(int csv_flag)  //std::string const & binary_name)
{
    using alphabet_type = <[alphabet_type]>; // to be replace by parser
    using inner_type = typename <[container_type]>; //std::vector<alphabet_type>;
    using size_type = typename ranges::v3::size_type_t<inner_type>;
    using time_type = long double;

    // setup
    std::vector<size_type> seq_lengths;
    for (short unsigned int p = POW1-1; p < POW2; ++p) seq_lengths.push_back(1 << (p+1));
    // container for gap-free sequence either std::vector<alphabet_type> or seqan3::bitcompressed_vector<alphabet_type>
    inner_type seq;

    // gapped sequence for measuring gap to dna ratio
    std::vector<gapped<alphabet_type>> gs;

    std::vector<size_type> gaps(0);
    // store accumulated statistics
    std::vector<stats<time_type, size_type>> results;
    std::array<time_type, NUM_OP*REPEAT> durations;
    time_type avg_duration, stddev_duration;
    for (auto seq_len : seq_lengths)
    {
        // reset durations
        std::fill(durations.begin(),durations.end(), 0);
        size_type pos = 0, pos_del;
        for (long unsigned int round = 0; round < REPEAT; ++round)
        {
            seq.resize(seq_len);

            //fill vector with A
            std::fill(seq.begin(), seq.end(), <[letter_A]>);
            // sample gap lengths
            gaps.resize(seq_len);
            sample<size_type>(&gaps, seq_len);

            // insert gaps
            size_type gap_acc = 0;
            if (GAP_FLAG)  // estimate new ungapped sequence length
            {
                // shorten gap vector s.t. final sequence length matches seq_len

                size_type letter_acc = 0;
                size_type gap_pos = 0;
                while (gap_pos < gaps.size() && gap_acc + letter_acc < seq_len)
                {
                    if (!gaps[gap_pos])
                        ++letter_acc;
                    else
                    {
                        if (letter_acc + gap_acc + gaps[gap_pos] > seq_len)
                        {
                            gaps[gap_pos] = seq_len - gap_acc - letter_acc;
                            gap_acc += gaps[gap_pos];
                            ++gap_pos;
                            break;
                        }
                        else
                            gap_acc += gaps[gap_pos];
                    }
                    ++gap_pos;
                }
                seq.resize(letter_acc); // resize ungapped sequence
                gaps.resize(gap_pos);
            }
            // initialize gap_decorator with (resized) sequence
            <[gap_decorator]><inner_type> gap_decorator(&seq);

            if (GAP_FLAG) // insert gaps
            {
                gap_acc = 0;
                for (size_type i = 0; i < gaps.size(); ++i)
                {
                    if (gaps[i]) if (LOG_LEVEL_<[LOG_LEVEL]>) std::cout << "insert gap (" << i + gap_acc << ", " << gaps[i] << ") into structure ...\n";
                    if (gaps[i])
                        gap_decorator.insert_gap(std::min(i + gap_acc, gap_decorator.size()), gaps[i]);
                    gap_acc += gaps[i];
                    //if (gap_decorator.size() >= (seq_len << 1)) break;
                }
            }

            int gap_ctr = 0;
            for (auto s : gs){ if (s == gap::GAP) ++gap_ctr;}

            std::chrono::high_resolution_clock::time_point t1, t2;

            // case gapped sequence: vector<alphabet_type>, else: vector<gapped<alphabet_type>>
            bool aux = true;

            for (size_type j = 0; j < std::min<size_type>(NUM_OP, seq_len); ++j)
            {
                // read from left to right
                pos = j;

                pos_del = (j%2) ? pos : pos+1;
                // insert & erase
                t1 = std::chrono::high_resolution_clock::now();
                aux &= gap_decorator.insert_gap(pos, 2);
                aux &= gap_decorator.erase_gap(pos_del, pos_del+1);
                t2 = std::chrono::high_resolution_clock::now();
                durations[round*REPEAT + j] = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/2;
                ++pos;
            } // REPEAT read ops
            //std::cerr << aux << std::endl;
        } // N experiment repetitions

        unsigned int quantile = 99;
        avg_duration = avg<REPEAT*NUM_OP>(&durations, quantile);
        stddev_duration = stddev<REPEAT*NUM_OP>(&durations, quantile);
        if (LOG_LEVEL_<[LOG_LEVEL]>) std::cout << "avg duration time of reading of approach: " << avg_duration << "\u00B1" << stddev_duration << " ns\n";
        results.push_back(stats<time_type, size_type>{seq_len, NUM_OP*REPEAT, avg_duration, stddev_duration});

    } // seq_len
    std::cout << "Benchmark 3: insert with subsequent erase at from left to right into gapped = " << GAP_FLAG << " sequence\nwith min(seq_len, " << REPEAT << ") repetitions and " << NUM_OP << " read operations per experiment\n";
    std::cout << "seq_len\t\tavg\u00B1stddev\t\tspace consumption\n---------------------------------------\n";

    for (auto it = results.begin(); it < results.end(); ++it)
    {
        std::cout << (*it).seq_len << "\t\t";
        std::cout << (*it).avg << "\u00B1" << (*it).stddev << "ns\t";
        std::cout << std::endl;
    }
    // write to csv file
    if (csv_flag)
    {
        std::ofstream outfile;
        outfile.open ("<[benchmark.csv]>");
        //outfile << "<[benchmark]>: " << "Read-only at random positions on GAP_FLAG= " << GAP_FLAG << " sequence\nwith " << REPEAT << " repetitions and " << NUM_OP << " read operations per experiment\n";
        outfile << "#seq_len,avg[ns],stddev[ns],space\n";
        for (auto it = results.begin(); it < results.end(); ++it)
            outfile << (*it).seq_len << "," << (*it).avg << "," << (*it).stddev << "\n";

        outfile.close();
    }
}

int main(int argc, char** argv)
{
    if (argc > 2)
    {
        std::cout << "Usage: <[benchmark]> [csv_flag]\n";
        return 2;
    }
    std::cout << "Start " << argv[0] << " ...\n";
    benchmark3((argc == 2) ? atoi(argv[1]) : 0);
    std::cout << "Finished " << argv[0] << " successfully.\n";
    return 1;
}
