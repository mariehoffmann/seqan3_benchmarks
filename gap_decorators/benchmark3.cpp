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
#include <seqan3/alphabet/nucleotide/dna15.hpp>

#include "anchor_list.hpp"
#include "anchor_set.hpp"
#include "gap_vector.hpp"
#include "utilities.hpp"

// g++ -std=c++17 -DNDEBUG -O3 -msse4.2 -I/<path_to>/include -L/<path_to>/lib -lsdsl -ldivsufsort -ldivsufsort64 -I/<path_to>/seqan3/include/ -I. -I/<path_to>/range-v3/include/ -fconcepts -Wall -Wextra benchmark3.cpp -o benchmark3

#define LOG_LEVEL 0
#define OVERFLOW_ALL 0
#define OVERFLOW_GV 6
#define OVERFLOW_AS 6
#define I 4   // number of methods

using namespace seqan3;

/*
* repeat N times, read R times, sequence lengths in [2^1, 2^2, ..., 2^P]
*/
template<long long int N, long long int R, short unsigned int P = 5, bool GAP_FLAG=1>
void benchmark3()
{
    std::cout << "Starting benchmark #3 with N = " << N << " repetitions ..." << std::endl;
    using inner_type = typename std::vector<dna15>;
    using size_type = typename ranges::v3::size_type_t<inner_type>;
    using time_t = long long int;
    using gap_vector = aligned_sequence_adaptor_constant_access<inner_type>;
    
    // setup
    std::vector<size_type> seq_lengths;
    for (short unsigned int p = 0; p < P; ++p) seq_lengths.push_back(1 << (p+1));
    std::vector<dna15> s_gv, s_as, s_al;    // gap-free sequences for all 3 implementations
    std::vector<gapped<dna15>> gs;          // gapped sequence
    std::vector<size_type> gaps(0);
    // store accumulated statistics
    std::unordered_map<unsigned short int, std::vector<stats<time_t, size_type>>> result_map = {
        {AS, std::vector<stats<time_t, size_type>>()},
        {AL, std::vector<stats<time_t, size_type>>()},
        {GV, std::vector<stats<time_t, size_type>>()},
        {GS, std::vector<stats<time_t, size_type>>()}
    };
    
    std::array<std::array<long long int, N*R>, I> durations;
    time_t avg_duration, stddev_duration;
    for (auto seq_len : seq_lengths)
    {
        if (LOG_LEVEL) std::cout << "seq_len = " << seq_len << std::endl;
        // reset durations
        if (LOG_LEVEL) std::cout << "reset durations ...\n";
        for (auto i = 0; i < I; ++i)
            std::fill(durations[i].begin(),durations[i].end(),0);
        if (LOG_LEVEL) std::cout << "done, continue benchmark with seq_len = " << seq_len << std::endl;
        for (long unsigned int round = 0; round < N; ++round)
        {
            if (LOG_LEVEL) std::cout << "start resizing ...\n";
            s_gv.resize(seq_len);
            gs.resize(seq_len);
            if (LOG_LEVEL) std::cout << "done resizing.\n";
            //fill vector with A
            if (LOG_LEVEL) std::cout << "start filling sequence ...\n";
            std::fill(s_gv.begin(), s_gv.end(), dna15::A);
            std::fill(gs.begin(), gs.end(), gapped<dna15>{dna15::A});
            if (LOG_LEVEL) std::cout << "done filling.\n";
            
            // copy for other 2
            if (LOG_LEVEL) std::cout << "start copying sequences ...\n";
            s_al = s_gv;
            s_as = s_gv;
            if (LOG_LEVEL) std::cout << "done.\n";
            
            if (LOG_LEVEL) std::cout << "init seq length is : " << s_gv.size() << std::endl;
            // sample gap lengths
            gaps.resize(seq_len);
            sample<size_type>(&gaps, seq_len);
            if (LOG_LEVEL) std::cout << "sampling done" << std::endl;
            
            // i) anchor list implementation
            anchor_list<inner_type> al(&s_al);
            // ii) anchor set implementation
            anchor_set<inner_type> as(&s_as);
            // iii) gap vector implementation
            gap_vector gv(&s_gv);
            
            // insert gaps
            if (GAP_FLAG)
            {
                size_type gap_acc = 0;
                for (size_type i = gaps.size()-1; i != 0; --i)
                {
                    if (LOG_LEVEL) std::cout << "gap_len = " << gaps[i] << std::endl;
                    if (gaps[i] > 0)
                    {
                        if (LOG_LEVEL) std::cout << "insert gap (" << i << ", " << gaps[i] << ") into AS ...\n";
                        as.insert_gap(i, gaps[i]);
                        if (LOG_LEVEL) print_sequence<anchor_set<std::vector<dna15>>>(as);
                        
                        if (LOG_LEVEL) std::cout << "and into AL ...\n";
                        al.insert_gap(i, gaps[i]);
                        if (LOG_LEVEL) print_sequence<anchor_list<std::vector<dna15>>>(al);
                        
                        if (LOG_LEVEL) std::cout << "and into GV ...\n";
                        gv.insert_gap(i, gaps[i]);
                        if (LOG_LEVEL) print_sequence<gap_vector>(gv);
                        
                        if (LOG_LEVEL) std::cout << "and into GS ...\n";
                        gs.insert(gs.begin()+i, gaps[i], gapped<dna15>{gap::GAP});
                        if (LOG_LEVEL) print_sequence<gap_vector>(gv);
                        
                        gap_acc += gaps[i];
                    }
                }
                
                if (LOG_LEVEL) std::cout << "num gaps inserted: " << gap_acc << std::endl;
                // shorten container for not exceeding target sequence length
                if (LOG_LEVEL) std::cout << "resize with final seq_len = " << seq_len << std::endl;
                al.resize(seq_len);
                if (LOG_LEVEL) std::cout << "... done for AL\n";
                as.resize(seq_len);
                if (LOG_LEVEL) std::cout << "... done for AS\n";
                gv.resize(seq_len);
                if (LOG_LEVEL) std::cout << "... done for GV\n";
                gs.resize(seq_len);
                if (LOG_LEVEL) std::cout << "... done for GS\n";
            }
            if (LOG_LEVEL) std::cout << "aligned size: " << as.size() << std::endl;
            
            if (LOG_LEVEL) print_sequence<anchor_set<inner_type>>(as);
            if (LOG_LEVEL) print_sequence<anchor_list<inner_type>>(al);
            if (LOG_LEVEL) print_sequence<gap_vector>(gv);
            
            // perform reading at random positions R times
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<> uni_dis(0.0, static_cast<double>(seq_len));
            
            size_type pos;
            std::chrono::high_resolution_clock::time_point t1, t2;
            
            for (size_type j = 0; j < R; ++j)
            {
                // insert/delete position from left to right
                pos = j % seq_len;
                
                if (LOG_LEVEL) std::cout << "seq_len = " << seq_len << ", insert gap at pos " << pos << std::endl;
                // read on anchor list (AS)
                t1 = std::chrono::high_resolution_clock::now();
                if (LOG_LEVEL) std::cout << "\tas.insert...\n";
                as.insert_gap(pos, 1);
                as.insert_gap(pos, 2);
                if (LOG_LEVEL) std::cout << "\tas.erase...\n";
                as.erase_gap(pos+1, pos+3);
                as.erase_gap(pos);
                t2 = std::chrono::high_resolution_clock::now();
                durations[AS][round*R + j] = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/4 >> OVERFLOW_AS;
                
                // read on anchor list (AL)
                t1 = std::chrono::high_resolution_clock::now();
                if (LOG_LEVEL) std::cout << "\tal.insert...\n";
                al.insert_gap(pos, 1);
                al.insert_gap(pos, 2);
                if (LOG_LEVEL) std::cout << "\tal.erase...\n";
                al.erase_gap(pos+1, pos+3);
                al.erase_gap(pos);
                t2 = std::chrono::high_resolution_clock::now();
                durations[AL][round*R + j] = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/4;
                
                // read on gap vector (GV)
                t1 = std::chrono::high_resolution_clock::now();
                if (LOG_LEVEL) std::cout << "\tgv.insert...\n";
                gv.insert_gap(pos, 1);
                gv.insert_gap(pos, 2);
                if (LOG_LEVEL) std::cout << "\tgv.erase...\n";
                gv.erase_gap(pos+1, pos+3);
                gv.erase_gap(pos);
                t2 = std::chrono::high_resolution_clock::now();
                durations[GV][round*R + j] = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/4 >> OVERFLOW_GV;
                
                // ins/del on gapped sequence vector (GS)
                t1 = std::chrono::high_resolution_clock::now();
                if (LOG_LEVEL) std::cout << "\tgs.insert...\n";
                gs.insert(gs.begin() + pos, 1, gapped<dna15>{gap::GAP});
                gs.insert(gs.begin() + pos, 2, gapped<dna15>{gap::GAP});
                if (LOG_LEVEL) std::cout << "\tgs.erase...\n";
                gs.erase(gs.begin()+pos+1, gs.begin()+pos+3);
                gs.erase(gs.begin()+pos);
                t2 = std::chrono::high_resolution_clock::now();
                durations[GS][round*R + j] = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/4;
                
                
                //std::cout << "read at " << pos << ": " << durations[round*R + j] << std::endl;
                
            } // R read ops
        } // N experiment repetitions
        
        time_t quantile = 99;
        std::cout << seq_len << " done\n";
        for (auto i = 0; i < I; ++i) // TODO: until I
        {
            avg_duration = avg<R*N>(&durations[i], quantile);
            stddev_duration = stddev<R*N>(&durations[i], quantile);
            if (LOG_LEVEL) std::cout << "avg duration time of reading of approach " << i << ": " << avg_duration << "\u00B1" << stddev_duration << " ns\n";
            result_map[i].push_back(stats<time_t, size_type>{seq_len, N*R, avg_duration, stddev_duration});
        }
    } // seq_len
    // reverse right shift for omitting overflow and print results
    
    std::cout << "Benchmark 3: insert/delete from left to right positions into GAP_FLAG=" << GAP_FLAG << " sequence\nwith " << N << " repetitions seq_len operations per experiment\n";
    std::cout << "seq_len\t\tavg_GS\t\tavg_GV\t\tavg_AL\t\tavg_AS\n-------------------------------------------------------------------------------------\n";
    auto it_AS = result_map[AS].begin(), it_AL = result_map[AL].begin(), it_GV = result_map[GV].begin(), it_GS = result_map[GS].begin();
    for (; it_AS < result_map[AS].end() & it_AL < result_map[AL].end(); ++it_AS, ++it_AL, ++it_GV, ++it_GS)
    {
        (*it_GV).avg <<= OVERFLOW_GV;
        (*it_GV).stddev <<= OVERFLOW_GV;
        (*it_AS).avg <<= OVERFLOW_AS;
        (*it_AS).stddev <<= OVERFLOW_AS;
        
        
        std::cout << (*it_AL).seq_len << "\t\t";
        std::cout << (*it_GS).avg << "\u00B1" << (*it_GS).stddev << "ns\t";
        std::cout << (*it_GV).avg << "\u00B1" << (*it_GV).stddev << "ns\t";
        std::cout << (*it_AL).avg << "\u00B1" << (*it_AL).stddev << "ns\t";
        std::cout << (*it_AS).avg << "\u00B1" << (*it_AS).stddev << "ns";
        std::cout << std::endl;
    }
    // write to csv file
    std::ofstream outfile;
    outfile.open ("benchmark3.csv");
    outfile << "Benchmark 3: insert/delete from left to right into GAP_FLAG=" << GAP_FLAG << " sequence\nwith " << N << " repetitions and seq_len operations per experiment\n";
    
    outfile << "#seq_len,avg_GV[ns],stddev,avg_GV[ns],stddev,tavg_AL,stddev,avg_AS,stddev\n";
    it_GS = result_map[GS].begin(), it_AS = result_map[AS].begin(), it_AL = result_map[AL].begin(), it_GV = result_map[GV].begin();
    for (; it_AS < result_map[AS].end() & it_AL < result_map[AL].end(); ++it_AS, ++it_AL, ++it_GV, ++it_GS)
    {
        outfile << (*it_AL).seq_len << ",";
        outfile << (*it_GS).avg << "," << (*it_GS).stddev << ",";
        outfile << (*it_GV).avg << "," << (*it_GV).stddev << ",";
        outfile << (*it_AL).avg << "," << (*it_AL).stddev << ",";
        outfile << (*it_AS).avg << "," << (*it_AS).stddev << "\n";
    }
    outfile.close();
    
}

int main(/*int argc, char** argv*/)
{
    benchmark3<1<<2, 1<<10, 18, 1>();
}
