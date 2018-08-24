#include <algorithm>
#include <array>
#include <cstdlib>
#include <chrono>
#include <climits>
#include <cmath>
#include <fstream>
#include <numeric>
#include <iostream>
#include <random>
#include <utility>

//#include <seqan3/alphabet/nucleotide/dna15.hpp>

//#include "anchor_list.hpp"
//#include "anchor_set.hpp"
//#include "gap_vector.hpp"

#define LOG_LEVEL_UTIL 0
#define AS 0    // anchor set
#define AL 1    // anchor list
#define GV 2    // gap vector
#define GS 3    // gapped sequence vector

typedef long double time_type;

//using namespace seqan3;

//############################################ HELPER START ###########################################
// sample gap lengths in range [0:9], in expectation 36% are gaps of length 1 or larger
template<typename SIZE_TYPE>
void sample(std::vector<SIZE_TYPE> * gap_vector, SIZE_TYPE size)
{
    std::default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::array<double,10> cumsum = {0.6395, 0.8263, 0.8871, 0.9257, 0.9544, 0.9709, 0.9813, 0.9890, 0.9955, 1.0000};
    for (SIZE_TYPE i = 0; i < size; ++i){
        double y = uni(generator);
        gap_vector->at(i) = y;
        auto it = std::find_if(cumsum.begin(), cumsum.end(), [y](double i){return y <= i;});
        if (LOG_LEVEL_UTIL) std::cout << it - cumsum.begin() << std::endl;
        gap_vector->at(i) = it - cumsum.begin();
    }
}

// sample single gap length
template<typename SIZE_TYPE>
SIZE_TYPE sample()
{
    std::default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> uni(0.6396, 1.0);
    std::array<double,10> cumsum = {0.6395, 0.8263, 0.8871, 0.9257, 0.9544, 0.9709, 0.9813, 0.9890, 0.9955, 1.0000};
    double y = uni(generator);
    std::cout << "sampled y from [0.6396, 1.0]: " << y << std::endl;
    auto it = std::find_if(cumsum.begin(), cumsum.end(), [y](double i){return y <= i;});
    return it - cumsum.begin();
}

// average of quantile
template<long unsigned int L>
time_type avg(std::array<time_type, L> * durations, time_type quantile)
{
    std::sort(durations->begin(), durations->end());
    
    time_type offset = (L - L*quantile/100)/2;
    //std::cout << "offset = " << offset << ", L = " << L << ", quantile = " << quantile << std::endl;
    time_type aux = 0;
    time_type overflow = 0;
    auto it = durations->begin();
    auto it_end = durations->end();
    for (std::advance(it, offset), std::prev(it, offset); it < it_end; ++it)
    {
        // prevent overflow by temporal right shift
        if (aux > (std::numeric_limits<time_type>::max()/2))
        {
            ++overflow;
            aux /= 2;
            std::cout << "aux above threshold\n";
        }
        aux += (*it/std::pow(2, overflow));
    }
    //return std::accumulate(durations->begin()+offset, durations->end()-offset, 0, std::plus<time_type>()) / (L-2*offset);
    std::cout << "overflow = " << overflow << std::endl;
    return (overflow) ? (aux/(L-2*offset))*(std::pow(2, overflow)) : (aux/(L-2*offset));
}

// standard deviation
template<long unsigned int L>
time_type stddev(std::array<time_type, L> * durations, long int quantile)
{
    time_type offset = (L - L*quantile/100)/2;
    std::cout << "offset = " << offset << std::endl;
    time_type avg_duration = avg<L>(durations, quantile);
    std::cout << "avg =\t" << avg_duration << std::endl;
    auto lambda = [avg_duration, offset](time_type acc, time_type t){return acc + ((t-avg_duration)/(L-1-2*offset)*(t-avg_duration)); };
    if (!std::is_sorted(durations->begin(), durations->end())) std::cout << "Error: durations are not sorted\n";
    time_type acc = 0;
    auto it = durations.begin();
    auto it_end = durations.end();
    for (std::advance(it, offset), std::prev(it_end, offset); it < it_end; ++it)
    {
        acc = lambda(acc, *it);
        std::cout << "intermediate acc = " << acc << std::endl;
    }
    //time_type d = std::sqrt(std::accumulate(durations->begin()+offset, durations->end()-offset, 0, lambda) / 1);
    return std::sqrt(acc);
}

template<typename time_type, typename SIZE_TYPE>
struct stats
{
    SIZE_TYPE seq_len;          // sequence length
    long unsigned int num_ops;  // number of read/write operations
    time_type avg;                 // mean exec time
    time_type stddev;              // std deviation
};

template <typename ALIGNED_SEQUENCE>
void print_sequence(ALIGNED_SEQUENCE &as)
{
    std::cout << "sequence length: " << as.size() << std::endl;
    for (unsigned int i = 0; i < as.size(); ++i)
        std::cout << as[i];
    std::cout << std::endl;
    for (unsigned int i = 0; i < as.size(); ++i)
        std::cout << (i%10);
    std::cout << std::endl;

}
