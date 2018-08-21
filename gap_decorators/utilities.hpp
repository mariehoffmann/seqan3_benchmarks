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

#include <seqan3/alphabet/nucleotide/dna15.hpp>

#include "anchor_list.hpp"
#include "anchor_set.hpp"
#include "gap_vector.hpp"

#define LOG_LEVEL_UTIL 0
#define AS 0    // anchor set
#define AL 1    // anchor list
#define GV 2    // gap vector
#define GS 3    // gapped sequence vector

using namespace seqan3;

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
template<long long int L>
time_t avg(std::array<long long int, L> * durations, long long int quantile)
{
    std::sort(durations->begin(), durations->end());
    
    time_t offset = (L - L*quantile/100)/2;
    //std::cout << "offset = " << offset << ", L = " << L << ", quantile = " << quantile << std::endl;
    return std::accumulate(durations->begin()+offset, durations->end()-offset, 0, std::plus<time_t>()) / (L-2*offset);
}

// standard deviation
template<long long int L>
time_t stddev(std::array<long long int, L> * durations, long int quantile)
{
    long int offset = (L - L*quantile/100)/2;
    long int avg_duration = avg<L>(durations, quantile);
    auto lambda = [avg_duration, offset](time_t acc, time_t t){return acc + ((t-avg_duration)*(t-avg_duration)/(L-1-2*offset)); };
    if (!std::is_sorted(durations->begin(), durations->end())) std::cout << "Error: durations are not sorted\n";
    time_t d = std::sqrt(std::accumulate(durations->begin()+offset, durations->end()-offset, 0, lambda) / 1);
    if (d >= 2384658)
    {
        for (long int i = 0; i < 10; ++i)
            std::cout << "duration[" << i << "] = " << durations->at(i) << std::endl;
        for (long int i = L-1; i >+ L-10; --i)
            std::cout << "duration[" << i << "] = " << durations->at(i) << std::endl;
        
    }
    
    return d;
}

template<typename TIME_T, typename SIZE_TYPE>
struct stats
{
    SIZE_TYPE seq_len;          // sequence length
    long unsigned int num_ops;  // number of read/write operations
    TIME_T avg;                 // mean exec time
    TIME_T stddev;              // std deviation
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
