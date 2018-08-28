# Gap Decorator Benchmarks for SeqAn3

## HOWTO
In order to
  1. Create `cpp` source files for a specific benchmark [1-3]
  2. Generate binaries for each gap decorator approach, for each each dna datatype (currently `dna15` and `bitcompressed_vector<gapped<dna4>>`) in parallel
  3. Execute the binaries in parallel with sequence lengths in range $[2^{POW1}, 2^{POW2}]$
  4. Produce `csv` result files of execution times and actual space consumption measurements

set environmental variables (or name directly) home and installation directories for `seqan3` and `range-v3`, and run:

`python main.py <benchmark_idx> $home_dir $seqan3_dir $rangev3_dir <REPEAT> <POW1> <POW2>`

E.g. called with benchmark index 1, this will create the following 20 benchmark source files and place them into the `./src/` subfolder:

  * `benchmark1_GS_dna{4,15}_{0,1}.cpp`
    * `std::vector<gapped<dna>>` as benchmark baseline
  * `benchmark1_GVsd_dna{4,15}_{0,1}.cpp`
    * compressed bit vector of type `sdsl::sd_vector` and length |sequence| + #gap symbols
  * `benchmark1_GVbit_dna{4,15}_{0,1}.cpp`
    * compressed bit vector of type `sdsl::bit_vector` same length as above
  * `benchmark1_AL_dna{4,15}_{0,1}.cpp`
    * contiguous gaps stored as anchors `[(gap_pos, gap_length)]` with positions relative to the underlying sequence
  * `benchmark1_AS_dna{4,15}_{0,1}.cpp`
    * gap anchor positions stored in `std::unordered_set`, i.e. as a red-black tree
    * positions are relative to the aligned sequence and gap lengths accumulated

where `dna4` uses a `bitcompressed_vector<dna4>` as underlying sequence container and `dna15` the `std::vector<dna15>` container.
The python script then compiles the source files and adds a suffix based on the parameter settings for `REPEAT` and power range.
Finally, the binaries are executed in parallel using Python's `Process` class from the `multiprocessing` library. Result files are written into `./results` subdirectory. 

## Benchmarks

Benchmark selection is done by setting the `benchmark_idx` to a value in [1:3] where

* `idx`=1 tests Read-only
* `idx`=2 tests gap erase and insert at random positions
* `idx`=3 tests gap erase and insert from left to right
