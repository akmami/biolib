#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "../sketch/blend.h"
#include "utils.h"


#define BLEND_BITS_HIFI 32
#define BLEND_NEIGHBOR_NUMBER_HIFI 5

void analyze(uint128_t *minimizers, uint64_t len,
             int &contiguous_counts,
             int (&distances)[DISTANCE_LENGTH],
             std::vector<uint32_t> &distancesXL,
			 gap &distance_gaps,
             int (&lengths)[DISTANCE_LENGTH],
             std::vector<uint32_t> &lengthsXL,
			 gap &length_gaps,
			 int (&overlap_lengths)[DISTANCE_LENGTH], 
			 std::vector<uint32_t> &overlap_lengthsXL, 
			 gap &overlap_length_gaps) {

	if (len > 0) {

		bool isOverlapped = false;

		if (BLEND_GET_LENGTH(minimizers[0]) < DISTANCE_LENGTH) {
			lengths[BLEND_GET_LENGTH(minimizers[0])] += 1;
		} else {
			lengthsXL.push_back(BLEND_GET_LENGTH(minimizers[0]));
		}

		for (uint64_t i = 1; i < len; i++) {

			uint32_t current_start = BLEND_GET_INDEX(minimizers[i]);
			uint32_t current_end = current_start + BLEND_GET_LENGTH(minimizers[i]);
			uint32_t previous_start = BLEND_GET_INDEX(minimizers[i-1]);
			uint32_t previous_end = previous_start + BLEND_GET_LENGTH(minimizers[i-1]);

			if (current_start <= previous_end) {
				contiguous_counts += 1;
				isOverlapped = true;
			}

			// Distance
			uint32_t distance = current_start - previous_start;
			if (distance < DISTANCE_LENGTH ) {
				distances[distance]++;
			} else {
				distancesXL.push_back(distance);
			}

			distance_gaps.minimum = min(distance_gaps.minimum, distance);
			distance_gaps.maximum = max(distance_gaps.maximum, distance);

			// Length
			uint32_t length = current_end - current_start;
			if (length < DISTANCE_LENGTH) {
				lengths[length] += 1;
			} else {
				lengthsXL.push_back(length);
			}

			length_gaps.minimum = min(length_gaps.minimum, length);
			length_gaps.maximum = max(length_gaps.maximum, length);

			// Overlap
			if (previous_end >= current_start) {
				uint32_t overlap = previous_end - current_start;
				if (overlap < DISTANCE_LENGTH) {
					overlap_lengths[overlap] += 1;
				} else {
					overlap_lengthsXL.push_back(overlap);
				}
				overlap_length_gaps.minimum = min(overlap_length_gaps.minimum, overlap);
				overlap_length_gaps.maximum = max(overlap_length_gaps.maximum, overlap);			
			}
		}

		if (isOverlapped) {
			contiguous_counts += 1;
		}
	}
};

void process(std::string &sequence,
             int window, int kmer_size,
             int &core_counts,
             int &contiguous_counts,
             std::unordered_map<uint32_t, size_t> &distinct_cores,
             std::chrono::milliseconds &durations,
             int (&distances)[DISTANCE_LENGTH],
             std::vector<uint32_t> &distancesXL,
			 gap &distance_gaps,
             int (&lengths)[DISTANCE_LENGTH],
             std::vector<uint32_t> &lengthsXL,
			 gap &length_gaps,
			 int (&overlap_lengths)[DISTANCE_LENGTH], 
			 std::vector<uint32_t> &overlap_lengthsXL, 
			 gap &overlap_length_gaps,
             double &sizes) {

	auto start = std::chrono::high_resolution_clock::now();

    uint128_t *minimizers;
    uint64_t minimizers_len = 0;
    
    minimizers_len = blend_sketch(sequence.c_str(), sequence.length(), window, kmer_size, BLEND_BITS_HIFI, BLEND_NEIGHBOR_NUMBER_HIFI, 0, &minimizers);

	auto extraction_end = std::chrono::high_resolution_clock::now();
	durations += std::chrono::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(extraction_end - start).count());
	core_counts += minimizers_len;
    sizes += sizeof(uint128_t) * minimizers_len;
	
	analyze(minimizers, minimizers_len, contiguous_counts, distances, distancesXL, distance_gaps, lengths, lengthsXL, length_gaps, overlap_lengths, overlap_lengthsXL, overlap_length_gaps);

    for (uint64_t index = 0; index < minimizers_len; index++) {
        distinct_cores[BLEND_GET_KMER(minimizers[index])]++;
    }

	std::cout << "Length of the processed sequence: " << format_int(sequence.size()) << std::endl;

	sequence.clear();

	if (minimizers_len) free(minimizers);
}

int main(int argc, char **argv) {

	if (argc < 4) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] [kmer-size] [window]" << std::endl;
		return -1;
	}

	std::ifstream input(argv[1]);
	if (!input.good()) {
		std::cerr << "Error opening: " << argv[1] << " . You have failed." << std::endl;
		return -1;
	}

    int kmer_size = atoi(argv[2]);
    int window = atoi(argv[3]);

	// variables
	std::string line;

	std::fstream genome;
	genome.open(argv[1], std::ios::in);

    // section 1
    int core_counts = 0;
	int contiguous_counts = 0;
    std::unordered_map<uint32_t, size_t> distinct_cores;
    std::chrono::milliseconds durations;
	// section 2
    int distances[DISTANCE_LENGTH] = {0};
    std::vector<uint32_t> distancesXL;
	gap distance_gaps;
    distance_gaps.minimum = UINT64_MAX;
    distance_gaps.maximum = 0;
	// section 3
	int lengths[DISTANCE_LENGTH] = {0};
    std::vector<uint32_t> lengthsXL;
	gap length_gaps;
    length_gaps.minimum = UINT64_MAX;
    length_gaps.maximum = 0;
	// section 4
	int overlap_lengths[DISTANCE_LENGTH] = {0};
	std::vector<uint32_t> overlap_lengthsXL;
	gap overlap_length_gaps;
    overlap_length_gaps.minimum = UINT64_MAX;
    overlap_length_gaps.maximum = 0;
    // section 5
    double sizes = 0;
    
    // other
	size_t genome_size = 0;
	size_t count_once = 0;
	size_t relaxed_unique = 0;

	// read file
	if (genome.is_open()) {

		// process genome
		std::string sequence, id;
		sequence.reserve(250000000);

		std::cout << "Program begins" << std::endl;

		while (getline(genome, line)) {

			if (line[0] == '>') {

				// process previous chromosome before moving into new one
				if (sequence.size() != 0) {
					genome_size += sequence.size();
					process(sequence, window, kmer_size, core_counts, contiguous_counts, distinct_cores, durations, distances, distancesXL, distance_gaps, lengths, lengthsXL, length_gaps, overlap_lengths, overlap_lengthsXL, overlap_length_gaps, sizes);
				}

				id = line.substr(1);
				std::cout << "Processing started for " << id << std::endl;
				continue;

			} else if (line[0] != '>') {
				sequence += line;
			}
		}

		if (sequence.size() != 0) {
			genome_size += sequence.size();
			process(sequence, window, kmer_size, core_counts, contiguous_counts, distinct_cores, durations, distances, distancesXL, distance_gaps, lengths, lengthsXL, length_gaps, overlap_lengths, overlap_lengthsXL, overlap_length_gaps, sizes);
		}

		for (const auto& [value, count] : distinct_cores) {
			if (count == 1) {
				++count_once;
			}
		}

		// test for uniqueness
		genome.clear();
		genome.seekg(0, std::ios::beg);

		while (getline(genome, line)) {

			if (line[0] == '>') {

				if (sequence.size() != 0) {
					uint128_t *minimizers;
					uint64_t minimizers_len = 0;
					
					minimizers_len = blend_sketch(sequence.c_str(), sequence.length(), window, kmer_size, BLEND_BITS_HIFI, BLEND_NEIGHBOR_NUMBER_HIFI, 0, &minimizers);

					if (minimizers_len && distinct_cores[BLEND_GET_KMER(minimizers[0])] == 1) {
						relaxed_unique++;
					}

					for (uint64_t i = 1; i < minimizers_len; i++) {
						if (distinct_cores[BLEND_GET_KMER(minimizers[i])] == 1 || distinct_cores[BLEND_GET_KMER(minimizers[i-1])] == 1) {
							relaxed_unique++;
						}
					}
					if (minimizers_len) free(minimizers);
				}

				continue;

			} else if (line[0] != '>') {
				sequence += line;
			}
		}

		if (sequence.size() != 0) {
			uint128_t *minimizers;
			uint64_t minimizers_len = 0;
			
			minimizers_len = blend_sketch(sequence.c_str(), sequence.length(), window, kmer_size, BLEND_BITS_HIFI, BLEND_NEIGHBOR_NUMBER_HIFI, 0, &minimizers);

			if (minimizers_len && distinct_cores[BLEND_GET_KMER(minimizers[0])] == 1) {
				relaxed_unique++;
			}

			for (uint64_t i = 1; i < minimizers_len; i++) {
				if (distinct_cores[BLEND_GET_KMER(minimizers[i])] == 1 || distinct_cores[BLEND_GET_KMER(minimizers[i-1])] == 1) {
					relaxed_unique++;
				}
			}
			if (minimizers_len) free(minimizers);
		}

		genome.close();
	}	

	// Total Cores
	std::cout << "Total \\# Cores: " << format_int(core_counts) << std::endl;

	// Contiguous Cores
	std::cout << "Contiguous Cores: " << format_int(contiguous_counts) << std::endl;

	// Distinct Cores
    std::cout << "Distinct Cores: " << format_int(distinct_cores.size()) << std::endl;

	// Unique Cores
    std::cout << "Unique Cores: " << format_int(count_once) << std::endl;

	// Uniqueness Ratio
    std::cout << "Uniqueness %: " << format_double(((double)count_once) / ((double)core_counts)) << std::endl;

	// Relaxed Uniqueness Ratio
    std::cout << "Relaxed-Uniqueness %: " << format_double(((double)relaxed_unique) / ((double)core_counts)) << std::endl;

	// Execution Time
	std::cout << "Exec. Time (sec): " << format_double(((double)durations.count()) / 1000) << std::endl;

	// Mean Core Distances
	std::cout << "Avg Distance: " << format_double(mean(distances, distancesXL)) << std::endl;

	// Std Dev of Distances
	std::cout << "StdDev Distance: " << format_double(stdev(distances, distancesXL)) << std::endl;

	// Min/Max Core Length
	std::cout << "Min/Max Distance: " << format_int(distance_gaps.minimum) << "/" << format_int(distance_gaps.maximum) << std::endl;

	// Mean Core Length
	std::cout << "Avg Length: " << format_double(mean(lengths, lengthsXL)) << std::endl;

	// Std Dev of Lengths
	std::cout << "StdDev Length: " << format_double(stdev(lengths, lengthsXL)) << std::endl;

	// Min/Max Core Length
	std::cout << "Min/Max Length: " << format_int(length_gaps.minimum) << "/" << format_int(length_gaps.maximum) << std::endl;

	// Mean Core Length
	std::cout << "Overlap Length: " << format_double(mean(overlap_lengths, overlap_lengthsXL)) << std::endl;

	// Std Dev of Lengths
	std::cout << "StdDev Overlap Length: " << format_double(stdev(overlap_lengths, overlap_lengthsXL)) << std::endl;

	// Min/Max Core Length
	std::cout << "Min/Max Overlap Length: " << format_int(overlap_length_gaps.minimum) << "/" << format_int(overlap_length_gaps.maximum) << std::endl;

	// Total Sizes
	std::cout << "Total Size (GB): " << format_double(sizes / (1024.0 * 1024.0 * 1024.0)) << std::endl << std::endl;

	return 0;
};