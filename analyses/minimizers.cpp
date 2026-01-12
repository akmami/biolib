#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <math.h>
#include "../sketch/hmin.h"
#include "utils.h"


#define LEVEL 5
#define ONLY_CANONICAL 0

void analyze(uint128_t *minimizers, uint64_t len,
             int level,
             int (&contiguous_counts)[LEVEL],
             int (&distances)[LEVEL][DISTANCE_LENGTH],
             std::vector<std::vector<uint32_t>> &distancesXL,
			 gap (&distance_gaps)[LEVEL],
             int (&lengths)[LEVEL][DISTANCE_LENGTH],
             std::vector<std::vector<uint32_t>> &lengthsXL,
			 gap (&length_gaps)[LEVEL],
			 int  (&overlap_lengths)[LEVEL][DISTANCE_LENGTH], 
			 std::vector<std::vector<uint32_t>> &overlap_lengthsXL, 
			 gap (&overlap_length_gaps)[LEVEL]) {

	if (len > 0) {

		bool isOverlapped = false;

		if (GET_128_LEN(minimizers[0]) < DISTANCE_LENGTH) {
			lengths[level][GET_128_LEN(minimizers[0])] += 1;
		} else {
			lengthsXL.at(level).push_back(GET_128_LEN(minimizers[0]));
		}

		for (uint64_t i = 1; i < len; i++) {

			uint32_t current_start = GET_128_INDEX(minimizers[i]);
			uint32_t current_end = current_start + GET_128_LEN(minimizers[i]);
			uint32_t previous_start = GET_128_INDEX(minimizers[i-1]);
			uint32_t previous_end = previous_start + GET_128_LEN(minimizers[i-1]);

			if (current_start <= previous_end) {
				contiguous_counts[level] += 1;
				isOverlapped = true;
			}

			// Distance
			uint32_t distance = current_start - previous_start;
			if (distance < DISTANCE_LENGTH ) {
				distances[level][distance]++;
			} else {
				distancesXL.at(level).push_back(distance);
			}

			distance_gaps[level].minimum = min(distance_gaps[level].minimum, distance);
			distance_gaps[level].maximum = max(distance_gaps[level].maximum, distance);

			// Length
			uint32_t length = current_end - current_start;
			if (length < DISTANCE_LENGTH) {
				lengths[level][length] += 1;
			} else {
				lengthsXL.at(level).push_back(length);
			}

			length_gaps[level].minimum = min(length_gaps[level].minimum, length);
			length_gaps[level].maximum = max(length_gaps[level].maximum, length);

			// Overlap
			if (previous_end >= current_start) {
				uint32_t overlap = previous_end - current_start;
				if (overlap < DISTANCE_LENGTH) {
					overlap_lengths[level][overlap] += 1;
				} else {
					overlap_lengthsXL.at(level).push_back(overlap);
				}
				overlap_length_gaps[level].minimum = min(overlap_length_gaps[level].minimum, overlap);
				overlap_length_gaps[level].maximum = max(overlap_length_gaps[level].maximum, overlap);			
			}
		}

		if (isOverlapped) {
			contiguous_counts[level] += 1;
		}
	}
};

void process(std::string &sequence,
             int window, int kmer_size, int only_symmetric,
             int (&core_counts)[LEVEL],
             int (&contiguous_counts)[LEVEL],
             std::set<uint32_t> (&distinct_cores)[LEVEL],
             std::vector<std::chrono::milliseconds> &durations,
             int (&distances)[LEVEL][DISTANCE_LENGTH],
             std::vector<std::vector<uint32_t>> &distancesXL,
			 gap (&distance_gaps)[LEVEL],
             int (&lengths)[LEVEL][DISTANCE_LENGTH],
             std::vector<std::vector<uint32_t>> &lengthsXL,
			 gap (&length_gaps)[LEVEL],
			 int  (&overlap_lengths)[LEVEL][DISTANCE_LENGTH], 
			 std::vector<std::vector<uint32_t>> &overlap_lengthsXL, 
			 gap (&overlap_length_gaps)[LEVEL],
             double (&sizes)[LEVEL]) {

	auto start = std::chrono::high_resolution_clock::now();

    uint128_t *minimizers;
    uint64_t minimizers_len = 0;
    
    minimizers_len = hmin_sketch(sequence.c_str(), sequence.length(), window, kmer_size, only_symmetric, 0, &minimizers);

	auto extraction_end = std::chrono::high_resolution_clock::now();
	durations[0] += std::chrono::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(extraction_end - start).count());
	core_counts[0] += minimizers_len;
    sizes[0] += sizeof(uint128_t) * minimizers_len;
	
	analyze(minimizers, minimizers_len, 0, contiguous_counts, distances, distancesXL, distance_gaps, lengths, lengthsXL, length_gaps, overlap_lengths, overlap_lengthsXL, overlap_length_gaps);

    for (uint64_t index = 0; index < minimizers_len; index++) {
        distinct_cores[0].insert(GET_128_KMER(minimizers[index]));
    }

	for (int i = 1; i < LEVEL; i++) {

		auto start_level = std::chrono::high_resolution_clock::now();

        uint128_t *hi_minimizers;
		uint64_t hi_minimizers_len = 0;
        hi_minimizers_len = hmin_hi_sketch(minimizers, minimizers_len, window, kmer_size, only_symmetric, i, &hi_minimizers);

		auto stop_level = std::chrono::high_resolution_clock::now();
		durations[i] += std::chrono::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(stop_level - start_level).count());
		core_counts[i] += hi_minimizers_len;
        sizes[i] += sizeof(uint128_t) * hi_minimizers_len;

        analyze(hi_minimizers, hi_minimizers_len, i, contiguous_counts, distances, distancesXL, distance_gaps, lengths, lengthsXL, length_gaps, overlap_lengths, overlap_lengthsXL, overlap_length_gaps);

        for (uint64_t index = 0; index < hi_minimizers_len; index++) {
            distinct_cores[i].insert(GET_128_KMER(hi_minimizers[index]));
        }

        if (minimizers_len) free(minimizers);
		minimizers = hi_minimizers;
		minimizers_len = hi_minimizers_len;
	}

	std::cout << "Length of the processed sequence: " << format_int(sequence.size()) << std::endl;

	sequence.clear();
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
    int core_counts[LEVEL] = {0};
	int contiguous_counts[LEVEL] = {0};
    std::set<uint32_t> distinct_cores[LEVEL];
    std::vector<std::chrono::milliseconds> durations(LEVEL);
	durations.resize(LEVEL);
	// section 2
    int distances[LEVEL][DISTANCE_LENGTH] = {0};
    std::vector<std::vector<uint32_t>> distancesXL(LEVEL);
	distancesXL.resize(LEVEL);
	gap distance_gaps[LEVEL];
	for (int i = 0; i < LEVEL; i++) {
		distance_gaps[i].minimum = UINT64_MAX;
		distance_gaps[i].maximum = 0;
	}
	// section 3
	int lengths[LEVEL][DISTANCE_LENGTH] = {0};
	std::vector<std::vector<uint32_t>> lengthsXL(LEVEL);
	lengthsXL.resize(LEVEL);
	gap length_gaps[LEVEL];
	for (int i = 0; i < LEVEL; i++) {
		length_gaps[i].minimum = UINT64_MAX;
		length_gaps[i].maximum = 0;
	}
	// section 4
	int overlap_lengths[LEVEL][DISTANCE_LENGTH] = {0};
	std::vector<std::vector<uint32_t>> overlap_lengthsXL(LEVEL);
	overlap_lengthsXL.resize(LEVEL);
	gap overlap_length_gaps[LEVEL];
	for (int i = 0; i < LEVEL; i++) {
		overlap_length_gaps[i].minimum = UINT64_MAX;
		overlap_length_gaps[i].maximum = 0;
	}	
    // section 5
    double sizes[LEVEL] = {0, 0, 0, 0, 0}; // 0, 0, 0};
    
    // other
	size_t genome_size = 0;

	// read file
	if (genome.is_open()) {

		std::string sequence, id;
		sequence.reserve(250000000);

		std::cout << "Program begins" << std::endl;

		while (getline(genome, line)) {

			if (line[0] == '>') {

				// process previous chromosome before moving into new one
				if (sequence.size() != 0) {
					genome_size += sequence.size();
					process(sequence, window, kmer_size, ONLY_CANONICAL, core_counts, contiguous_counts, distinct_cores, durations, distances, distancesXL, distance_gaps, lengths, lengthsXL, length_gaps, overlap_lengths, overlap_lengthsXL, overlap_length_gaps, sizes);
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
			process(sequence, window, kmer_size, ONLY_CANONICAL, core_counts, contiguous_counts, distinct_cores, durations, distances, distancesXL, distance_gaps, lengths, lengthsXL, length_gaps, overlap_lengths, overlap_lengthsXL, overlap_length_gaps, sizes);
		}

		genome.close();
	}

	std::string sep = " & ";
	double previous, current;

	std::cout << "Hierarchy level";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << i + 1;
	}
	std::cout << std::endl;

	// Total Cores
	std::cout << "Total \\# Cores";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_int(core_counts[i]);
	}
	std::cout << " \\\\" << std::endl;

	// Contiguous Cores
	std::cout << "Contiguous Cores";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_int(contiguous_counts[i]);
	}
	std::cout << " \\\\" << std::endl;

	// Distinct Cores
    std::cout << "Unique Cores";
    for (int i = 0; i < LEVEL; i++) {
        std::cout << sep << format_int(distinct_cores[i].size());
    }
    std::cout << " \\\\" << std::endl;

	// Distinctness Ratio
    std::cout << "Distinctness (%)";
    for (int i = 0; i < LEVEL; i++) {
        std::cout << sep << format_double(((double)distinct_cores[i].size()) / ((double)contiguous_counts[i]));
    }
    std::cout << " \\\\" << std::endl;

	// Execution Time
	std::cout << "Exec. Time (sec)";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_double(((double)durations[i].count()) / 1000);
	}
	std::cout << " \\\\" << std::endl;
    std::cout << "\\midrule" << std::endl;

	// Mean Core Distances
	std::cout << "Avg Distance";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_double(mean(distances[i], distancesXL[i]));
	}
	std::cout << " \\\\" << std::endl;

	// Std Dev of Distances
	std::cout << "StdDev Distance";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_double(stdev(distances[i], distancesXL[i]));
	}
	std::cout << " \\\\" << std::endl;

	// Min/Max Core Length
	std::cout << "Min/Max Distance";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_int(distance_gaps[i].minimum) << "/" << format_int(distance_gaps[i].maximum);
	}
	std::cout << " \\\\" << std::endl;
    std::cout << "\\midrule" << std::endl;

	// Mean Core Length
	std::cout << "Avg Length";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_double(mean(lengths[i], lengthsXL[i]));
	}
	std::cout << " \\\\" << std::endl;

	// Std Dev of Lengths
	std::cout << "StdDev Length";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_double(stdev(lengths[i], lengthsXL[i]));
	}
	std::cout << " \\\\" << std::endl;

	// Min/Max Core Length
	std::cout << "Min/Max Length";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_int(length_gaps[i].minimum) << "/" << format_int(length_gaps[i].maximum);
	}
	std::cout << " \\\\" << std::endl;
    std::cout << "\\midrule" << std::endl;

	// Mean Core Length
	std::cout << "Overlap Length";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_double(mean(overlap_lengths[i], overlap_lengthsXL[i]));
	}
	std::cout << " \\\\" << std::endl;

	// Std Dev of Lengths
	std::cout << "StdDev Overlap Length";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_double(stdev(overlap_lengths[i], overlap_lengthsXL[i]));
	}
	std::cout << " \\\\" << std::endl;

	// Min/Max Core Length
	std::cout << "Min/Max Overlap Length";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_int(overlap_length_gaps[i].minimum) << "/" << format_int(overlap_length_gaps[i].maximum);
	}
	std::cout << " \\\\" << std::endl;
    std::cout << "\\midrule" << std::endl;

	// Decrease in Total Counts
	previous = genome_size;
	std::cout << "Decrease in Core Count";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_double(static_cast<double>(core_counts[i]) / previous);
		previous = static_cast<double>(core_counts[i]);
	}
	std::cout << " \\\\" << std::endl;

	// Increase in Mean Lengths
	previous = 1;
	std::cout << "Increase in Avg Length";
	for (int i = 0; i < LEVEL; i++) {
		current = mean(lengths[i], lengthsXL[i]);
		std::cout << sep << format_double(current / previous);
		previous = current;
	}
	std::cout << " \\\\" << std::endl;

	// Increase in Mean Overlap Lengths
	previous = 1;
	std::cout << "Increase in Overlap Length";
	for (int i = 0; i < LEVEL; i++) {
		current = mean(overlap_lengths[i], overlap_lengthsXL[i]);
		std::cout << sep << format_double(current / previous);
		previous = current;
	}
	std::cout << " \\\\" << std::endl;

	// Increase in Mean Distances
	previous = 1;
	std::cout << "Increase in Avg Distance";
	for (int i = 0; i < LEVEL; i++) {
		current = mean(distances[i], distancesXL[i]);
		std::cout << sep << format_double(current / previous);
		previous = current;
	}
	std::cout << " \\\\" << std::endl;
    std::cout << "\\midrule" << std::endl;

	// Total Sizes
	std::cout << "Total Size (GB)";
	for (int i = 0; i < LEVEL; i++) {
		std::cout << sep << format_double(sizes[i] / (1024.0 * 1024.0 * 1024.0));
	}
	std::cout << " \\\\" << std::endl;
    std::cout << "\\bottomrule" << std::endl << std::endl;

	return 0;
};