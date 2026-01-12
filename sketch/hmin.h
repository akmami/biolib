#ifndef __HIERARCHICAL_MINIMIZERS__
#define __HIERARCHICAL_MINIMIZERS__

#ifdef __cplusplus
extern "C" {
#endif

#include "struct_def.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>


uint64_t hmin_sketch(const char *str, int len, int window, int kmer_size, int only_symmetric, uint32_t reference_id, uint128_t **minimizers) {

	if (3 * len / window == 0) return 0;

	uint64_t shift = 2 * (kmer_size - 1); 					// helper shift for to skip right k - 1 chars in kmer
	uint64_t kmer_bit_mask = (1ULL << (2 * kmer_size)) - 1; // helper mask for removing bits outside of 2 * k right most bits
	uint64_t current_kmer[2];	
	int char_index, helper_index, minimizers_buffer_pos, min_pos, kmer_len, kmer_span = 0;
	uint128_t minimizers_buffer[256];
	uint128_t min = { UINT64_MAX, UINT64_MAX };

	assert(len > 0 && (window > 0 && window < 256) && kmer_size > 0); 
	memset(minimizers_buffer, 0xff, 256 * sizeof(uint128_t));
	memset(current_kmer, 0, 16);

	kmer_len = 0;
	minimizers_buffer_pos = 0;
	min_pos = 0;

	*minimizers = (uint128_t *)malloc(sizeof(uint128_t) * 3 * len / window);

	uint128_t *kmers = *minimizers;
	uint64_t kmers_len = 0;


	for (char_index = 0; char_index < len; ++char_index) {
		int current_char = seq_nt4_table[(uint8_t)str[char_index]];
		uint128_t current_info = { UINT64_MAX, UINT64_MAX };

		// not an ambiguous base
		if (current_char < 4) { 
			uint64_t strand;
			kmer_span = kmer_len + 1 < kmer_size ? kmer_len + 1 : kmer_size;
			current_kmer[0] = (current_kmer[0] << 2 | current_char) & kmer_bit_mask;    // forward k-mer
			current_kmer[1] = (current_kmer[1] >> 2) | (3ULL^current_char) << shift; 	// reverse k-mer
			if (only_symmetric && current_kmer[0] == current_kmer[1]) continue; 		// skip "symmetric k-mers"
			strand = current_kmer[0] < current_kmer[1] ? 0 : 1; 	// strand
			kmer_len++;
			if (kmer_size <= kmer_len && kmer_span < 256) {
				SET_128_X(current_info, MurmurHash3_32(current_kmer + strand, 4, 32), reference_id);
				SET_128_Y(current_info, kmer_span, strand, char_index - kmer_span + 1);
			}
		} else {
			kmer_len = 0, kmer_span = 0;
		}

		minimizers_buffer[minimizers_buffer_pos] = current_info; // need to do this here as appropriate minimizers_buffer_pos and minimizers_buffer[minimizers_buffer_pos] are needed below
		
        // special case for the first window - because identical k-mers are not stored yet
		if (kmer_len == window + kmer_size - 1 && min.x != UINT64_MAX) { 
			for (helper_index = minimizers_buffer_pos + 1; helper_index < window; ++helper_index) {
				if (min.x == minimizers_buffer[helper_index].x && minimizers_buffer[helper_index].y != min.y) {
					kmers[kmers_len++] = minimizers_buffer[helper_index];
				}
			}
			for (helper_index = 0; helper_index < minimizers_buffer_pos; ++helper_index) {
				if (min.x == minimizers_buffer[helper_index].x && minimizers_buffer[helper_index].y != min.y) {
					kmers[kmers_len++] = min;
				}
			}
		}
		
		if (current_info.x <= min.x) { // a new minimum; then write the old min
			if (kmer_len >= window + kmer_size && min.x != UINT64_MAX) {
				kmers[kmers_len++] = min;
			}
			min = current_info, min_pos = minimizers_buffer_pos;
		} else if (minimizers_buffer_pos == min_pos) { // old min has moved outside the window
			if (kmer_len >= window + kmer_size - 1 && min.x != UINT64_MAX) {
				kmers[kmers_len++] = min;
			}
            // the two loops are necessary when there are identical k-mers
			for (helper_index = minimizers_buffer_pos + 1, min.x = UINT64_MAX; helper_index < window; ++helper_index) { 
				if (min.x >= minimizers_buffer[helper_index].x) { // >= is important s.t. min is always the closest k-mer
					min = minimizers_buffer[helper_index];
					min_pos = helper_index; 
				}
			}
			for (helper_index = 0; helper_index <= minimizers_buffer_pos; ++helper_index) {
				if (min.x >= minimizers_buffer[helper_index].x) {
					min = minimizers_buffer[helper_index];
					min_pos = helper_index;
				}
			}
			if (kmer_len >= window + kmer_size - 1 && min.x != UINT64_MAX) { // write identical k-mers
                // these two loops make sure the output is sorted
				for (helper_index = minimizers_buffer_pos + 1; helper_index < window; ++helper_index) { 
					if (min.x == minimizers_buffer[helper_index].x && min.y != minimizers_buffer[helper_index].y) {
						kmers[kmers_len++] = minimizers_buffer[helper_index];
					}
				}
				for (helper_index = 0; helper_index <= minimizers_buffer_pos; ++helper_index) {
					if (min.x == minimizers_buffer[helper_index].x && min.y != minimizers_buffer[helper_index].y) {
						kmers[kmers_len++] = minimizers_buffer[helper_index];
					}
				}
			}
		}

		if (++minimizers_buffer_pos == window) {
			minimizers_buffer_pos = 0;
		}
	}
	if (min.x != UINT64_MAX) {
		kmers[kmers_len++] = min;
	}

	if (kmers_len == 0) free(kmers);
	
	return kmers_len;
}

uint64_t hmin_hi_sketch(uint128_t *minimizers, int len, int window, int kmer_size, int only_symmetric, uint32_t reference_id, uint128_t **h_minimizers) {
	
	if (((7 * len / window) >> 1) == 0) return 0;
	
	uint64_t shift = 2 * (kmer_size - 1);
	uint64_t kmer_bit_mask = (1ULL << (2 * kmer_size)) - 1;
	uint64_t current_kmer[2];
	int minimizer_index, helper_index, minimizers_buffer_pos, min_pos, kmer_len, kmer_span;
	uint128_t minimizers_buffer[256];
	uint128_t min = { UINT64_MAX, UINT64_MAX };

	memset(minimizers_buffer, 0xff, 256 * sizeof(uint128_t));
	memset(current_kmer, 0, 16);

	kmer_len = 0;
	minimizers_buffer_pos = 0;
	min_pos = 0;

	*h_minimizers = (uint128_t *)malloc(sizeof(uint128_t) * ((7 * len / window) >> 1));

	uint128_t *kmers = *h_minimizers;
	uint64_t kmers_len = 0;

	for (minimizer_index = 0; minimizer_index < len; ++minimizer_index) {

		uint64_t current_char = GET_128_KMER(minimizers[minimizer_index]) & 3;
		uint128_t current_info = { UINT64_MAX, UINT64_MAX };

		uint64_t strand;
		kmer_span = kmer_len + 1 < kmer_size ? kmer_len + 1 : kmer_size;
		current_kmer[0] = (current_kmer[0] << 2 | current_char) & kmer_bit_mask;
		current_kmer[1] = (current_kmer[1] >> 2) | (current_char << shift);
		if (only_symmetric && current_kmer[0] == current_kmer[1]) continue;
		strand = current_kmer[0] < current_kmer[1] ? 0 : 1;
		kmer_len++;
		if (kmer_size <= kmer_len && kmer_span < 256) {
			uint64_t hkmer_start = GET_128_INDEX(minimizers[minimizer_index-kmer_size+1]);
			uint64_t hkmer_end = GET_128_INDEX(minimizers[minimizer_index]) + GET_128_LEN(minimizers[minimizer_index]);
			uint64_t hkmer_len = hkmer_end - hkmer_start;

			SET_128_X(current_info, MurmurHash3_32(current_kmer + strand, 4, 32), reference_id);
			SET_128_Y(current_info, hkmer_len, strand, hkmer_start);
		}

		minimizers_buffer[minimizers_buffer_pos] = current_info;
		
		if (kmer_len == window + kmer_size - 1 && min.x != UINT64_MAX) {
			for (helper_index = minimizers_buffer_pos + 1; helper_index < window; ++helper_index) {
				if (min.x == minimizers_buffer[helper_index].x && minimizers_buffer[helper_index].y != min.y) {
					kmers[kmers_len++] = minimizers_buffer[helper_index];
				}
			}
			for (helper_index = 0; helper_index < minimizers_buffer_pos; ++helper_index) {
				if (min.x == minimizers_buffer[helper_index].x && minimizers_buffer[helper_index].y != min.y) {
					kmers[kmers_len++] = min;
				}
			}
		}
		
		if (current_info.x <= min.x) {
			if (kmer_len >= window + kmer_size && min.x != UINT64_MAX) {
				kmers[kmers_len++] = min;
			}
			min = current_info, min_pos = minimizers_buffer_pos;
		} else if (minimizers_buffer_pos == min_pos) {
			if (kmer_len >= window + kmer_size - 1 && min.x != UINT64_MAX) {
				kmers[kmers_len++] = min;
			}
			for (helper_index = minimizers_buffer_pos + 1, min.x = UINT64_MAX; helper_index < window; ++helper_index) {
				if (min.x >= minimizers_buffer[helper_index].x) {
					min = minimizers_buffer[helper_index];
					min_pos = helper_index; 
				}
			}
			for (helper_index = 0; helper_index <= minimizers_buffer_pos; ++helper_index) {
				if (min.x >= minimizers_buffer[helper_index].x) {
					min = minimizers_buffer[helper_index];
					min_pos = helper_index;
				}
			}
			if (kmer_len >= window + kmer_size - 1 && min.x != UINT64_MAX) {
				for (helper_index = minimizers_buffer_pos + 1; helper_index < window; ++helper_index) {
					if (min.x == minimizers_buffer[helper_index].x && min.y != minimizers_buffer[helper_index].y) {
						kmers[kmers_len++] = minimizers_buffer[helper_index];
					}
				}
				for (helper_index = 0; helper_index <= minimizers_buffer_pos; ++helper_index) {
					if (min.x == minimizers_buffer[helper_index].x && min.y != minimizers_buffer[helper_index].y) {
						kmers[kmers_len++] = minimizers_buffer[helper_index];
					}
				}
			}
		}

		if (++minimizers_buffer_pos == window) {
			minimizers_buffer_pos = 0;
		}
	}
	if (min.x != UINT64_MAX) {
		kmers[kmers_len++] = min;
	}

	if (kmers_len == 0) free(kmers);

	return kmers_len;
}

#ifdef __cplusplus
}
#endif

#endif