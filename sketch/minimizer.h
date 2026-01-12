#ifndef __MINIMIZERS__
#define __MINIMIZERS__

#ifdef __cplusplus
extern "C" {
#endif

#include "struct_def.h"
#include <assert.h>
#include <string.h>

/**
 * Find symmetric (w, k)-minimizers on a DNA sequence
 *
 * @param str    		 sequence
 * @param len    		 length of sequence
 * @param window      	 find a minimizer for every `window` consecutive k-mers
 * @param kmer_size      k-mer size
 * @param only_symmetric include only symmetric kmers
 * @param reference_id   reference ID; will be embeded to minimizer struct
 * 
 * @note kmer is encoded as:
 * 		minimizer.x = hash
 * 		minimizer.y = pos << 32 | len 
 * 		minimizer.z = ref << 1 | strand
 */
uint64_t sketch_minimizers(const char *str, int len, int window, int kmer_size, int only_symmetric, uint32_t reference_id, uint128_t **minimizers) {

	if (len <= 0) return 0;

	uint64_t shift = 2 * (kmer_size - 1);
	uint64_t kmer_bit_mask = (1ULL << (2 * kmer_size)) - 1;
	uint64_t current_kmer[2] = {0, 0};
	int char_index, helper_index, buf_pos, min_pos, kmer_len, kmer_span = 0;
	uint128_t buf[256];
	uint128_t min = { UINT64_MAX, UINT64_MAX };

	// legacy code
	assert(len > 0 && (window > 0 && window < 256) && (kmer_size > 0 && kmer_size <= 28)); 
	memset(buf, 0xff, window * 16);

	kmer_len = 0;
	buf_pos = 0;
	min_pos = 0;

	*minimizers = (uint128_t *)malloc(sizeof(uint128_t) * 3 * len / window);

	uint128_t *kmers = *minimizers;
	uint64_t kmers_len = 0;


	for (char_index = 0; char_index < len; ++char_index) {
		int current_char = seq_nt4_table[(uint8_t)str[char_index]];
		uint128_t current_info = { UINT64_MAX, UINT64_MAX };

		if (current_char < 4) { 
			int strand;
			kmer_span = kmer_len + 1 < kmer_size? kmer_len + 1 : kmer_size;
			current_kmer[0] = (current_kmer[0] << 2 | current_char) & kmer_bit_mask;    // forward k-mer
			current_kmer[1] = (current_kmer[1] >> 2) | (3ULL^current_char) << shift; 	// reverse k-mer
			if (only_symmetric && current_kmer[0] == current_kmer[1]) continue;
			strand = current_kmer[0] < current_kmer[1] ? 0 : 1;
			kmer_len++;
			if (kmer_len >= kmer_size && kmer_span < 256) {
				current_info.x = hash64(current_kmer[strand], kmer_bit_mask) << 8 | kmer_span;
				current_info.y = (uint64_t)reference_id<<32 | (uint32_t)(char_index-kmer_len+1) << 1 | strand;
			}
		} else {
			kmer_len = 0, kmer_span = 0;
		}

		buf[buf_pos] = current_info;
		
		if (kmer_len == window + kmer_size - 1 && min.x != UINT64_MAX) { 
			for (helper_index = buf_pos + 1; helper_index < window; ++helper_index) {
				if (min.x == buf[helper_index].x && buf[helper_index].y != min.y) {
					kmers[kmers_len++] = buf[helper_index];
				}
			}
			for (helper_index = 0; helper_index < buf_pos; ++helper_index) {
				if (min.x == buf[helper_index].x && buf[helper_index].y != min.y) {
					kmers[kmers_len++] = min;
				}
			}
		}
		
		if (current_info.x <= min.x) {
			if (kmer_len >= window + kmer_size && min.x != UINT64_MAX) {
				kmers[kmers_len++] = min;
			}
			min = current_info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) {
			if (kmer_len >= window + kmer_size - 1 && min.x != UINT64_MAX) {
				kmers[kmers_len++] = min;
			}
			for (helper_index = buf_pos + 1, min.x = UINT64_MAX; helper_index < window; ++helper_index) { 
				if (min.x >= buf[helper_index].x) {
					min = buf[helper_index];
					min_pos = helper_index; 
				}
			}
			for (helper_index = 0; helper_index <= buf_pos; ++helper_index) {
				if (min.x >= buf[helper_index].x) {
					min = buf[helper_index];
					min_pos = helper_index;
				}
			}
			if (kmer_len >= window + kmer_size - 1 && min.x != UINT64_MAX) {
				for (helper_index = buf_pos + 1; helper_index < window; ++helper_index) { 
					if (min.x == buf[helper_index].x && min.y != buf[helper_index].y) {
						kmers[kmers_len++] = buf[helper_index];
					}
				}
				for (helper_index = 0; helper_index <= buf_pos; ++helper_index) {
					if (min.x == buf[helper_index].x && min.y != buf[helper_index].y) {
						kmers[kmers_len++] = buf[helper_index];
					}
				}
			}
		}

		if (++buf_pos == window) {
			buf_pos = 0;
		}
	}
	if (min.x != UINT64_MAX) {
		kmers[kmers_len++] = min;
	}
	
	return kmers_len;
}

#ifdef __cplusplus
}
#endif

#endif