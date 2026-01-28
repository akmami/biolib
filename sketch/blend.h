#ifndef __BLEND_SKETCH__
#define __BLEND_SKETCH__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <x86intrin.h>

#define __blend_get_kmer(minimizer) 		((minimizer).x >> 14)
#define __blend_get_length(minimizer) 		((minimizer).x & 0x3FFF)
#define __blend_get_reference_id(minimizer) (((minimizer).y >> 32))
#define __blend_get_index(minimizer) 		(((minimizer).y & 0xFFFFFFFF) >> 1)
#define __blend_get_strand(minimizer) 		(((minimizer).y & 1))

#ifndef __UINT128_T__
#define __UINT128_T__
typedef struct {
	uint64_t x; // kmer(50) + span(14)
	uint64_t y; // reference_id(32) + index(31) + strand(1)
} uint128_t;
#endif

/**
 * Find (w,k)-minimizers on a DNA sequence and BLEND neighbor_window many consecutive minimizers
 *
 * @param str               sequence to be processed
 * @param str_len           length of `str`
 * @param window_size       find a BLEND value for every `window` consecutive k-mers
 * @param kmer_size         k-mer size
 * @param blend_bits        use blend_bits many bits when generating the hash values of seeds
 * @param neighbor_window   How many neighbors consecutive minimizer k-mers should be combined to generate a strobemer seed
 * @param str_id            string ID, will be encoded to output `fuzzy_seeds` array
 * @param fuzzy_seeds       the array of fuzzy seeds being processed and found
 * 
 * @return                  Number of fuzzy_seeds being detected
 */
uint64_t blend_sb_sketch(const char *str, int str_len, int window_size, int kmer_size, int blend_bits, uint64_t neighbor_window, uint32_t str_id, uint128_t **fuzzy_seeds);

/**
 * Find (w,k)-minimizers on a DNA sequence and BLEND neighbor_window many consecutive minimizers
 *
 * @param str               sequence to be processed
 * @param str_len           length of `str`
 * @param window_size       find a BLEND value for every `window` consecutive k-mers
 * @param kmer_size         k-mer size
 * @param blend_bits        use blend_bits many bits when generating the hash values of seeds
 * @param neighbor_window   How many neighbors consecutive minimizer k-mers should be combined to generate a strobemer seed
 * @param str_id            string ID, will be encoded to output `fuzzy_seeds` array
 * @param fuzzy_seeds       the array of fuzzy seeds being processed and found
 * 
 * @return                  Number of fuzzy_seeds being detected
 */
uint64_t blend_sketch(const char *str, int str_len, int window_size, int kmer_size, int blend_bits, int neighbor_window, uint32_t str_id, uint128_t **fuzzy_seeds);

#ifdef __cplusplus
}
#endif

#endif