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

#define BLEND_GET_KMER(minimizer) ((minimizer).x >> 14)
#define BLEND_GET_LENGTH(minimizer) ((minimizer).x & 0x3FFF)
#define BLEND_GET_INDEX(minimizer) (((minimizer).y & 0xFFFFFFFF) >> 1)
#define BLEND_GET_REFERENCE_IDX(minimizer) (((minimizer).y >> 32))
#define BLEND_GET_STRAND(minimizer) (((minimizer).y & 1))

typedef struct {
	uint64_t x; // kmer(32) + reference(32)
	uint64_t y; // length(31) + strand(1) + index(32)
} uint128_t;

/**
 * Find (w,k)-minimizers on a DNA sequence and BLEND n_neighbors many consecutive minimizers
 *
 * @param str           DNA sequence
 * @param len           length of $str
 * @param window        find a BLEND value for every $w consecutive k-mers
 * @param kmer_size     k-mer size
 * @param blend_bits    use blend_bits many bits when generating the hash values of seeds
 * @param n_neighbors   How many neighbors consecutive minimizer k-mers should be combined to generate a strobemer seed
 * @param rid           reference ID; will be copied to the output $p array
 * @param fuzzy_seeds
 */
uint64_t blend_sb_sketch(const char *str, int len, int window, int kmer_size, int blend_bits, uint64_t n_neighbors, uint32_t rid, uint128_t **fuzzy_seeds);

uint64_t blend_sketch(const char *str, int len, int window, int kmer_size, int blend_bits, int n_neighbors, uint32_t rid, uint128_t **fuzzy_seeds);

#endif