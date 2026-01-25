#include "blend.h"


#ifdef __cplusplus
extern "C" {
#endif

static const uint64_t mash_lsb_32     = (1ULL << 32) - 1;   // = 0xFFFFFFFF, get right-most 32 bits
static const uint64_t kmer_span_mask  = (1ULL << 14) - 1;   // = 0x3FFF, get kmer span
static const uint64_t kmer_index_mask = (1ULL << 31) - 1;   // for iMask, the actual value should be shifted right first        

static unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask) {
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

// Source: https://stackoverflow.com/a/21673221
static inline __m256i movemask_inverse(const uint32_t hash_value) {
    __m256i vmask = _mm256_set1_epi32(hash_value);
    
    const __m256i shuffle = _mm256_setr_epi64x(0x0000000000000000, 0x0101010101010101, 0x0202020202020202, 0x0303030303030303);
    vmask = _mm256_shuffle_epi8(vmask, shuffle);
    
    const __m256i bit_mask = _mm256_set1_epi64x(0x7fbfdfeff7fbfdfe);
    vmask = _mm256_or_si256(vmask, bit_mask);
    
    return _mm256_cmpeq_epi8(vmask, _mm256_set1_epi64x(-1));
}

static inline void calc_blend_simd(__m256i* blndcnt_lsb, __m256i* blndcnt_msb, __m256i* ma, __m256i* mb, uint64_t val, const uint64_t mask, const int bits) {

    (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*ma), (*mb), movemask_inverse(val&mask)));
    uint64_t blendVal = (uint64_t)_mm256_movemask_epi8((*blndcnt_lsb)) & mask;
    
    if(bits > 32) {
        (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb), _mm256_blendv_epi8((*ma), (*mb), movemask_inverse((val >> 32) & mask)));
        blendVal |= ((uint64_t)_mm256_movemask_epi8((*blndcnt_msb))) << 32;
    }
}

static inline uint64_t calc_blend_rm_simd(__m256i* blndcnt_lsb, __m256i* blndcnt_msb, __m256i* ma, __m256i* mb, uint64_t val, uint64_t remval, const uint64_t mask, const int bits) {

    (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*ma), (*mb), movemask_inverse(val & mask)));
    uint64_t blendVal = (uint64_t)_mm256_movemask_epi8((*blndcnt_lsb))&mask;
    
    // Removal of the oldest item
    (*blndcnt_lsb) = _mm256_adds_epi8((*blndcnt_lsb), _mm256_blendv_epi8((*mb), (*ma), movemask_inverse(remval & mask)));

    if(bits > 32) {
        (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb), _mm256_blendv_epi8((*ma), (*mb),movemask_inverse((val >> 32) & mask)));
        blendVal |= ((uint64_t)_mm256_movemask_epi8((*blndcnt_msb)))<<32;

        (*blndcnt_msb) = _mm256_adds_epi8((*blndcnt_msb),_mm256_blendv_epi8((*mb), (*ma), movemask_inverse((remval >> 32) & mask)));
    }
    
    return blendVal;
}

#define process_buffer(blend_buffer, blend_buffer_pos, current_seed, valid_latest_seed, n_neighbors, processed_neighbors, blend_count_lsb, blend_count_msb, ma, mb, blend_value, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, debug) \
\
    blend_buffer[blend_buffer_pos].x = current_seed.x; \
    blend_buffer[blend_buffer_pos].y = (current_seed.y & kmer_index_mask) >> 1; \
    uint32_t tot_span = blend_buffer[blend_buffer_pos].y; \
\
    if (++blend_buffer_pos == n_neighbors) { \
        blend_buffer_pos = 0; \
    } \
\
    if (++processed_neighbors >= n_neighbors) { \
        blendVal = calc_blend_rm_simd(&blend_count_lsb, &blend_count_msb, &ma, &mb, current_seed.x >> 14, blend_buffer[blend_buffer_pos].x >> 14, mash_lsb_32, blend_bits_count); \
        tot_span = tot_span - blend_buffer[blend_buffer_pos].y + (blend_buffer[blend_buffer_pos].x & kmer_span_mask); \
        valid_latest_seed.x = blend_value << 14 | tot_span; \
        valid_latest_seed.y = current_seed.y; \
\
        temp_fuzzy_seeds[seed_count++] = valid_latest_seed; \
        if (seed_count == seed_capacity) { \
            seed_capacity *= 2; \
            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, seed_capacity * sizeof(uint128_t)); \
            temp_fuzzy_seeds = *fuzzy_seeds; \
            if (!temp_fuzzy_seeds) return 0; \
        } \
    } else { /* only addition of current_seed */ \
        calc_blend_simd(&blend_count_lsb, &blend_count_msb, &ma, &mb, current_seed.x >> 14, mash_lsb_32, blend_bits_count); \
    } \

/**
 * Find (w,k)-minimizers on a DNA sequence and BLEND n_neighbors many consecutive minimizers
 *
 * @param str           sequence to be processed
 * @param len           length of `str`
 * @param window        find a BLEND value for every `window` consecutive k-mers
 * @param kmer_size     k-mer size
 * @param blend_bits    use blend_bits many bits when generating the hash values of seeds
 * @param n_neighbors   How many neighbors consecutive minimizer k-mers should be combined to generate a strobemer seed
 * @param str_id        string ID, will be encoded to output `fuzzy_seeds` array
 * @param fuzzy_seeds   the array of fuzzy seeds being processed and found
 * 
 * @return              Number of fuzzy_seeds being detected
 */
uint64_t blend_sb_sketch(const char *str, int len, int window, int kmer_size, int blend_bits, uint64_t n_neighbors, uint32_t str_id, uint128_t **fuzzy_seeds) {
    
    assert(len > 0 && (window > 0 && window + kmer_size < 8192) && (kmer_size > 0 && kmer_size <= 28) && (n_neighbors > 0 && n_neighbors + kmer_size < 8192) && (blend_bits <= 56));
    
    const int blend_bits_count = (blend_bits > 0) ? blend_bits : 2 * kmer_size;
    const uint64_t blend_bits_mask = (1ULL << blend_bits_count) - 1; // mash to get BLEND bits

    const uint64_t rc_shift_bits = 2 * (kmer_size- 1);
    const uint64_t kmer_mask = (1ULL << 2 * kmer_size) - 1;          // mask to get kmer (2 * kmer-size number of ones)

    uint64_t kmer[2] = {0, 0};
    int char_index = 0, temp_seed_buffer_pos = 0, current_kmer_span = 0, seed_buffer_pos = 0, min_pos = 0, kmer_span = 0;

    const static int seed_buffer_len = 256;
    uint128_t *seed_buffer = (uint128_t *)malloc(sizeof(uint128_t) * seed_buffer_len);
    uint128_t min = { UINT64_MAX, UINT64_MAX };
    
    // BLEND Variables
    uint64_t blendVal = 0;
    uint128_t *forward_blend_buffer = (uint128_t*)malloc(sizeof(uint128_t) * n_neighbors);
    uint128_t *reverse_blend_buffer = (uint128_t*)malloc(sizeof(uint128_t) * n_neighbors);
    uint64_t forward_blend_pos = 0, reverse_blend_pos = 0;
    uint64_t f_nm = 0, r_nm = 0; // number of minimizers processed
    
    // SSE4-related variables
    __m256i ma = _mm256_set1_epi8(1);
    __m256i mb = _mm256_set1_epi8(-1);
    __m256i f_blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i f_blndcnt_msb = _mm256_set1_epi8(0);
    __m256i r_blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i r_blndcnt_msb = _mm256_set1_epi8(0);
    
    memset(seed_buffer, 0xff, sizeof(uint128_t) * seed_buffer_len);
    memset(forward_blend_buffer, 0, sizeof(uint128_t) * n_neighbors);
    memset(reverse_blend_buffer, 0, sizeof(uint128_t) * n_neighbors);

    uint64_t seed_capacity = 3 * len / window;
    *fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * seed_capacity);
    uint128_t *temp_fuzzy_seeds = *fuzzy_seeds;
    uint64_t seed_count = 0;

    if (!seed_buffer || !forward_blend_buffer || !reverse_blend_buffer || !fuzzy_seeds) {
        fprintf(stderr, "Poor programming skills. Please curse akmami...\n");
        abort();
    }

    for ( ; char_index < len; char_index++ ) {
        
        int ch = seq_nt4_table[(uint8_t)str[char_index]];
        
        uint128_t info = { UINT64_MAX, UINT64_MAX };
        uint128_t f_info = { UINT64_MAX, UINT64_MAX };
        
        if (ch < 4) { // skip an ambiguous base
            int strand;
            kmer_span = current_kmer_span + 1 < kmer_size ? current_kmer_span + 1 : kmer_size;
            
            kmer[0] = ((kmer[0] << 2 | ch) & kmer_mask); // forward k-mer
            kmer[1] = ((kmer[1] >> 2) | (3ULL ^ ch) << rc_shift_bits); // reverse k-mer k-mer
            
            if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know its strand
            
            strand = kmer[0] < kmer[1] ? 0 : 1; // strand
            current_kmer_span++;
            if (current_kmer_span >= kmer_size && kmer_span < 256) {
                f_info.x = hash64(kmer[strand], blend_bits_mask) << 14 | kmer_span;
                f_info.y = (uint64_t)str_id << 32 | (uint32_t)char_index << 1 | strand;
            }
        } else {
            current_kmer_span = 0;
            kmer_span = 0;
        }
        seed_buffer[seed_buffer_pos] = f_info; // need to do this here as appropriate seed_buffer[seed_buffer_pos] are needed below
        
        // special case for the first window - because identical k-mers are not stored yet
        if (current_kmer_span == window + kmer_size - 1 && min.x != UINT64_MAX) {
            for (temp_seed_buffer_pos = seed_buffer_pos + 1; temp_seed_buffer_pos < seed_buffer_len; temp_seed_buffer_pos++) {
                if (min.x == seed_buffer[temp_seed_buffer_pos].x && seed_buffer[temp_seed_buffer_pos].y != min.y) {
                    if ((seed_buffer[temp_seed_buffer_pos].y & 1) == 1) { // reverse strand
                        process_buffer(reverse_blend_buffer, reverse_blend_pos, seed_buffer[temp_seed_buffer_pos], info, n_neighbors, r_nm, r_blndcnt_lsb, r_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 218);
                    } else {
                        process_buffer(forward_blend_buffer, forward_blend_pos, seed_buffer[temp_seed_buffer_pos], info, n_neighbors, f_nm, f_blndcnt_lsb, f_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 220);
                    }
                }
            }

            for (temp_seed_buffer_pos = 0; temp_seed_buffer_pos < seed_buffer_pos; temp_seed_buffer_pos++)
                if (min.x == seed_buffer[temp_seed_buffer_pos].x && seed_buffer[temp_seed_buffer_pos].y != min.y) {
                    if ((seed_buffer[temp_seed_buffer_pos].y & 1) == 1) { // reverse strand
                        process_buffer(reverse_blend_buffer, reverse_blend_pos, seed_buffer[temp_seed_buffer_pos], info, n_neighbors, r_nm, r_blndcnt_lsb, r_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 228);
                    } else {
                        process_buffer(forward_blend_buffer, forward_blend_pos, seed_buffer[temp_seed_buffer_pos], info, n_neighbors, f_nm, f_blndcnt_lsb, f_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 230);
                    }
                }
        }

        if (f_info.x <= min.x) { // a new minimum; then write the old min
            if (current_kmer_span >= window + kmer_size && min.x != UINT64_MAX) {
                if ((min.y & 1) == 1) { // reverse strand
                    process_buffer(reverse_blend_buffer, reverse_blend_pos, min, info, n_neighbors, r_nm, r_blndcnt_lsb, r_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 238);
                } else {
                    process_buffer(forward_blend_buffer, forward_blend_pos, min, info, n_neighbors, f_nm, f_blndcnt_lsb, f_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 240);
                }
            }
            min = f_info, min_pos = seed_buffer_pos;
        } else if (seed_buffer_pos == min_pos) { // old min has moved outside the window
            if (current_kmer_span >= window + kmer_size - 1 && min.x != UINT64_MAX) {
                if ((min.y & 1) == 1) { // reverse strand
                    process_buffer(reverse_blend_buffer, reverse_blend_pos, min, info, n_neighbors, r_nm, r_blndcnt_lsb, r_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 247);
                } else {
                    process_buffer(forward_blend_buffer, forward_blend_pos, min, info, n_neighbors, f_nm, f_blndcnt_lsb, f_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 249);
                }
            }
            for (temp_seed_buffer_pos = seed_buffer_pos + 1, min.x = UINT64_MAX; temp_seed_buffer_pos < seed_buffer_len; temp_seed_buffer_pos++) { // the two loops are necessary when there are identical k-mers. This one starts from the oldest (first in: seed_buffer_pos + 1) 
                if (min.x >= seed_buffer[temp_seed_buffer_pos].x) {
                    min = seed_buffer[temp_seed_buffer_pos], min_pos = temp_seed_buffer_pos; // >= is important s.t. min is always the closest k-mer... chosing the last one because we will store *all* identical minimizers and continue with the last one once all identical ones are stored. @IMPORTANT: should we really store all identical k-mers? it is true that they are the minimizers for their own window but not sure how effective it is
                }
            }
            for (temp_seed_buffer_pos = 0; temp_seed_buffer_pos <= seed_buffer_pos; temp_seed_buffer_pos++) {
                if (min.x >= seed_buffer[temp_seed_buffer_pos].x) {
                    min = seed_buffer[temp_seed_buffer_pos], min_pos = temp_seed_buffer_pos;
                }
            }
            if (current_kmer_span >= window + kmer_size - 1 && min.x != UINT64_MAX) { // write identical k-mers
                for (temp_seed_buffer_pos = seed_buffer_pos + 1; temp_seed_buffer_pos < seed_buffer_len; temp_seed_buffer_pos++) { // these two loops make sure the output is sorted
                    if (min.x == seed_buffer[temp_seed_buffer_pos].x && min.y != seed_buffer[temp_seed_buffer_pos].y) {
                        if ((seed_buffer[temp_seed_buffer_pos].y & 1) == 1) { // reverse strand
                            process_buffer(reverse_blend_buffer, reverse_blend_pos, seed_buffer[temp_seed_buffer_pos], info, n_neighbors, r_nm, r_blndcnt_lsb, r_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 266);
                        } else {
                            process_buffer(forward_blend_buffer, forward_blend_pos, seed_buffer[temp_seed_buffer_pos], info, n_neighbors, f_nm, f_blndcnt_lsb, f_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 268);
                        }
                    }
                }
                for (temp_seed_buffer_pos = 0; temp_seed_buffer_pos <= seed_buffer_pos; temp_seed_buffer_pos++)
                    if (min.x == seed_buffer[temp_seed_buffer_pos].x && min.y != seed_buffer[temp_seed_buffer_pos].y) {
                        if ((seed_buffer[temp_seed_buffer_pos].y & 1) == 1) { // reverse strand
                            process_buffer(reverse_blend_buffer, reverse_blend_pos, seed_buffer[temp_seed_buffer_pos], info, n_neighbors, r_nm, r_blndcnt_lsb, r_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 275);
                        } else {
                            process_buffer(forward_blend_buffer, forward_blend_pos, seed_buffer[temp_seed_buffer_pos], info, n_neighbors, f_nm, f_blndcnt_lsb, f_blndcnt_msb, ma, mb, blendVal, blend_bits_count, temp_fuzzy_seeds, seed_count, seed_capacity, len, 277);
                        }
                    }
            }
        }

        if (++seed_buffer_pos == window) {
            seed_buffer_pos = 0;
        }
    }

    if (min.x != UINT64_MAX) {
        if ((min.y & 1) == 1) { // Reverse strand
            reverse_blend_buffer[reverse_blend_pos].x = min.x;
            reverse_blend_buffer[reverse_blend_pos].y = (min.y & kmer_index_mask) >> 1;
            uint32_t tot_span = reverse_blend_buffer[reverse_blend_pos].y;

            if (++reverse_blend_pos == n_neighbors) {
                reverse_blend_pos = 0;
            }

            if (++r_nm >= n_neighbors) {
                blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, reverse_blend_buffer[reverse_blend_pos].x >> 14, mash_lsb_32, blend_bits_count);
                tot_span = tot_span - reverse_blend_buffer[reverse_blend_pos].y + (reverse_blend_buffer[reverse_blend_pos].x & kmer_span_mask);
                min.x = blendVal << 14 | tot_span;

                temp_fuzzy_seeds[seed_count++] = min;
                if (seed_count == seed_capacity) {
                    seed_capacity *= 2;
                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, seed_capacity * sizeof(uint128_t));
                    temp_fuzzy_seeds = *fuzzy_seeds;
                    if (!temp_fuzzy_seeds) return 0;
                }
            } else { // only addition of min
                calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, mash_lsb_32, blend_bits_count);
            }
        } else {
            forward_blend_buffer[forward_blend_pos].x = min.x;
            forward_blend_buffer[forward_blend_pos].y = (min.y & kmer_index_mask) >> 1;
            uint32_t tot_span = forward_blend_buffer[forward_blend_pos].y;

            if(++forward_blend_pos == n_neighbors) forward_blend_pos = 0;
            if(++f_nm >= n_neighbors) {
                blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, forward_blend_buffer[forward_blend_pos].x >> 14, mash_lsb_32, blend_bits_count);
                tot_span = tot_span - forward_blend_buffer[forward_blend_pos].y + (forward_blend_buffer[forward_blend_pos].x & kmer_span_mask);
                min.x = blendVal << 14 | tot_span;

                temp_fuzzy_seeds[seed_count++] = min;
                if (seed_count == seed_capacity) {
                    seed_capacity *= 2;
                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, seed_capacity * sizeof(uint128_t));
                    temp_fuzzy_seeds = *fuzzy_seeds;
                    if (!temp_fuzzy_seeds) return 0;
                }
            } else { // only addition of min
                calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, mash_lsb_32, blend_bits_count);
            }
        }
    }

    free(seed_buffer);
    free(forward_blend_buffer);
    free(reverse_blend_buffer);

    return seed_count;
}

uint64_t blend_sketch(const char *str, int len, int window, int kmer_size, int blend_bits, int n_neighbors, uint32_t str_id, uint128_t **fuzzy_seeds) {
	return blend_sb_sketch(str, len, window, kmer_size, blend_bits, (uint64_t)n_neighbors, str_id, fuzzy_seeds);
}

#ifdef __cplusplus
}
#endif