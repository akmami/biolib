#include "blend.h"


#ifdef __cplusplus
extern "C" {
#endif

static const uint64_t MASK_LS_32_BITS = (1ULL << 32) - 1;   // = 0xFFFFFFFF, get right-most 32 bits
static const uint64_t MASK_KMER_SPAN  = (1ULL << 14) - 1;   // = 0x3FFF, get kmer span
static const uint64_t MASK_KMER_INDEX = ((1ULL << 31) - 1) << 1;   // for iMask, the actual value should be shifted right first        

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

static inline void calc_blend_simd(__m256i *blend_count_lsb, __m256i *blend_count_msb, __m256i *vote_plus, __m256i *vote_minus, uint64_t new_val, const uint64_t mask, const int bits) {

    (*blend_count_lsb) = _mm256_adds_epi8((*blend_count_lsb), _mm256_blendv_epi8((*vote_plus), (*vote_minus), movemask_inverse(new_val & mask)));
    uint64_t blend_value = (uint64_t)_mm256_movemask_epi8((*blend_count_lsb)) & mask;
    
    if(bits > 32) {
        (*blend_count_msb) = _mm256_adds_epi8((*blend_count_msb), _mm256_blendv_epi8((*vote_plus), (*vote_minus), movemask_inverse((new_val >> 32) & mask)));
        blend_value |= ((uint64_t)_mm256_movemask_epi8((*blend_count_msb))) << 32;
    }
}

static inline uint64_t calc_blend_rm_simd(__m256i *blend_count_lsb, __m256i *blend_count_msb, __m256i *vote_plus, __m256i *vote_minus, uint64_t new_val, uint64_t old_val, const uint64_t mask, const int bits) {

    (*blend_count_lsb) = _mm256_adds_epi8((*blend_count_lsb), _mm256_blendv_epi8((*vote_plus), (*vote_minus), movemask_inverse(new_val & mask)));
    uint64_t blend_value = (uint64_t)_mm256_movemask_epi8((*blend_count_lsb))&mask;
    
    // Removal of the oldest item
    (*blend_count_lsb) = _mm256_adds_epi8((*blend_count_lsb), _mm256_blendv_epi8((*vote_minus), (*vote_plus), movemask_inverse(old_val & mask)));

    if(bits > 32) {
        (*blend_count_msb) = _mm256_adds_epi8((*blend_count_msb), _mm256_blendv_epi8((*vote_plus), (*vote_minus), movemask_inverse((new_val >> 32) & mask)));
        blend_value |= ((uint64_t)_mm256_movemask_epi8((*blend_count_msb)))<<32;

        (*blend_count_msb) = _mm256_adds_epi8((*blend_count_msb),_mm256_blendv_epi8((*vote_minus), (*vote_plus), movemask_inverse((old_val >> 32) & mask)));
    }
    
    return blend_value;
}

#define process_buffer(prefix, current_seed, valid_latest_seed, debug)                              \
                                                                                                    \
    prefix##_blend_buffer[prefix##_blend_buffer_pos].x = current_seed.x;                            \
    prefix##_blend_buffer[prefix##_blend_buffer_pos].y = (current_seed.y & MASK_KMER_INDEX) >> 1;   \
                                                                                                    \
    uint64_t newest_span = prefix##_blend_buffer[prefix##_blend_buffer_pos].x & MASK_KMER_SPAN;     \
    uint64_t current_end = prefix##_blend_buffer[prefix##_blend_buffer_pos].y + newest_span;        \
                                                                                                    \
    if (++prefix##_blend_buffer_pos == neighbor_window) {                                           \
        prefix##_blend_buffer_pos = 0;                                                              \
    }                                                                                               \
                                                                                                    \
    if (++prefix##_neighbors_processed >= neighbor_window) {                                        \
        blend_value = calc_blend_rm_simd(&prefix##_blend_count_lsb, &prefix##_blend_count_msb,      \
                                         &vote_plus, &vote_minus, current_seed.x >> 14,             \
                                         prefix##_blend_buffer[prefix##_blend_buffer_pos].x >> 14,  \
                                         MASK_LS_32_BITS, blend_bits_count);                        \
        valid_latest_seed.x = blend_value << 14 | current_end - prefix##_blend_buffer[prefix##_blend_buffer_pos].y;          \
        valid_latest_seed.y = prefix##_blend_buffer[prefix##_blend_buffer_pos].y;                   \
                                                                                                    \
        if (str_len < current_end) {                                                                \
            printf("[DEBUG] Overflow detected info %d (seed: %lu, extra: %lu).\n",                  \
                    debug, seed_count, current_end - str_len);                                      \
        }                                                                                           \
                                                                                                    \
        temp_fuzzy_seeds[seed_count++] = valid_latest_seed;                                         \
        if (seed_count == seed_capacity) {                                                          \
            seed_capacity *= 2;                                                                     \
            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, seed_capacity * sizeof(uint128_t)); \
            temp_fuzzy_seeds = *fuzzy_seeds;                                                        \
            if (!temp_fuzzy_seeds) return 0;                                                        \
        }                                                                                           \
        valid_latest_seed.x = blend_value << 14 | newest_span;                                      \
        valid_latest_seed.y = current_seed.y;                                                       \
    } else { /* only addition of current_seed */                                                    \
        calc_blend_simd(&prefix##_blend_count_lsb, &prefix##_blend_count_msb,                       \
                        &vote_plus, &vote_minus,                                                    \
                        current_seed.x >> 14,                                                       \
                        MASK_LS_32_BITS, blend_bits_count);                                         \
    }


uint64_t blend_sb_sketch(const char *str, int str_len, int window_size, int kmer_size, int blend_bits, uint64_t neighbor_window, uint32_t str_id, uint128_t **fuzzy_seeds) {
    
    assert(str_len > 0 && (window_size > 0 && window_size + kmer_size < 8192) && (kmer_size > 0 && kmer_size <= 28) && (neighbor_window > 0 && neighbor_window + kmer_size < 8192) && (blend_bits <= 56));
    
    const int blend_bits_count = (blend_bits > 0) ? blend_bits : 2 * kmer_size;
    const uint64_t mask_blend_bits = (1ULL << blend_bits_count) - 1; // mask to get BLEND bits

    const uint64_t rc_shift_bits = 2 * (kmer_size- 1);
    const uint64_t mask_kmer = (1ULL << 2 * kmer_size) - 1;          // mask to get kmer (2 * kmer-size number of ones)

    uint64_t current_kmer[2] = {0, 0};
    int char_index = 0, temp_seed_buffer_pos = 0, current_kmer_span = 0, seed_buffer_pos = 0, min_pos = 0, kmer_span = 0;

    const static int seed_buffer_len = 256;
    uint128_t *seed_buffer = (uint128_t *)malloc(sizeof(uint128_t) * seed_buffer_len);
    uint128_t min = { UINT64_MAX, UINT64_MAX };
    
    // BLEND Variables
    uint64_t blend_value = 0;
    uint128_t *forward_blend_buffer = (uint128_t*)malloc(sizeof(uint128_t) * neighbor_window);
    uint128_t *reverse_blend_buffer = (uint128_t*)malloc(sizeof(uint128_t) * neighbor_window);
    uint64_t forward_blend_buffer_pos = 0, reverse_blend_buffer_pos = 0;
    uint64_t forward_neighbors_processed = 0, reverse_neighbors_processed = 0; // number of minimizers processed
    uint64_t forward_window_start_index = 0, reverse_window_start_index = 0;
    
    // SSE4-related variables
    __m256i vote_plus = _mm256_set1_epi8(1);
    __m256i vote_minus = _mm256_set1_epi8(-1);
    __m256i forward_blend_count_lsb = _mm256_set1_epi8(0);
    __m256i forward_blend_count_msb = _mm256_set1_epi8(0);
    __m256i reverse_blend_count_lsb = _mm256_set1_epi8(0);
    __m256i reverse_blend_count_msb = _mm256_set1_epi8(0);
    
    memset(seed_buffer, 0xff, sizeof(uint128_t) * seed_buffer_len);
    memset(forward_blend_buffer, 0, sizeof(uint128_t) * neighbor_window);
    memset(reverse_blend_buffer, 0, sizeof(uint128_t) * neighbor_window);

    uint64_t seed_capacity = 3 * str_len / window_size;
    *fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * seed_capacity);
    uint128_t *temp_fuzzy_seeds = *fuzzy_seeds;
    uint64_t seed_count = 0;

    if (!seed_buffer || !forward_blend_buffer || !reverse_blend_buffer || !fuzzy_seeds) {
        fprintf(stderr, "Poor programming skills. Please curse akmami...\n");
        abort();
    }

    for ( ; char_index < str_len; char_index++ ) {
        
        int ch = seq_nt4_table[(uint8_t)str[char_index]];
        
        uint128_t info = { UINT64_MAX, UINT64_MAX };
        uint128_t f_info = { UINT64_MAX, UINT64_MAX };
        
        if (ch < 4) { // skip an ambiguous base
            int strand;
            kmer_span = current_kmer_span + 1 < kmer_size ? current_kmer_span + 1 : kmer_size;
            
            current_kmer[0] = ((current_kmer[0] << 2 | ch) & mask_kmer); // forward k-mer
            current_kmer[1] = ((current_kmer[1] >> 2) | (3ULL ^ ch) << rc_shift_bits); // reverse k-mer k-mer
            
            if (current_kmer[0] == current_kmer[1]) continue; // skip symmetric k-mers
            
            strand = current_kmer[0] < current_kmer[1] ? 0 : 1; // strand
            current_kmer_span++;
            if (current_kmer_span >= kmer_size && kmer_span < 256) {
                f_info.x = hash64(current_kmer[strand], mask_blend_bits) << 14 | kmer_span;
                f_info.y = (uint64_t)str_id << 32 | (uint32_t)char_index << 1 | strand;
            }
        } else {
            current_kmer_span = 0;
            kmer_span = 0;
        }
        seed_buffer[seed_buffer_pos] = f_info; // need to do this here as appropriate seed_buffer[seed_buffer_pos] are needed below
        
        // special case for the first window_size - because identical k-mers are not stored yet
        if (current_kmer_span == window_size + kmer_size - 1 && min.x != UINT64_MAX) {
            for (temp_seed_buffer_pos = seed_buffer_pos + 1; temp_seed_buffer_pos < seed_buffer_len; temp_seed_buffer_pos++) {
                if (min.x == seed_buffer[temp_seed_buffer_pos].x && seed_buffer[temp_seed_buffer_pos].y != min.y) {
                    if ((seed_buffer[temp_seed_buffer_pos].y & 1) == 1) { // reverse strand
                        process_buffer(reverse, seed_buffer[temp_seed_buffer_pos], info, 198)
                    } else {
                        process_buffer(forward, seed_buffer[temp_seed_buffer_pos], info, 200)
                    }
                }
            }
            for (temp_seed_buffer_pos = 0; temp_seed_buffer_pos < seed_buffer_pos; temp_seed_buffer_pos++) {
                if (min.x == seed_buffer[temp_seed_buffer_pos].x && seed_buffer[temp_seed_buffer_pos].y != min.y) {
                    if ((seed_buffer[temp_seed_buffer_pos].y & 1) == 1) { // reverse strand
                        process_buffer(reverse, seed_buffer[temp_seed_buffer_pos], info, 207)
                    } else {
                        process_buffer(forward, seed_buffer[temp_seed_buffer_pos], info, 209)
                    }
                }
            }
        }

        if (f_info.x <= min.x) { // a new minimum; then write the old min
            if (current_kmer_span >= window_size + kmer_size && min.x != UINT64_MAX) {
                if ((min.y & 1) == 1) { // reverse strand
                    process_buffer(reverse, min, info, 218)
                } else {
                    process_buffer(forward, min, info, 220)
                }
            }
            min = f_info, min_pos = seed_buffer_pos;
        } else if (seed_buffer_pos == min_pos) { // old min has moved outside the window_size
            if (current_kmer_span >= window_size + kmer_size - 1 && min.x != UINT64_MAX) {
                if ((min.y & 1) == 1) { // reverse strand
                    process_buffer(reverse, min, info, 227)
                } else {
                    process_buffer(forward, min, info, 229)
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
            if (current_kmer_span >= window_size + kmer_size - 1 && min.x != UINT64_MAX) { // write identical k-mers
                for (temp_seed_buffer_pos = seed_buffer_pos + 1; temp_seed_buffer_pos < seed_buffer_len; temp_seed_buffer_pos++) { // these two loops make sure the output is sorted
                    if (min.x == seed_buffer[temp_seed_buffer_pos].x && min.y != seed_buffer[temp_seed_buffer_pos].y) {
                        if ((seed_buffer[temp_seed_buffer_pos].y & 1) == 1) { // reverse strand
                            process_buffer(reverse, seed_buffer[temp_seed_buffer_pos], info, 246)
                        } else {
                            process_buffer(forward, seed_buffer[temp_seed_buffer_pos], info, 248)
                        }
                    }
                }
                for (temp_seed_buffer_pos = 0; temp_seed_buffer_pos <= seed_buffer_pos; temp_seed_buffer_pos++)
                    if (min.x == seed_buffer[temp_seed_buffer_pos].x && min.y != seed_buffer[temp_seed_buffer_pos].y) {
                        if ((seed_buffer[temp_seed_buffer_pos].y & 1) == 1) { // reverse strand
                            process_buffer(reverse, seed_buffer[temp_seed_buffer_pos], info, 255)
                        } else {
                            process_buffer(forward, seed_buffer[temp_seed_buffer_pos], info, 257)
                        }
                    }
            }
        }

        if (++seed_buffer_pos == window_size) {
            seed_buffer_pos = 0;
        }
    }

    if (min.x != UINT64_MAX) {
        if ((min.y & 1) == 1) { // Reverse strand
            process_buffer(reverse, min, min, 270)
        } else {
            process_buffer(forward, min, min, 272)
        }
    }

    free(seed_buffer);
    free(forward_blend_buffer);
    free(reverse_blend_buffer);

    return seed_count;
}

uint64_t blend_sketch(const char *str, int str_len, int window_size, int kmer_size, int blend_bits, int neighbor_window, uint32_t str_id, uint128_t **fuzzy_seeds) {
	return blend_sb_sketch(str, str_len, window_size, kmer_size, blend_bits, (uint64_t)neighbor_window, str_id, fuzzy_seeds);
}


#ifdef __cplusplus
}
#endif