#include "blend.h"


#ifdef __cplusplus
extern "C" {
#endif

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
    const uint64_t rc_shift_bits = 2 * (kmer_size- 1);
    const uint64_t kmer_mask = (1ULL << 2 * kmer_size) - 1;          // mask to get kmer (2 * kmer-size number of ones)
    const uint64_t mask32 = (1ULL << 32) - 1;                        // = 0xFFFFFFFF, get right-most 32 bits
    const uint64_t kmer_span_mask = (1ULL << 14) - 1;                // = 0x3FFF, get kmer span
    const uint64_t blend_bits_mask = (1ULL << blend_bits_count) - 1; // mash to get BLEND bits
    const uint64_t kmer_index_mask = (1ULL << 31) - 1;               // for iMask, the actual value should be shifted right first
    uint64_t kmer[2] = {0, 0};
    int char_index = 0, temp_buffer_index = 0, current_kmer_span = 0, buf_pos = 0, min_pos = 0, kmer_span = 0;

    uint32_t tot_span;

    const static int buf_len = 256;
    uint128_t *buf = (uint128_t *)malloc(sizeof(uint128_t) * buf_len);
    uint128_t min = { UINT64_MAX, UINT64_MAX };
    
    // BLEND Variables
    uint64_t blendVal = 0;
    uint128_t *forward_blend_buf = (uint128_t*)malloc(sizeof(uint128_t) * n_neighbors);
    uint128_t *reverse_blend_buf = (uint128_t*)malloc(sizeof(uint128_t) * n_neighbors);
    uint64_t forward_blend_pos = 0, reverse_blend_pos = 0;
    uint64_t f_nm = 0, r_nm = 0; // number of minimizers processed
    
    // SSE4-related variables
    __m256i ma = _mm256_set1_epi8(1);
    __m256i mb = _mm256_set1_epi8(-1);
    __m256i f_blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i f_blndcnt_msb = _mm256_set1_epi8(0);
    __m256i r_blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i r_blndcnt_msb = _mm256_set1_epi8(0);
    
    memset(buf, 0xff, sizeof(uint128_t) * buf_len);
    memset(forward_blend_buf, 0, sizeof(uint128_t) * n_neighbors);
    memset(reverse_blend_buf, 0, sizeof(uint128_t) * n_neighbors);

    uint64_t fuzzy_seeds_capacity = 3 * len / window;
    *fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * fuzzy_seeds_capacity);
    uint128_t *temp_fuzzy_seeds = *fuzzy_seeds;
    uint64_t seed_count = 0;

    if (!buf || !forward_blend_buf || !reverse_blend_buf || !fuzzy_seeds) {
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
        buf[buf_pos] = f_info; // need to do this here as appropriate buf[buf_pos] are needed below
        
        // special case for the first window - because identical k-mers are not stored yet
        if (current_kmer_span == window + kmer_size - 1 && min.x != UINT64_MAX) {
            for (temp_buffer_index = buf_pos + 1; temp_buffer_index < buf_len; temp_buffer_index++) {
                if (min.x == buf[temp_buffer_index].x && buf[temp_buffer_index].y != min.y) {
                    if ((buf[temp_buffer_index].y & 1) == 1) { // reverse strand
                        reverse_blend_buf[reverse_blend_pos].x = buf[temp_buffer_index].x;
                        reverse_blend_buf[reverse_blend_pos].y = (buf[temp_buffer_index].y & kmer_index_mask) >> 1;
                        tot_span = reverse_blend_buf[reverse_blend_pos].y;

                        if (++reverse_blend_pos == n_neighbors) {
                            reverse_blend_pos = 0;
                        }

                        if (++r_nm >= n_neighbors) {
                            blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, reverse_blend_buf[reverse_blend_pos].x >> 14, mask32, blend_bits_count);
                            tot_span = tot_span - reverse_blend_buf[reverse_blend_pos].y + (reverse_blend_buf[reverse_blend_pos].x & kmer_span_mask);
                            info.x = blendVal << 14 | tot_span;
                            info.y = buf[temp_buffer_index].y;
                            temp_fuzzy_seeds[seed_count++] = info;
                            if (seed_count == fuzzy_seeds_capacity) {
                                fuzzy_seeds_capacity *= 2;
                                *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                temp_fuzzy_seeds = *fuzzy_seeds;
                                if (!temp_fuzzy_seeds) return 0;
                            }
                        } else { // only addition of buf[temp_buffer_index]
                            calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x>>14, mask32, blend_bits_count);
                        }
                    } else {
                        forward_blend_buf[forward_blend_pos].x = buf[temp_buffer_index].x;
                        forward_blend_buf[forward_blend_pos].y = (buf[temp_buffer_index].y & kmer_index_mask) >> 1;
                        tot_span = forward_blend_buf[forward_blend_pos].y;

                        if (++forward_blend_pos == n_neighbors) {
                            forward_blend_pos = 0;
                        }

                        if (++f_nm >= n_neighbors) {
                            blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, forward_blend_buf[forward_blend_pos].x >> 14, mask32, blend_bits_count);
                            tot_span = tot_span - forward_blend_buf[forward_blend_pos].y + (forward_blend_buf[forward_blend_pos].x & kmer_span_mask);
                            info.x = blendVal << 14 | tot_span;
                            info.y = buf[temp_buffer_index].y;
                            temp_fuzzy_seeds[seed_count++] = info;
                            if (seed_count == fuzzy_seeds_capacity) {
                                fuzzy_seeds_capacity *= 2;
                                *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                temp_fuzzy_seeds = *fuzzy_seeds;
                                if (!temp_fuzzy_seeds) return 0;
                            }
                        } else { // only addition of buf[temp_buffer_index]
                            calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, mask32, blend_bits_count);
                        }
                    }
                }
            }

            for (temp_buffer_index = 0; temp_buffer_index < buf_pos; temp_buffer_index++)
                if (min.x == buf[temp_buffer_index].x && buf[temp_buffer_index].y != min.y) {
                    if ((buf[temp_buffer_index].y & 1) == 1) { // reverse strand
                        reverse_blend_buf[reverse_blend_pos].x = buf[temp_buffer_index].x;
                        reverse_blend_buf[reverse_blend_pos].y = (buf[temp_buffer_index].y & kmer_index_mask) >> 1;
                        tot_span = reverse_blend_buf[reverse_blend_pos].y;

                        if (++reverse_blend_pos == n_neighbors) {
                            reverse_blend_pos = 0;
                        }

                        if (++r_nm >= n_neighbors) {
                            blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, reverse_blend_buf[reverse_blend_pos].x >> 14, mask32, blend_bits_count);
                            tot_span = tot_span - reverse_blend_buf[reverse_blend_pos].y + (reverse_blend_buf[reverse_blend_pos].x & kmer_span_mask);
                            info.x = blendVal << 14 | tot_span;
                            info.y = buf[temp_buffer_index].y;
                            temp_fuzzy_seeds[seed_count++] = info;
                            if (seed_count == fuzzy_seeds_capacity) {
                                fuzzy_seeds_capacity *= 2;
                                *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                temp_fuzzy_seeds = *fuzzy_seeds;
                                if (!temp_fuzzy_seeds) return 0;
                            }
                        } else { // only addition of buf[temp_buffer_index]
                            calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, mask32, blend_bits_count);
                        }
                    } else {
                        forward_blend_buf[forward_blend_pos].x = buf[temp_buffer_index].x;
                        forward_blend_buf[forward_blend_pos].y = (buf[temp_buffer_index].y & kmer_index_mask) >> 1;
                        tot_span = forward_blend_buf[forward_blend_pos].y;

                        if (++forward_blend_pos == n_neighbors) {
                            forward_blend_pos = 0;
                        }

                        if (++f_nm >= n_neighbors) {
                            blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, forward_blend_buf[forward_blend_pos].x >> 14, mask32, blend_bits_count);
                            tot_span = tot_span - forward_blend_buf[forward_blend_pos].y + (forward_blend_buf[forward_blend_pos].x & kmer_span_mask);
                            info.x = blendVal << 14 | tot_span;
                            info.y = buf[temp_buffer_index].y;
                            temp_fuzzy_seeds[seed_count++] = info;
                            if (seed_count == fuzzy_seeds_capacity) {
                                fuzzy_seeds_capacity *= 2;
                                *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                temp_fuzzy_seeds = *fuzzy_seeds;
                                if (!temp_fuzzy_seeds) return 0;
                            }
                        } else { // only addition of buf[temp_buffer_index]
                            calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, mask32, blend_bits_count);
                        }
                    }
                }
        }

        if (f_info.x <= min.x) { // a new minimum; then write the old min
            if (current_kmer_span >= window + kmer_size && min.x != UINT64_MAX) {
                if ((min.y & 1) == 1) { // reverse strand
                    reverse_blend_buf[reverse_blend_pos].x = min.x;
                    reverse_blend_buf[reverse_blend_pos].y = (min.y & kmer_index_mask) >> 1;
                    tot_span = reverse_blend_buf[reverse_blend_pos].y;

                    if (++reverse_blend_pos == n_neighbors) {
                        reverse_blend_pos = 0;
                    }

                    if (++r_nm >= n_neighbors) {
                        blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, reverse_blend_buf[reverse_blend_pos].x >> 14, mask32, blend_bits_count);
                        tot_span = tot_span - reverse_blend_buf[reverse_blend_pos].y + (reverse_blend_buf[reverse_blend_pos].x & kmer_span_mask);
                        info.x = blendVal << 14 | tot_span;
                        info.y = min.y;
                        temp_fuzzy_seeds[seed_count++] = info;
                        if (seed_count == fuzzy_seeds_capacity) {
                            fuzzy_seeds_capacity *= 2;
                            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                            temp_fuzzy_seeds = *fuzzy_seeds;
                            if (!temp_fuzzy_seeds) return 0;
                        }
                    } else { // only addition of min
                        calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blend_bits_count);
                    }
                } else {
                    forward_blend_buf[forward_blend_pos].x = min.x;
                    forward_blend_buf[forward_blend_pos].y = (min.y & kmer_index_mask) >> 1;
                    tot_span = forward_blend_buf[forward_blend_pos].y;

                    if (++forward_blend_pos == n_neighbors) {
                        forward_blend_pos = 0;
                    }

                    if (++f_nm >= n_neighbors) {
                        blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, forward_blend_buf[forward_blend_pos].x >> 14, mask32, blend_bits_count);
                        tot_span = tot_span - forward_blend_buf[forward_blend_pos].y + (forward_blend_buf[forward_blend_pos].x & kmer_span_mask);
                        info.x = blendVal << 14 | tot_span;
                        info.y = min.y;
                        temp_fuzzy_seeds[seed_count++] = info;
                        if (seed_count == fuzzy_seeds_capacity) {
                            fuzzy_seeds_capacity *= 2;
                            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                            temp_fuzzy_seeds = *fuzzy_seeds;
                            if (!temp_fuzzy_seeds) return 0;
                        }
                    } else { // only addition of min
                        calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blend_bits_count);
                    }
                }
            }
            min = f_info, min_pos = buf_pos;
        } else if (buf_pos == min_pos) { // old min has moved outside the window
            if (current_kmer_span >= window + kmer_size - 1 && min.x != UINT64_MAX) {
                if ((min.y & 1) == 1) { // reverse strand
                    reverse_blend_buf[reverse_blend_pos].x = min.x;
                    reverse_blend_buf[reverse_blend_pos].y = (min.y & kmer_index_mask) >> 1;
                    tot_span = reverse_blend_buf[reverse_blend_pos].y;

                    if (++reverse_blend_pos == n_neighbors) {
                        reverse_blend_pos = 0;
                    }

                    if (++r_nm >= n_neighbors) {
                        blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, reverse_blend_buf[reverse_blend_pos].x >> 14, mask32, blend_bits_count);
                        tot_span = tot_span - reverse_blend_buf[reverse_blend_pos].y + (reverse_blend_buf[reverse_blend_pos].x & kmer_span_mask);
                        info.x = blendVal << 14 | tot_span;
                        info.y = min.y;
                        temp_fuzzy_seeds[seed_count++] = info;
                        if (seed_count == fuzzy_seeds_capacity) {
                            fuzzy_seeds_capacity *= 2;
                            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                            temp_fuzzy_seeds = *fuzzy_seeds;
                            if (!temp_fuzzy_seeds) return 0;
                        }
                    } else { // only addition of min
                        calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blend_bits_count);
                    }
                } else {
                    forward_blend_buf[forward_blend_pos].x = min.x;
                    forward_blend_buf[forward_blend_pos].y = (min.y & kmer_index_mask) >> 1;
                    tot_span = forward_blend_buf[forward_blend_pos].y;

                    if(++forward_blend_pos == n_neighbors) {
                        forward_blend_pos = 0;
                    }

                    if(++f_nm >= n_neighbors) {
                        blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, forward_blend_buf[forward_blend_pos].x >> 14, mask32, blend_bits_count);
                        tot_span = tot_span - forward_blend_buf[forward_blend_pos].y + (forward_blend_buf[forward_blend_pos].x & kmer_span_mask);
                        info.x = blendVal << 14 | tot_span;
                        info.y = min.y;
                        temp_fuzzy_seeds[seed_count++] = info;
                        if (seed_count == fuzzy_seeds_capacity) {
                            fuzzy_seeds_capacity *= 2;
                            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                            temp_fuzzy_seeds = *fuzzy_seeds;
                            if (!temp_fuzzy_seeds) return 0;
                        }
                    } else { // only addition of min
                        calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blend_bits_count);
                    }
                }
            }
            for (temp_buffer_index = buf_pos + 1, min.x = UINT64_MAX; temp_buffer_index < buf_len; temp_buffer_index++) { // the two loops are necessary when there are identical k-mers. This one starts from the oldest (first in: buf_pos + 1) 
                if (min.x >= buf[temp_buffer_index].x) {
                    min = buf[temp_buffer_index], min_pos = temp_buffer_index; // >= is important s.t. min is always the closest k-mer... chosing the last one because we will store *all* identical minimizers and continue with the last one once all identical ones are stored. @IMPORTANT: should we really store all identical k-mers? it is true that they are the minimizers for their own window but not sure how effective it is
                }
            }
            for (temp_buffer_index = 0; temp_buffer_index <= buf_pos; temp_buffer_index++) {
                if (min.x >= buf[temp_buffer_index].x) {
                    min = buf[temp_buffer_index], min_pos = temp_buffer_index;
                }
            }
            if (current_kmer_span >= window + kmer_size - 1 && min.x != UINT64_MAX) { // write identical k-mers
                for (temp_buffer_index = buf_pos + 1; temp_buffer_index < buf_len; temp_buffer_index++) { // these two loops make sure the output is sorted
                    if (min.x == buf[temp_buffer_index].x && min.y != buf[temp_buffer_index].y) {
                        if ((buf[temp_buffer_index].y & 1) == 1) { // reverse strand
                            reverse_blend_buf[reverse_blend_pos].x = buf[temp_buffer_index].x;
                            reverse_blend_buf[reverse_blend_pos].y = (buf[temp_buffer_index].y & kmer_index_mask) >> 1;
                            tot_span = reverse_blend_buf[reverse_blend_pos].y;

                            if (++reverse_blend_pos == n_neighbors) {
                                reverse_blend_pos = 0;
                            }

                            if (++r_nm >= n_neighbors) {
                                blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, reverse_blend_buf[reverse_blend_pos].x >> 14, mask32, blend_bits_count);
                                tot_span = tot_span - reverse_blend_buf[reverse_blend_pos].y + (reverse_blend_buf[reverse_blend_pos].x & kmer_span_mask);
                                info.x = blendVal << 14 | tot_span;
                                info.y = buf[temp_buffer_index].y;
                                temp_fuzzy_seeds[seed_count++] = info;
                                if (seed_count == fuzzy_seeds_capacity) {
                                    fuzzy_seeds_capacity *= 2;
                                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                    temp_fuzzy_seeds = *fuzzy_seeds;
                                    if (!temp_fuzzy_seeds) return 0;
                                }
                            } else { // only addition of buf[temp_buffer_index]
                                calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, mask32, blend_bits_count);
                            }
                        } else {
                            forward_blend_buf[forward_blend_pos].x = buf[temp_buffer_index].x;
                            forward_blend_buf[forward_blend_pos].y = (buf[temp_buffer_index].y & kmer_index_mask) >> 1;
                            tot_span = forward_blend_buf[forward_blend_pos].y;

                            if (++forward_blend_pos == n_neighbors) {
                                forward_blend_pos = 0;
                            }

                            if (++f_nm >= n_neighbors) {
                                blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, forward_blend_buf[forward_blend_pos].x >> 14, mask32, blend_bits_count);
                                tot_span = tot_span - forward_blend_buf[forward_blend_pos].y + (forward_blend_buf[forward_blend_pos].x & kmer_span_mask);
                                info.x = blendVal << 14 | tot_span;
                                info.y = buf[temp_buffer_index].y;
                                temp_fuzzy_seeds[seed_count++] = info;
                                if (seed_count == fuzzy_seeds_capacity) {
                                    fuzzy_seeds_capacity *= 2;
                                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                    temp_fuzzy_seeds = *fuzzy_seeds;
                                    if (!temp_fuzzy_seeds) return 0;
                                }
                            } else { // only addition of buf[temp_buffer_index]
                                calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, mask32, blend_bits_count);
                            }
                        }
                    }
                }
                for (temp_buffer_index = 0; temp_buffer_index <= buf_pos; temp_buffer_index++)
                    if (min.x == buf[temp_buffer_index].x && min.y != buf[temp_buffer_index].y) {
                        if ((buf[temp_buffer_index].y & 1) == 1) { // reverse strand
                            reverse_blend_buf[reverse_blend_pos].x = buf[temp_buffer_index].x;
                            reverse_blend_buf[reverse_blend_pos].y = (buf[temp_buffer_index].y & kmer_index_mask) >> 1;
                            tot_span = reverse_blend_buf[reverse_blend_pos].y;

                            if (++reverse_blend_pos == n_neighbors) { 
                                reverse_blend_pos = 0;
                            }

                            if (++r_nm >= n_neighbors) {
                                blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, reverse_blend_buf[reverse_blend_pos].x >> 14, mask32, blend_bits_count);
                                tot_span = tot_span - reverse_blend_buf[reverse_blend_pos].y + (reverse_blend_buf[reverse_blend_pos].x & kmer_span_mask);
                                info.x = blendVal << 14 | tot_span;
                                info.y = buf[temp_buffer_index].y;
                                temp_fuzzy_seeds[seed_count++] = info;
                                if (seed_count == fuzzy_seeds_capacity) {
                                    fuzzy_seeds_capacity *= 2;
                                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                    temp_fuzzy_seeds = *fuzzy_seeds;
                                    if (!temp_fuzzy_seeds) return 0;
                                }
                            } else { // only addition of buf[temp_buffer_index]
                                calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, mask32, blend_bits_count);
                            }
                        } else {
                            forward_blend_buf[forward_blend_pos].x = buf[temp_buffer_index].x;
                            forward_blend_buf[forward_blend_pos].y = (buf[temp_buffer_index].y & kmer_index_mask) >> 1;
                            tot_span = forward_blend_buf[forward_blend_pos].y;

                            if (++forward_blend_pos == n_neighbors) {
                                forward_blend_pos = 0;
                            }

                            if (++f_nm >= n_neighbors) {
                                blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x >> 14, forward_blend_buf[forward_blend_pos].x >> 14, mask32, blend_bits_count);
                                tot_span = tot_span - forward_blend_buf[forward_blend_pos].y + (forward_blend_buf[forward_blend_pos].x & kmer_span_mask);
                                info.x = blendVal << 14 | tot_span;
                                info.y = buf[temp_buffer_index].y;
                                temp_fuzzy_seeds[seed_count++] = info;
                                if (seed_count == fuzzy_seeds_capacity) {
                                    fuzzy_seeds_capacity *= 2;
                                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                    temp_fuzzy_seeds = *fuzzy_seeds;
                                    if (!temp_fuzzy_seeds) return 0;
                                }
                            } else { // only addition of buf[temp_buffer_index]
                                calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[temp_buffer_index].x>>14, mask32, blend_bits_count);
                            }
                        }
                    }
            }
        }

        if (++buf_pos == window) {
            buf_pos = 0;
        }
    }

    if (min.x != UINT64_MAX) {
        if ((min.y & 1) == 1) { // Reverse strand
            reverse_blend_buf[reverse_blend_pos].x = min.x;
            reverse_blend_buf[reverse_blend_pos].y = (min.y & kmer_index_mask) >> 1;
            tot_span = reverse_blend_buf[reverse_blend_pos].y;

            if (++reverse_blend_pos == n_neighbors) {
                reverse_blend_pos = 0;
            }

            if (++r_nm >= n_neighbors) {
                blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, reverse_blend_buf[reverse_blend_pos].x >> 14, mask32, blend_bits_count);
                tot_span = tot_span - reverse_blend_buf[reverse_blend_pos].y + (reverse_blend_buf[reverse_blend_pos].x & kmer_span_mask);
                min.x = blendVal << 14 | tot_span;
                temp_fuzzy_seeds[seed_count++] = min;
                if (seed_count == fuzzy_seeds_capacity) {
                    fuzzy_seeds_capacity *= 2;
                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                    temp_fuzzy_seeds = *fuzzy_seeds;
                    if (!temp_fuzzy_seeds) return 0;
                }
            } else { // only addition of min
                calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blend_bits_count);
            }
        } else {
            forward_blend_buf[forward_blend_pos].x = min.x;
            forward_blend_buf[forward_blend_pos].y = (min.y & kmer_index_mask) >> 1;
            tot_span = forward_blend_buf[forward_blend_pos].y;

            if(++forward_blend_pos == n_neighbors) forward_blend_pos = 0;
            if(++f_nm >= n_neighbors) {
                blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, forward_blend_buf[forward_blend_pos].x >> 14, mask32, blend_bits_count);
                tot_span = tot_span - forward_blend_buf[forward_blend_pos].y + (forward_blend_buf[forward_blend_pos].x & kmer_span_mask);
                min.x = blendVal << 14 | tot_span;
                temp_fuzzy_seeds[seed_count++] = min;
                if (seed_count == fuzzy_seeds_capacity) {
                    fuzzy_seeds_capacity *= 2;
                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                    temp_fuzzy_seeds = *fuzzy_seeds;
                    if (!temp_fuzzy_seeds) return 0;
                }
            } else { // only addition of min
                calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blend_bits_count);
            }
        }
    }

    free(buf);
    free(forward_blend_buf);
    free(reverse_blend_buf);

    return seed_count;
}

uint64_t blend_sketch(const char *str, int len, int window, int kmer_size, int blend_bits, int n_neighbors, uint32_t str_id, uint128_t **fuzzy_seeds) {
	return blend_sb_sketch(str, len, window, kmer_size, blend_bits, (uint64_t)n_neighbors, str_id, fuzzy_seeds);
}

#ifdef __cplusplus
}
#endif