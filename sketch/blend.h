#ifndef __BLEND_SKETCH__
#define __BLEND_SKETCH__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "struct_def.h"
#include <x86intrin.h>


#define BLEND_GET_KMER(minimizer) ((minimizer).x >> 14)
#define BLEND_GET_LENGTH(minimizer) ((minimizer).x & 0x3FFF)
#define BLEND_GET_INDEX(minimizer) (((minimizer).y & 0xFFFFFFFF) >> 1)
#define BLEND_GET_REFERENCE_IDX(minimizer) (((minimizer).y >> 32))
#define BLEND_GET_STRAND(minimizer) (((minimizer).y & 1))

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
 * @param str           DNA sequence
 * @param len           length of $str
 * @param window        find a BLEND value for every $w consecutive k-mers
 * @param kmer_size     k-mer size
 * @param blend_bits    use blend_bits many bits when generating the hash values of seeds
 * @param n_neighbors   How many neighbors consecutive minimizer k-mers should be combined to generate a strobemer seed
 * @param rid           reference ID; will be copied to the output $p array
 * @param fuzzy_seeds
 */
uint64_t blend_sb_sketch(const char *str, int len, int window, int kmer_size, int blend_bits, uint64_t n_neighbors, uint32_t rid, uint128_t **fuzzy_seeds) {
    
    assert(len > 0 && (window > 0 && window + kmer_size < 8192) && (kmer_size > 0 && kmer_size <= 28) && (n_neighbors > 0 && n_neighbors + kmer_size < 8192) && (blend_bits <= 56));
    
    const int blndK = (blend_bits > 0) ? blend_bits : 2 * kmer_size;
    const uint64_t shift1 = 2 * (kmer_size- 1);
    const uint64_t mask = (1ULL << 2 * kmer_size) - 1;
    const uint64_t mask32 = (1ULL << 32) - 1;
    const uint64_t sMask = (1ULL << 14) - 1;
    const uint64_t blndMask = (1ULL << blndK) - 1;
    const uint64_t iMask = (1ULL << 31) - 1; // for iMask, the actual value should be shifted right first
    uint64_t kmer[2] = {0, 0};
    int i, j, l, buf_pos, min_pos, kmer_span = 0;

    uint32_t tot_span;

    // uint128_t buf[window];
    uint128_t* buf = (uint128_t*)malloc(sizeof(uint128_t) * window);
    uint128_t min = { UINT64_MAX, UINT64_MAX };
    
    // BLEND Variables
    uint64_t blendVal = 0;
    // uint128_t f_blndBuf[n_neighbors];
    // uint128_t r_blndBuf[n_neighbors];
    uint128_t* f_blndBuf = (uint128_t*)malloc(sizeof(uint128_t) * n_neighbors);
    uint128_t* r_blndBuf = (uint128_t*)malloc(sizeof(uint128_t) * n_neighbors);
    uint64_t f_blendpos = 0, r_blendpos = 0;
    uint64_t f_nm = 0, r_nm = 0; // number of minimizers processed
    
    // SSE4-related variables
    __m256i ma = _mm256_set1_epi8(1);
    __m256i mb = _mm256_set1_epi8(-1);
    __m256i f_blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i f_blndcnt_msb = _mm256_set1_epi8(0);
    __m256i r_blndcnt_lsb = _mm256_set1_epi8(0);
    __m256i r_blndcnt_msb = _mm256_set1_epi8(0);
    
    memset(buf, 0xff, sizeof(uint128_t) * window);
    memset(f_blndBuf, 0, sizeof(uint128_t) * n_neighbors);
    memset(r_blndBuf, 0, sizeof(uint128_t) * n_neighbors);

    uint64_t fuzzy_seeds_capacity = 3 * len / window;
    *fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * fuzzy_seeds_capacity);
    uint128_t *temp_fuzzy_seeds = *fuzzy_seeds;
    uint64_t seed_count = 0;

    if (!buf || !f_blndBuf || !r_blndBuf || !fuzzy_seeds) {
        fprintf(stderr, "Poor programming skills. Please curse akmami...\n");
        abort();
    }

    for (i = l = f_blendpos = r_blendpos = buf_pos = min_pos = 0; i < len; ++i) {
        
        int c = seq_nt4_table[(uint8_t)str[i]];
        
        uint128_t info = { UINT64_MAX, UINT64_MAX };
        uint128_t f_info = { UINT64_MAX, UINT64_MAX };
        
        if (c < 4) { // not an ambiguous base
            int z;
            kmer_span = l + 1 < kmer_size ? l + 1 : kmer_size;
            
            kmer[0] = ((kmer[0] << 2 | c) & mask); // forward k-mer
            kmer[1] = ((kmer[1] >> 2) | (3ULL ^ c) << shift1); //reverse k-mer k-mer
            
            if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
            
            z = kmer[0] < kmer[1] ? 0 : 1; // strand
            ++l;
            if (l >= kmer_size && kmer_span < 256) {
                f_info.x = hash64(kmer[z], blndMask) << 14 | kmer_span;
                f_info.y = (uint64_t)rid << 32 | (uint32_t)i << 1 | z;
            }
        } else {
            l = 0, kmer_span = 0;
        }
        buf[buf_pos] = f_info; // need to do this here as appropriate buf[buf_pos] are needed below
        
        // special case for the first window - because identical k-mers are not stored yet
        if (l == window + kmer_size - 1 && min.x != UINT64_MAX) {
            for (j = buf_pos + 1; j < window; ++j) {
                if (min.x == buf[j].x && buf[j].y != min.y) {
                    if ((buf[j].y & 1) == 1) { // reverse strand
                        r_blndBuf[r_blendpos].x = buf[j].x;
                        r_blndBuf[r_blendpos].y = (buf[j].y & iMask) >> 1;
                        tot_span = r_blndBuf[r_blendpos].y;

                        if (++r_blendpos == n_neighbors) {
                            r_blendpos = 0;
                        }

                        if (++r_nm >= n_neighbors) {
                            blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[j].x >> 14, r_blndBuf[r_blendpos].x >> 14, mask32, blndK);
                            tot_span = tot_span - r_blndBuf[r_blendpos].y + (r_blndBuf[r_blendpos].x & sMask);
                            info.x = blendVal << 14 | tot_span;
                            info.y = buf[j].y;
                            // kv_push(mm128_t, km, *p, info);
                            temp_fuzzy_seeds[seed_count++] = info;
                            if (seed_count == fuzzy_seeds_capacity) {
                                fuzzy_seeds_capacity *= 2;
                                *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                temp_fuzzy_seeds = *fuzzy_seeds;
                                if (!temp_fuzzy_seeds) return 0;
                            }
                        } else { // only addition of buf[j]
                            calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[j].x>>14, mask32, blndK);
                        }
                    } else {
                        f_blndBuf[f_blendpos].x = buf[j].x;
                        f_blndBuf[f_blendpos].y = (buf[j].y & iMask) >> 1;
                        tot_span = f_blndBuf[f_blendpos].y;

                        if (++f_blendpos == n_neighbors) {
                            f_blendpos = 0;
                        }

                        if (++f_nm >= n_neighbors) {
                            blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[j].x >> 14, f_blndBuf[f_blendpos].x >> 14, mask32, blndK);
                            tot_span = tot_span - f_blndBuf[f_blendpos].y + (f_blndBuf[f_blendpos].x & sMask);
                            info.x = blendVal << 14 | tot_span;
                            info.y = buf[j].y;
                            // kv_push(mm128_t, km, *p, info);
                            temp_fuzzy_seeds[seed_count++] = info;
                            if (seed_count == fuzzy_seeds_capacity) {
                                fuzzy_seeds_capacity *= 2;
                                *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                temp_fuzzy_seeds = *fuzzy_seeds;
                                if (!temp_fuzzy_seeds) return 0;
                            }
                        } else { // only addition of buf[j]
                            calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[j].x >> 14, mask32, blndK);
                        }
                    }
                }
            }

            for (j = 0; j < buf_pos; ++j)
                if (min.x == buf[j].x && buf[j].y != min.y) {
                    if ((buf[j].y & 1) == 1) { // reverse strand
                        r_blndBuf[r_blendpos].x = buf[j].x;
                        r_blndBuf[r_blendpos].y = (buf[j].y & iMask) >> 1;
                        tot_span = r_blndBuf[r_blendpos].y;

                        if (++r_blendpos == n_neighbors) {
                            r_blendpos = 0;
                        }

                        if (++r_nm >= n_neighbors) {
                            blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[j].x >> 14, r_blndBuf[r_blendpos].x >> 14, mask32, blndK);
                            tot_span = tot_span - r_blndBuf[r_blendpos].y + (r_blndBuf[r_blendpos].x & sMask);
                            info.x = blendVal << 14 | tot_span;
                            info.y = buf[j].y;
                            // kv_push(mm128_t, km, *p, info);
                            temp_fuzzy_seeds[seed_count++] = info;
                            if (seed_count == fuzzy_seeds_capacity) {
                                fuzzy_seeds_capacity *= 2;
                                *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                temp_fuzzy_seeds = *fuzzy_seeds;
                                if (!temp_fuzzy_seeds) return 0;
                            }
                        } else { // only addition of buf[j]
                            calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[j].x >> 14, mask32, blndK);
                        }
                    } else {
                        f_blndBuf[f_blendpos].x = buf[j].x;
                        f_blndBuf[f_blendpos].y = (buf[j].y & iMask) >> 1;
                        tot_span = f_blndBuf[f_blendpos].y;

                        if (++f_blendpos == n_neighbors) {
                            f_blendpos = 0;
                        }

                        if (++f_nm >= n_neighbors) {
                            blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[j].x >> 14, f_blndBuf[f_blendpos].x >> 14, mask32, blndK);
                            tot_span = tot_span - f_blndBuf[f_blendpos].y + (f_blndBuf[f_blendpos].x & sMask);
                            info.x = blendVal << 14 | tot_span;
                            info.y = buf[j].y;
                            // kv_push(mm128_t, km, *p, info);
                            temp_fuzzy_seeds[seed_count++] = info;
                            if (seed_count == fuzzy_seeds_capacity) {
                                fuzzy_seeds_capacity *= 2;
                                *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                temp_fuzzy_seeds = *fuzzy_seeds;
                                if (!temp_fuzzy_seeds) return 0;
                            }
                        } else { // only addition of buf[j]
                            calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[j].x >> 14, mask32, blndK);
                        }
                    }
                }
        }

        if (f_info.x <= min.x) { // a new minimum; then write the old min
            if (l >= window + kmer_size && min.x != UINT64_MAX) {
                if ((min.y & 1) == 1) { // reverse strand
                    r_blndBuf[r_blendpos].x = min.x;
                    r_blndBuf[r_blendpos].y = (min.y & iMask) >> 1;
                    tot_span = r_blndBuf[r_blendpos].y;

                    if (++r_blendpos == n_neighbors) {
                        r_blendpos = 0;
                    }

                    if (++r_nm >= n_neighbors) {
                        blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, r_blndBuf[r_blendpos].x >> 14, mask32, blndK);
                        tot_span = tot_span - r_blndBuf[r_blendpos].y + (r_blndBuf[r_blendpos].x & sMask);
                        info.x = blendVal << 14 | tot_span;
                        info.y = min.y;
                        // kv_push(mm128_t, km, *p, info);
                        temp_fuzzy_seeds[seed_count++] = info;
                        if (seed_count == fuzzy_seeds_capacity) {
                            fuzzy_seeds_capacity *= 2;
                            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                            temp_fuzzy_seeds = *fuzzy_seeds;
                            if (!temp_fuzzy_seeds) return 0;
                        }
                    } else { //only addition of min
                        calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blndK);
                    }
                } else {
                    f_blndBuf[f_blendpos].x = min.x;
                    f_blndBuf[f_blendpos].y = (min.y & iMask) >> 1;
                    tot_span = f_blndBuf[f_blendpos].y;

                    if (++f_blendpos == n_neighbors) {
                        f_blendpos = 0;
                    }

                    if (++f_nm >= n_neighbors) {
                        blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, f_blndBuf[f_blendpos].x >> 14, mask32, blndK);
                        tot_span = tot_span - f_blndBuf[f_blendpos].y + (f_blndBuf[f_blendpos].x & sMask);
                        info.x = blendVal << 14 | tot_span;
                        info.y = min.y;
                        // kv_push(mm128_t, km, *p, info);
                        temp_fuzzy_seeds[seed_count++] = info;
                        if (seed_count == fuzzy_seeds_capacity) {
                            fuzzy_seeds_capacity *= 2;
                            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                            temp_fuzzy_seeds = *fuzzy_seeds;
                            if (!temp_fuzzy_seeds) return 0;
                        }
                    } else { //only addition of min
                        calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blndK);
                    }
                }
            }
            min = f_info, min_pos = buf_pos;
        } else if (buf_pos == min_pos) { // old min has moved outside the window
            if (l >= window + kmer_size - 1 && min.x != UINT64_MAX) {
                if ((min.y & 1) == 1) { //reverse strand
                    r_blndBuf[r_blendpos].x = min.x;
                    r_blndBuf[r_blendpos].y = (min.y & iMask) >> 1;
                    tot_span = r_blndBuf[r_blendpos].y;

                    if (++r_blendpos == n_neighbors) {
                        r_blendpos = 0;
                    }

                    if (++r_nm >= n_neighbors) {
                        blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, r_blndBuf[r_blendpos].x >> 14, mask32, blndK);
                        tot_span = tot_span - r_blndBuf[r_blendpos].y + (r_blndBuf[r_blendpos].x & sMask);
                        info.x = blendVal << 14 | tot_span;
                        info.y = min.y;
                        // kv_push(mm128_t, km, *p, info);
                        temp_fuzzy_seeds[seed_count++] = info;
                        if (seed_count == fuzzy_seeds_capacity) {
                            fuzzy_seeds_capacity *= 2;
                            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                            temp_fuzzy_seeds = *fuzzy_seeds;
                            if (!temp_fuzzy_seeds) return 0;
                        }
                    } else { //only addition of min
                        calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blndK);
                    }
                } else {
                    f_blndBuf[f_blendpos].x = min.x;
                    f_blndBuf[f_blendpos].y = (min.y & iMask) >> 1;
                    tot_span = f_blndBuf[f_blendpos].y;

                    if(++f_blendpos == n_neighbors) {
                        f_blendpos = 0;
                    }

                    if(++f_nm >= n_neighbors) {
                        blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, f_blndBuf[f_blendpos].x >> 14, mask32, blndK);
                        tot_span = tot_span - f_blndBuf[f_blendpos].y + (f_blndBuf[f_blendpos].x & sMask);
                        info.x = blendVal << 14 | tot_span;
                        info.y = min.y;
                        // kv_push(mm128_t, km, *p, info);
                        temp_fuzzy_seeds[seed_count++] = info;
                        if (seed_count == fuzzy_seeds_capacity) {
                            fuzzy_seeds_capacity *= 2;
                            *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                            temp_fuzzy_seeds = *fuzzy_seeds;
                            if (!temp_fuzzy_seeds) return 0;
                        }
                    } else { //only addition of min
                        calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blndK);
                    }
                }
            }
            for (j = buf_pos + 1, min.x = UINT64_MAX; j < window; ++j) { // the two loops are necessary when there are identical k-mers. This one starts from the oldest (first in: buf_pos + 1) 
                if (min.x >= buf[j].x) {
                    min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer... chosing the last one because we will store *all* identical minimizers and continue with the last one once all identical ones are stored. @IMPORTANT: should we really store all identical k-mers? it is true that they are the minimizers for their own window but not sure how effective it is
                }
            }
            for (j = 0; j <= buf_pos; ++j) {
                if (min.x >= buf[j].x) {
                    min = buf[j], min_pos = j;
                }
            }
            if (l >= window + kmer_size - 1 && min.x != UINT64_MAX) { // write identical k-mers
                for (j = buf_pos + 1; j < window; ++j) { // these two loops make sure the output is sorted
                    if (min.x == buf[j].x && min.y != buf[j].y) {
                        if ((buf[j].y & 1) == 1) { // reverse strand
                            r_blndBuf[r_blendpos].x = buf[j].x;
                            r_blndBuf[r_blendpos].y = (buf[j].y & iMask) >> 1;
                            tot_span = r_blndBuf[r_blendpos].y;

                            if (++r_blendpos == n_neighbors) {
                                r_blendpos = 0;
                            }

                            if (++r_nm >= n_neighbors) {
                                blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[j].x >> 14, r_blndBuf[r_blendpos].x >> 14, mask32, blndK);
                                tot_span = tot_span - r_blndBuf[r_blendpos].y + (r_blndBuf[r_blendpos].x & sMask);
                                info.x = blendVal << 14 | tot_span;
                                info.y = buf[j].y;
                                // kv_push(mm128_t, km, *p, info);
                                temp_fuzzy_seeds[seed_count++] = info;
                                if (seed_count == fuzzy_seeds_capacity) {
                                    fuzzy_seeds_capacity *= 2;
                                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                    temp_fuzzy_seeds = *fuzzy_seeds;
                                    if (!temp_fuzzy_seeds) return 0;
                                }
                            } else { //only addition of buf[j]
                                calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[j].x >> 14, mask32, blndK);
                            }
                        } else {
                            f_blndBuf[f_blendpos].x = buf[j].x;
                            f_blndBuf[f_blendpos].y = (buf[j].y & iMask) >> 1;
                            tot_span = f_blndBuf[f_blendpos].y;

                            if (++f_blendpos == n_neighbors) {
                                f_blendpos = 0;
                            }

                            if (++f_nm >= n_neighbors) {
                                blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[j].x >> 14, f_blndBuf[f_blendpos].x >> 14, mask32, blndK);
                                tot_span = tot_span - f_blndBuf[f_blendpos].y + (f_blndBuf[f_blendpos].x & sMask);
                                info.x = blendVal << 14 | tot_span;
                                info.y = buf[j].y;
                                // kv_push(mm128_t, km, *p, info);
                                temp_fuzzy_seeds[seed_count++] = info;
                                if (seed_count == fuzzy_seeds_capacity) {
                                    fuzzy_seeds_capacity *= 2;
                                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                    temp_fuzzy_seeds = *fuzzy_seeds;
                                    if (!temp_fuzzy_seeds) return 0;
                                }
                            } else { // only addition of buf[j]
                                calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[j].x >> 14, mask32, blndK);
                            }
                        }
                    }
                }
                for (j = 0; j <= buf_pos; ++j)
                    if (min.x == buf[j].x && min.y != buf[j].y) {
                        if ((buf[j].y & 1) == 1) { //reverse strand
                            r_blndBuf[r_blendpos].x = buf[j].x;
                            r_blndBuf[r_blendpos].y = (buf[j].y & iMask) >> 1;
                            tot_span = r_blndBuf[r_blendpos].y;

                            if (++r_blendpos == n_neighbors) { 
                                r_blendpos = 0;
                            }

                            if (++r_nm >= n_neighbors) {
                                blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[j].x >> 14, r_blndBuf[r_blendpos].x >> 14, mask32, blndK);
                                tot_span = tot_span - r_blndBuf[r_blendpos].y + (r_blndBuf[r_blendpos].x & sMask);
                                info.x = blendVal << 14 | tot_span;
                                info.y = buf[j].y;
                                // kv_push(mm128_t, km, *p, info);
                                temp_fuzzy_seeds[seed_count++] = info;
                                if (seed_count == fuzzy_seeds_capacity) {
                                    fuzzy_seeds_capacity *= 2;
                                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                    temp_fuzzy_seeds = *fuzzy_seeds;
                                    if (!temp_fuzzy_seeds) return 0;
                                }
                            } else { //only addition of buf[j]
                                calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, buf[j].x >> 14, mask32, blndK);
                            }
                        } else {
                            f_blndBuf[f_blendpos].x = buf[j].x;
                            f_blndBuf[f_blendpos].y = (buf[j].y & iMask) >> 1;
                            tot_span = f_blndBuf[f_blendpos].y;

                            if (++f_blendpos == n_neighbors) {
                                f_blendpos = 0;
                            }

                            if (++f_nm >= n_neighbors) {
                                blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[j].x >> 14, f_blndBuf[f_blendpos].x >> 14, mask32, blndK);
                                tot_span = tot_span - f_blndBuf[f_blendpos].y + (f_blndBuf[f_blendpos].x & sMask);
                                info.x = blendVal << 14 | tot_span;
                                info.y = buf[j].y;
                                // kv_push(mm128_t, km, *p, info);
                                temp_fuzzy_seeds[seed_count++] = info;
                                if (seed_count == fuzzy_seeds_capacity) {
                                    fuzzy_seeds_capacity *= 2;
                                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                                    temp_fuzzy_seeds = *fuzzy_seeds;
                                    if (!temp_fuzzy_seeds) return 0;
                                }
                            } else { //only addition of buf[j]
                                calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, buf[j].x>>14, mask32, blndK);
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
            r_blndBuf[r_blendpos].x = min.x;
            r_blndBuf[r_blendpos].y = (min.y & iMask) >> 1;
            tot_span = r_blndBuf[r_blendpos].y;

            if (++r_blendpos == n_neighbors) {
                r_blendpos = 0;
            }

            if (++r_nm >= n_neighbors) {
                blendVal = calc_blend_rm_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, r_blndBuf[r_blendpos].x >> 14, mask32, blndK);
                tot_span = tot_span - r_blndBuf[r_blendpos].y + (r_blndBuf[r_blendpos].x & sMask);
                min.x = blendVal << 14 | tot_span;
                // kv_push(mm128_t, km, *p, min);
                temp_fuzzy_seeds[seed_count++] = min;
                if (seed_count == fuzzy_seeds_capacity) {
                    fuzzy_seeds_capacity *= 2;
                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                    temp_fuzzy_seeds = *fuzzy_seeds;
                    if (!temp_fuzzy_seeds) return 0;
                }
            } else { // only addition of min
                calc_blend_simd(&r_blndcnt_lsb, &r_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blndK);
            }
        } else {
            f_blndBuf[f_blendpos].x = min.x;
            f_blndBuf[f_blendpos].y = (min.y & iMask) >> 1;
            tot_span = f_blndBuf[f_blendpos].y;

            if(++f_blendpos == n_neighbors) f_blendpos = 0;
            if(++f_nm >= n_neighbors) {
                blendVal = calc_blend_rm_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, f_blndBuf[f_blendpos].x >> 14, mask32, blndK);
                tot_span = tot_span - f_blndBuf[f_blendpos].y + (f_blndBuf[f_blendpos].x & sMask);
                min.x = blendVal << 14 | tot_span;
                // kv_push(mm128_t, km, *p, min);
                temp_fuzzy_seeds[seed_count++] = min;
                if (seed_count == fuzzy_seeds_capacity) {
                    fuzzy_seeds_capacity *= 2;
                    *fuzzy_seeds = (uint128_t *)realloc(temp_fuzzy_seeds, fuzzy_seeds_capacity * sizeof(uint128_t));
                    temp_fuzzy_seeds = *fuzzy_seeds;
                    if (!temp_fuzzy_seeds) return 0;
                }
            } else { //only addition of min
                calc_blend_simd(&f_blndcnt_lsb, &f_blndcnt_msb, &ma, &mb, min.x >> 14, mask32, blndK);
            }
        }
    }

    free(buf);
    free(f_blndBuf);
    free(r_blndBuf);

    return seed_count;
}

uint64_t blend_sketch(const char *str, int len, int window, int kmer_size, int blend_bits, int n_neighbors, uint32_t rid, uint128_t **fuzzy_seeds) {
	return blend_sb_sketch(str, len, window, kmer_size, blend_bits, (uint64_t)n_neighbors, rid, fuzzy_seeds);
}

}

#endif