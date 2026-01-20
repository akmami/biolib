#ifndef __STRUCT_DEF__
#define __STRUCT_DEF__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

unsigned char seq_nt4_table[256] = {
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

#ifndef __UINT128_T__
#define __UINT128_T__
typedef struct {
	uint64_t x; // kmer(32) + reference(32)
	uint64_t y; // length(31) + strand(1) + index(32)
} uint128_t;
#endif

#ifndef __UINT160_T__
#define __UINT160_T__
typedef struct {
	uint64_t x;
	uint64_t y;
	uint32_t z;
} uint160_t;
#endif

#define MASK_HEAD_64    0xFFFFFFFF00000000ULL
#define MASK_TAIL_64    0x00000000FFFFFFFFULL
#define MASK_STRAND_64  0xFFFFFFFEFFFFFFFFULL
#define MASK_LENGTH_64  0x00000001FFFFFFFFULL

#define SET_128_KMER(minimizer, id) ((minimizer).x = ((minimizer).x & MASK_TAIL_64) | ((uint64_t)(id) << 32))
#define SET_128_REFERENCE(minimizer, reference) ((minimizer).x = ((minimizer).x & MASK_HEAD_64) | ((uint32_t)(reference)))
#define SET_128_LEN(minimizer, length) ((minimizer).y = ((minimizer).y & MASK_LENGTH_64) | ((uint64_t)(length) << 33))
#define SET_128_STRAND(minimizer, strand) ((minimizer).y = ((minimizer).y & MASK_STRAND_64) | (((uint64_t)(strand & 1)) << 32))
#define SET_128_INDEX(minimizer, index) ((minimizer).y = ((minimizer).y & MASK_HEAD_64) | ((uint32_t)(index)))

#define SET_128_X(minimizer, id, reference) ((minimizer).x = (((uint64_t)(id) << 32) | (uint32_t)(reference)))
#define SET_128_Y(minimizer, length, strand, index) ((minimizer).y = (((uint64_t)(length) << 33) | ((uint64_t)(strand & 1) << 32) | (uint32_t)(index)))

#define GET_128_KMER(minimizer) ((minimizer).x >> 32)
#define GET_128_REFERENCE(minimizer) ((minimizer).x & MASK_TAIL_64)
#define GET_128_LEN(minimizer) ((minimizer).y >> 33) 
#define GET_128_STRAND(minimizer) (((minimizer).y >> 32) & 1)
#define GET_128_INDEX(minimizer) ((minimizer).y & MASK_TAIL_64)

#define SET_160_KMER(minimizer, hash) ((minimizer).x = hash)
#define SET_160_INDEX(minimizer, index) ((minimizer).y |= (index << 32))
#define SET_160_LEN(minimizer, len) ((minimizer).y |= 0xFFFFFFFF)
#define SET_160_REFERENCE(minimizer, ref) ((minimizer).z |= (ref << 1))
#define SET_160_STRAND(minimizer, strand) ((minimizer).z |= (strand & 1))

#define GET_160_KMER(minimizer) ((minimizer).x)
#define GET_160_INDEX(minimizer) ((minimizer).y >> 32)
#define GET_160_LEN(minimizer) ((minimizer).y & 0xFFFFFFFF)
#define GET_160_REFERENCE(minimizer) ((minimizer).z >> 1)
#define GET_160_STRAND(minimizer) ((minimizer).z & 1)

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

static inline uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed) {
    const uint8_t *data = (const uint8_t *)key;
    const int nblocks = len / 4;

    uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    // Body: Process blocks of 4 bytes at a time
    const uint32_t *blocks = (const uint32_t *)(data + nblocks * 4);

    for (int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;

        h1 ^= k1;
        h1 = (h1 << 13) | (h1 >> (32 - 13));
        h1 = h1 * 5 + 0xe6546b64;
    }

    // Tail: Process remaining bytes
    const uint8_t *tail = (const uint8_t *)(data + nblocks * 4);

    uint32_t k1 = 0;

    switch (len & 3) {
    case 3:
        k1 ^= tail[2] << 16;
        break;
    case 2:
        k1 ^= tail[1] << 8;
        break;
    case 1:
        k1 ^= tail[0];
        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;
        h1 ^= k1;
    }

    // Finalization: Mix the hash to ensure the last few bits are fully mixed
    h1 ^= len;

    /* fmix32 */
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    return h1;
}

#ifdef __cplusplus
}
#endif

#endif