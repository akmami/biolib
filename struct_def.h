#ifndef __STRUCT_DEF_H__
#define __STRUCT_DEF_H__

#include <stdlib.h>
#include <stdint.h>

#define BLEND_K_SHORT 21
#define BLEND_W_SHORT 11
#define BLEND_BITS_SHORT 32
#define BLEND_NEIGHBOR_NUMBER_SHORT 5
#define BLEND_K_HIFI 19
#define BLEND_W_HIFI 50
#define BLEND_BITS_HIFI 38
#define BLEND_NEIGHBOR_NUMBER_HIFI 5
#define SKETCH_CAPACITY 40000000

#define PROGRESS 0

#define MATCH     2
#define MISMATCH -5
#define GAP      -3
#define BAND 20        // ±8 bp band
#define MAX_LEN 2048   // max extension length
#define NEG_INF -100000000


typedef struct {
    char *fasta;
    char *fastq;
    int k;
    int w;
    int blend_bits;
    int n_neighbors;
    int progress;
} params;

typedef struct {
    char *header;
    char *chrom;
    uint64_t len;
} ref_seq;

#endif