#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <zlib.h>
#include "struct_def.h"
#include "sketch/blend.h"
#include "utils/kseq.h"
#include "utils/khashl.h"


KHASHL_MAP_INIT(KH_LOCAL, map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)
KSEQ_INIT(gzFile, gzread)


int store_seqs(const char *path, ref_seq **seqs);

uint64_t process_fasta(params *p, uint128_t **fuzzy_seeds, uint64_t *fuzzy_seeds_len, map32_t **index_table, ref_seq **seqs, int *chrom_count);

int banded_align_and_report(const char *ref, uint64_t ref_span, const char *read, uint64_t read_span, uint64_t ref_pos, uint64_t ref_id, uint64_t read_id, int ref_strand, int read_strand);

#endif