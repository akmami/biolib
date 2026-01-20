#include <stdlib.h>
#include "struct_def.h"
#include "opt_parser.h"
#include "utils.h"


int main(int argc, char **argv) {

    params p;
    init_params(&p);
    parse_args(argc, argv, &p);

    fprintf(stderr,
        "Params:\n"
        "  fasta       = %s\n"
        "  fastq       = %s\n"
        "  k           = %d\n"
        "  w           = %d\n"
        "  blend-bits  = %d\n"
        "  n-neighbors = %d\n",
        p.fasta, p.fastq, p.k, p.w,
        p.blend_bits, p.n_neighbors
    );

    map32_t *index_table = map32_init();
    uint128_t *fuzzy_seeds;
    uint64_t relaxed_fuzzy_seeds_len;
    ref_seq *seqs;
    int seq_count = 0;
    uint64_t unique_fuzzy_seeds_len = process_fasta(&p, &fuzzy_seeds, &relaxed_fuzzy_seeds_len, &index_table, &seqs, &seq_count);

    gzFile fp = gzopen(p.fastq, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open reads %s\n", p.fastq);
        exit(1);
    }

    kseq_t *seq = kseq_init(fp);

    int found = 0, not_found = 0, mismatches = 0;
    uint64_t total_len = 0, read_id = 0;

    while (kseq_read(seq) >= 0) {
        const char *bases = seq->seq.s;
        int len = seq->seq.l;

        uint128_t *temp_fuzzy_seeds;
        uint64_t temp_fuzzy_seeds_len = 0;
        temp_fuzzy_seeds_len = blend_sketch(bases, len, p.w, p.k, p.blend_bits, p.n_neighbors, 0, &temp_fuzzy_seeds);

        // store in set
        for (uint64_t i = 0; i < temp_fuzzy_seeds_len; i++) {
            khint_t k = map32_get(index_table, BLEND_GET_KMER(temp_fuzzy_seeds[i]));
            
            if (k < kh_end(index_table)) {

                found++;

                uint32_t seed_index = kh_val(index_table, k);

                if (seed_index < unique_fuzzy_seeds_len) { // check if unique
                    // // uint64_t ref_kmer = BLEND_GET_KMER(fuzzy_seeds[seed_index]);
                    // uint64_t ref_span = BLEND_GET_LENGTH(fuzzy_seeds[seed_index]);
                    // uint64_t ref_index = BLEND_GET_INDEX(fuzzy_seeds[seed_index]);
                    // uint64_t ref_id = BLEND_GET_REFERENCE_IDX(fuzzy_seeds[seed_index]);
                    // int ref_strand = BLEND_GET_STRAND(fuzzy_seeds[seed_index]);

                    // uint64_t read_span = BLEND_GET_LENGTH(temp_fuzzy_seeds[i]);
                    // uint64_t read_index = BLEND_GET_INDEX(temp_fuzzy_seeds[i]);
                    // int read_strand = BLEND_GET_STRAND(temp_fuzzy_seeds[i]);

                    // // get substrings
                    // const char *ref = seqs[ref_id].chrom + ref_index;
                    // const char *read = bases + read_index;

                    // // align and report variants in alignment
                    // int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, ref_index, ref_id, read_id, ref_strand, read_strand);
                    // // if (alignment_mismatches && ref_span != read_span) printf("ref_span: %lu, read_span: %lu\n", ref_span, read_span);
                    // mismatches += alignment_mismatches;

                    continue;
                }

                if (i) { // check left
                    k = map32_get(index_table, BLEND_GET_KMER(temp_fuzzy_seeds[i-1]));

                    if (k < kh_end(index_table)) {

                        uint32_t left_seed_index = kh_val(index_table, k);

                        if (left_seed_index < unique_fuzzy_seeds_len) {

                            for (uint64_t temp_seed_index = seed_index; seed_index < relaxed_fuzzy_seeds_len; temp_seed_index++) {
                                if (BLEND_GET_KMER(fuzzy_seeds[temp_seed_index]) != BLEND_GET_KMER(temp_fuzzy_seeds[i])) {
                                    break;
                                }
                                if (BLEND_GET_REFERENCE_IDX(fuzzy_seeds[temp_seed_index]) == BLEND_GET_REFERENCE_IDX(fuzzy_seeds[left_seed_index]) &&
                                    abs_diff(BLEND_GET_INDEX(fuzzy_seeds[temp_seed_index]), BLEND_GET_INDEX(fuzzy_seeds[left_seed_index])) < DISTANCE_THRESHOLD) {
                                    // found
                                    uint64_t ref_span = BLEND_GET_LENGTH(fuzzy_seeds[temp_seed_index]);
                                    uint64_t ref_index = BLEND_GET_INDEX(fuzzy_seeds[temp_seed_index]);
                                    uint64_t ref_id = BLEND_GET_REFERENCE_IDX(fuzzy_seeds[temp_seed_index]);
                                    int ref_strand = BLEND_GET_STRAND(fuzzy_seeds[temp_seed_index]);

                                    uint64_t read_span = BLEND_GET_LENGTH(temp_fuzzy_seeds[i]);
                                    uint64_t read_index = BLEND_GET_INDEX(temp_fuzzy_seeds[i]);
                                    int read_strand = BLEND_GET_STRAND(temp_fuzzy_seeds[i]);

                                    // get substrings
                                    const char *ref = seqs[ref_id].chrom + ref_index;
                                    const char *read = bases + read_index;

                                    // align and report variants in alignment
                                    int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, ref_index, ref_id, read_id, ref_strand, read_strand);
                                    // if (alignment_mismatches && ref_span != read_span) printf("ref_span: %lu, read_span: %lu\n", ref_span, read_span);
                                    mismatches += alignment_mismatches;
                                }
                            }
                        }
                    }
                } else if (i < temp_fuzzy_seeds_len -1) { // check right
                    k = map32_get(index_table, BLEND_GET_KMER(temp_fuzzy_seeds[i+1]));

                    if (k < kh_end(index_table)) {

                        uint32_t right_seed_index = kh_val(index_table, k);
                        
                        if (right_seed_index < unique_fuzzy_seeds_len) {

                            for (uint64_t temp_seed_index = seed_index; seed_index < relaxed_fuzzy_seeds_len; temp_seed_index++) {
                                if (BLEND_GET_REFERENCE_IDX(fuzzy_seeds[temp_seed_index]) != BLEND_GET_KMER(temp_fuzzy_seeds[i])) {
                                    break;
                                }
                                if (BLEND_GET_REFERENCE_IDX(fuzzy_seeds[temp_seed_index]) == BLEND_GET_REFERENCE_IDX(fuzzy_seeds[right_seed_index]) &&
                                    abs_diff(BLEND_GET_INDEX(fuzzy_seeds[temp_seed_index]), BLEND_GET_INDEX(fuzzy_seeds[right_seed_index])) < DISTANCE_THRESHOLD) {
                                    // found
                                    uint64_t ref_span = BLEND_GET_LENGTH(fuzzy_seeds[temp_seed_index]);
                                    uint64_t ref_index = BLEND_GET_INDEX(fuzzy_seeds[temp_seed_index]);
                                    uint64_t ref_id = BLEND_GET_REFERENCE_IDX(fuzzy_seeds[temp_seed_index]);
                                    int ref_strand = BLEND_GET_STRAND(fuzzy_seeds[temp_seed_index]);

                                    uint64_t read_span = BLEND_GET_LENGTH(temp_fuzzy_seeds[i]);
                                    uint64_t read_index = BLEND_GET_INDEX(temp_fuzzy_seeds[i]);
                                    int read_strand = BLEND_GET_STRAND(temp_fuzzy_seeds[i]);

                                    // get substrings
                                    const char *ref = seqs[ref_id].chrom + ref_index;
                                    const char *read = bases + read_index;

                                    // align and report variants in alignment
                                    int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, ref_index, ref_id, read_id, ref_strand, read_strand);
                                    // if (alignment_mismatches && ref_span != read_span) printf("ref_span: %lu, read_span: %lu\n", ref_span, read_span);
                                    mismatches += alignment_mismatches;
                                }
                            }
                        }
                    }
                }
            } else {
                not_found++;
            }
        }

        if (temp_fuzzy_seeds_len) free(temp_fuzzy_seeds);

        total_len += temp_fuzzy_seeds_len;

        read_id++;
    }

    kseq_destroy(seq);
    gzclose(fp);

    printf("Found: %d, not found: %d, mismatches: %d, total: %lu\n", found, not_found, mismatches, total_len);

    // sort reported variatns and take consensus and finalize report

    // cleanup
    map32_destroy(index_table);
    free_seqs(&seqs, seq_count);

    if (relaxed_fuzzy_seeds_len) free(fuzzy_seeds);

    return 0;
}
