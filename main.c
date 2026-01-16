#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <sys/stat.h>
#include "sketch/blend.h"
#include "utils/kseq.h"
#include "utils/khashl.h"
#include <zlib.h>


KHASHL_MAP_INIT(KH_LOCAL, map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)
KSEQ_INIT(gzFile, gzread)


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

int file_exists(char *filename) {
    struct stat buffer;   
    return (stat(filename, &buffer) == 0);
}

void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options]\n\n"
        "Options:\n"
        "  -f, --fasta <file>        FASTA file\n"
        "  -q, --fastq <file>        FASTQ file\n"
        "  -k, --kmer <int>          k-mer size [%d]\n"
        "  -w, --window <int>        window size [%d]\n"
        "  -b, --blend-bits <int>    number of bits for hash [%d]\n"
        "  -n, --n-neighbors <int>   number of neighbour [%d]\n"
        "  -p, --progress            Display progress\n"
        "  -h, --help                Show this help\n",
        prog, BLEND_K_SHORT, BLEND_W_SHORT, BLEND_BITS_SHORT, BLEND_NEIGHBOR_NUMBER_SHORT
    );
}

void init_params(params *p) {
    p->k = BLEND_K_SHORT;
    p->w = BLEND_W_SHORT;
    p->blend_bits = BLEND_BITS_SHORT;
    p->n_neighbors = BLEND_NEIGHBOR_NUMBER_SHORT;
    p->progress = PROGRESS;
}

void parse_args(int argc, char **argv, params *p) {
    static struct option long_opts[] = {
        {"fasta",       required_argument, 0, 'f'},
        {"fastq",       required_argument, 0, 'q'},
        {"kmer",        required_argument, 0, 'k'},
        {"window",      required_argument, 0, 'w'},
        {"blend-bits",  required_argument, 0, 'b'},
        {"n-neighbors", required_argument, 0, 'n'},
        {"progress",    no_argument,       0, 'p'},
        {"help",        no_argument,       0, 'h' },
        {0, 0, 0, 0}
    };

    int is_f_set = 0; 
    int is_q_set = 0;

    int opt, idx;
    while ((opt = getopt_long(argc, argv, "f:q:k:w:b:n:ph", long_opts, &idx)) != -1) {
        switch (opt) {
        case 'f':
            p->fasta = optarg;
            is_f_set = 1;
            break;
        case 'q':
            p->fastq = optarg;
            is_q_set = 1;
            break;
        case 'k':
            p->k = atoi(optarg);
            break;
        case 'w':
            p->w = atoi(optarg);
            break;
        case 'b':
            p->blend_bits = atoi(optarg);
            break;
        case 'n':
            p->n_neighbors = atoi(optarg);
            break;
        case 'p':
            p->progress = 1;
            break;
        case 0:
            print_usage(argv[0]);
            exit(0);
        default:
            print_usage(argv[0]);
            exit(1);
        }
    }

    /* sanity checks */
    if (p->k <= 0 || p->w <= 0) {
        fprintf(stderr, "Error: k and w must be > 0\n");
        exit(1);
    }

    if (!is_f_set) {
        fprintf(stderr, "Error: fasta is not provided\n");
        print_usage(argv[0]);
        exit(1);
    }

    if (!is_q_set) {
        fprintf(stderr, "Error: query is not provided\n");
        print_usage(argv[0]);
        exit(1);
    }

    if (!file_exists(p->fasta)) {
        fprintf(stderr, "Error: fasta file does not exists\n");
        exit(1);
    }

    if (!file_exists(p->fastq)) {
        fprintf(stderr, "Error: fastq file does not exists\n");
        exit(1);
    }
}

int cmp_fuzzy_seeds(const void *a, const void *b) {
    const uint128_t *x = (uint128_t *)a;
    const uint128_t *y = (uint128_t *)b;
    return (x->x >> 14) > (y->x >> 14) ? 1 : -1;
}


int store_seqs(const char *path, ref_seq **seqs) {

    char *fai_path = (char *)malloc(strlen(path) + 5);
    if (fai_path == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        exit(-1);
    }
    sprintf(fai_path, "%s.fai", path);

    size_t line_cap = 1024;
    char *line = (char *)malloc(line_cap);
    int chrom_index = 0;

    FILE *fai = fopen(fai_path, "r");
    if (!fai) {
        fprintf(stderr, "Error: Couldn't open %s\n", fai_path);
        exit(-1);
    }

    while (fgets(line, line_cap, fai)) {
        chrom_index++;
    }

    if (chrom_index == 0) {
        fprintf(stderr, "Error: Index file is empty.\n");
        exit(-1);
    }
    
    *seqs = (ref_seq *)malloc(sizeof(ref_seq) * chrom_index);
    if (!(*seqs)) {
        fprintf(stderr, "Error: Couldn't allocate memory to ref sequences\n");
        exit(-1);
    }

    rewind(fai);
    chrom_index = 0;

    while (fgets(line, line_cap, fai)) {
        char *name, *length;

        // assign name
        char *saveptr;
        name = strtok_r(line, "\t", &saveptr);
        uint64_t name_len = strlen(name);
        (*seqs)[chrom_index].header = (char *)malloc(name_len+1);
        memcpy((*seqs)[chrom_index].header, name, name_len);
        (*seqs)[chrom_index].header[name_len] = '\0';

        // assign size and allocate in memory
        length = strtok_r(NULL, "\t", &saveptr);
        (*seqs)[chrom_index].len = strtol(length, NULL, 10);
        chrom_index++;
    }

    fclose(fai);
    free(fai_path);
    free(line);

    return chrom_index;
}

uint64_t process_fasta(params *p, uint128_t **fuzzy_seeds, map32_t **index_table, ref_seq **seqs, int *chrom_count) {

    *chrom_count = store_seqs(p->fasta, seqs);

    uint128_t *all_fuzzy_seeds;
    uint64_t all_fuzzy_seeds_len = 0;
    uint64_t all_fuzzy_seeds_cap = SKETCH_CAPACITY;

    all_fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * all_fuzzy_seeds_cap);
    if (!all_fuzzy_seeds) {
        fprintf(stderr, "Error: couldn't allocate array\n");
        exit(1);
    }

    gzFile fp = gzopen(p->fasta, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open reference %s\n", p->fasta);
        exit(1);
    }

    kseq_t *seq = kseq_init(fp);

    int chrom_index = 0;

    while (kseq_read(seq) >= 0) {
        const char *name = seq->name.s;
        const char *bases = seq->seq.s;
        int len = seq->seq.l;

        (*seqs[chrom_index]).chrom = (char *)malloc(len);
        memcpy((*seqs[chrom_index]).chrom, bases, len);

        uint128_t *chr_fuzzy_seeds;
        uint64_t chr_fuzzy_seeds_len = 0;
        chr_fuzzy_seeds_len = blend_sketch(bases, len, p->w, p->k, p->blend_bits, p->n_neighbors, chrom_index, &chr_fuzzy_seeds);

        if (chr_fuzzy_seeds_len) {
            if (all_fuzzy_seeds_cap <= all_fuzzy_seeds_len + chr_fuzzy_seeds_len) {
                while (all_fuzzy_seeds_cap <= all_fuzzy_seeds_len + chr_fuzzy_seeds_len) {
                    all_fuzzy_seeds_cap *= 2;
                }
                uint128_t *temp = (uint128_t *)realloc(all_fuzzy_seeds, sizeof(uint128_t) * all_fuzzy_seeds_cap);
                if (!temp) {
                    fprintf(stderr, "Error: couldn't reallocate array\n");
                    exit(1);
                }
                all_fuzzy_seeds = temp;
            }
            
            memcpy(all_fuzzy_seeds + all_fuzzy_seeds_len, chr_fuzzy_seeds, sizeof(uint128_t) * chr_fuzzy_seeds_len);
            all_fuzzy_seeds_len += chr_fuzzy_seeds_len;

            free(chr_fuzzy_seeds);
        }

        chrom_index++;
    }
    // clean-up fasta reading 
    kseq_destroy(seq);
    gzclose(fp);

    // get unique fuzzy_seeds
    qsort(all_fuzzy_seeds, all_fuzzy_seeds_len, sizeof(uint128_t), cmp_fuzzy_seeds);
    
    uint64_t unique_fuzzy_seeds_len = 0;

    for (uint64_t i = 0; i < all_fuzzy_seeds_len; ) {
        uint64_t j = i + 1;

        while (j < all_fuzzy_seeds_len && (all_fuzzy_seeds[j].x >> 14) == (all_fuzzy_seeds[i].x >> 14)) {
            j++;
        }

        if (j == i + 1) {
            if (i != unique_fuzzy_seeds_len) {
                uint128_t tmp = all_fuzzy_seeds[unique_fuzzy_seeds_len];
                all_fuzzy_seeds[unique_fuzzy_seeds_len] = all_fuzzy_seeds[i];
                all_fuzzy_seeds[i] = tmp;
            }
            unique_fuzzy_seeds_len++;
        }

        i = j;
    }
    
    printf("Unique %lu / %lu (%.2f)\n", unique_fuzzy_seeds_len, all_fuzzy_seeds_len, (double)unique_fuzzy_seeds_len / (double)all_fuzzy_seeds_len);

    all_fuzzy_seeds_len = unique_fuzzy_seeds_len;

    map32_t *temp_index_table = *index_table;
    for (uint64_t i = 0; i < all_fuzzy_seeds_len; i++) {
        khint_t k; int absent;
        k = map32_put(temp_index_table, (all_fuzzy_seeds[i].x >> 14), &absent);
        kh_val(temp_index_table, k) = i;
    }

    printf("Log: Indexed reference\n");

    if (all_fuzzy_seeds) {
        *fuzzy_seeds = all_fuzzy_seeds;
    }

    return all_fuzzy_seeds_len;
}

static inline char rc_base(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

void reverse_complement(const char *in, char *out, int len) {
    for (int i = 0; i < len; i++)
        out[len - i - 1] = rc_base(in[i]);
}

int banded_align_and_report(const char *ref, uint64_t ref_span, const char *read, uint64_t read_span, uint64_t ref_pos, uint64_t ref_id, uint64_t read_id, int ref_strand, int read_strand) {
    
    int ref_len = ref_span;
    int read_len = read_span;

    if (abs(ref_len - read_len) > BAND) return 0; // why bother?

    char ref_buf[MAX_LEN];
    char read_buf[MAX_LEN];

    // Strand normalization
    if (ref_strand) reverse_complement(ref, ref_buf, ref_len);
    else memcpy(ref_buf, ref, ref_len);
    
    if (read_strand) reverse_complement(read, read_buf, read_len);
    else memcpy(read_buf, read, read_len);

    static int dp[MAX_LEN + 1][2 * BAND + 1];
    static int bt[MAX_LEN + 1][2 * BAND + 1];

    // init
    for (int i = 0; i <= ref_len; i++)
        for (int k = 0; k <= 2 * BAND; k++)
            dp[i][k] = NEG_INF;

    dp[0][BAND] = 0;

    // DP
    for (int i = 1; i <= ref_len; i++) {

        int j_start = 0 < i - BAND ? i - BAND : 0;
        int j_end = read_len < i + BAND ? read_len : i + BAND;

        for (int j = j_start; j <= j_end; j++) {

            int k = j - i + BAND;
            if (k < 0 || 2 * BAND < k) continue;

            int best = NEG_INF, op = -1;

            // diagonal
            if (0 < j && dp[i-1][k] != NEG_INF) {
                int d = dp[i-1][k] + (ref_buf[i-1] == read_buf[j-1] ? MATCH : MISMATCH);
                if (d > best) { best = d; op = 0; }
            }

            // deletion
            if (k + 1 <= 2 * BAND && dp[i-1][k+1] != NEG_INF) {
                int d = dp[i-1][k+1] + GAP;
                if (d > best) { best = d; op = 1; }
            }

            // insertion
            if (0 < j && k > 0 && dp[i][k-1] != NEG_INF) {
                int d = dp[i][k-1] + GAP;
                if (d > best) { best = d; op = 2; }
            }

            dp[i][k] = best;
            bt[i][k] = op;
        }
    }

    int end_k = read_len - ref_len + BAND;
    if (end_k < 0 || end_k > 2 * BAND) return 0;

    // alignment acceptance
    if (dp[ref_len][end_k] < 0) return 0;

    // traceback
    int i = ref_len, j = read_len, k = end_k;
    int mismatches = 0;
#ifdef __DEBUG
    char aln_ref[2*MAX_LEN];
    char aln_read[2*MAX_LEN];
    char aln_mid[2*MAX_LEN];
    int aln_len = 0;
#endif
    while (i > 0 && j > 0) {
        int op = bt[i][k];

        if (op == 0) {
#ifdef __DEBUG
            aln_ref[aln_len]  = ref_buf[i-1];
            aln_read[aln_len] = read_buf[j-1];
            aln_mid[aln_len]  = (ref_buf[i-1] == read_buf[j-1]) ? '|' : '*';
#endif
            if (ref_buf[i-1] != read_buf[j-1]) {
                mismatches++;
                // printf("SNP\tREFID=%lu\tREADID=%lu\tPOS=%lu\tREF=%c\tALT=%c\n", ref_id, read_id, ref_pos + i - 1, ref_buf[i-1], read_buf[j-1]);
            }
            i--; j--;
        }
        else if (op == 1) {
#ifdef __DEBUG
            aln_ref[aln_len]  = ref_buf[i-1];
            aln_read[aln_len] = '-';
            aln_mid[aln_len]  = ' ';
#endif
            // printf("DEL\tREFID=%lu\tREADID=%lu\tPOS=%lu\tBASE=%c\n", ref_id, read_id, ref_pos + i - 1, ref_buf[i-1]);
            i--; k++;
        }
        else {
#ifdef __DEBUG
            aln_ref[aln_len]  = '-';
            aln_read[aln_len] = read_buf[j-1];
            aln_mid[aln_len]  = ' ';
#endif
            // printf("INS\tREFID=%lu\tREADID=%lu\tPOS=%lu\tBASE=%c\n", ref_id, read_id, ref_pos + i, read_buf[j-1]);
            j--; k--;
        }

        // if (mismatches > 20) return 0;
#ifdef __DEBUG
        aln_len++;
#endif
    }
#ifdef __DEBUG

    // reverse alignment strings
    for (int x = 0; x < aln_len / 2; x++) {
        char t;

        t = aln_ref[x];
        aln_ref[x] = aln_ref[aln_len - 1 - x];
        aln_ref[aln_len - 1 - x] = t;

        t = aln_read[x];
        aln_read[x] = aln_read[aln_len - 1 - x];
        aln_read[aln_len - 1 - x] = t;

        t = aln_mid[x];
        aln_mid[x] = aln_mid[aln_len - 1 - x];
        aln_mid[aln_len - 1 - x] = t;
    }

    aln_ref[aln_len]  = '\0';
    aln_read[aln_len] = '\0';
    aln_mid[aln_len]  = '\0';

    // print alignment
    printf("\n");
    printf("REF : %s (%d)\n", aln_ref, ref_len);
    printf("      %s\n", aln_mid);
    printf("READ: %s (%d)\n", aln_read, read_len);

    printf("--------------------------------------------------------------------------\n");

    // print seqs
    printf("REF:  ");
    fwrite(ref, 1, ref_len, stdout);
    printf("\n");

    printf("READ: ");
    fwrite(read, 1, read_len, stdout);
    printf("\n");
#endif
    return mismatches;
}


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
    ref_seq *seqs;
    int seq_count = 0;
    uint64_t fuzzy_seeds_len = process_fasta(&p, &fuzzy_seeds, &index_table, &seqs, &seq_count);

    gzFile fp = gzopen(p.fastq, "r");
    if (!fp) {
        fprintf(stderr, "Error: cannot open reads %s\n", p.fastq);
        exit(1);
    }

    kseq_t *seq = kseq_init(fp);

    int found = 0, not_found = 0, mismatches = 0;
    uint64_t total_len = 0, read_id = 0;

    while (kseq_read(seq) >= 0) {
        const char *name = seq->name.s;
        const char *bases = seq->seq.s;
        int len = seq->seq.l;

        uint128_t *temp_fuzzy_seeds;
        uint64_t temp_fuzzy_seeds_len = 0;
        temp_fuzzy_seeds_len = blend_sketch(bases, len, p.w, p.k, p.blend_bits, p.n_neighbors, 0, &temp_fuzzy_seeds);

        // store in set
        for (uint64_t i = 0; i < temp_fuzzy_seeds_len; i++) {
            khint_t k = map32_get(index_table, temp_fuzzy_seeds[i].x >> 14);
            
            if (k < kh_end(index_table)) {
                found++;

                uint32_t seed_index = kh_val(index_table, k);
                // uint64_t ref_kmer = BLEND_GET_KMER(fuzzy_seeds[seed_index]);
                uint64_t ref_span = BLEND_GET_LENGTH(fuzzy_seeds[seed_index]);
                uint64_t ref_index = BLEND_GET_INDEX(fuzzy_seeds[seed_index]);
                uint64_t rid = BLEND_GET_REFERENCE_IDX(fuzzy_seeds[seed_index]);
                int ref_strand = BLEND_GET_STRAND(fuzzy_seeds[seed_index]);
                
                uint64_t read_span = BLEND_GET_LENGTH(temp_fuzzy_seeds[i]);
                uint64_t read_index = BLEND_GET_INDEX(temp_fuzzy_seeds[i]);
                int read_strand = BLEND_GET_STRAND(temp_fuzzy_seeds[i]);

                // get substrings
                const char *ref = seqs[rid].chrom + ref_index;
                const char *read = bases + read_index;

                // align and report variants in alignment
                int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, ref_index, rid, read_id, ref_strand, read_strand);
#ifdef __DEBUG
                if (alignment_mismatches && ref_span != read_span) printf("ref_span: %lu, read_span: %lu\n", ref_span, read_span);
#endif
                mismatches += alignment_mismatches;
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

    for (int i = 0; i < seq_count; i++) {
        free(seqs[i].header);
        free(seqs[i].chrom);
    }
    if (seq_count) free(seqs);

    return 0;
}
