#include "utils.h"


int cmp_fuzzy_seeds(const void *a, const void *b) {
    const uint128_t *x = (uint128_t *)a;
    const uint128_t *y = (uint128_t *)b;
    return (x->x >> 14) > (y->x >> 14) ? 1 : -1;
}

static inline char rc_base(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 'T';
        case 'c': return 'G';
        case 'g': return 'C';
        case 't': return 'A';
        default:  return 'N';
    }
}

void reverse_complement(const char *in, char *out, int len) {
    for (int i = 0; i < len; i++)
        out[len - i - 1] = rc_base(in[i]);
}

void free_seqs(ref_seq **seqs, int seq_len) {
    if (seq_len) {
        ref_seq *temp = *seqs;
        for (int i = 0; i < seq_len; i++) {
            if (temp[i].len) {
                free(temp[i].chrom);
            }
            if (temp[i].header) free(temp[i].header);
        }

        free(temp);
    }
}
int store_seqs(const char *path, ref_seq **seqs) {

    char *fai_path = (char *)malloc(strlen(path) + 5);
    if (fai_path == NULL) {
        fprintf(stderr, "[ERROR] Memory allocation failed\n");
        exit(-1);
    }
    sprintf(fai_path, "%s.fai", path);

    size_t line_cap = 1024;
    char *line = (char *)malloc(line_cap);
    int chrom_index = 0;

    FILE *fai = fopen(fai_path, "r");
    if (!fai) {
        fprintf(stderr, "[ERROR] Couldn't open %s\n", fai_path);
        exit(-1);
    }

    while (fgets(line, line_cap, fai)) {
        chrom_index++;
    }

    if (chrom_index == 0) {
        fprintf(stderr, "[ERROR] Index file is empty.\n");
        exit(-1);
    }
    
    *seqs = (ref_seq *)malloc(sizeof(ref_seq) * chrom_index);
    if (!(*seqs)) {
        fprintf(stderr, "[ERROR] Couldn't allocate memory to ref sequences\n");
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

uint64_t process_fasta(params *p, uint128_t **fuzzy_seeds, uint64_t *fuzzy_seeds_len, map32_t **index_table, ref_seq **seqs, int *chrom_count) {

    *chrom_count = store_seqs(p->fasta, seqs);
    *fuzzy_seeds_len = 0;

    uint128_t *all_fuzzy_seeds;
    uint64_t all_fuzzy_seeds_len = 0;
    uint64_t all_fuzzy_seeds_cap = SKETCH_CAPACITY;

    all_fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * all_fuzzy_seeds_cap);
    if (!all_fuzzy_seeds) {
        fprintf(stderr, "[ERROR] couldn't allocate array\n");
        exit(1);
    }

    gzFile fp = gzopen(p->fasta, "r");
    if (!fp) {
        fprintf(stderr, "[ERROR] cannot open reference %s\n", p->fasta);
        exit(1);
    }

    kseq_t *seq = kseq_init(fp);

    int chrom_index = 0;

    while (kseq_read(seq) >= 0) {
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
                    fprintf(stderr, "[ERROR] couldn't reallocate array\n");
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

    // if no seeds are found, then no need to proceed
    if (!all_fuzzy_seeds_len) {
        free(all_fuzzy_seeds);
        return 0;
    }

    uint128_t *temp_all_fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * all_fuzzy_seeds_len);
    if (!temp_all_fuzzy_seeds) {
        fprintf(stderr, "[ERROR] couldn't allocate array\n");
    }
    memcpy(temp_all_fuzzy_seeds, all_fuzzy_seeds, sizeof(uint128_t) * all_fuzzy_seeds_len);

    // get unique fuzzy_seeds
    qsort(all_fuzzy_seeds, all_fuzzy_seeds_len, sizeof(uint128_t), cmp_fuzzy_seeds);
    
    uint64_t unique_fuzzy_seeds_len = 0;

    for (uint64_t i = 0; i < all_fuzzy_seeds_len; ) {
        uint64_t j = i + 1;

        while (j < all_fuzzy_seeds_len && __blend_get_kmer(all_fuzzy_seeds[j]) == __blend_get_kmer(all_fuzzy_seeds[i])) {
            j++;
        }

        if (j == i + 1) {
            all_fuzzy_seeds[unique_fuzzy_seeds_len++] = all_fuzzy_seeds[i];
        }

        i = j;
    }
    
    map32_t *temp_index_table = *index_table;
    for (uint64_t i = 0; i < unique_fuzzy_seeds_len; i++) {
        khint_t k; int absent;
        k = map32_put(temp_index_table, __blend_get_kmer(all_fuzzy_seeds[i]), &absent);
        kh_val(temp_index_table, k) = i;
    }

    if (all_fuzzy_seeds) {
        *fuzzy_seeds = all_fuzzy_seeds;
    }

    uint64_t relaxed_fuzzy_seeds_len = unique_fuzzy_seeds_len;

    // todo: handle i=0 and i=all_fuzzy_seeds_len-1
    for (uint64_t i = 1; i < all_fuzzy_seeds_len - 1; i++) {
        khint_t k = map32_get(temp_index_table, __blend_get_kmer(temp_all_fuzzy_seeds[i]));
        
        if (k < kh_end(temp_index_table)) { // if unique
            continue;
        }

        k = map32_get(temp_index_table, __blend_get_kmer(temp_all_fuzzy_seeds[i-1]));

        if (k < kh_end(temp_index_table)) { // if left is unique
            all_fuzzy_seeds[relaxed_fuzzy_seeds_len++] = temp_all_fuzzy_seeds[i];
            continue;
        }

        k = map32_get(temp_index_table, __blend_get_kmer(temp_all_fuzzy_seeds[i+1]));

        if (k < kh_end(temp_index_table)) { // if right is unique
            all_fuzzy_seeds[relaxed_fuzzy_seeds_len++] = temp_all_fuzzy_seeds[i];
        }
    }

    // no need for temp_all_fuzzy_seeds any more
    free(temp_all_fuzzy_seeds);

    if (relaxed_fuzzy_seeds_len != unique_fuzzy_seeds_len) {
        qsort(all_fuzzy_seeds + unique_fuzzy_seeds_len, relaxed_fuzzy_seeds_len - unique_fuzzy_seeds_len, sizeof(uint128_t), cmp_fuzzy_seeds);
    }

    for (uint64_t i = unique_fuzzy_seeds_len; i < relaxed_fuzzy_seeds_len; ) {
        uint64_t j = i + 1;

        while (j < relaxed_fuzzy_seeds_len && __blend_get_kmer(all_fuzzy_seeds[j]) == __blend_get_kmer(all_fuzzy_seeds[i])) {
            j++;
        }

        khint_t k; int absent;
        k = map32_put(temp_index_table, __blend_get_kmer(all_fuzzy_seeds[i]), &absent);
        kh_val(temp_index_table, k) = i;

        i = j;
    }

    *fuzzy_seeds_len = relaxed_fuzzy_seeds_len;

    printf("[INFO] Unique %lu / %lu (%.2f), Relaxted %lu (%.2f)\n", unique_fuzzy_seeds_len, all_fuzzy_seeds_len, (double)unique_fuzzy_seeds_len / (double)all_fuzzy_seeds_len, relaxed_fuzzy_seeds_len, (double)relaxed_fuzzy_seeds_len / (double)all_fuzzy_seeds_len);

    return unique_fuzzy_seeds_len;
}

int banded_align_and_report(const char *ref, uint64_t ref_span, const char *read, uint64_t read_span, uint64_t ref_pos, uint64_t ref_id, uint64_t read_id, int ref_strand, int read_strand, uint64_t read_pos, int len) {
    
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

    static int dp[MAX_LEN + 1][MAX_LEN + 1];
    static int bt[MAX_LEN + 1][MAX_LEN + 1];

    // init
    for (int i = 0; i <= ref_len; i++)
        for (int j = 0; j <= read_len; j++)
            dp[i][j] = NEG_INF;

    // free leading gaps within band
    for (int j = 0; j <= BAND && j <= read_len; j++)
        dp[0][j] = 0;

    for (int i = 0; i <= BAND && i <= ref_len; i++)
        dp[i][0] = 0;

    // DP
    for (int i = 1; i <= ref_len; i++) {

        int j_start = (i - BAND > 1) ? i - BAND : 1;
        int j_end   = (i + BAND < read_len) ? i + BAND : read_len;

        for (int j = j_start; j <= j_end; j++) {

            int best = NEG_INF, op = -1;

            // diagonal
            int diag = dp[i-1][j-1] + (ref_buf[i-1] == read_buf[j-1] ? MATCH : MISMATCH);
            best = diag; 
            op = 0; // diagonal movement

            int del = dp[i-1][j] + GAP;
            if (del > best) {
                best = del;
                op = 1; // downward movement
            }

            // insertion
            int ins = dp[i][j-1] + GAP;
            if (ins > best) {
                best = ins;
                op = 2; // right movement
            }

            dp[i][j] = best;
            bt[i][j] = op;
        }
    }

    // traceback
    int best = NEG_INF;
    int bi = ref_len, bj = read_len;

    // last row
    for (int j = 0; j <= read_len; j++) {
        if (dp[ref_len][j] > best) {
            best = dp[ref_len][j];
            bi = ref_len;
            bj = j;
        }
    }

    // last column
    for (int i = 0; i <= ref_len; i++) {
        if (dp[i][read_len] > best) {
            best = dp[i][read_len];
            bi = i;
            bj = read_len;
        }
    }

    if (best < MIN_SCORE) return 0;

    int i = bi, j = bj;

    int mismatches = 0;
#ifdef __DEBUG__
    char aln_ref[2*MAX_LEN];
    char aln_read[2*MAX_LEN];
    char aln_mid[2*MAX_LEN];
    int aln_len = 0;
#endif
    while (i > 0 && j > 0) {
        int op = bt[i][j];

        if (op == 0) {
#ifdef __DEBUG__
            aln_ref[aln_len]  = ref_buf[i-1];
            aln_read[aln_len] = read_buf[j-1];
            aln_mid[aln_len]  = (ref_buf[i-1] == read_buf[j-1]) ? '|' : '*';
#endif
            if (ref_buf[i-1] != read_buf[j-1]) {
                mismatches++;
                if (i-1 != 0 && j-1 != 0 && i != ref_len && j != read_len) {
                    // printf("SNP\tREFID=%lu\tREADID=%lu\tPOS=%lu\tREF=%c\tALT=%c\n", ref_id, read_id, ref_pos + i - 1, ref_buf[i-1], read_buf[j-1]);
                }
            }
            i--; j--;
        }
        else if (op == 1) {
#ifdef __DEBUG__
            aln_ref[aln_len]  = ref_buf[i-1];
            aln_read[aln_len] = '-';
            aln_mid[aln_len]  = ' ';
            mismatches++;
#endif
            if (i-1 != 0 && j-1 != 0 && i != ref_len && j != read_len) {
                // printf("DEL\tREFID=%lu\tREADID=%lu\tPOS=%lu\tREF=%c\tALT=%c\n", ref_id, read_id, ref_pos + i - 1, ref_buf[i-1], '-');
            }
            i--;
        }
        else {
#ifdef __DEBUG__
            aln_ref[aln_len]  = '-';
            aln_read[aln_len] = read_buf[j-1];
            aln_mid[aln_len]  = ' ';
            mismatches++;
#endif
            if (i-1 != 0 && j-1 != 0 && i != ref_len && j != read_len) {
                // printf("INS\tREFID=%lu\tREADID=%lu\tPOS=%lu\tREF=%c\tALT=%c\n", ref_id, read_id, ref_pos + i - 1, '-', read_buf[j-1]);
            }
            j--;
        }

        // if (mismatches > 20) return 0;
#ifdef __DEBUG__
        aln_len++;
#endif
    }
#ifdef __DEBUG__
    if (mismatches) {
        // reverse alignment strings
        for (int x = 0; x < aln_len / 2; x++) {
            char t;

            t = aln_ref[x];
            aln_ref[x] = aln_ref[aln_len - 1 - x];
            aln_ref[aln_len - 1 - x] = t;

            t = aln_mid[x];
            aln_mid[x] = aln_mid[aln_len - 1 - x];
            aln_mid[aln_len - 1 - x] = t;

            t = aln_read[x];
            aln_read[x] = aln_read[aln_len - 1 - x];
            aln_read[aln_len - 1 - x] = t;
        }

        aln_ref[aln_len]  = '\0';
        aln_mid[aln_len]  = '\0';
        aln_read[aln_len] = '\0';

        // print alignment
        printf("\n");
        printf("REF : %s\n", aln_ref);
        printf("      %s\n", aln_mid);
        printf("READ: %s\n", aln_read);

        printf("--------------------------------------------------------------------------\n");

        // print seqs
        printf("REF:  %.*s (%d)\n", ref_len, ref, ref_len);
        printf("READ: %.*s (%d)\n\n", read_len, read, read_len);

        for (int i = 0; i < read_len; i++) {
            printf("'%c' ", read[i]);
        }
        printf(" read: %lu, Pos: %lu, Len: %d\n", read_id, read_pos, len);
    }
#endif
    return mismatches;
}