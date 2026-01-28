#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include "../sketch/hmin.h"
#include "../sketch/blend.h"

#define MAX_SEQ_LEN 10000
#define READ_LEN 10000
#define MAX_SKETCH 10000
#define HMIN_K 15
#define HMIN_W 4
#define HMIN_LEVEL 2
#define BLEND_K_SHORT 21
#define BLEND_W_SHORT 11
#define BLEND_BITS_SHORT 32
#define BLEND_NEIGHBOR_NUMBER_SHORT 5
#define BLEND_K_HIFI 19
#define BLEND_W_HIFI 50
#define BLEND_BITS_HIFI 38
#define BLEND_NEIGHBOR_NUMBER_HIFI 5

#define READ_COUNT 10000
#define READ_COUNT_BREAKPOINT 1000
#define PROGRESS 0


/* ============================
   Arguments representation
   ============================ */
typedef enum {
    GENERATE,
    ANALIZE,
} program_e;

typedef enum {
    BLEND,
    HMIN
} sketch_e; 

typedef enum {
    HIFI,
    SHORT
} read_e; 

typedef struct {
    uint32_t max_seq_len;
    uint32_t read_len;
    uint32_t max_sketch;
    int k;
    int w;
    int blend_bits;
    int n_neighbors;
    int h_level;   
    sketch_e sketch_type;
    read_e read_type;
    program_e program_type;
    const char *dir;
    int progress;
} params;

/* ============================
   Sketch representation
   ============================ */

typedef struct {
    uint128_t *anchors;
    uint64_t count;
} minimizer_sketch;

/* ============================
   Matched Anchor representation
   ============================ */

typedef struct {
    uint64_t ref_pos;
    uint64_t read_pos;
} match;

/* ============================
   Utility
   ============================ */
void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options]\n\n"
        "Options:\n"
        "  -L, --max-seq-len <int>   Max reference length [%d]\n"
        "  -r, --read-len <int>      Read length [%d]\n"
        "  -s, --max-sketch <int>    Max sketch size [%d]\n"
        "  -k, --kmer <int>          k-mer size [NUM]\n"
        "  -w, --window <int>        window size [NUM]\n"
        "  -b, --blend-bits <int>    number of bits for hash [NUM]\n"
        "  -n, --n-neighbors <int>   number of neighbour [NUM]\n"
        "  -h, --h-level <int>       hierarchy level [%d]\n"
        "  --blend                   process BLEND sketching [DEFAULT]\n"
        "  --h-min                   process Hierarchical minimizers sketching\n"
        "  --hifi                    PacBio-HiFi long reads\n"
        "  --short                   Illumina short reads [DEFAULT]\n"
        "  -a, --analize             Analize sketching method\n"
        "  -g, --simulate            Simulate reads and write for file\n"
        "  -d, --dir <directory>     Directory to put simulated reads [sim]\n"
        "  -p, --progress            Display progress\n"
        "      --help                Show this help\n",
        prog, MAX_SEQ_LEN, READ_LEN, MAX_SKETCH, HMIN_LEVEL
    );
}

int ensure_dir(const char *path) {
    struct stat st;

    if (stat(path, &st) == 0) {
        if (S_ISDIR(st.st_mode)) {
            return 0;
        } else {
            return -1;
        }
    }

    // Does not exist -> creating one
    if (mkdir(path, 0755) == 0) {
        return 0;
    }

    return -1;
}

void init_params(params *p) {
    p->max_seq_len = MAX_SEQ_LEN;
    p->read_len = READ_LEN;
    p->max_sketch = MAX_SKETCH;
    p->sketch_type = BLEND;
    p->read_type = SHORT;
    p->program_type = ANALIZE;
    p->h_level = HMIN_LEVEL;
    p->progress = PROGRESS;
}

void parse_args(int argc, char **argv, params *p) {
    static struct option long_opts[] = {
        {"max-seq-len", required_argument, 0, 'L'},
        {"read-len",    required_argument, 0, 'r'},
        {"max-sketch",  required_argument, 0, 's'},
        {"kmer",        required_argument, 0, 'k'},
        {"window",      required_argument, 0, 'w'},
        {"blend-bits",  required_argument, 0, 'b'},
        {"n-neighbors", required_argument, 0, 'n'},
        {"h-level",     required_argument, 0, 'h'},
        {"h-min",       no_argument,       0,  1 },
        {"blend",       no_argument,       0,  2 },
        {"hifi",        no_argument,       0,  3 },
        {"short",       no_argument,       0,  4 },
        {"analize",     no_argument,       0, 'a'},
        {"simulate",    no_argument,       0, 'g'},
        {"dir",         required_argument, 0, 'd'},
        {"progress",    no_argument,       0, 'p'},
        {"help",        no_argument,       0,  0 },
        {0, 0, 0, 0}
    };

    int is_k_set = 0;
    int is_w_set = 0;
    int is_b_set = 0;
    int is_n_set = 0;
    int is_d_set = 0;

    int opt, idx;
    while ((opt = getopt_long(argc, argv, "L:r:s:k:w:h:d:agv", long_opts, &idx)) != -1) {
        switch (opt) {
        case 'L':
            p->max_seq_len = (uint32_t)atoi(optarg);
            break;
        case 'r':
            p->read_len = (uint32_t)atoi(optarg);
            break;
        case 's':
            p->max_sketch = (uint32_t)atoi(optarg);
            break;
        case 'k':
            p->k = atoi(optarg);
            is_k_set = 1;
            break;
        case 'w':
            p->w = atoi(optarg);
            is_w_set = 1;
            break;
        case 'b':
            p->blend_bits = atoi(optarg);
            is_b_set = 1;
            break;
        case 'n':
            p->n_neighbors = atoi(optarg);
            is_n_set = 1;
            break;
        case 'h':
            p->h_level = atoi(optarg);
            break;
        case 1:
            p->sketch_type = HMIN;
            break;
        case 2:
            p->sketch_type = BLEND;
            break;
        case 3:
            p->read_type = HIFI;
            break;
        case 4:
            p->read_type = SHORT;
            break;
        case 'a':
            p->program_type = ANALIZE;
            break;
        case 'g':
            p->program_type = GENERATE;
            break;
        case 'd':
            p->dir = optarg;
            is_d_set = 1;
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

    if (p->program_type == GENERATE) {
        return;
    }

    if (!is_k_set) {
        if (p->sketch_type == BLEND) {
            if (p->read_type == HIFI) {
                p->k = BLEND_K_HIFI;
            } else if (p->read_type == SHORT) {
                p->k = BLEND_K_SHORT;
            }
        } else if (p->sketch_type == HMIN) {
            p->k = HMIN_K;
        }
    }
    if (!is_w_set) {
        if (p->sketch_type == BLEND) {
            if (p->read_type == HIFI) {
                p->w = BLEND_W_HIFI;
            } else if (p->read_type == SHORT) {
                p->w = BLEND_W_SHORT;
            }
        } else if (p->sketch_type == HMIN) {
            p->w = HMIN_W;
        }
    }
    if (!is_b_set) {
        if (p->sketch_type == BLEND) {
            if (p->read_type == HIFI) {
                p->blend_bits = BLEND_BITS_HIFI;
            } else if (p->read_type == SHORT) {
                p->blend_bits = BLEND_BITS_SHORT;
            }
        }
    }
    if (!is_n_set) {
        if (p->sketch_type == BLEND) {
            if (p->read_type == HIFI) {
                p->n_neighbors = BLEND_NEIGHBOR_NUMBER_HIFI;
            } else if (p->read_type == SHORT) {
                p->n_neighbors = BLEND_NEIGHBOR_NUMBER_SHORT;
            }
        }
    }

    /* sanity checks */
    if (p->k <= 0 || p->w <= 0) {
        fprintf(stderr, "Error: k and w must be > 0\n");
        exit(1);
    }
    
    if (p->read_len > p->max_seq_len) {
        fprintf(stderr, "Error: read_len must <= max_seq_len\n");
        exit(1);
    }

    if (!is_d_set) {
        p->dir = "sim";
    }
}

char rand_base() {
    char bases[] = {'A','C','G','T'};
    return bases[rand() % 4];
}

void random_sequence(char *seq, int len) {
    for (int i = 0; i < len; i++) seq[i] = rand_base();
    seq[len] = '\0';
}

/* ============================
   Error simulation
   ============================ */
int sample_geometric(double p_extend) {
    int len = 1;
    while (((double)rand() / RAND_MAX) < p_extend)
        len++;
    return len;
}

int simulate_errors(const char *ref, char *read, int len, double snp_rate, int *edit_distance) {
    int count = 0;
    int r = 0;
    for (int i = 0; i < len; i++) {
        double p = (double)rand() / RAND_MAX;
        
        if (p < 0.33 * snp_rate) { /* Del */
            count++;
            continue;
        } else if (p < 0.67 * snp_rate) { /* Ins */
            char b = rand_base();
            read[r++] = b;
            count++;
            i--;
        } else if (p < snp_rate) { /* SNP */
            char b;
            do { b = rand_base(); } while (b == ref[i]);
            read[r++] = b;
            count++;
        } else {
            read[r++] = ref[i];
        }
    }
    read[r] = '\0';
    *edit_distance = count;

    return r;
}

/* ============================
   SKETCH FUNCTION
   ============================ */

void sketch_blend_sequence(const char *seq, int len, params *p, minimizer_sketch *sk) {

    uint128_t *minimizers;
    uint64_t minimizers_len = 0;

    minimizers_len = blend_sketch(seq, len, p->w, p->k, p->blend_bits, p->n_neighbors, 0, &minimizers);

    sk->anchors = minimizers;
    sk->count = minimizers_len;
}

void sketch_hmin_sequence(const char *seq, int len, params *p, minimizer_sketch *sk) {

    uint128_t *minimizers;
    uint64_t minimizers_len = 0;

    minimizers_len = hmin_sketch(seq, len, p->w, p->k, 0, 0, &minimizers);
    
    for (int i = 1; i < p->h_level; i++) {
        uint128_t *hi_minimizers;
        uint64_t hi_minimizers_len = 0;
        hi_minimizers_len = hmin_hi_sketch(minimizers, minimizers_len, p->w, p->k, 0, i, &hi_minimizers);

        if (minimizers_len) free(minimizers);
        minimizers = hi_minimizers;
        minimizers_len = hi_minimizers_len;
    }

    sk->anchors = minimizers;
    sk->count = minimizers_len;
}

/* ============================
   Evaluation
   ============================ */

int cmp_match(const void *a, const void *b) {
    const match *x = (match *)a;
    const match *y = (match *)b;
    if (x->ref_pos != y->ref_pos) return x->ref_pos - y->ref_pos;
    return x->read_pos - y->read_pos;
}

void evaluate(const minimizer_sketch *ref, const minimizer_sketch *read, int *TP, int *FP, int *FN, params *p) {

    /* ----------------------------
       1. Collect candidate matches
       ---------------------------- */

    if (ref->count == 0) {
        *TP = 0;
        *FP = 0;
        *FN = ref->count;
        return;
    }

    uint64_t capacity = 2 * ref->count;
    match *matches = (match *)malloc(sizeof(match) * 2 * ref->count);
    uint64_t m = 0;

    for (uint64_t i = 0; i < ref->count; i++) {
        for (uint64_t j = 0; j < read->count; j++) {
            if (p->sketch_type == HMIN) {
                if (GET_128_KMER(ref->anchors[i]) == GET_128_KMER(read->anchors[j])) {
                    matches[m].ref_pos = GET_128_INDEX(ref->anchors[i]);
                    matches[m++].read_pos = GET_128_INDEX(read->anchors[j]);

                    if (m == capacity) {
                        capacity *= 2;
                        matches = (match *)realloc(matches, capacity * sizeof(match));

                        if (!matches) {
                            fprintf(stderr, "realloc failed\n");
                            exit(-1);
                        }
                    }
                }
            } else if (p->sketch_type == BLEND) {
                if (__blend_get_kmer(ref->anchors[i]) == __blend_get_kmer(read->anchors[j])) {
                    matches[m].ref_pos = __blend_get_index(ref->anchors[i]);
                    matches[m++].read_pos = __blend_get_index(read->anchors[j]);

                    if (m == capacity) {
                        capacity *= 2;
                        matches = (match *)realloc(matches, capacity * sizeof(match));

                        if (!matches) {
                            fprintf(stderr, "realloc failed\n");
                            exit(-1);
                        }
                    }
                }
            }
        }
    }

    if (m == 0) {
        *TP = 0;
        *FP = 0;
        *FN = ref->count;
        free(matches);
        return;
    }

    /* ----------------------------
       2. Sort by reference position
       ---------------------------- */

    qsort(matches, m, sizeof(match), cmp_match);

    /* ----------------------------
       3. DP chaining
       ---------------------------- */

    int *dp = (int *)malloc(sizeof(int) * m);
    int best = 0;

    for (uint64_t i = 0; i < m; i++) {
        dp[i] = 1;
        for (uint64_t j = 0; j < i; j++) {
            if (matches[j].ref_pos < matches[i].ref_pos && matches[j].read_pos < matches[i].read_pos) {

                if (dp[j] + 1 > dp[i])
                    dp[i] = dp[j] + 1;
            }
        }
        if (dp[i] > best)
            best = dp[i];
    }

    /* ----------------------------
       4. Metrics
       ---------------------------- */

    *TP = best;
    *FP = ref->count - best;
    *FN = read->count - best;

    free(dp);
    free(matches);
}

void print_metrics(int read_len, int TP, int FP, int FN) {
    double precision = TP / (double)(TP + FP + 1e-9);
    double recall = TP / (double)(TP + FN + 1e-9);
    double f1 = 2 * precision * recall / (precision + recall + 1e-9);

    printf("Read Len=%d TP=%d FP=%d FN=%d | Precision=%.3f Recall=%.3f F1=%.3f\n", read_len, TP, FP, FN, precision, recall, f1);
}

/* ============================
   Main experiment
   ============================ */

int main(int argc, char **argv) {
    srand(time(NULL));

    params p;
    init_params(&p);
    parse_args(argc, argv, &p);

    double error_rates[] = {0.001, 0.005, 0.01, 0.02};
    int n_rates = sizeof(error_rates) / sizeof(double);

    if (p.program_type == ANALIZE) {
        
        const char *sketch_mode = p.sketch_type == HMIN ? "HMIN" : "BLEND";
        /* debug print */
        if (p.sketch_type == HMIN) {
            fprintf(stderr,
                "Params:\n"
                "  max_seq_len = %u\n"
                "  read_len    = %u\n"
                "  max_sketch  = %u\n"
                "  k           = %d\n"
                "  w           = %d\n"
                "  h_level     = %d\n"
                "  mode        = %s\n",
                p.max_seq_len, p.read_len, p.max_sketch,
                p.k, p.w, p.h_level, sketch_mode
            );
        } else if (p.sketch_type == BLEND) {
            fprintf(stderr,
                "Params:\n"
                "  max_seq_len = %u\n"
                "  read_len    = %u\n"
                "  max_sketch  = %u\n"
                "  k           = %d\n"
                "  w           = %d\n"
                "  blend_bits  = %d\n"
                "  n_neighbors = %d\n"
                "  mode        = %s\n",
                p.max_seq_len, p.read_len, p.max_sketch,
                p.k, p.w, p.blend_bits, p.n_neighbors, sketch_mode
            );
        }

        printf("\n");

        for (int e = 0; e < n_rates; e++) {

            printf("Error rate %.1f%c\n", error_rates[e] * 100, '%');

            char maf_name[256];
            snprintf(maf_name, sizeof(maf_name), "%s/sim.%d.maf", p.dir, (int)(error_rates[e] * 1000));

            FILE *maf = fopen(maf_name, "r");
            if (!maf) {
                fprintf(stderr, "fopen maf %s\n", maf_name);
                continue;
            }

            char *line = (char *)malloc(p.max_seq_len * 2 + 1);
            char *ref = (char *)malloc(p.max_seq_len + 1);
            char *read = (char *)malloc(p.max_seq_len * 2 + 1);

            int line_size = p.max_seq_len * 2 + 1;

            double TP = 0, FP = 0, FN = 0, precision = 0, recall = 0, f1 = 0, read_len = 0;    
            int read_count = 0;

            if (p.progress) {
                printf("\rProgress: %.2f%%", 0.0);
                fflush(stdout);
            }

            while (fgets(line, line_size, maf)) {
                
                /* header */
                if (line[0] != '@') continue;

                /* read ref line */
                if (!fgets(line, line_size, maf)) break;
                line[strcspn(line, "\n")] = '\0';
                strncpy(ref, line, p.max_seq_len);
                ref[p.max_seq_len] = '\0';
                /* '+' separator */
                if (!fgets(line, line_size, maf)) break;
                /* read sequence */
                if (!fgets(line, line_size, maf)) break;
                line[strcspn(line, "\n")] = '\0';
                strncpy(read, line, p.max_seq_len * 2);
                read[p.max_seq_len * 2] = '\0';

                int temp_read_len = strlen(read);

                minimizer_sketch ref_sk, read_sk;
                ref_sk.count = read_sk.count = 0;
                if (p.sketch_type == HMIN) {
                    sketch_hmin_sequence(ref, p.max_seq_len, &p, &ref_sk);
                    sketch_hmin_sequence(read, temp_read_len, &p, &read_sk);
                } else if (p.sketch_type == BLEND) {
                    sketch_blend_sequence(ref, p.max_seq_len, &p, &ref_sk);
                    sketch_blend_sequence(read, temp_read_len, &p, &read_sk);
                }

                int temp_TP, temp_FP, temp_FN;
                evaluate(&ref_sk, &read_sk, &temp_TP, &temp_FP, &temp_FN, &p);

                double temp_precision = temp_TP / (double)(temp_TP + temp_FP + 1e-9);
                double temp_recall = temp_TP / (double)(temp_TP + temp_FN + 1e-9);
                double temp_f1 = 2 * temp_precision * temp_recall / (temp_precision + temp_recall + 1e-9);
                
                TP += temp_TP;
                FP += temp_FP;
                FN += temp_FN;
                precision += temp_precision;
                recall += temp_recall;
                f1 += temp_f1;
                read_len += temp_read_len;
                read_count++;

                // print_metrics(read_len, TP, FP, FN); // legacy :')

                if (ref_sk.count) free(ref_sk.anchors);
                if (read_sk.count) free(read_sk.anchors);

                if (p.progress && read_count % READ_COUNT_BREAKPOINT == 0) {
                    double percent = (100.0 * read_count) / READ_COUNT;
                    printf("\rProgress: %.2f%%", percent);
                    fflush(stdout);
                }
            }

            TP = TP / read_count;
            FP = FP / read_count;
            FN = FN / read_count;
            precision = precision / read_count;
            recall = recall / read_count;
            f1 = f1 / read_count;
            read_len = read_len / read_count;

            printf("Read Cnt=%d Read Len=%.2f TP=%.2f FP=%.2f FN=%.2f | Precision=%.3f Recall=%.3f F1=%.3f\n\n", read_count, read_len, TP, FP, FN, precision, recall, f1);

            fclose(maf);

            free(line);
            free(ref);
            free(read);
        }
    } else if (p.program_type == GENERATE) {

        if (ensure_dir(p.dir) != 0) {
            fprintf(stderr, "Failed to create directory\n");
            exit(1);
        }

        for (int e = 0; e < n_rates; e++) {

            printf("Simulating reads for error rate %.1f%c\n", error_rates[e] * 100, '%');
            
            char maf_name[256];
            snprintf(maf_name, sizeof(maf_name), "%s/sim.%d.maf", p.dir, (int)(error_rates[e] * 1000));

            FILE *maf = fopen(maf_name, "w");
            if (!maf) {
                fprintf(stderr, "fopen maf %s\n", maf_name);
                exit(1);
            }

            double total_edit_dist = 0, total_read_len = 0;

            for (int i = 0; i < READ_COUNT; i++) {
                char *ref = (char *)malloc(p.max_seq_len + 1);
                random_sequence(ref, p.max_seq_len);

                int edit_distance = 0;
                char *read = (char *)malloc(p.max_seq_len * 2);
                total_read_len += simulate_errors(ref, read, p.read_len, error_rates[e], &edit_distance);

                fprintf(maf, "@read_%d\n%s\n+\n%s\n", i, ref, read);

                free(ref);
                free(read);

                total_edit_dist += edit_distance;
            }

            fclose(maf);

            printf("Avg edit distance: %0.2f, avg read len: %0.2f\n\n", total_edit_dist / READ_COUNT, total_read_len / READ_COUNT);
        }
    }

    return 0;
}
