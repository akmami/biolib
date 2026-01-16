#include "opt_parser.h"


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