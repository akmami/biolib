#ifndef __OPT_PARSER_H__
#define __OPT_PARSER_H__

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <sys/stat.h>
#include "struct_def.h"


int file_exists(char *filename);

void print_usage(const char *prog);

void init_params(params *p);

void parse_args(int argc, char **argv, params *p);

#endif