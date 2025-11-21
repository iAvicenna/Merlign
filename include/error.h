#ifndef ERROR_H
#define ERROR_H
#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>

#ifndef VERBOSE
#define fprintf(...) ((void)0)
#define vfprintf(...) ((void)0)
#endif

void assertm(int condition, char *file, const char *function, int line, 
    const char *msg, ...);

void check_read(char* id, char* seq, char* qual, int seq_len, int qual_len);

void check_file(const char* file_name, const char* arg_name);

void check_parameters_mer(int match, int mismatch, int gap_open, 
    int gap_extension, int end_free1, int end_free2, 
    int begin_free1, int begin_free2, int mask_threshold,
    int merge_debug);

void check_parameters_ign(int match, int mismatch, int gap_open,
    int gap_extension, int read_end_free, int read_begin_free, 
    int ref_end_free, int ref_begin_free, int mid_start, int mid_end,
    int align_debug);
#endif
