#ifndef UTILS_H
#define UTILS_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "structures.h"
#include "error.h"

#ifdef UNIT_TESTING
  #define UTILITY_STATIC(DECLARATION) extern DECLARATION; DECLARATION
  #define UTILITY_STATIC_HEADER(DECLARATION) extern DECLARATION; DECLARATION
#else
  #define UTILITY_STATIC(DECLARATION) static inline DECLARATION
  #define UTILITY_STATIC_HEADER(DECLARATION)
#endif

UTILITY_STATIC_HEADER(double *decode_to_error_probabilities(char *quals, 
                                                            char* qtype));

UTILITY_STATIC_HEADER(int count_operations(char *cigar));

UTILITY_STATIC_HEADER(char *get_cigar_array(char *cigar, 
                                            int num_operations));

UTILITY_STATIC_HEADER(char encode_from_error_probability(double error_prob, 
                                                         char *qtype));

UTILITY_STATIC_HEADER(int decode_char(char qual, char *qtype));


void *safe_malloc(size_t n, char* name, char* file, const char *function, 
    int line);

int get_offset(char* qtype);

double decode_to_error_probability(char qual, int offset);

void extend_read_and_reference_with_gaps(char *seq1, char *seq2, char *quals1,
																				 char *quals2, char *cigar, 
                                         struct Seq *extended_seqs);

void soft_trim_cigar(char *cigar, int read_begin_free,int read_end_free,
                     int ref_begin_free, int ref_end_free);
                                                           
int count_matches(char *cigar);

void get_cigar_stats(char *cigar, int soft_clip, int convert_cigar, 
                     int begin_free1, int begin_free2, int end_free1,
                     int end_free2, int cigar_stats[6]);

void merge_seqs(char *seq1, char *quals1, char *seq2, char *quals2, 
    char *qtype, struct Seq *merged_seq, int mask_threshold);

#endif
