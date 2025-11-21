#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include "../ext/WFA2/wavefront/wfa.h"
#include "../include/structures.h"
#include "../include/utils.h"
#include "../include/io.h"
#include "../ext/edlib/edlib/include/edlib.h"

UTILITY_STATIC_HEADER(int argextremum(int *nums, int nums_len, char *type));
UTILITY_STATIC_HEADER(int choose_loc(int *starts, int starts_len, char *type));
UTILITY_STATIC_HEADER(int rescmp(int edist1, int edist2, int start1, int start2, 
                      char *best_loc));
UTILITY_STATIC_HEADER(int get_start(char *extended_ref, int ref_len));
UTILITY_STATIC_HEADER(int get_end(char *extended_ref, int ref_len));


/*
  Extended_read, extended_ref, extended_qual are read, reference and quality 
  strings extended with space characters for indels after alignment of read 
  and reference (qual string has a gap whenever the read has a gap).
  cigar_array is the array of operations that map the read to reference
  after alignment with WFA2 alignment. qtype is the type of quality score
  to be used for instance Phred+33 (see utils.h for more details).
  dif_threshold is the edit threshold after which alignments are only
  reported as OTHER. The insert alignment *iar is a struct defined in structures.c 
  which contains the following information: seq_mods, qual_mods, condensed_cigar, 
  aligned_err_prob, index, s (for size), l (for length). See structures.h
  more explanation on these. Given the first 6 inputs, this functions
  fills an instance of iar. iar should have been initialized apriori however
  if the allocated size for members of iar dont have sufficient size, 
  then they will be rellocated inside this function. 
  
  iar should be allocated before usage.
  iar should be freed after usage.
*/
void get_insert_alignment_result(char *extended_read, char *extended_ref, 
												         char *extended_quals, char *cigar_array,
                                 char *qtype, int dif_threshold,
                                 struct InsertAlignmentResult *iar,
                                 int min_len);


/*
  Given read, qual and list of references, it goes through the references one by
  one, aligns them to the read. If the edit_distance is smaller then the best_distance 
  obtained so far, it will record the information in an insert alignment result
  by first constructing extended_read, extended_ref, extended_qual and then calling
  get_insert_alignment_result above. If the edit_distance is smaller then stop_when
  then the process will terminate and will not look at the rest of the references.
  wf_aligner is the WFA2 aligner instance which should be initialized before 
  (see ign.c for an example). extended_seqs and cigar_buffer are preallocated memmory
  buffers for objects used in this function. extended_seqs should be an struct Seq
  array of length 2. These Seq's members should be preallocated and they maybe
  reallocated in this function based on need. cigar_buffer is a pointer to a string
  of preallocated size which holds the cigar_string produced from alignment. This
  may also be reallocated during the process. buffer_size is the size allocated to
  cigar_string which chances if it is reallocated in this function and is used
  to keep track of it. 

  extended_seqs, iar, cigar_buffer should be allocated before usage.
  if the supplied iar already has its distance initialized to something other than
  -1 then the starting best_distance (see best_distance above) otherwise the initial 
  best_distance will be set to length of the read.

  extended_seqs, iar, cigar_buffer should be freed after usage.
  if iar is going to be reused then minimally its iar.distance should be set
  to -1 (align_long reallocates other members depending on need).
*/
void align_long(char *read, char *quals, struct SeqList references, 
 								wavefront_aligner_t *wf_aligner, float stop_when,
 								char *qtype, int dif_threshold, struct Seq *extended_seqs,
 								struct InsertAlignmentResult *iar, char **cigar_buffer,
								int *buffer_size
 							);


/*
  Similiar to above however uses EDLIB instance of WFA2 for alignment of a read to 
  list of short references (adapters and MIDs). max_edit_distance is the upper threshold
  of the edit_distance for alignment to be considered succesful. stop_when
  functions as above, if the edit_distance is less than or equal to threshold,
  algorithm terminates and does not look at the rest. Results are recorded in 
  ar (see structures.h for details).

  ar should be initialized before usage.
  ar should be freed after usage.
*/
void align_short(char *read, struct SeqList adapters, int max_edit_distance,
								 int stop_when, char *best_loc, int start, int end, 
                 struct AlignmentResult *ar);


void align_pairs(kseq_t *seq_record1, kseq_t *seq_record2, 
                 struct Seq *rc_seq_buffer, char **cigar_buffer,
                 int *cigar_size, wavefront_aligner_t *wf_aligner, 
                 struct Seq *extended_seqs);
     

#endif
