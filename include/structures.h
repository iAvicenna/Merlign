#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <stdbool.h> 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>


/*
  One of the main structs in merlign, iar contains information about
  alignment of read to a long reference (such a reference viral genome).

  seq_mods is the parts of the extended read which does not match
  with the extended ref.

  qual_mods is the read qualities when read and reference does not
  match. In compliance with RFC-4180, Section 2 [1], quality strings are enclosed
  in double quotation marks and if there are actual double quotations then
  its escaped via adding another double quotation. So a quality string 
  AAC"A would be written as "AAC""A. This is because iars are later
  written into csv files. An exception to this rule is when the dif_threshold
  is passed in which case it will be NA.

  1-https://www.rfc-editor.org/rfc/rfc4180#section-2

  condensed_cigar is the condensed cigar string (for the part where the
  reference aligns to the read, so gaps at the start and end of read are
  not shown)

  aligned_err_prob is the sum of base error probabilities of the part of the
  read that aligns to the reference.

  Since alignment to multiple references is allowed, the references are tracked
  by their unique index number and the index note here is the index of the
  reference for which the alignment information is given in an instance of 
  this struct.

  s is respectively the sizes of seq_mods, qual_mods, condensed_cigar. 
  l is respectively the lengths of seq_mods, qual_mods, condensed_cigar.
  l <= s-1 holds true.  
 
 */
typedef struct InsertAlignmentResult{
	char *seq_mods;
	char *qual_mods;
	char *condensed_cigar;
	int distance;
  int s[3];
	int l[3];
	int index;
  int is_rc;
	double aligned_err_prob;
}InsertAlignmentResult;

void write_to_ins_alig_res(struct InsertAlignmentResult *iar, char *seq_mods,
																	 char *qual_mods, char *condensed_cigar, 
																	 double aligned_err_prob, int index
																	);

void initialize_ins_alig_res(struct InsertAlignmentResult *iar, int length1, 
																		 int length2, int length3);

void reinitialize_ins_alig_res(struct InsertAlignmentResult *iar, int length1,
																			 int length2, int length3);

void free_ins_alig_res(struct InsertAlignmentResult *iar);


typedef struct Seq
{
  char *s;
  char *q;
	char *id;
  int l;
	int size;
}Seq;

void initialize_seq(struct Seq *seq, int length);

void reinitialize_seq(struct Seq *seq, int length);

void write_to_seq(struct Seq* seq, char* s, char* q, char* id);

void free_seq(struct Seq *seq);


typedef struct SeqList{
  struct Seq* seqs;
  int nseqs;
  int min_length;
  int max_length;
}SeqList;

void initialize_seq_list(struct SeqList *seq_list, int nseqs);

void add_seq_to_seqlist(struct SeqList* seq_list, int index, char* s, char* q, char* id);

void print_seq_list(struct SeqList seq_list);

void free_seq_list(struct SeqList *seq_list);


typedef struct AlignmentResult{
  int index;
  int distance;
  int start;
  int end;
  char *cigar;  
	int cigar_size;
}AlignmentResult;

void initialize_alignment_result(struct AlignmentResult *res, int cigar_len);

void free_alignment_result(struct AlignmentResult *res);

void invert_cigar(char *cigar);

void reverse_complement(struct Seq *rc_seq, char* sequence, char* qual);

void reverse_complement_list(struct SeqList *rev_seq_list, struct SeqList seq_list);

#endif
