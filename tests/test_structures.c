#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include <math.h>
#include <float.h>
#include <stddef.h>
#include <check_stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include "../include/utils.h"
#include "../include/structures.h"



START_TEST(test_iar_initialization){
	
	struct InsertAlignmentResult iar;
	int lengths1[3] = {3, 4, 5};
	int lengths2[3] = {5, 6, 8};
	int i;
	
	//initialize
	initialize_ins_alig_res(&iar, 3, 4, 5);
	strcpy(iar.seq_mods, "AAA");
	strcpy(iar.qual_mods, "BBBB");
	strcpy(iar.condensed_cigar, "CCCCC");
	iar.aligned_err_prob = 0.4;
	
	ck_assert_str_eq(iar.seq_mods, "AAA");
	ck_assert_str_eq(iar.qual_mods, "BBBB");
	ck_assert_str_eq(iar.condensed_cigar, "CCCCC");
	ck_assert_int_eq(iar.aligned_err_prob, 0.4);
	
	for (i=0; i<3; i++){
		ck_assert_int_eq(iar.l[i], lengths1[i]);
		ck_assert_int_eq(iar.s[i], lengths1[i]+1);
	}
	
	//reinitialize to longer lengths
	reinitialize_ins_alig_res(&iar, 5, 6, 8);
	strcpy(iar.seq_mods, "AAAAA");
	strcpy(iar.qual_mods, "BBBBBB");
	strcpy(iar.condensed_cigar, "CCCCCCCC");
	iar.aligned_err_prob = 0.1;
	
	ck_assert_str_eq(iar.seq_mods, "AAAAA");
	ck_assert_str_eq(iar.qual_mods, "BBBBBB");
	ck_assert_str_eq(iar.condensed_cigar, "CCCCCCCC");
	ck_assert_int_eq(iar.aligned_err_prob, 0.1);
	
	for (i=0; i<3; i++){
		ck_assert_int_eq(iar.l[i], lengths2[i]);
		ck_assert_int_eq(iar.s[i], lengths2[i]+1);
	}
	
	//reinitiliaze to shorter lengths
	reinitialize_ins_alig_res(&iar, 3, 4, 5);
	strcpy(iar.seq_mods, "AAA");
	strcpy(iar.qual_mods, "BBBB");
	strcpy(iar.condensed_cigar, "CCCCC");
	iar.aligned_err_prob = 0.2;
	
	ck_assert_str_eq(iar.seq_mods, "AAA");
	ck_assert_str_eq(iar.qual_mods, "BBBB");
	ck_assert_str_eq(iar.condensed_cigar, "CCCCC");
	ck_assert_int_eq(iar.aligned_err_prob, 0.2);
	
	for (i=0; i<3; i++){
		ck_assert_int_eq(iar.l[i], lengths1[i]);
		ck_assert_int_eq(iar.s[i], lengths2[i]+1);
	}
	
	
	//reinitialize to longer lengths
	reinitialize_ins_alig_res(&iar, 5, 6, 8 );
	strcpy(iar.seq_mods, "AAAAA");
	strcpy(iar.qual_mods, "BBBBBB");
	strcpy(iar.condensed_cigar, "CCCCCCCC");
	iar.aligned_err_prob = 0.1;
	
	ck_assert_str_eq(iar.seq_mods, "AAAAA");
	ck_assert_str_eq(iar.qual_mods, "BBBBBB");
	ck_assert_str_eq(iar.condensed_cigar, "CCCCCCCC");
	ck_assert_int_eq(iar.aligned_err_prob, 0.1);
	
	for (i=0; i<3; i++){
		ck_assert_int_eq(iar.l[i], lengths2[i]);
		ck_assert_int_eq(iar.s[i], lengths2[i]+1);
	}
	
	//reinitialize to same lengths
	reinitialize_ins_alig_res(&iar, 5, 6, 8 );
	strcpy(iar.seq_mods, "AAAAA");
	strcpy(iar.qual_mods, "BBBBBB");
	strcpy(iar.condensed_cigar, "CCCCCCCC");
	iar.aligned_err_prob = 0.1;
	
	ck_assert_str_eq(iar.seq_mods, "AAAAA");
	ck_assert_str_eq(iar.qual_mods, "BBBBBB");
	ck_assert_str_eq(iar.condensed_cigar, "CCCCCCCC");
	ck_assert_int_eq(iar.aligned_err_prob, 0.1);
	
	for (i=0; i<3; i++){
		ck_assert_int_eq(iar.l[i], lengths2[i]);
		ck_assert_int_eq(iar.s[i], lengths2[i]+1);
	}
	
	free_ins_alig_res(&iar);
	
}END_TEST


START_TEST(test_seq)
{
  /*
  various ways of initializing seqs and then freeing 
  them.
  */

  struct Seq seq;
  char *s = "ABCDEF";
  char *q = "AAAAAA";
	char *snew1 = "GHI";
	char *qnew1 = "BBB";
	char *snew2 = "GHIJKLMON";
	char *qnew2 = "BBBBBBBBB";

  seq.s = (char*)calloc(strlen(s)+1, sizeof(char));
  seq.q = (char*)calloc(strlen(q)+1, sizeof(char));
	seq.id = NULL; 
	
  strncpy(seq.s, s, strlen(s));
  strncpy(seq.q, q, strlen(q));
  seq.l = strlen(s);

  ck_assert_str_eq(seq.s, s);
  ck_assert_str_eq(seq.q, q);
  ck_assert_int_eq(seq.l, strlen(s));
  ck_assert_int_eq((int)strlen(seq.s), (int)strlen(s));
  ck_assert_int_eq((int)strlen(seq.q), (int)strlen(q));

  free_seq(&seq);

  ck_assert(!seq.s);
  ck_assert(!seq.q);

  initialize_seq(&seq, strlen(s));
  write_to_seq(&seq, s, q, NULL);
  ck_assert_str_eq(seq.s, s);
  ck_assert_str_eq(seq.q, q);
  ck_assert_int_eq(seq.l, strlen(s));
  ck_assert_int_eq((int)strlen(seq.s), (int)strlen(s));
  ck_assert_int_eq((int)strlen(seq.q), (int)strlen(q));
	
	reinitialize_seq(&seq, strlen(snew1));
	write_to_seq(&seq, snew1, qnew1, NULL);
	ck_assert_str_eq(seq.s, snew1);
	ck_assert_str_eq(seq.q, qnew1);
	ck_assert_int_eq((int)seq.l, (int)strlen(snew1));
	
	reinitialize_seq(&seq, strlen(snew2));
	write_to_seq(&seq, snew2, qnew2, NULL);
	ck_assert_str_eq(seq.s, snew2);
	ck_assert_str_eq(seq.q, qnew2);
	ck_assert_int_eq((int)seq.l, (int)strlen(snew2));
	
	reinitialize_seq(&seq, strlen(snew2));
	write_to_seq(&seq, snew2, NULL, NULL);
	ck_assert_str_eq(seq.s, snew2);
	ck_assert_int_eq((int)seq.l, (int)strlen(snew2));
	ck_assert(seq.q == NULL);

  free_seq(&seq);
  ck_assert(!seq.s);
  ck_assert(!seq.q);

  initialize_seq(&seq, strlen(s));

  strncpy(seq.s, s, strlen(s));
  strncpy(seq.q, q, strlen(q));
  ck_assert_str_eq(seq.s, s);
  ck_assert_str_eq(seq.q, q);
  ck_assert_int_eq((int)strlen(seq.s), (int)strlen(s));
  ck_assert_int_eq((int)strlen(seq.q), (int)strlen(q));

  free_seq(&seq);
  ck_assert(!seq.s);
  ck_assert(!seq.q);
	
	
	

}
END_TEST

START_TEST(test_seq_list)
{
  struct SeqList seq_list;
  char *s = "AAAAAAAAAAAAAAAAAAAA";
  char *q = "AAAAAAAAAAAAAAAAAAAA";
  int i;

  initialize_seq_list(&seq_list, 100);
    
  for (i=0; i<100; i++){
    add_seq_to_seqlist(&seq_list, i, s, q, NULL);
    ck_assert_str_eq(seq_list.seqs[i].s, s);
    ck_assert_str_eq(seq_list.seqs[i].q, q); 
  }

  free_seq_list(&seq_list);
  ck_assert(seq_list.seqs==NULL);
}
END_TEST


START_TEST(test_alignment_result)
{
  struct AlignmentResult ar;
  char *cigar = "10M4I";

  initialize_alignment_result(&ar, strlen(cigar));
	ar.index = 2;
	ar.start = 0;
	ar.end = 10;
	strcpy(ar.cigar, cigar);
	
  ck_assert_int_eq(ar.index, 2);
  ck_assert_int_eq(ar.start, 0);
  ck_assert_int_eq(ar.end, 10);
  ck_assert_str_eq(ar.cigar, cigar);

  free_alignment_result(&ar);

  ck_assert(ar.cigar == NULL);

}
END_TEST


START_TEST(test_misc)
{

  char *cigar = (char*)calloc(14, sizeof(char));
  strcpy(cigar,"IIIX===DDDMMI");
  invert_cigar(cigar);
  invert_cigar(cigar);
  ck_assert_str_eq(cigar,  "IIIX===DDDMMI");
	free(cigar);
	cigar=NULL;
	invert_cigar(cigar);
	ck_assert(cigar == NULL);

  struct Seq seq;
	initialize_seq(&seq, 9);
  write_to_seq(&seq, "AGCTGCTGC", "ABCABCBCD", NULL);
  free_seq(&seq);
  
  initialize_seq(&seq, 1);
  reverse_complement(&seq, "AGCTGCTGC", "ABCABCBCD");
  ck_assert_str_eq(seq.s, "GCAGCAGCT");
  ck_assert_str_eq(seq.q, "DCBCBACBA");
	
  reverse_complement(&seq, "GCAGCAGCT", "DCBCBACBA");
	ck_assert_str_eq(seq.s, "AGCTGCTGC");
  ck_assert_str_eq(seq.q, "ABCABCBCD");
  
	
  reverse_complement(&seq, "MRWSYKVHDBXNP", "ABCABCBCDDBCB");
  ck_assert_str_eq(seq.s, "PNXVHDBMRSWYK");
  ck_assert_str_eq(seq.q, "BCBDDCBCBACBA");
  free_seq(&seq);

		
  int i; 
  struct SeqList seqlist;
  struct SeqList* seqlist_ptr=NULL;
  struct SeqList rev_complement_seqlist;
  char *seqs[] = {"AGCTGCTGC", "MRWSYKVHDBXNP"};
  char *quals[] = {"ABCABCBCD", "ABCABCBCDDBCB"};
  char *seq_answers[] = {"GCAGCAGCT", "PNXVHDBMRSWYK"};
  char *qual_answers[] = {"DCBCBACBA", "BCBDDCBCBACBA"};
  initialize_seq_list(&seqlist, 2);
	for (i=0; i<2; i++){
		add_seq_to_seqlist(&seqlist, i, seqs[i], quals[i], NULL);
	}
	
  reverse_complement_list(&rev_complement_seqlist, seqlist);
	
	for (i=0; i<2; i++){
		ck_assert_str_eq(rev_complement_seqlist.seqs[i].s, seq_answers[i]);
    ck_assert_str_eq(rev_complement_seqlist.seqs[i].q, qual_answers[i]);
		ck_assert_int_eq((int)rev_complement_seqlist.seqs[i].l, (int)strlen(qual_answers[i]));
	}
	free_seq_list(&seqlist);
	free_seq_list(seqlist_ptr);
	free_seq_list(&rev_complement_seqlist);
 	
 	struct AlignmentResult *ar_ptr = NULL;
 	free_alignment_result(ar_ptr);
	initialize_seq(&seq, 1);
 	reverse_complement(&seq, "AGCGTC", NULL);
 	ck_assert(seq.q == NULL);
	ck_assert_str_eq(seq.s, "GACGCT");
	ck_assert_int_eq((int)seq.l, 6);
 	free_seq(&seq);
 	
  
}


Suite *structures_suite(void)
{
  Suite *s;
  TCase *tc_core;
    
  s = suite_create("STRUCTURES");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_seq);
  tcase_add_test(tc_core, test_iar_initialization);
  tcase_add_test(tc_core, test_seq_list); 
  tcase_add_test(tc_core, test_alignment_result);
  tcase_add_test(tc_core, test_misc);
  suite_add_tcase(s, tc_core);
  return s;
}


int main(void)
{

  int number_failed=0;
  Suite *s;
  SRunner *sr;

  s = structures_suite();
  sr = srunner_create(s);
  srunner_set_fork_status(sr, CK_NOFORK);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
