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
#include "../include/structures.h"
#include "../include/utils.h"
#include "../include/io.h"
#include "../include/alignment.h"

START_TEST(test_bad_alignadapter){

	char *adapter_seqs[] = {"AGCTGTGCTC", "GCTAGTGTCGTC"};
	char *seq = {"AGCTGCTAGCTCGATCGA"};
	int starts[] = {-1, 0};
	int ends[] = {10, -1};
  int i; 
	struct AlignmentResult ar;
	initialize_alignment_result(&ar, 1000);
 
	struct SeqList adapters;
	initialize_seq_list(&adapters, 3);
	
	for (i=0; i<2; i++){
		add_seq_to_seqlist(&adapters, i, adapter_seqs[i], NULL, NULL);
	}
	
	
	for (i=0; i<2; i++){
  	pid_t pid = fork();
		if (pid == -1) {
		    ck_abort_msg("fork() failed");
		} else if (pid == 0) {
		    align_short(seq, adapters, 1, 0, "start", starts[i], ends[i], &ar);
		    exit(EXIT_SUCCESS);
		} else {
		    int status;
		    waitpid(pid, &status, 0);
		    ck_assert(WIFEXITED(status));
		    ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
		}
  }	

}
END_TEST

START_TEST(test_bad_cigar){
	/*
	Testing some bad inputs to count_operations
	*/
	char *cigars[] = {"2M2", "A21M", "2M431I31DD", "7I8II325D"};
	int i=0;	
	
  for (i=0; i<4; i++){
		pid_t pid = fork();
		if (pid == -1) {
		    ck_abort_msg("fork() failed");
		} else if (pid == 0) {
		    count_operations(cigars[i]);
		    exit(EXIT_SUCCESS);
		} else {
		    int status;
		    waitpid(pid, &status, 0);
		    ck_assert(WIFEXITED(status));
		    ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
		}
  }
  
}
END_TEST


START_TEST(test_bad_merge_input){
	/*
	Testing some inputs to merge_seqs which should exit properly. 
	*/
	
	char *seqs[] = {"AGCTG", "AGCTG",};
	char *quals[] ={"AAAAA", "AAAAA"};
	char *quals2_wrong = "AAAA";
	char *seq2_wrong = "AGCT";
	struct Seq merged_seq;
	pid_t pid;
	
	initialize_seq(&merged_seq, 1);
	
	pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
      merge_seqs(seqs[0], quals[0], seqs[1], quals2_wrong, "Phred+33", 
      					 &merged_seq, 100);
      free_seq(&merged_seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
      merge_seqs(seqs[0], quals[0], seq2_wrong, quals[1], "Phred+33", 
      				   &merged_seq, 100);
      free_seq(&merged_seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
      merge_seqs(seqs[0], quals[0], seqs[1], quals[1], "Phred+34", 
      					 &merged_seq, 100);
      free_seq(&merged_seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
      merge_seqs(seqs[0], quals2_wrong, seqs[1], quals[1], "Phred+33", 
      					 &merged_seq, 100);
      free_seq(&merged_seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
} 
END_TEST


START_TEST(test_bad_structures){

  struct Seq *seq=NULL;
	pid_t pid;
  	
	pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
      initialize_seq(seq, -2); 
      free_seq(seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
      struct Seq seq;
			initialize_seq(&seq, 4);
      write_to_seq(&seq, "ABC", "ABCD", NULL);
      free_seq(&seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
  		struct SeqList seqlist;
  		initialize_seq_list(&seqlist, -1);
  		free_seq_list(&seqlist);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
  		struct AlignmentResult ar;
  		initialize_alignment_result(&ar,-1);
  		free_alignment_result(&ar);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
			struct SeqList seqlist;
			initialize_seq_list(&seqlist, 10);
			add_seq_to_seqlist(&seqlist, 10, "ABC", "ABC", NULL);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
			struct Seq seq;
			initialize_seq(&seq, 1);
			write_to_seq(&seq, NULL, NULL, NULL);
			free_seq(&seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
			struct Seq *seq = NULL;
			initialize_seq(seq, 10);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
			struct AlignmentResult *ar = NULL;
			initialize_alignment_result(ar, 0);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
			struct SeqList *seq_list = NULL;
			initialize_seq_list(seq_list, 10);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
			struct Seq *seq = NULL;
			reinitialize_seq(seq, 4);
			free_seq(seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
			struct Seq seq;
			seq.s = NULL;
			reinitialize_seq(&seq, 4);
			free_seq(&seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
			struct Seq *seq = NULL;
			write_to_seq(seq, "ABCD", "ABCD", NULL);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
  
  pid = fork();
  if (pid == -1) {
      ck_abort_msg("fork() failed");
  } else if (pid == 0) {
			struct Seq seq;
			initialize_seq(&seq, 1);
			reverse_complement(&seq, NULL, "ABC");
			free_seq(&seq);
      exit(EXIT_SUCCESS);
  } else {
      int status;
      waitpid(pid, &status, 0);
      ck_assert(WIFEXITED(status));
      ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
  }
  
}
END_TEST


START_TEST(test_bad_io_inputs){

		pid_t pid1 = fork();
		if (pid1 == -1) {
		    ck_abort_msg("fork() failed");
		} else if (pid1 == 0) {
				struct SeqList seq_list;
		    seq_list = load_sequences("fakefile.fastq", 0);
		    free_seq_list(&seq_list);
		    exit(EXIT_SUCCESS);
		} else {
		    int status;
		    waitpid(pid1, &status, 0);
		    ck_assert(WIFEXITED(status));
		    ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
		}
		
		pid_t pid2 = fork();
		if (pid2 == -1) {
		    ck_abort_msg("fork() failed");
		} else if (pid2 == 0) {
				struct SeqList seq_list;
		    seq_list = load_sequences("./tests_files/empty.fasta", 0);
		    free_seq_list(&seq_list);
		    exit(EXIT_SUCCESS);
		} else {
		    int status;
		    waitpid(pid2, &status, 0);
		    ck_assert(WIFEXITED(status));
		    ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
		}

}
END_TEST


START_TEST(test_bad_cigar_for_extension){
	/*
	Testing some bad CIGAR strings incompatible with sequences to be extended
	*/
	char *cigars[] = {"DDDD", "IIII", "DDDDIIIIMMMM", "DIMMMMM", "DIMMMMM", 
										"DIMMMMM", "DIMMMMM"};
	
	char *seqs1[] = {"AGCT", "AGCT", "AGCT", "AGCT", NULL, "AGCT", "AGCT"};
	char *seqs2[] = {"CGTA", "CGTA", "CGTA", "CGTA", "CGTA", "CGTA", "CGTA"};
	char *quals1[] = {"ABCD", "ABCD", "ABCD", "ABCD", "ABCD", NULL, "ABC"};
	char *quals2[] = {"EFGH", "EFGH", "EFGH", "EFGH", "EFGH", "EFGH", "EFG"};
	
	int i=0;	
  struct Seq extended_seqs[2];
	
	for (i=0; i<2; i++)initialize_seq(&extended_seqs[i], 10);
			

	for (i=0; i<7; i++){
		pid_t pid = fork();
		if (pid == -1) {
		    ck_abort_msg("fork() failed");
		} else if (pid == 0) {
		    extend_read_and_reference_with_gaps(seqs1[i], seqs2[i],
																					  quals1[i], quals2[i],
																					  cigars[i], extended_seqs);
		    exit(EXIT_SUCCESS);
		} else {
		    int status;
		    waitpid(pid, &status, 0);
		    ck_assert(WIFEXITED(status));
		    ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
		}
  }
  
  
}
END_TEST


START_TEST(test_bad_get_insert_alignment_result_input){

  char *extended_reads[] = {
		"AGCTGCTAG",
		"AGCTGCTAG",
		"AGCTGCTAG",
		"AGCTGCTAG",
		"AGCTGCTAG ",
		"AGCTGCTAGAA"
	};
	
	char *extended_refs[] = {
		"AGCTGCTA",
		"AGCTGCTAG",
		"AGCTGCTAG",
		"AGCTGCTAG",
		"AGCTGCTAGA",
		"AGCTGCTAG A"
	};
	
	char *extended_quals[] = {
		"ABCDEFGHI",
		"ABCDEFGH",
		"ABCDEFGHI",
		"ABCDEFGHI",
		"ABCDEFGHIJ",
		"ABCDEFGHI J"
	};
	
	char *cigar_arrays[] = {
		"MMMMMMMMM",
		"MMMMMMMMM",
		"MMMMMMMM",
		"MMMMMMMMX",
		"MMMMMMMMMI",
		"MMMMMMMMMDM",
	};
	
	struct InsertAlignmentResult iar;
	initialize_ins_alig_res(&iar, 1, 1, 1);
	int i=0;
	
	for (i=0; i<6; i++){
		pid_t pid = fork();
		if (pid == -1) {
		    ck_abort_msg("fork() failed");
		} else if (pid == 0) {
				get_insert_alignment_result(extended_reads[i], extended_refs[i], 
														 extended_quals[i], cigar_arrays[i],
														"Phred+33", 10, &iar, 100);
				
		    exit(EXIT_SUCCESS);
		} else {
		    int status;
		    waitpid(pid, &status, 0);
		    ck_assert(WIFEXITED(status));
		    ck_assert_int_eq(WEXITSTATUS(status), EXIT_FAILURE);
		}
  }	
  
  free_ins_alig_res(&iar);

}
END_TEST


Suite *badinputs_suite(void)
{
    Suite *s;
    TCase *tc_core;
    
    s = suite_create("BAD INPUTS");
    tc_core = tcase_create("Core");
		tcase_add_test(tc_core, test_bad_merge_input);
		tcase_add_test(tc_core, test_bad_cigar);
		tcase_add_test(tc_core, test_bad_cigar_for_extension);
		tcase_add_test(tc_core, test_bad_alignadapter);
		tcase_add_test(tc_core, test_bad_structures);
		tcase_add_test(tc_core, test_bad_io_inputs);
		tcase_add_test(tc_core, test_bad_get_insert_alignment_result_input);
    suite_add_tcase(s, tc_core);
    
    return s;
}


int main(void)
{

  int number_failed=0;
  Suite *s;
  SRunner *sr;

  s = badinputs_suite();
  sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
