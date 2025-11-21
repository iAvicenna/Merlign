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
#include "../include/io.h"


START_TEST(test_load_sequences){

	int i;
	struct SeqList seq_list1;
	struct SeqList seq_list2;
	
	seq_list1 = load_sequences("./test_files/example_file1_R1.fastq.gz", 0);
	seq_list2 = load_sequences("./test_files/example_file1_R1.fastq.gz", 0);
	
	ck_assert_int_eq(seq_list1.nseqs, seq_list2.nseqs);
	
	for (i=0; i<seq_list1.nseqs; i++){
		ck_assert_int_eq((int)strlen(seq_list1.seqs[i].s), (int)strlen(seq_list1.seqs[i].q));
	  ck_assert_int_eq((int)strlen(seq_list2.seqs[i].s), (int)strlen(seq_list2.seqs[i].q));
		ck_assert_int_eq((int)strlen(seq_list1.seqs[i].s), (int)strlen(seq_list2.seqs[i].q));
	}
	
	free_seq_list(&seq_list1);
	free_seq_list(&seq_list2);	

}
END_TEST


START_TEST(test_load_adapters){

	struct SeqList seq_list1;
	
	seq_list1 = load_sequences("./test_files/adapters.fasta", 0);
	
	ck_assert_int_eq(seq_list1.nseqs, 1);
	ck_assert_str_eq(seq_list1.seqs[0].s, "AGATCGGAAGAG");
	free_seq_list(&seq_list1);

	seq_list1 = load_sequences("./test_files/adapters.fasta", 1);
	
	ck_assert_int_eq(seq_list1.nseqs, 2); 
	ck_assert_str_eq(seq_list1.seqs[0].s, "AGATCGGAAGAG");
	ck_assert_str_eq(seq_list1.seqs[1].s, "CTCTTCCGATCT");
	
	free_seq_list(&seq_list1);

}
END_TEST


int check_output(FILE *fp, char *text){
  
  char line[1000];
  int counter=0; 
	
  while (fgets(line, 1000, fp)){
		
		counter += 1;
		line[strcspn(line, "\n\r")]=0;
		
		if (counter==2 && strcmp(text, line)==0) return 1;
	}
	
  return 0;

}


START_TEST(test_get_seq){

		FILE *fp2 = NULL;
		char *text = "Seq1,0,0,9,10,0,0,0,1,222,233,11,0,0,1,0, 16M1X66M1I11M1D39M1I81M,T C ,\"J J \",0.059267";

		system("cat ./test_files/ign_stdin.fasta | ../src/ign  --references ./test_files/BE92.fasta "
					 "--mids ./test_files/MIDs.fasta --mid_start 0 --mid_end 20 --maxN 100 --output ./test_files/get_seq_output1.csv 2>/dev/null"); 
		
		// the case where get_seq gets input from stdin is not possible to test by calling the function
		// so just doing a system call to test the case of fp=NULL for get_seq. 
		
		fp2 = fopen("./test_files/get_seq_output1.csv", "r");

		if (fp2==NULL){
			fprintf(stderr, "./test_files/get_seq_output1.csv");
			exit(EXIT_FAILURE);
		}
		
		ck_assert(check_output(fp2, text)==1);
		fclose(fp2);
		
		system("../src/ign --reads ./test_files/ign_stdin.fasta  --references ./test_files/BE92.fasta "
					 "--mids ./test_files/MIDs.fasta --mid_start 0 --mid_end 20 --maxN 100 --output ./test_files/get_seq_output2.csv 2>/dev/null"); 
		
		fp2 = fopen("./test_files/get_seq_output2.csv", "r");

		if (fp2==NULL){
			fprintf(stderr, "./test_files/get_seq_output2.csv");
			exit(EXIT_FAILURE);
		}
		ck_assert(check_output(fp2, text)==1);
		
		fclose(fp2);
	
}


Suite *io_suite(void)
{
    Suite *s;
    TCase *tc_core;
    
    s = suite_create("IO");
    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_load_sequences);
    tcase_add_test(tc_core, test_load_adapters);
		tcase_add_test(tc_core, test_get_seq);
    suite_add_tcase(s, tc_core);
    return s;
}


int main(void)
{

  int number_failed=0;
  Suite *s;
  SRunner *sr;

  s = io_suite();
  sr = srunner_create(s);
  srunner_set_fork_status(sr, CK_NOFORK);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
