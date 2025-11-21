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

char *encode_from_error_probabilities(double *error_probs, int len, 
                                      char *qtype){
    
  int i;
  char *quals = (char*)safe_malloc((len+1)*sizeof(char), "quals", 
                                   __FILE__, __FUNCTION__, __LINE__);
  quals[len] = '\0';
  
  for (i=0; i<len; i++){
    quals[i] = encode_from_error_probability(error_probs[i], qtype);
  }
  
  quals[len] = '\0';
  
  return quals;
}


START_TEST(test_merge_seqs)
{
	/*
		This tests the case given in the page 3478, paper Fig1
		Robert C. Edgar, and Henrik Flyvbjerg

	  "Error filtering, pair assembly and error
	  correction for next-generation sequencing reads"

	  Bioinformatics, 31(21), 2015, 3476–3482
	  doi: 10.1093/bioinformatics/btv401
	*/
	
	int i;
  char *seq1 = "CATTGACA  ";
	char *seq2 = "  TAGACATT";
	char *expected_seq = "CATTGACATT";
	char *expected_seq_small = "CAttGACATT";
	double expected_num_quals[10] = {32,34,22,16,35,28,30,34,38,40};
	double fake_num_quals[10] = {32,34,22,15,35,28,30,34,38,40};
	char *expected_quals = calloc(11, sizeof(char));
	char *fake_quals = calloc(11, sizeof(char));
	const int num_quals1[10] = {32,34,20,20,28,16,14,10,-1,-1};
	const int num_quals2[10] = {-1,-1,2,5,4,8,12,20,38,40};
	char *char_quals1 = calloc(11, sizeof(char));
	char *char_quals2 = calloc(11, sizeof(char));
	struct Seq merged_seq1;
	
	for (i=0; i<10; i++){
		expected_quals[i] = (char)(expected_num_quals[i]+33);
		fake_quals[i] = (char)(fake_num_quals[i]+33);
	}
	
	for (i=0; i<10; i++){
		char_quals1[i] = (char)(num_quals1[i]+33);
		char_quals2[i] = (char)(num_quals2[i]+33);
	}

	initialize_seq(&merged_seq1, 1);
	merge_seqs(seq1, char_quals1, seq2, char_quals2, "Phred+33", &merged_seq1, 23);

	ck_assert_str_eq(merged_seq1.s, expected_seq_small);
	ck_assert_str_eq(merged_seq1.q, expected_quals);
	
	/*should give the same in reverse */
	merge_seqs(seq2, char_quals2, seq1, char_quals1, "Phred+33", &merged_seq1, 14);

	ck_assert_str_eq(merged_seq1.s, expected_seq);
	ck_assert_str_eq(merged_seq1.q, expected_quals);
	
	/*
	We also do some simpler tests
	*/
	
	char *seqs[] = {"AGCTAG", "AGCTAG", "GAAGCT", "TGCGAA"};
	char *quals[] = {"BBBBBB","AAAAAA","AAAAAA","AAAAAC"};
	char *seq_answers[] = {"AGCTAG", "AGCTAG", "AGCTAA"};
	
	for (i=0; i<3; i++){
		
		merge_seqs(seqs[0], quals[0], seqs[i+1], quals[i+1], "Phred+33", &merged_seq1, 1);
		ck_assert_str_eq(merged_seq1.s, seq_answers[i]);
	}	
	
	for (i=0; i<3; i++){
		
		merge_seqs(seqs[i+1], quals[i+1], seqs[0], quals[0],  "Phred+33", &merged_seq1, 1);
		ck_assert_str_eq(merged_seq1.s, seq_answers[i]);
	}	
	
	
	char *seq6 = "AGCTG";
	char *seq7 = "AGCTG";
	char *quals6 = "AAAAA";
	char *quals7 = "AAAAA";
	
	merge_seqs(seq6, quals6, seq7, quals7, "Phred+33", &merged_seq1, 10);
	ck_assert_str_eq(merged_seq1.s, seq6);
	
	free_seq(&merged_seq1);
  free(expected_quals);
  free(fake_quals);
  free(char_quals1);
  free(char_quals2);
}
END_TEST


START_TEST(test_decode_encode){
	int i;
	char *quals = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
	double *error_probs_33 = decode_to_error_probabilities(quals, "Phred+33");
	char *decoded_quals_33 = encode_from_error_probabilities(error_probs_33, strlen(quals), "Phred+33");
	double *error_probs_64 = decode_to_error_probabilities(quals, "Phred+64");
	char *decoded_quals_64 = encode_from_error_probabilities(error_probs_64, strlen(quals), "Phred+64");

	ck_assert_str_eq(quals, decoded_quals_33);
	ck_assert_str_eq(quals, decoded_quals_64);
	
	for (i=0; i < (int)strlen(quals); i++){

		ck_assert_double_ne(error_probs_33[i], error_probs_64[i]);
		ck_assert_int_le(error_probs_33[i], 1);
		if (i>=33){ck_assert_int_le(error_probs_64[i], 1);}
	}
	
	free(error_probs_33);
	free(decoded_quals_33);
	free(error_probs_64);
	free(decoded_quals_64);
	
}
END_TEST


START_TEST(test_cigar_count){
	
	char *cigars[] = {"2M", "21M", "2M431I31D", "7I8I323D", ""};
	int expected_counts[] = {2, 21, 464, 338, 0};
	int i,count;
	
	for (i=0; i<5; i++){
		count = count_operations(cigars[i]);
		ck_assert_int_eq(count, expected_counts[i]);
	}	
}
END_TEST


START_TEST(test_get_cigar_array){
	
	char *cigars[] = {"2M", "2M5D", "1I3D4M", "4X4X"};
	char *expected_cigar_arrays[] = {"MM", "MMDDDDD", "IDDDMMMM", "XXXXXXXX"};
	int i, num_operations;
	char *cigar_array;
	
	for (i=0; i<4; i++){
		num_operations = count_operations(cigars[i]);
		cigar_array = get_cigar_array(cigars[i], num_operations);
		ck_assert_str_eq(expected_cigar_arrays[i], cigar_array);
		free(cigar_array);
	}	
}
END_TEST


START_TEST(test_extend_read_and_reference_with_gaps){

	int i;
  char *cigar1 = "4D2M1D1I1M5I";
  char *query1[] = {"CGATAGCG", "AGTGCTGAG", "AAAAAAAA", "AAAAAAAAA"};
  char *seq_answer1[] = {"CGATAGC G     ", "    AG TGCTGAG"};
  char *qual_answer1[] = {"AAAAAAA A     ", "    AA AAAAAAA"};
	int num_operations = count_operations(cigar1);
	char *cigar_array1 = get_cigar_array(cigar1, num_operations);

	struct Seq extended_seqs1[2];
	for (i=0; i<2; i++)
	{
		initialize_seq(&extended_seqs1[i], strlen(cigar_array1));
	}
  extend_read_and_reference_with_gaps(query1[0], query1[1], query1[2], query1[3],
																		  cigar_array1, extended_seqs1);
  
  char *cigar2 = "4D4I";
  char *query2[] = {"AGCT", "CGTA", "ABCD", "DEFG"};
  char *seq_answer2[] = {"AGCT    ", "    CGTA"};
  char *qual_answer2[] = {"ABCD    ", "    DEFG"};
	num_operations = count_operations(cigar2);
	char *cigar_array2 = get_cigar_array(cigar2, num_operations);

  struct Seq extended_seqs2[2];
	for (i=0; i<2; i++)
	{
		initialize_seq(&extended_seqs2[i], strlen(cigar_array1));
	}
  extend_read_and_reference_with_gaps(query2[0], query2[1], query2[2], query2[3],
																		  cigar_array2, extended_seqs2);
  																													 
  for (i=0; i<2; i++){
  	ck_assert_str_eq(extended_seqs1[i].s, seq_answer1[i]);
  	ck_assert_str_eq(extended_seqs2[i].s, seq_answer2[i]);
    ck_assert_str_eq(extended_seqs1[i].q, qual_answer1[i]);
  	ck_assert_str_eq(extended_seqs2[i].q, qual_answer2[i]);
    free_seq(&extended_seqs1[i]);
    free_seq(&extended_seqs2[i]);
  }

  free(cigar_array1);
  free(cigar_array2);
  
}
END_TEST


START_TEST(test_get_cigar_stats)
{

		int i,j, convert;
		int cigar_stats1[6] = {0,0,0,0,0,0};
		int cigar_stats2[6] = {0,0,0,0,0,0};
		char *cigar_arrays[14] = {
		"DMMMMMI",
		"MMMMMMD",
		"DMMMMMM",
		"DDXMMXDIXMM",
		"IDX===ID",
		"DIMDD",
		"1D1I1M2D",
		"1I1D1X3=1I1D",
    "IDMMMDII",
    "DDMMMMIII",
    "DDMMMMIII",
    "DDDDIIIII",
    "DDDD",
    "DDDD"
		};

    int begin_frees1[14] = {1, 1, 1, 2, 2, 2, 2, 2, 2, 
                            2, 0, 4, 4, 4};
    int begin_frees2[14] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 
                            0, 0, 4, 4, 0};

    int end_frees1[14] = {1, 1, 1, 1, 2, 3, 2, 2, 2, 0, 
                          0, 4, 0, 4};
	  int end_frees2[14] = {0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 
                          0, 4, 4, 4}; 

		int answers1[14][6] = {{1,0,5,0,0,1},
								 				{0,1,6,0,0,0},
								 				{1,0,6,0,0,0},
								 				{2,0,4,3,1,1},
								 				{0,1,3,1,1,2},
								 				{1,2,1,0,0,1},
								 				{1,2,1,0,0,1},
								 				{0,1,3,1,1,2},
                        {1,2,3,0,2,0},
                        {2,2,4,0,0,1},    
                        {0,0,4,0,2,3},       
                        {4,4,0,0,0,1},
                        {4,0,0,0,0,0},
                        {4,0,0,0,0,0},
                        };
								 				
    int answers2[14][6] = {{0,0,5,0,1,1},
								 				{0,0,6,0,1,0},
								 				{0,0,6,0,1,0},
								 				{0,0,4,3,3,1},
								 				{0,0,3,1,2,2},
								 				{0,0,1,0,3,1},
								 				{0,0,1,0,3,1},
								 				{0,0,3,1,2,2},
                        {0,0,3,0,2,3},
                        {0,0,4,0,2,3},
                        {0,0,4,0,2,3},
                        {0,0,0,0,4,5},
                        {0,0,0,0,4,0},
                        {0,0,0,0,4,0},                      
                        };
	
    //nsoftstart;
    //nsoftend;
    //nmatches;
    //nmismatches;
    //ndeletions;
    //ninsertions;

		for (i=0; i<14; i++)
		{
		  
			if (i<=5 || i>7) convert=0;
			else convert=1;
		
			get_cigar_stats(cigar_arrays[i], 1, convert, begin_frees1[i], 
          begin_frees2[i], end_frees1[i], end_frees2[i], cigar_stats1);
			get_cigar_stats(cigar_arrays[i], 0, convert, begin_frees1[i], 
          begin_frees2[i], end_frees1[i], end_frees2[i], cigar_stats2);
			
			for (j=0; j<6; j++)
			{

				ck_assert_int_eq(cigar_stats1[j], answers1[i][j]);
				ck_assert_int_eq(cigar_stats2[j], answers2[i][j]);
				
			}
			
			
		}
		
		get_cigar_stats(NULL, 0, 0, 0, 0, 0, 0, cigar_stats1);
		
		for (j=0; j<6; j++){
			ck_assert(cigar_stats1[j] == -1);
		}
		
		
		
}


Suite *utils_suite(void)
{
    Suite *s;
    TCase *tc_core;
    
    s = suite_create("UTILS");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_merge_seqs);
    tcase_add_test(tc_core, test_decode_encode);
    tcase_add_test(tc_core, test_cigar_count);
    tcase_add_test(tc_core, test_get_cigar_array);
    tcase_add_test(tc_core, test_extend_read_and_reference_with_gaps);
    tcase_add_test(tc_core, test_get_cigar_stats);
    suite_add_tcase(s, tc_core);
    return s;
}


int main(void)
{

  int number_failed=0;
  Suite *s;
  SRunner *sr;

  s = utils_suite();
  sr = srunner_create(s);
  srunner_set_fork_status(sr, CK_NOFORK);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
