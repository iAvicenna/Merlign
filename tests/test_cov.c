#include <string.h>
#include <stdio.h>
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
#include <unistd.h>
#include <dirent.h>
#include <time.h>
#include <regex.h>        


int check_log(FILE *fp, double threshold){
  
  char line[1000];
  regmatch_t match;
  regex_t regex;
	int reti;
	double percentage;
	char matched[100];
	char msgbuf[100];
	char pattern[] = "[0-9]+\\.[0-9]+";
	reti = regcomp(&regex, pattern, REG_EXTENDED);
	
	if (reti) {
		  fprintf(stderr, "Could not compile regex\n");
		  exit(EXIT_FAILURE);
	}

 
  while (fgets(line, 1000, fp)){
      
      if (strstr(line, "Lines executed:"))
      {
		    reti = regexec(&regex, line, 1, &match, 0);
				if (!reti) {
					strncpy(matched, line + match.rm_so, match.rm_eo - match.rm_so);
					percentage = atof(matched);
        	if (percentage>=threshold) return 1;
        	
				}
				else if (reti == REG_NOMATCH) {
						continue;
				}
				else {
						regerror(reti, &regex, msgbuf, sizeof(msgbuf));
						fprintf(stderr, "Regex match failed: %s\n", msgbuf);
						exit(EXIT_FAILURE);
				}

      }
      
    }
	
	regfree(&regex);
  return 0;

}


START_TEST(test_structures_cov)
{
  FILE *fp;
  fp = fopen("./logs/structures.cov.log", "r");

  if (fp==NULL){
    fprintf(stderr,"./log/structures.cov.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp, 94.0));
}
END_TEST


START_TEST(test_utils_cov)
{
  FILE *fp;
  fp = fopen("./logs/utils.cov.log", "r");

  if (fp==NULL){
    fprintf(stderr,"./log/utils.cov.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp, 99.0));
}
END_TEST


START_TEST(test_io_cov)
{
  FILE *fp;
  fp = fopen("./logs/io.cov.log", "r");

  if (fp==NULL){
    fprintf(stderr,"./log/io.cov.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp, 95.0));  
																	 
																	 
}																	 
END_TEST


START_TEST(test_alignment_cov)
{
  FILE *fp;
  fp = fopen("./logs/alignment.cov.log", "r");

  if (fp==NULL){
    fprintf(stderr,"./log/alignment.cov.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp, 98.0));
}
END_TEST


Suite *mem_suite(void)
{
    Suite *s;
    TCase *tc_core;
    
    s = suite_create("COVERAGE");
    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_structures_cov);
    tcase_add_test(tc_core, test_utils_cov);
    tcase_add_test(tc_core, test_io_cov);
		tcase_add_test(tc_core, test_alignment_cov);
    suite_add_tcase(s, tc_core);
    
    return s;
}


int main(void)
{

  int number_failed=0;
  Suite *s;
  SRunner *sr;

  s = mem_suite();
  sr = srunner_create(s);
  srunner_set_fork_status(sr, CK_NOFORK);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
