#include <string.h>
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

int check_log(FILE *fp){
  
  char line[1000];

  int found1 = 0;
  int found2 = 0;
  int found3 = 0;
  int counter = 0;
  int max_lines = 20;
  
 
  while (fgets(line, 1000, fp)){
      
      if (strstr(line, "in use at exit: 0 bytes in 0 blocks"))
      {
        found1 = 1;
      }
      
      if (strstr(line, "All heap blocks were freed -- no leaks are possible"))
      {
        found2 = 1;
      }
      
      if (strstr(line, "ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)"))
      {
      	found3 = 1;
      }
      
      counter += 1;
      }

  return found1 && found2 && found3 && counter<max_lines;

}


START_TEST(test_structures_mem)
{
  FILE *fp;
  fp = fopen("./logs/test_structures.mem.log", "r");

  if (fp==NULL){
    fprintf(stderr,"./log/test_structures.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp));
	fclose(fp);
}
END_TEST


START_TEST(test_utils_mem)
{
  FILE *fp;
  fp = fopen("./logs/test_utils.mem.log", "r");

  if (fp==NULL){
    fprintf(stderr,"./log/test_utils.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp));
	fclose(fp);
}
END_TEST


START_TEST(test_io_mem)
{
  FILE *fp;
  fp = fopen("./logs/test_io.mem.log", "r");

  if (fp==NULL){
    fprintf(stderr,"./log/test_io.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp));
	fclose(fp);
}
END_TEST


START_TEST(test_alignment_mem)
{
  FILE *fp;
  fp = fopen("./logs/test_alignment.mem.log", "r");

  if (fp==NULL){
    fprintf(stderr,"./log/test_alignment.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp));
	fclose(fp);
}
END_TEST


START_TEST(test_mer_mem)
{

	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	FILE *fp4;
	FILE *fp5; 

	fp1 = fopen("./logs/mer1.mem.log", "r");
  if (fp1==NULL){
    fprintf(stderr,"./log/mer1.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp1));

	fp2 = fopen("./logs/mer2.mem.log", "r");
  if (fp2==NULL){
    fprintf(stderr,"./log/mer2.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp2));
  
	fp3 = fopen("./logs/mer3.mem.log", "r");
  if (fp3==NULL){
    fprintf(stderr,"./log/mer3.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp3));
  
	fp4 = fopen("./logs/mer4.mem.log", "r");
  if (fp4==NULL){
    fprintf(stderr,"./log/mer4.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp4));

	fp5 = fopen("./logs/mer5.mem.log", "r");
  if (fp5==NULL){
    fprintf(stderr,"./log/mer5.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp5));
	
    
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
  fclose(fp5);

}
END_TEST


START_TEST(test_ign_mem)
{

	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	FILE *fp4;

	fp1 = fopen("./logs/ign1.mem.log", "r");
  if (fp1==NULL){
    fprintf(stderr,"./log/ign1.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp1));
	
	
	fp2 = fopen("./logs/ign2.mem.log", "r");
  if (fp2==NULL){
    fprintf(stderr,"./log/ign2.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp2));

  fp3 = fopen("./logs/ign3.mem.log", "r");
  if (fp3==NULL){
    fprintf(stderr,"./log/ign3.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp3));

	
	fp4 = fopen("./logs/merlign.mem.log", "r");
  if (fp4==NULL){
    fprintf(stderr,"./log/ign1.mem.log not found");
    exit(EXIT_FAILURE);
  }
  ck_assert(check_log(fp4));
	
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
  fclose(fp4);
	
}
END_TEST


Suite *mem_suite(void)
{
    Suite *s;
    TCase *tc_core;
    
    s = suite_create("MEMORY");
    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_structures_mem);
    tcase_add_test(tc_core, test_utils_mem);
    tcase_add_test(tc_core, test_io_mem);
		tcase_add_test(tc_core, test_alignment_mem);
    tcase_add_test(tc_core, test_mer_mem);
    tcase_add_test(tc_core, test_ign_mem);
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
