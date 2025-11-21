/*
the code for pasing the csv file was obtained from the following source
and modified to fit purpose
source: https://stackoverflow.com/questions/12911299/read-csv-file-in-c
author: Gus Gator

*/

#include <string.h>
#include <errno.h>
#include <libgen.h>
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
#include "../include/io.h"
#include "../include/alignment.h" 

// adjust BUFFER_SIZE to suit longest line
#define BUFFER_SIZE 1024 * 1024
#define NUM_FIELDS 30
#define RET_OK 0
#define RET_FAIL 1
#define FALSE 0
#define TRUE 1
#define MAXERRS 50
// char* array will point to fields
char *pFields1[NUM_FIELDS];
char *pFields2[NUM_FIELDS];
// field offsets into pFields array:



long loadFile(FILE *pFile, long *errcount);
static int  loadValues(char *line, char **pFields);
static char delim;

char * removeSpacesFromStr(char *string)
{
    // non_space_count to keep the frequency of non space characters
    int non_space_count = 0;
 
    //Traverse a string and if it is non space character then, place it at index non_space_count
    for (int i = 0; string[i] != '\0'; i++)
    {
        if (string[i] != ' ')
        {
            string[non_space_count] = string[i];
            non_space_count++;//non_space_count incremented
        }    
    }
    
    //Finally placing final character at the string end
    string[non_space_count] = '\0';
    return string;
}
 
 
long compareOutputs(FILE *pFile1, FILE *pFile2, long *errcount, char *test_type){

    char sInputBuf1 [BUFFER_SIZE];
    char sInputBuf2 [BUFFER_SIZE];
    int *indices;
    int i,length;
    int lineno=0;

    if (strcmp(test_type, "ign")==0) {
      indices = malloc(sizeof(int)*5);
      indices[0]=1;
      indices[1]=8;
      indices[2]=15;
      indices[3]=16;
      indices[4]=18;
      length=5;
    }
    else if (strcmp(test_type, "mer")==0) {
      indices = malloc(sizeof(int)*2);
      indices[0]=9;
      indices[1]=16;
      length=2;
    }
    else
    {
      return RET_FAIL;
    }
    
    if(pFile1 == NULL || pFile2 == NULL)
    { 
      free(indices); 
      return RET_FAIL;
    }
    
    while (!feof(pFile1) && !feof(pFile2)) {
      // load line into static buffer
      if(fgets(sInputBuf1, BUFFER_SIZE-1, pFile1)==NULL || fgets(sInputBuf2, BUFFER_SIZE-1, pFile2)==NULL)
          break;

      // skip first line (headers)
      if(++lineno==1)
          continue;

      // jump over empty lines
      if(strlen(sInputBuf1)==0 && strlen(sInputBuf2)==0)
          continue;
      // set pFields array pointers to null-terminated string fields in sInputBuf
      if(loadValues(sInputBuf1, pFields1)==RET_FAIL || loadValues(sInputBuf2, pFields2)==RET_FAIL){
         (*errcount)++;
          if(*errcount > MAXERRS)
              break;
      } else {
          // On return pFields array pointers point to loaded fields ready for load into DB or whatever
          // Fields can be accessed via pFields, e.g.
          
          for (i=0; i<length; i++){
          
            if (i<4){
              removeSpacesFromStr(pFields1[indices[i]]);
              removeSpacesFromStr(pFields2[indices[i]]);
              }
          
            if (strcmp(pFields1[indices[i]], pFields2[indices[i]])!=0)
            {
              
              if (strcmp(pFields1[indices[i]], pFields2[indices[i]])!=0){
                *errcount += 1;
                break;
                
              }
              
            }
            
          }
                
        }
      
    }
    
    free(indices);
    return RET_OK;


}


static int  loadValues(char *line, char** pFields){
    if(line == NULL)
        return RET_FAIL;

    // chop of last char of input if it is a CR or LF (e.g.Windows file loading in Unix env.)
    // can be removed if sure fgets has removed both CR and LF from end of line
    if(*(line + strlen(line)-1) == '\r' || *(line + strlen(line)-1) == '\n')
        *(line + strlen(line)-1) = '\0';
    if(*(line + strlen(line)-1) == '\r' || *(line + strlen(line)-1 )== '\n')
        *(line + strlen(line)-1) = '\0';

    char *cptr = line;
    int fld = 0;
    int inquote = FALSE;
    char ch;

    pFields[fld]=cptr;

    while((ch=*cptr) != '\0' && fld < NUM_FIELDS){

        if(ch == '"') {
            if(! inquote)
                pFields[fld]=cptr+1;
            else {
                *cptr = '\0';               // zero out " and jump over it
            }
            inquote = ! inquote;
        } else if(ch == delim && ! inquote){
            *cptr = '\0';                   // end of field, null terminate it
            pFields[++fld]=cptr+1;
        }
        cptr++;
    }

    

    return RET_OK;
}

START_TEST(test_ign)
{
   FILE *fp1,  *fp2;
   long errcount = 0L;
   int status;

   char *path1, *path2;
   path1 = "./test_files/simulated_ign_output.csv";
   path2 = "./test_files/simulated_ign_output_test.csv";
   
   delim = ',';

   fp1 = fopen(path1 , "r");
   if(fp1 == NULL) {
      fprintf(stderr,"Error opening file: %d\n",errno);
      ck_assert(1==0);
   }
   
   fp2 = fopen(path2 , "r");
   if(fp2 == NULL) {
      fprintf(stderr,"Error opening file: %d\n",errno);
      ck_assert(1==0);
   }
   
   status = compareOutputs(fp1, fp2, &errcount, "ign");
 
   fclose(fp1);
   fclose(fp2);

   ck_assert_int_le(errcount, 30);
   ck_assert(status==0);
}


START_TEST(test_mer)
{
   FILE *fp1,  *fp2;
   long errcount = 0L;
   int status;

   char *path1, *path2;
   path1 = "./test_files/simulated_mer_output.csv";
   path2 = "./test_files/simulated_mer_output_test.csv";
   
   delim = ',';

   fp1 = fopen(path1 , "r");
   if(fp1 == NULL) {
      fprintf(stderr,"Error opening file: %d\n",errno);
      ck_assert(1==0);
   }
   
   fp2 = fopen(path2 , "r");
   if(fp2 == NULL) {
      fprintf(stderr,"Error opening file: %d\n",errno);
      ck_assert(1==0);
   }
   
   status = compareOutputs(fp1, fp2, &errcount, "mer");
 
   fclose(fp1);
   fclose(fp2);

   ck_assert_int_le(errcount, 30);
   ck_assert(status==0);
}


Suite *structures_suite(void)
{
  Suite *s;
  TCase *tc_core;
    
  s = suite_create("MERLIGN");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_ign);
  tcase_add_test(tc_core, test_mer);
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


