/*
source: https://stackoverflow.com/questions/12911299/read-csv-file-in-c
author: Gus Gator
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <libgen.h>
// adjust BUFFER_SIZE to suit longest line
#define BUFFER_SIZE 1024 * 1024
#define NUM_FIELDS 22
#define MAXERRS 5000
#define RET_OK 0
#define RET_FAIL 1
#define FALSE 0
#define TRUE 1

// char* array will point to fields
char *pFields1[NUM_FIELDS];
char *pFields2[NUM_FIELDS];
// field offsets into pFields array:
#define M1_index 1
#define M2_index 8
#define align_index 15
#define is_rc 16
#define align_mods 18


long loadFile(FILE *pFile, long *errcount);
static int  loadValues(char *line, char **pFields, long lineno);
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
 
 
long compareOutputs(FILE *pFile1, FILE *pFile2, long *errcount){

    char sInputBuf1 [BUFFER_SIZE];
    char sInputBuf2 [BUFFER_SIZE];
    long lineno = 0L;
    int i;
    int indices[5] = {1, 8, 15, 16, 18};
    int nerrors = 0;
    if(pFile1 == NULL || pFile2 == NULL)
        return RET_FAIL;

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
        if(loadValues(sInputBuf1, pFields1, lineno)==RET_FAIL || loadValues(sInputBuf2, pFields2, lineno)==RET_FAIL){
           (*errcount)++;
            if(*errcount > MAXERRS)
                break;
        } else {
            // On return pFields array pointers point to loaded fields ready for load into DB or whatever
            // Fields can be accessed via pFields, e.g.
            
            for (i=0; i<5; i++){
            
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
    
    return RET_OK;


}


static int  loadValues(char *line, char** pFields, long lineno){
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

int main(int argc, char **argv)
{
   FILE *fp1,  *fp2;
   long errcount = 0L;
   int status;

   if(argc!=4){
       printf("Usage: %s csvfilepath1 csvfilepath2 delimiter\n", basename(argv[0]));
       return (RET_FAIL);
   }
   if((delim=argv[3][0])=='\0'){
       fprintf(stderr,"delimiter must be specified\n");
       return (RET_FAIL);
   }
   
   fp1 = fopen(argv[1] , "r");
   if(fp1 == NULL) {
      fprintf(stderr,"Error opening file: %d\n",errno);
      return(RET_FAIL);
   }
   
   fp2 = fopen(argv[2] , "r");
   if(fp2 == NULL) {
      fprintf(stderr,"Error opening file: %d\n",errno);
      return(RET_FAIL);
   }
   
   status = compareOutputs(fp1, fp2, &errcount);
 

   fclose(fp1);
   fclose(fp2);
   printf("Encountered %ld error(s)\n", errcount);
   if(errcount>0)
     return(RET_FAIL);
   return status;
}
