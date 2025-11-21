#include "../include/error.h"

void print_trace()
{
  void *array[10];
  char **strings;
  int size, i;

  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);
  if (strings != NULL)
  {

    for (i = 2; i < size-2; i++)
      fprintf (stderr, "%s\n", strings[i]);
  }

  free (strings);
}

void crash() //only because I wanted to call this crash
{
  exit(EXIT_FAILURE);
}

void assertm(int condition, char *file, const char *function, int line,  
            const char *msg, ...){
    
  if (condition != 1){
    va_list args;
    va_start(args, msg);
    if (file != NULL){
        fprintf(stderr, "%s:%s:%d: ", file, function, line);
    }
    vfprintf(stderr, msg, args);
    fprintf(stderr, ", exiting.\n");
    va_end(args);
    
#ifdef DEBUG
#ifdef VERBOSE
    if (file != NULL)print_trace();
#endif
#endif
    
    crash();
  }
}

void check_read(char *id, char* seq, char* qual, int seq_len, int qual_len){
	
	assertm(seq != NULL, __FILE__, __FUNCTION__, __LINE__,
  			 "sequence for read %s is NULL", id);
  assertm(qual != NULL, __FILE__, __FUNCTION__, __LINE__,
  			 "quality for read %s is NULL", id);
  assertm(seq_len==qual_len && seq_len > 0,
  			 __FILE__, __FUNCTION__, __LINE__,
    	   "read %s sequence has length %d but quality has length %d " 
     	   "which should be the same and positive", id, (int)strlen(seq), 
     	   (int)strlen(qual));

}

void check_file(const char* file_name, const char* arg_name){
    struct stat buffer;
    int exist = stat(file_name,&buffer);
    assertm(exist == 0, NULL, NULL, 0, "file %s, supplied as %s does not exist", 
    			file_name, arg_name); 
}


void check_parameters_mer(int match, int mismatch, int gap_open, 
                          int gap_extension, int end_free1, int end_free2, 
                          int begin_free1, int begin_free2,
                          int mask_threshold, int merge_debug)
                          
{
    
  assertm(match<=0, NULL, NULL, -1, "match is %d, but must be non-postive", 
         match);
  assertm(mismatch>0, NULL, NULL, -1, "mismatch is %d, but must be positive", 
         mismatch);
  assertm(gap_open>0, NULL, NULL, -1, "gap_open is %d, but must be positive", 
         gap_open);
  assertm(gap_extension>0, NULL, NULL, -1, 
         "gap_extension is %d, but must be positive", gap_extension);
  assertm(end_free1>=0, NULL, NULL, -1, 
         "end_free1 is %d, but must be non-negative", end_free1);
  assertm(end_free2>=0, NULL, NULL, -1, 
         "end_free2 is %d, but must be non-negative", end_free2);
  assertm(begin_free1>=0, NULL, NULL, -1, 
         "begin_free1 is %d, but must be non-negative", begin_free1);
  assertm(begin_free2>=0, NULL, NULL, -1, 
         "begin_free2 is %d, but must be non-negative", begin_free2);
  assertm(mask_threshold>=0, NULL, NULL, -1, 
         "mask_threshold is %d, but must be non-negative", mask_threshold);
  assertm(merge_debug == 0 || merge_debug ==1, NULL, NULL, -1, 
         "merge_debug is %d, but must be 0 or 1", merge_debug);       
}

void check_parameters_ign(int match, int mismatch, int gap_open, 
												  int gap_extension, int read_end_free, 
												  int read_begin_free, int ref_end_free, 
													int ref_begin_free, int mid_start, int mid_end,
													int alignment_debug)
{
  assertm(match<=0, NULL, NULL, -1, "match is %d, but must be non-postive", 
  			 match);
  assertm(mismatch>0, NULL, NULL, -1, "mismatch is %d, but must be positive", 
         mismatch);
  assertm(gap_open>0, NULL, NULL, -1, "gap_open is %d, but must be positive", 
         gap_open);
  assertm(gap_extension>0, NULL, NULL, -1, 
         "gap_extension is %d, but must be positive", gap_extension);
  assertm(read_end_free>=-1, NULL, NULL, -1, 
         "read_end_free is %d, but must be positive or -1", read_end_free);
	assertm(read_begin_free>=-1, NULL, NULL, -1, 
	       "read_begin_free is %d, but must be positive or -1", read_begin_free);
	assertm(ref_end_free>=-1, NULL, NULL, -1, 
	       "ref_end_free is %d, but must be positive or -1", ref_end_free);
	assertm(ref_begin_free>=-1, NULL, NULL, -1, 
	       "ref_begin_free is %d, but must be positive or -1", ref_begin_free);
  assertm(mid_start>=-1, NULL, NULL, -1, 
         "adapter_start is %d, but must be positive or -1", mid_start);
  assertm(mid_end>=-1, NULL, NULL, -1, 
         "adapter_end is %d, but must be positive or -1", mid_end);
  assertm(alignment_debug == 0 || alignment_debug ==1, NULL, NULL, -1, 
         "alignment_debug is %d, but must be 0 or 1", alignment_debug);
}
