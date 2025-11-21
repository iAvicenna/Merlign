#include "../include/alignment.h"
#include "../include/error.h"

UTILITY_STATIC(int argextremum(int *nums, int nums_len, char *type)){

  int i;
  int ext_index=0;

  for (i=0; i<nums_len; i++){
    if (strcmp(type, "min")==0 && nums[i]<nums[ext_index]){
      ext_index=i;
    }
    else if (strcmp(type, "max")==0 && nums[i]>nums[ext_index]){
      ext_index=i;
    }

  }

  return ext_index;

}

UTILITY_STATIC(int choose_loc(int *starts, int starts_len, char *type)){

  if (strcmp(type, "start")==0){
    return argextremum(starts, starts_len, "min");
  }
  else if (strcmp(type, "end")==0){
    return argextremum(starts, starts_len, "max");
  }
  else{
    return 0;
  }

}

UTILITY_STATIC(int rescmp(int edist1, int edist2, int start1, int start2, 
               char *best_loc)){
  
  if (edist2 < edist1) return 1;
  else if(edist2==edist1){
    if (strcmp(best_loc, "start")==0 && start2 < start1) return 1;
    else if (strcmp(best_loc, "end")==0 && start1 < start2) return 1;
  }

  return 0;
  
}

UTILITY_STATIC(int get_start(char *extended_ref, int ref_len))
{

  int i;
  int start = 0;
  
  for (i=0; i<ref_len; i++){
    if (extended_ref[i]!=' ')
    {
      return start;
    }
    else start += 1;
  }
  
  return start;

}


UTILITY_STATIC(int get_end(char *extended_ref, int ref_len))
{
  int i;
  int end = ref_len - 1;
  
  for (i=ref_len-1; i>-1; i--){
    if (extended_ref[i]!=' ')
    {
      return end;
    }
    else end -= 1;
  }
  
  return end;

}


static inline int min(int x, int y) {
  return (x < y) ? x : y;
}


void get_length(int start, int end, char *extended_quals, char *cigar_array, 
                int mcounters[3])
{
  int i, j, num_len;
  char num[100];
  int op_counter=1;
  
  mcounters[1] += 2;
  
  for (i=start; i<=end; i++){
  
    if (cigar_array[i] == 'X' || cigar_array[i] == 'I' || cigar_array[i] == 'D')
    {
      
      mcounters[0] += 1;

      if (extended_quals[i] != '"'){
        mcounters[1] += 1;
      }
      else{
        mcounters[1] += 2;
      }
      
    }
    
    if (i<end && cigar_array[i] == cigar_array[i+1])
    {
      op_counter += 1;
    }
    else{
      sprintf(num, "%d", op_counter);
      num_len = strlen(num);
      
      for (j=0; j<num_len; j++){
        mcounters[2] += 1;
        }
        
      mcounters[2] += 1;
      op_counter=1;  
    }
    
  }  
  
}


void get_insert_alignment_result(char *extended_read, char *extended_ref, 
			                           char *extended_quals, char *cigar_array,
       				                   char *qtype, int dif_threshold, 
       				                   struct InsertAlignmentResult *iar,
       				                   int min_len)
{

  int i,j;
  int mcounter0=0;
  int mcounter1=0;
  int mcounter2=0;
  char num[100];
  double err_prob=0;
  int offset = get_offset(qtype);
  int read_len = (int)strlen(extended_read);
  int ref_len = (int)strlen(extended_ref);
  int qual_len = (int)strlen(extended_quals);
  int cigar_len = (int)strlen(cigar_array);
  int num_len = 0;
  
  assertm(read_len == ref_len, __FILE__, __FUNCTION__, __LINE__,
         "Extended read has length %d but extended ref has length %d " 
         "which should be the same", read_len, ref_len);
  
  assertm(read_len == cigar_len, __FILE__, __FUNCTION__, __LINE__,
         "Extended read has length %d but cigar_array has length %d " 
         "which should be the same", read_len, cigar_len);
  
  assertm(read_len == qual_len, __FILE__, __FUNCTION__, __LINE__,
         "Extended read has length %d but extended qual has length %d " 
         "which should be the same", read_len, qual_len);

  int start = 0;//get_start(extended_ref, ref_len);
  int end = ref_len-1;//get_end(extended_ref, ref_len);
  int op_counter = 1;
  int lengths[3]={0,0,0};
  int nmatches = 0;
  
  get_length(start, end, extended_quals, cigar_array, lengths);

  reinitialize_ins_alig_res(iar, lengths[0], lengths[1], lengths[2]);
  
  iar->qual_mods[0] = '"';
  mcounter1 += 1;                           
  
  for (i=start; i<=end; i++){
    if (extended_quals[i] != ' '){
      err_prob += decode_to_error_probability(extended_quals[i], offset);
    }
  
    if (cigar_array[i] == 'X' || cigar_array[i] == 'I' || cigar_array[i] == 'D')
    {
      iar->seq_mods[mcounter0] = extended_read[i];
      mcounter0 += 1;
      
      assertm(extended_read[i] != extended_ref[i], __FILE__, __FUNCTION__, __LINE__, 
             "Extended read and ref is equal while cigar string is %c", 
              cigar_array[i]);
      
      assertm(((extended_read[i] != ' ' && extended_quals[i] != ' ') ||
             (extended_read[i] == ' ' && extended_quals[i] == ' ')), 
             __FILE__, __FUNCTION__, __LINE__,
             "%d^th position  extended_read is \"%c\" and extended_quals" 
             "is \"%c\", gaps do not match",
             i, extended_read[i], extended_quals[i]);
    
      if (extended_quals[i] != '"'){
        iar->qual_mods[mcounter1] = extended_quals[i];
        mcounter1 += 1;
      }
      else{
        iar->qual_mods[mcounter1] = '"';
        mcounter1 += 1;
        iar->qual_mods[mcounter1] = '"';
        mcounter1 += 1;
      }
      
    }
    else if(cigar_array[i] == 'M'){nmatches += 1;}
    
    if (i<end && cigar_array[i] == cigar_array[i+1])
    {
      op_counter += 1;
    }
    else{
      sprintf(num, "%d", op_counter);
      num_len = strlen(num);
      
      for (j=0; j<num_len; j++){
        iar->condensed_cigar[mcounter2] = num[j];
        mcounter2 += 1;
        }
        
      iar->condensed_cigar[mcounter2] = cigar_array[i];
      mcounter2 += 1;
      op_counter=1;  
    }
    
  }

  assertm(mcounter0==lengths[0], __FILE__, __FUNCTION__, __LINE__,
         "post and precalculated lengths for seq mods not equal"
         "(%d/%d)", mcounter0, lengths[0]);
  assertm(mcounter1+1==lengths[1], __FILE__, __FUNCTION__, __LINE__,
         "post and precalculated lengths for qual mods not equal"
         "(%d/%d)", mcounter1, lengths[1]);
  assertm(mcounter2==lengths[2], __FILE__, __FUNCTION__, __LINE__,
         "post and precalculated lengths for qual mods not equal"
         "(%d/%d)", mcounter2, lengths[2]);

  iar->qual_mods[mcounter1] = '"';
  iar->aligned_err_prob = err_prob;
  iar->distance = (min_len - nmatches>0) ? (min_len - nmatches):0;
  
  if (iar->distance > dif_threshold){
    reinitialize_ins_alig_res(iar, 5, 2, 2);
    write_to_ins_alig_res(iar, "OTHER", "NA", "NA", -1.0, -1);
    return;
  }
  
  
}


void align_long(char *read, char *quals, struct SeqList references, 
                wavefront_aligner_t *wf_aligner, float stop_when,
                char* qtype, int dif_threshold, struct Seq *extended_seqs,
                struct InsertAlignmentResult *iar, char** cigar_buffer,
                int *buffer_size){
  
  assertm(cigar_buffer != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "cigar_buffer is NULL");
  assertm(iar != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "iar is NULL");

  int i = 0;
  int cigar_len;
  int read_len = (int)strlen(read);
  int best_nmatches=-1;
  int current_nmatches;
  int min_len = read_len;
  
  for (i=0; i<references.nseqs; i++){

    wavefront_align(wf_aligner, references.seqs[i].s, references.seqs[i].l,
    								read, read_len);

                    
    if (references.seqs[i].l < min_len){min_len = references.seqs[i].l;}
    
    cigar_len =  wf_aligner->cigar->end_offset - wf_aligner->cigar->begin_offset;

    if (cigar_len+1 > *buffer_size){
      *cigar_buffer = realloc(*cigar_buffer, sizeof(char)*(cigar_len+1));
      *buffer_size = cigar_len+1;
      (*cigar_buffer)[cigar_len] = '\0';
    }
    else{
      (*cigar_buffer)[cigar_len] = '\0';
    }
    
    strcpy(*cigar_buffer,  wf_aligner->cigar->operations + 
           wf_aligner->cigar->begin_offset);

    current_nmatches = count_matches(*cigar_buffer);
    
    
    
    if (current_nmatches > best_nmatches){
        
        best_nmatches = current_nmatches;
        
        extend_read_and_reference_with_gaps(references.seqs[i].s, read,
                                            NULL, quals, *cigar_buffer,
                                            extended_seqs);
                                            
        get_insert_alignment_result(extended_seqs[1].s, extended_seqs[0].s,
                                    extended_seqs[1].q, *cigar_buffer,
                                    qtype, dif_threshold, iar, min_len);
        iar->index = i;
        
      }

    if ((float)current_nmatches/(float)min(read_len, references.seqs[i].l)
          >stop_when){
      break;
    }
  }
    
} 


void align_short(char *read, struct SeqList references, 
                 int max_edit_distance, int stop_when, char *best_loc,
                 int start, int buffer, struct AlignmentResult *ar){
                                      
  assertm(read != NULL, __FILE__, __FUNCTION__, __LINE__,
         "read is NULL");
  assertm(ar != NULL, __FILE__, __FUNCTION__, __LINE__,
         "AlignmentResult ar is NULL");
  
  int read_len = (int)strlen(read);
  
  assertm(buffer>=0 && start>=0, __FILE__, __FUNCTION__, __LINE__,
         "buffer and start must be non-negative but they are %d %d", 
         buffer, start);
  
  
  assertm(stop_when>=0 && max_edit_distance>=0,
         __FILE__, __FUNCTION__, __LINE__,
         "max_edit_distance (%d) or stop_when (%d) is negative",
         max_edit_distance, stop_when); 
  
  int i = 0;
  char* best_cigar = NULL;
  int cigar_len;
  int best_distance = read_len;
  int best_start=-1, best_end=-1, best_index=-1; 
  int newlen, end, loc_ind;
  
  EdlibAlignResult result;
 
  for (i=0; i<references.nseqs; i++){
    newlen = references.seqs[i].l + buffer;
    end = min(start + newlen, read_len);

    assertm(start<=end,
           __FILE__, __FUNCTION__, __LINE__,
           "read range start %d should be smaller than end %d",
           start, end);
  
    result = edlibAlign(references.seqs[i].s, references.seqs[i].l, 
                        read + start, end-start, 
                        edlibNewAlignConfig(max_edit_distance,
                        EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                                  
    if (result.status == EDLIB_STATUS_OK) {
      loc_ind = choose_loc(result.startLocations, result.numLocations,
                           best_loc);

      if (result.editDistance>-1 &&
          rescmp(best_distance, result.editDistance, 
                 best_start, result.startLocations[loc_ind], 
                 best_loc) ==1)
      {
        best_index = i;
        best_distance = result.editDistance;

        best_start = start + result.startLocations[loc_ind];
        best_end = start + result.endLocations[loc_ind];
        
        if(best_cigar != NULL){free(best_cigar);}

        best_cigar = edlibAlignmentToCigar(result.alignment, 
                                           result.alignmentLength, 
                                           EDLIB_CIGAR_EXTENDED);
      }
      
    if (result.editDistance<=stop_when && result.editDistance>-1){
      edlibFreeAlignResult(result);
      break;
      }
      
    }

    edlibFreeAlignResult(result);
    
  }


  ar->index = best_index;
  ar->distance = best_distance;
  ar->start = best_start;
  ar->end = best_end;
  
  if (best_cigar != NULL){
    invert_cigar(best_cigar); // for implementation reasons, 
                              // one has to put shorter sequence
                              // as query which means if there 
                              // is a deletion in the read
                              // in the cigar it will appear as I. 
                              // That can lead to confusions
                              // down the line so I am preemptively
                              // inverting cigar strings instead,
                              // changing I to D and D to I.
    
    cigar_len = (int)strlen(best_cigar);
    
    if (1 + cigar_len > ar->cigar_size){
        ar->cigar = realloc(ar->cigar, sizeof(char)*(1+cigar_len));
        ar->cigar_size = 1+cigar_len;
    }
    ar->cigar[cigar_len] = '\0';
    strcpy(ar->cigar, best_cigar);
    free(best_cigar); 
  }
    
  else{ ar->cigar[0] = '\0';}
    
} 


void align_pairs(kseq_t *seq_record1, kseq_t *seq_record2, 
                 struct Seq *rc_seq_buffer, char **cigar_buffer,
                 int *cigar_size, wavefront_aligner_t *wf_aligner, 
                 struct Seq *extended_seqs){
     
  int cigar_len;

  reverse_complement(rc_seq_buffer, seq_record2->seq.s, seq_record2->qual.s);
  wavefront_align(wf_aligner,seq_record1->seq.s, seq_record1->seq.l, 
                  rc_seq_buffer->s, seq_record2->seq.l);
 
  
  cigar_len =  wf_aligner->cigar->end_offset-wf_aligner->cigar->begin_offset;

  if (cigar_len+1 > *cigar_size){
    *cigar_buffer = realloc(*cigar_buffer, cigar_len+1);
    (*cigar_buffer)[cigar_len] = '\0';
    *cigar_size = cigar_len + 1;
  }
  else{
    (*cigar_buffer)[cigar_len] = '\0';
  }
  
  strcpy(*cigar_buffer,wf_aligner->cigar->operations+wf_aligner->cigar->begin_offset);

  extend_read_and_reference_with_gaps(seq_record1->seq.s, 
                                      rc_seq_buffer->s, seq_record1->qual.s,
                                      rc_seq_buffer->q, *cigar_buffer, extended_seqs);

}
