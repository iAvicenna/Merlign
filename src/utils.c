#include "../include/utils.h"
#include "../include/error.h"

UTILITY_STATIC(char encode_from_error_probability(double error_prob, char *qtype)){
  int offset = get_offset(qtype);
  int Q = (int)round(-10*log10(error_prob)) + offset;
  
  /*merged quality scores can be out of the usual boundaries 
  so they should be dealt seperately in order not to produce 
  weird decoded letters down the line */
  
  if (Q>126) return '-'; 
  if (Q<33) return '!';  

  return (char)((int)round(-10*log10(error_prob)) + offset);
}

UTILITY_STATIC(int decode_char(char qual, char *qtype)){
  
  int offset = get_offset(qtype);
  return (int)qual - offset;

}


UTILITY_STATIC(double *decode_to_error_probabilities(char *quals, char *qtype)){
  assertm(quals != NULL, __FILE__, __FUNCTION__, __LINE__, "quals is NULL");
    
  int i;
  int offset = get_offset(qtype);
  int len = (int)strlen(quals);
  double *probs = (double*)safe_malloc(len*sizeof(double), "probs", 
                                       __FILE__, __FUNCTION__, __LINE__);
  
  for (i=0; i < len; i++){
    probs[i] = decode_to_error_probability(quals[i], offset);
  }
  
  return probs;

}


UTILITY_STATIC(int count_operations(char *cigar)){
  assertm(cigar != NULL, __FILE__, __FUNCTION__, __LINE__, "cigar is NULL");
    
  int i,j;
  int sum = 0;
  int prev_index = -1;
  int len = (int)strlen(cigar);
  
  if (strcmp(cigar, "")==0){
    return 0;
  }
 
  assertm(isdigit(cigar[0])>0, __FILE__, __FUNCTION__, __LINE__, 
         "Invalid CIGAR starting with non-digit: %s", cigar);
  assertm(isdigit(cigar[len-1])==0, __FILE__, __FUNCTION__, __LINE__, 
         "Invalid CIGAR ending with digit: %s", cigar);
 
  for (i=0; i < len; i++){
    if (!isdigit(cigar[i])){
      for (j=i-1; j>prev_index; j--){
        sum = sum + pow(10.0,(double)(i-1-j))*(cigar[j]-'0');
      }
      
      prev_index = i;
      assertm(i >= len-1 || isdigit(cigar[i+1])>0, __FILE__, __FUNCTION__, __LINE__, 
             "Invalid CIGAR starting after index %d: %s", i, cigar);

    }
  } 
  
  return sum;
}


UTILITY_STATIC(char *get_cigar_array(char *cigar, int num_operations)){

  assertm(cigar != NULL, __FILE__, __FUNCTION__, __LINE__, "cigar is NULL");

  char *cigar_array = (char*)safe_malloc((num_operations+1)*sizeof(char), 
                                         "cigar_array", __FILE__, __FUNCTION__, __LINE__);
  cigar_array[num_operations] = '\0';
  
  int i,j;
  int sum=0;
  int counter=-1;
  int prev_index=-1;
  int len = (int)strlen(cigar);
  
  for(i=0; i < len; i++){
    if (!isdigit(cigar[i])){
      for (j=i-1; j>prev_index; j--){
        sum = sum + pow(10.0,(double)(i-1-j))*(cigar[j]-'0');
      }
      prev_index = i;

      for (j=0; j<sum; j++){
        counter += 1;        
        cigar_array[counter] = cigar[i];
      }
      
      sum=0;
    }
  }
  
  cigar_array[num_operations] = '\0';
  
  return cigar_array;
  
}


void* safe_malloc(size_t n, char* name, char* file, const char *function, int line)
{
    assertm(n>0, file, function, line, 
           "%s with negative size %ld can not be allocated", name, n);
           
    void* p = malloc(n);
    
    assertm(p!=NULL, file, function, line, 
           "%s of size %ld could not be allocated", name, n);
    
    return p;
}


int get_offset(char* qtype){
    
  assertm(strcmp(qtype,"Solexa+64")==0 || strcmp(qtype,"Phred+33")==0
         || strcmp(qtype,"Phred+64")==0, __FILE__, __FUNCTION__, __LINE__, 
         "Unknown quality type %s", qtype);

  if (strcmp(qtype,"Phred+33")==0){
    return 33;
  }
  
  else if (strcmp(qtype,"Solexa+64")==0 || strcmp(qtype, "Phred+64")==0){
    return 64;
  }
  
  exit(EXIT_FAILURE); // will never reach here due to assert but here for disabling warnings
}


double decode_to_error_probability(char qual, int offset){

  double num_qual = (int)qual - offset;

  return pow(10.0, -(double)num_qual/10);

}


void merge_seqs(char *seq1, char *quals1, char *seq2, char *quals2, 
                char *qtype, struct Seq *merged_seq, int mask_threshold)
{
    /*
      This merges two sequences and quals based on the 
      Bayesian framework implemented here is from

      Robert C. Edgar, and Henrik Flyvbjerg

      "Error filtering, pair assembly and error
      correction for next-generation sequencing reads"

      Bioinformatics, 31(21), 2015, 3476–3482
      doi: 10.1093/bioinformatics/btv401
      
      Note that sequences and qualities must be of equal length and where
      one of the sequences has a deletion, corresponding sequence and quality
      letter should be empty space. Example
      
      Seq1  = 'AGCT G'
      Qual1 = 'AAAC D'
      Seq2  = ' GCTTG'
      Qual2 = ' CCCDA
      
      The paper above uses letter px,py for error probabilities and
      X,Y for read letters. We will follow the same convention in the loop
      that computes the merge.
      
      
    */
    int slen1, slen2, qlen1, qlen2, i;
    double px,py,merged_prob;
    double *err_probs1;
    double *err_probs2;
    
    assertm(seq1 != NULL,  __FILE__, __FUNCTION__, __LINE__, "seq1 is NULL");
    assertm(seq2 != NULL,  __FILE__, __FUNCTION__, __LINE__, "seq2 is NULL");
    assertm(quals1 != NULL,  __FILE__, __FUNCTION__, __LINE__, "quals1 is NULL");
    assertm(quals2 != NULL,  __FILE__, __FUNCTION__, __LINE__, "quals2 is NULL");
    assertm(merged_seq != NULL,  __FILE__, __FUNCTION__, __LINE__, 
           "merged_seq is NULL");
    
    slen1 = (int)strlen(seq1);
    slen2 = (int)strlen(seq2);
    qlen1 = (int)strlen(quals1);
    qlen2 = (int)strlen(quals2);
    
    assertm(slen1 == slen2, __FILE__, __FUNCTION__, __LINE__, 
           "Length of sequences is not equal: len1=%d, len2=%d", slen1, slen2);
    assertm(slen1 == qlen1, __FILE__, __FUNCTION__, __LINE__, 
           "Length of qual1 %d and length of seq1 %d is not equal", qlen1, 
           slen1);
    assertm(slen2 == qlen2, __FILE__, __FUNCTION__, __LINE__, 
           "Length of qual2 %d and length of seq2 %d is not equal", qlen2, 
           slen2);
    
    
    if (slen1 > merged_seq->size)
    {
        reinitialize_seq(merged_seq, slen1);
    }
    else{
        merged_seq->s[slen1] = '\0';
        merged_seq->q[slen1] = '\0';
        merged_seq->l = slen1;
    }
              
    err_probs1 = decode_to_error_probabilities(quals1, qtype);
    err_probs2 = decode_to_error_probabilities(quals2, qtype);

    for (i=0; i<slen1; i++){
    
      if (seq1[i]==' '){
        merged_seq->s[i] = seq2[i];
        merged_seq->q[i] = quals2[i];
      }
      else if (seq2[i]==' '){
        merged_seq->s[i] = seq1[i];
        merged_seq->q[i] = quals1[i];
      }
      
      else if (seq1[i]==seq2[i]){
        px = err_probs1[i];
        py = err_probs2[i];
        
        merged_seq->s[i] = seq1[i];
        
        merged_prob = (px*py/3.0)/(1.0-px-py+4.0*px*py/3.0);
        merged_seq->q[i] = encode_from_error_probability(merged_prob, qtype);
        
      }
      else{ /*seq1[i]!=seq2[i] and non-gap*/
        if (err_probs1[i]<=err_probs2[i]){
          px = err_probs1[i];
          py = err_probs2[i];
          merged_seq->s[i] = seq1[i];
        }
        else{
          px = err_probs2[i];
          py = err_probs1[i];
          merged_seq->s[i] = seq2[i];
        }

        merged_prob = px*(1-py/3)/(px+py-4*px*py/3);
        merged_seq->q[i] = encode_from_error_probability(merged_prob, qtype);
        
      }
      
      /*if threshold not satisfied mask the letter by making it small case*/
    	if (decode_char(merged_seq->q[i], qtype)<mask_threshold){
    		merged_seq->s[i] = tolower(merged_seq->s[i]);
    		}
    }
    
    free(err_probs1);
    free(err_probs2);
}  


void extend_read_and_reference_with_gaps(char *seq1, char *seq2, char *quals1,
                                         char *quals2, char *cigar_array,
                                         struct Seq extended_seqs[2])
{

    
  assertm(seq1 != NULL && seq2 != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "seq1 or seq2 is NULL");
  
  int slen1 = (int)strlen(seq1);
  int slen2 = (int)strlen(seq2);
  int qlen1;
  int qlen2;
  
	if (quals1 != NULL){
		qlen1 = strlen(quals1);
		assertm(slen1 == qlen1, __FILE__, __FUNCTION__, __LINE__, 
					"seq1 and quals1 do not have the same lengths");
	}
	if (quals2 != NULL){
		qlen2 = strlen(quals2);
		assertm(slen2 == qlen2, __FILE__, __FUNCTION__, __LINE__, 
					"seq2 and quals2 do not have the same lengths");
	}
		
  int i;
  int counter1=-1;
  int counter2=-1;
  int counter3=0;
  int cigar_len = (int)strlen(cigar_array);
  
  for (i=0; i<2; i++){
    reinitialize_seq(&extended_seqs[i], cigar_len);
  }
  
  for(i=0; i < cigar_len; i++){
  
    assertm(counter1<slen1 && counter2<slen2, __FILE__, __FUNCTION__, __LINE__, 
           "CIGAR string exceeds either seq1 or seq2 length");
    
    if (cigar_array[i]=='I' || cigar_array[i]=='S'){
      counter2 += 1;
      extended_seqs[0].s[counter3] = ' ';
      extended_seqs[1].s[counter3] = seq2[counter2];
      if (quals1 != NULL) extended_seqs[0].q[counter3] = ' ';
      if (quals2 != NULL) extended_seqs[1].q[counter3] = quals2[counter2];
      
    }
    else if (cigar_array[i]=='D'){
      counter1 += 1;
      extended_seqs[0].s[counter3] = seq1[counter1];
      extended_seqs[1].s[counter3] = ' ';
      if (quals1 != NULL)extended_seqs[0].q[counter3] = quals1[counter1];
      if (quals2 != NULL)extended_seqs[1].q[counter3] = ' ';
      
    }
    else if (cigar_array[i]=='M' || cigar_array[i]=='=' || cigar_array[i]=='X'){
      counter1 += 1;
      counter2 += 1;
      extended_seqs[0].s[counter3] = seq1[counter1];
      extended_seqs[1].s[counter3] = seq2[counter2];
      if (quals1 != NULL)extended_seqs[0].q[counter3] = quals1[counter1];
      if (quals2 != NULL)extended_seqs[1].q[counter3] = quals2[counter2];
    }
    
    counter3 += 1;

  }
  
  assertm(counter1 == slen1-1, __FILE__, __FUNCTION__, __LINE__,  
         "CIGAR %s does not exhaust seq1 %s", cigar_array, seq1);
  assertm(counter2 == slen2-1, __FILE__, __FUNCTION__, __LINE__,  
         "CIGAR %s does not exhaust seq2 %s", cigar_array, seq2);
  
}


int count_matches(char *cigar){

  assertm(cigar != NULL, __FILE__, __FUNCTION__, __LINE__,
         "cigar cant be NULL");

  int slen = strlen(cigar);
  int i = 0;
  int mcounter=0; 

  for (i=0; i<slen; i++){
    if (cigar[i] == 'M') mcounter++;
  }

  return mcounter;

}

void get_cigar_stats(char *cigar, int soft_clip, int convert_cigar, 
                     int begin_free1, int begin_free2, int end_free1, 
                     int end_free2, int result[6]){

  assertm(soft_clip == 0 || soft_clip ==1, __FILE__, __FUNCTION__, __LINE__, 
         "soft_clip can be 0 or 1 but is %d", soft_clip);
  assertm(convert_cigar == 0 || convert_cigar ==1, __FILE__, __FUNCTION__, 
         __LINE__, "convert_cigar can be 0 or 1 but is %d", convert_cigar);  
  
  char *cigar_array; 
  
  int i, start, end, nmismatches, nmatches, ninsertions, 
      ndeletions, nsoftstart, nsoftend;

  nmismatches = 0;
  nmatches = 0;
  ninsertions = 0;
  ndeletions = 0;
  nsoftstart = 0;
  nsoftend = 0;
  
  result[0] = nsoftstart;
  result[1] = nsoftend;
  result[2] = nmatches;
  result[3] = nmismatches;
  result[4] = ndeletions;
  result[5] = ninsertions;
  
  if (cigar == NULL){
    for (i=0; i<6; i++) result[i] = -1;
    return;
  }
  
  if (convert_cigar==1){
    int num_operations = count_operations(cigar);
    cigar_array = get_cigar_array(cigar, num_operations);
  }
  else{
    cigar_array = cigar;
  }
  
  int cigar_len = strlen(cigar_array);
  
  start  = 0;
  end = cigar_len; 
  
  if (soft_clip==1){
  
    if(begin_free1>0 && cigar_array[0] == 'D')
    {
      nsoftstart += 1;
      start += 1;
    }
    
    else if(begin_free2>0 && cigar_array[0] == 'I')
    {
      nsoftstart += 1;
      start += 1;
    }

    for (i=1; i<cigar_len; i++)
    {
      if(
         cigar_array[i] == cigar_array[i-1] &&
         ((cigar_array[i] == 'D' && nsoftstart < begin_free1) ||
          (cigar_array[i] == 'I' && nsoftstart < begin_free2))  
        )
      {   
        start += 1;
        nsoftstart += 1; 
      } 
      else{
        break;
      }
    
    }

    if(end_free1>0 && cigar_array[end-1] == 'D' &&
       nsoftstart != cigar_len)
    {
      nsoftend += 1;
      end -= 1;
    }
  
    else if(end_free2>0 && cigar_array[end-1] == 'I' &&
            nsoftstart != cigar_len){
      nsoftend += 1;
      end -= 1;
    }
   
    for (i=cigar_len-2; i>=nsoftstart; i--)
    {
      if(
         cigar_array[i] == cigar_array[i+1] &&
         ((cigar_array[i] == 'D' && nsoftend < end_free1) ||
          (cigar_array[i] == 'I' && nsoftend < end_free2))
        )  
      {
          end -= 1;
          nsoftend += 1;
      }
      else{
        break;
      } 
    
    }
    
  }


  for (i=start; i<end; i++)
  {
    if (cigar_array[i]=='M' || cigar_array[i]=='=')nmatches +=1;
    if (cigar_array[i]=='X')nmismatches +=1;  
    if (cigar_array[i]=='D')ndeletions += 1;
    if (cigar_array[i]=='I')ninsertions += 1;
  }
  
  result[0] = nsoftstart;
  result[1] = nsoftend;
  result[2] = nmatches;
  result[3] = nmismatches;
  result[4] = ndeletions;
  result[5] = ninsertions;
  
  if (convert_cigar){free(cigar_array);}

}

static inline int min(int x, int y) {
  return (x < y) ? x : y;
}

static inline int max(int x, int y) {
  return (x > y) ? x : y;
}

