#include "../include/utils.h"
#include "../include/error.h"
#include "../include/structures.h"

inline static char* ambiguous_dna_complement(char *aa){

  switch(*aa){
    case 'A':
      return "T";
    case 'C':
      return "G";
    case 'G':
      return "C";
    case 'T':
      return "A";
    case 'M':
      return "K";
    case 'R':
      return "Y";
    case 'W':
      return "W";
    case 'S':
      return "S";
    case 'Y':
      return "R";
    case 'K':
      return "M";
    case 'V':
      return "B";
    case 'H':
      return "D";
    case 'D':
      return "H";
    case 'B':
      return "V";
    case 'X':
      return "X";
    case 'N':
      return "N";
    default:
      return aa;

  }

}



void write_to_ins_alig_res(struct InsertAlignmentResult *iar, char *seq_mods,
                           char *qual_mods, char *condensed_cigar, 
                           double aligned_err_prob, int index
                          )
{

  int len1 = (int)strlen(seq_mods);
  int len2 = (int)strlen(qual_mods);
  int len3 = (int)strlen(condensed_cigar);
  
  assertm(len1 <= iar->l[0],  __FILE__, __FUNCTION__, __LINE__, 
         "seq_mods length %d is bigger than iar->l[0] %d",
         len1, iar->l[0]);
  assertm(len2 <= iar->l[1],  __FILE__, __FUNCTION__, __LINE__, 
         "seq_mods length %d is bigger than iar->l[1] %d",
         len2, iar->l[1]);
  assertm(len3 <= iar->l[2],  __FILE__, __FUNCTION__, __LINE__, 
         "seq_mods length %d is bigger than iar->l[2] %d",
         len3, iar->l[2]);
  
  strcpy(iar->seq_mods, seq_mods);
  iar->l[0] = len1;
  iar->seq_mods[len1] = '\0';
  
  strcpy(iar->qual_mods, qual_mods);
  iar->qual_mods[len2] = '\0';
  iar->l[1] = '\0';
  
  strcpy(iar->condensed_cigar, condensed_cigar);
  iar->qual_mods[len3] = '\0';
  iar->l[2] = '\0';
  
  iar->aligned_err_prob = aligned_err_prob;
  iar->index = index;
}

void initialize_ins_alig_res(struct InsertAlignmentResult *iar, int length1,
                                     int length2, int length3){
  int i;

  assertm(iar != NULL, __FILE__, __FUNCTION__, __LINE__,
         "iar is NULL");
  
  iar->seq_mods = (char*)safe_malloc(sizeof(char)*(length1+1), "iar->seq_mods",
                                          __FILE__, __FUNCTION__, __LINE__);
  iar->qual_mods = (char*)safe_malloc(sizeof(char)*(length2+1), "iar->qual_mods",
                                          __FILE__, __FUNCTION__, __LINE__);
  iar->condensed_cigar = (char*)safe_malloc(sizeof(char)*(length3+1), 
                                                "iar->condensed_cigar",
                                               __FILE__, __FUNCTION__, __LINE__);
  iar->aligned_err_prob = -1;
  iar->index = -1;
  iar->distance = -1; 
  iar->l[0] = length1;
  iar->l[1] = length2;
  iar->l[2] = length3;
  
  for (i=0; i<3; i++){
    iar->s[i] = iar->l[i] + 1;
  }
}


void reinitialize_ins_alig_res(struct InsertAlignmentResult *iar, int length1,
                               int length2, int length3){
    
    assertm(length1>=0 && length2>=0 && length3>=0, __FILE__, __FUNCTION__, __LINE__,
           "lengths should be all non-negative but are %d %d %d", length1,
           length2, length3);
  
    if (length1+1>iar->s[0]){
      assertm(iar->seq_mods != NULL, __FILE__, __FUNCTION__, __LINE__,
              "iar->seq_mods is NULL");
      
      iar->seq_mods = (char*)realloc(iar->seq_mods, sizeof(char)*(length1+1));
      iar->l[0] = length1;
      iar->s[0] = length1+1;
      iar->seq_mods[length1] = '\0';
    }
    else if (length1+1 <= iar->s[0]){
      iar->seq_mods[length1] = '\0';
      iar->l[0] = length1;
    }
    
    
    if (length2+1>iar->s[1]){
      assertm(iar->qual_mods != NULL, __FILE__, __FUNCTION__, __LINE__,
              "iar->qual_mods is NULL");
      
      iar->qual_mods = (char*)realloc(iar->qual_mods, sizeof(char)*(length2+1));
      iar->l[1] = length2;
      iar->s[1] = length2+1;
      iar->qual_mods[length2] = '\0';
    }
    else if (length2+1 <= iar->s[1]){
      iar->qual_mods[length2] = '\0';
      iar->l[1] = length2;
    }
      
    
    if (length3+1>iar->s[2]){
      assertm(iar->condensed_cigar != NULL, __FILE__, __FUNCTION__, __LINE__,
              "iar->condensed_cigar is NULL");
      
      iar->condensed_cigar = (char*)realloc(iar->condensed_cigar, sizeof(char)*(length3+1));
      iar->l[2] = length3;
      iar->s[2] = length3+1;
      iar->condensed_cigar[length3] = '\0';
    }
    else if (length3+1 <= iar->s[2]){
      iar->condensed_cigar[length3] = '\0';
      iar->l[2] = length3;
    }
    
    iar->aligned_err_prob = -1;
    iar->index = -1;
  
}

void free_ins_alig_res(struct InsertAlignmentResult *iar)
{
  int i = 0;
  if (iar == NULL)return;
  if (iar->seq_mods != NULL) free(iar->seq_mods);
  if (iar->qual_mods != NULL) free(iar->qual_mods);
  if (iar->condensed_cigar != NULL) free(iar->condensed_cigar);

  for (i=0; i<3; i++){
      iar->s[i] = 0;
      iar->l[i] = -1;
  }
}



void initialize_seq(struct Seq *seq, int length){

  assertm(seq != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "seq is NULL");
 
  seq->s = (char*)safe_malloc((length+1)*sizeof(char), 
                              "seq->s", __FILE__, __FUNCTION__, 
                              __LINE__);
  seq->s[length] = '\0';
  
  seq->q = (char*)safe_malloc((length+1)*sizeof(char), 
                              "seq->q", __FILE__, __FUNCTION__, 
                              __LINE__);
  seq->q[length] = '\0';

  seq->q[length] = '\0';
  seq->l = length;
  seq->size = length;
  
  seq->id = NULL;

}


void reinitialize_seq(struct Seq* seq, int length){

  assertm(seq != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "seq is NULL");
  assertm(seq->s != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "seq->s is NULL");
  assertm(length>=0, __FILE__, __FUNCTION__, __LINE__, 
         "length must be non-negative but it is %d", length);
  
  if (length > seq->l) seq->s = realloc(seq->s, length+1);
  assertm(seq->s != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "seq->s is NULL after realloc");
  seq->s[length] = '\0';
  
  if (seq->q != NULL){
    if (length > seq->l)
    {
      seq->q = realloc(seq->q, length+1);
      assertm(seq->q != NULL, __FILE__, __FUNCTION__, __LINE__, 
             "seq->q is NULL after realloc");
    }  
    seq->q[length] = '\0';
  }
  
  if (seq->id != NULL){
      free(seq->id);
      seq->id = NULL;
  }
  
  seq->l = length;
  seq->size = length;
  
  
  
}


void write_to_seq(struct Seq* seq, char* s, char* q, char *id){

  assertm(seq != NULL, __FILE__, __FUNCTION__, __LINE__, "seq is NULL");
  assertm(s != NULL, __FILE__, __FUNCTION__, __LINE__, "s is NULL");
  
  int slen = strlen(s);
  int idlen; 
  
  if (id != NULL) idlen = strlen(id);
  
  if( q!= NULL){ 
      int qlen = strlen(q);
      assertm(qlen == slen, __FILE__, __FUNCTION__, __LINE__, 
             "length of q (%ld) not equal to length of seq (%d)",
              qlen, slen);
  }
  
  
  if (seq->size < slen){
    reinitialize_seq(seq, slen);
  }
  else{
    seq->s[slen] = '\0';
    seq->q[slen] = '\0';
    seq->l = slen;
  }

  strcpy(seq->s, s);
  
  if (q != NULL){
    strcpy(seq->q, q);
  }
  else{
    if (seq->q != NULL)free(seq->q);
    seq->q = NULL;
    }
        
  if (id != NULL && (seq->id == NULL || idlen > (int)strlen(seq->id))){
      seq->id = realloc(seq->id, idlen+1);
      assertm(seq->id != NULL, __FILE__, __FUNCTION__, __LINE__,
             "seq->id could not be reallocated");
      strcpy(seq->id, id);
      seq->id[idlen] = '\0';
  }
  else if (id == NULL && seq->id != NULL){
      free(seq->id);
      seq->id = NULL;
  }
  
}


void free_seq(struct Seq *seq){
  
  if (seq == NULL){
    return;
  }

  if (seq->s != NULL){free(seq->s); seq->s=NULL;}
  if (seq->q != NULL){free(seq->q); seq->q=NULL;}
  if (seq->id != NULL){free(seq->id); seq->id=NULL;}

}


void initialize_alignment_result(struct AlignmentResult *res, 
                                 int cigar_len)
{
  assertm(res != NULL, __FILE__, __FUNCTION__, __LINE__,
         "AlignmentResult res is NULL");
  
  res->cigar = (char*)safe_malloc((cigar_len+1)*sizeof(char), 
                                  "res->cigar", __FILE__, __FUNCTION__, 
                                  __LINE__);
  res->cigar[cigar_len] = '\0';
  
  res->cigar_size = cigar_len+1;
  res->index = -1;
  res->start = -1;
  res->end = -1;
  res->distance = -1;
  

}


void free_alignment_result(struct AlignmentResult *res){
  
  if (res == NULL){
    return;
  }

  if (res->cigar != NULL){
    free(res->cigar);
    res->cigar = NULL;
  }
}


void invert_cigar(char *cigar){
  
  int i;
  
  if (cigar == NULL){
    return;
  }

  int len = (int)strlen(cigar);
  
  for (i=0; i<len; i++){
    if (cigar[i] == 'I'){
      cigar[i] = 'D';
    }
    else if (cigar[i] == 'D'){
      cigar[i] = 'I';
    }
  }


}


void initialize_seq_list(struct SeqList *seq_list, int nseqs){

  assertm(seq_list != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "SeqList is NULL");

  seq_list->seqs = (struct Seq*)safe_malloc(nseqs*sizeof(struct Seq), 
                                            "seq_list->seqs", 
                                            __FILE__, __FUNCTION__, 
                                            __LINE__);
  seq_list->nseqs = nseqs;
  seq_list->min_length = -1;
  seq_list->max_length = -1;
}


void add_seq_to_seqlist(struct SeqList* seq_list, int index, char* s, 
                        char* q, char* id){

  assertm(index < seq_list->nseqs, __FILE__, __FUNCTION__, __LINE__, 
         "index %d is out of bounds for SeqList of size %d", 
         index, seq_list->nseqs);
  initialize_seq(&(seq_list->seqs[index]), strlen(s));
  write_to_seq(&(seq_list->seqs[index]), s, q, id);
  
  if ((int)strlen(s) < seq_list->min_length || seq_list->min_length==-1){
  	seq_list->min_length = (int)strlen(s);
  }
  
  if ((int)strlen(s) > seq_list->max_length){
  	seq_list->max_length = (int)strlen(s);
  }
    
}


void free_seq_list(struct SeqList *seq_list){
  int i;

  if (seq_list == NULL){
    return;
  }

  for (i=0; i<seq_list->nseqs; i++)
  {
    free_seq(&(seq_list->seqs[i]));
  }

  free(seq_list->seqs);
  seq_list->seqs = NULL;
}


void reverse_complement(struct Seq *rc_seq, char* sequence, 
                        char* qual) {

  int i, j;
  
  assertm(sequence != NULL, __FILE__, __FUNCTION__, __LINE__,  
         "sequence is NULL");
  assertm(rc_seq != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "rc_seq is NULL");
  
  int slen = (int)strlen(sequence);
  
  if (qual != NULL)
  {
    int qlen = strlen(qual);
    assertm(qlen == slen,  __FILE__, __FUNCTION__, __LINE__, 
           "sequence length (%d) should be equal to qual length (%ld)", 
           slen, qlen);
  }
  
  if (rc_seq->size < slen)
  {
    reinitialize_seq(rc_seq, slen); 
  }
  else{
    rc_seq->s[slen] = '\0';
    rc_seq->q[slen] = '\0';
    rc_seq->l = slen;
  }
  
  for (i = slen - 1, j = 0; i >= 0; i--, j++) {
      rc_seq->s[i] = *ambiguous_dna_complement(&sequence[j]);
      
      if (qual != NULL){
        rc_seq->q[i] = qual[j];
      }
  }
  
  if (qual == NULL && rc_seq->q != NULL){
    free(rc_seq->q);
    rc_seq->q = NULL;
  }
  

}


void reverse_complement_list(struct SeqList *seq_list_reverse,
                             struct SeqList seq_list){

  struct Seq rev_seq;
  int i;  
  
  initialize_seq_list(seq_list_reverse, seq_list.nseqs);
  
  for (i=0; i<seq_list.nseqs; i++){
    initialize_seq(&rev_seq, seq_list.seqs[i].l);
    reverse_complement(&rev_seq, seq_list.seqs[i].s, seq_list.seqs[i].q);
    
    initialize_seq(&seq_list_reverse->seqs[i], rev_seq.l);
    write_to_seq(&seq_list_reverse->seqs[i], rev_seq.s, rev_seq.q, seq_list.seqs[i].id);
    free_seq(&rev_seq);
  }  
  
  
}

