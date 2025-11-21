#include "../include/io.h"
#include "../include/error.h"
#define MAX_NAME_LEN 1000
#define BUFFERSIZE 10000

UTILITY_STATIC(int count_sequences(char *file_name))
{
  gzFile fp;
  kseq_t *seq;
  int l,counter;
  fp = gzopen(file_name, "r"); // STEP 2: open the file handler
  assertm(fp != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "Can not open or no sequences in file path %s", file_name);
  
  seq = kseq_init(fp); // STEP 3: initialize seq
  counter = 0;
  while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
      counter += 1;
  }
  kseq_destroy(seq); // STEP 5: destroy seq
  gzclose(fp); // STEP 6: close the file handler

  return counter;
}


struct SeqList load_sequences(char *file_name, int add_rc){

  int l,counter,nseqs;
  struct SeqList seq_list;
  struct Seq rc_seq;
  
  gzFile fp;
  kseq_t *seq_record;
  fp = gzopen(file_name, "r");
  assertm(fp != NULL, __FILE__, __FUNCTION__, __LINE__, 
         "Can not open or no sequences in file path %s", file_name);
  
  nseqs = count_sequences(file_name);
  
  if (add_rc==1){
    nseqs = 2*nseqs;
  }
  
  /*seq_list.seqs = (struct Seq*)malloc(nseqs*sizeof(struct Seq));
  seq_list.nseqs = nseqs;*/
  initialize_seq_list(&seq_list, nseqs);

  seq_record = kseq_init(fp);
  counter = 0;

  while ((l = kseq_read(seq_record)) >= 0) { // STEP 4: read sequence

			add_seq_to_seqlist(&seq_list, counter, seq_record->seq.s,
												 seq_record->qual.s, seq_record->name.s);

      /*initialize_seq(&seq_list.seqs[counter], seq_record->seq.l);
      write_to_seq(&seq_list.seqs[counter], seq_record->seq.s, 
                   seq_record->qual.s, seq_record->name.s);*/
      counter += 1;
      
      if (add_rc==1){
          initialize_seq(&rc_seq, 1);
          reverse_complement(&rc_seq, seq_record->seq.s, 
                             seq_record->qual.s);
          /*initialize_seq(&seq_list.seqs[counter], rc_seq.l);
          write_to_seq(&seq_list.seqs[counter], rc_seq.s, rc_seq.q, seq_record->name.s);*/
          add_seq_to_seqlist(&seq_list, counter, rc_seq.s, rc_seq.q, 
          									 seq_record->name.s);
          
          counter += 1;
          free_seq(&rc_seq);
      }
      
  }

  kseq_destroy(seq_record);
  gzclose(fp);

  return seq_list;

}


int get_seq(gzFile* fp, char* buffer, struct Seq *seq, kseq_t *seq_record,
						int debug){
  
    int i=0;
    int l1=0;
    char seq_id[MAX_NAME_LEN+1];
    
    
    if (*fp == NULL){
    
      for (i=0; i<4; i++){
        if (fgets(buffer, BUFFERSIZE , stdin)==NULL) return 0;
        buffer[strcspn(buffer, "\n\r")] = 0;
        
        if (i==0){
          assertm(buffer[0]=='@', __FILE__, __FUNCTION__, __LINE__,
                 "seq_id should start with @ but is %s", buffer);
          buffer[strcspn(buffer, " ")] = 0; //spaces in ids are interpreted as start of comments
          strncpy(seq_id, buffer+1, MAX_NAME_LEN);
          seq_id[MAX_NAME_LEN] = '\0';
          
          if (debug==1){
						fprintf(stderr, "%s\n", seq_id);
					}	

          
        }
        
        if (i==1){
          reinitialize_seq(seq, strlen(buffer));
          seq->id = malloc(strlen(seq_id)+1);
          assertm(seq->id != NULL, __FILE__, __FUNCTION__, __LINE__, 
                 "seq->id is NULL");
                 
          strcpy(seq->s, buffer);
          strcpy(seq->id, seq_id);

        }
        if (i==2){
          assertm(strcmp(buffer,"+")==0, __FILE__, __FUNCTION__, __LINE__, 
                 "Line %d is %s but it should be +", i, buffer);
        }
        if (i==3)
        {
          strcpy(seq->q, buffer);
          return 1;
        }
        
      }
    }
    else{
			
			assertm(seq_record != NULL, __FILE__, __FUNCTION__, __LINE__, 
						 "seq_record is NULL");
			
      l1 = kseq_read(seq_record);
      
      if (debug==1){
    		fprintf(stderr, "%s\n", seq_record->name.s);

    	}	

      
      if (l1<0) return 0;
      
      reinitialize_seq(seq, seq_record->seq.l);
      write_to_seq(seq, seq_record->seq.s, seq_record->qual.s, 
                   seq_record->name.s);
      
      return 1;
    }
  
  return 0;
}

