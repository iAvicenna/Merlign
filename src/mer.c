#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include <getopt.h>
#include <ctype.h>
#include "../ext/WFA2/wavefront/wfa.h"
#include "../include/structures.h"
#include "../include/io.h"
#include "../include/alignment.h"
#include "../include/error.h"


static inline int min(int x, int y) {
  return (x < y) ? x : y;
}


void parse_args(int argc, char **argv, char **seq_file1_path, char **seq_file2_path,
                char **output_merged_path, char** output_path, int *mismatch, int *match,
                int *gap_opening, int *gap_extension, int *begin_free1, int *begin_free2, 
                int *end_free1, int *end_free2, int *N0, int *N1, char **qtype, int *trim,
                int *mask_threshold, int *merge_debug)
{
  int c;

  while (1)
  {
    static struct option long_options[] =
    {
      {"file1",    required_argument, 0, 'F'},
      {"file2",    required_argument, 0, 'G'},
      {"merged_path",    required_argument, 0, 'M'},
      {"output_path", required_argument, 0, 'O'},
      {"match", required_argument, 0, 'q'},
      {"mismatch", required_argument, 0, 'w'},
      {"gap_open", required_argument, 0, 'e'},
      {"gap_extend", required_argument, 0, 'r'},
      {"end_free1", required_argument, 0, 't'},
      {"end_free2", required_argument, 0, 'y'},
      {"begin_free1", required_argument, 0, 'u'},
      {"begin_free2", required_argument, 0, 'i'},
      {"N0", required_argument, 0, 'B'},
      {"N1", required_argument, 0, 'N'},
      {"qtype", required_argument, 0, 'Q'},
      {"trim", required_argument, 0, 'T'},
      {"mask_threshold", required_argument, 0, 'K'},
      {"merge_debug", required_argument, 0, 'D'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "F:G:M:O:q:w:e:r:t:y:u:i:B:N:Q:T:K:D:",
                      long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
      {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
          break;
        printf ("option %s", long_options[option_index].name);
        if (optarg)
          printf (" with arg %s", optarg);
        printf ("\n");
        break;
      case 'B':
        *N0 = atoi(optarg);
        break;
      case 'N':
        *N1 = atoi(optarg);
        break;
      case 'q':
        *match = atoi(optarg);
        break;
      case 'w':
        *mismatch = atoi(optarg);
        break;
      case 'e':
        *gap_opening = atoi(optarg);
        break;
      case 'r':
        *gap_extension = atoi(optarg);
        break;
      case 't':
        *end_free1 = atoi(optarg);
        break;
      case 'y':
        *end_free2 = atoi(optarg);
        break;
      case 'u':
        *begin_free1 = atoi(optarg);
        break;
      case 'i':
        *begin_free2 = atoi(optarg);
        break;    
      case 'F':
        *seq_file1_path = optarg;
        break;
      case 'G':
        *seq_file2_path = optarg;
        break;
      case 'M':
        *output_merged_path = optarg;
        break;
      case 'O':
        *output_path = optarg;
        break;
      case 'Q':
        *qtype = optarg;
        break;
      case 'T':
        *trim = atoi(optarg);
        break;
      case 'K':
        *mask_threshold = atoi(optarg);
        break;
      case 'D':
      	*merge_debug = atoi(optarg);
      	break;
      case '?':
        if (optopt == 'c')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                    "Unknown option character `\\x%x'.\n",
                    optopt);
        return;

      default:
        abort ();
      }
  }
  
  check_parameters_mer(*match, *mismatch, *gap_opening, *gap_extension, 
      *end_free1, *end_free2, *begin_free1, *begin_free2, *mask_threshold,
      *merge_debug);
  
  assertm(*seq_file1_path != NULL && seq_file2_path != NULL, NULL, NULL, -1,  
         "Minimal usage us mer --file1 seq_file1_path --file2 seq_file2_path");
  check_file(*seq_file1_path, "--file1");
  check_file(*seq_file2_path, "--file2");
  
  char *str_time;
  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  str_time = (char*)calloc(strlen(asctime(timeinfo)), sizeof(char));
  assertm(str_time != NULL, NULL, NULL, -1, "Out of storage for size %ld", 
      strlen(asctime(timeinfo))*sizeof(char));
      strncpy(str_time, asctime(timeinfo), strlen(asctime(timeinfo))-1);

  fprintf(stderr, "\nStart time: %s (UT: %u).\n", str_time, (unsigned)time(NULL));
  fprintf(stderr, "sequences file1 path=%s, sequences file2 path=%s\n"
          "N0=%d, N1=%d, qtype=%s, trim=%d, mask_threshold=%d, merge_debug=%d", 
          *seq_file1_path, *seq_file2_path, *N0, *N1, *qtype, *trim, 
          *mask_threshold, *merge_debug);

  if (*output_path != NULL){
    fprintf(stderr, "output path = %s\n", *output_path);
  }
  if (*output_merged_path != NULL){
    fprintf(stderr, "merged sequences path = %s\n", *output_merged_path);
  }

  fprintf(stderr, "\npenalties: \n");
  fprintf(stderr, "match=%d, mismatch=%d, ", *match, *mismatch);
  fprintf(stderr, "gap open=%d, gap extension=%d, ", *gap_opening, *gap_extension);
  fprintf(stderr, "begin_free1=%d, begin_free2=%d, end_free1=%d, end_free2=%d\n\n", 
                  *begin_free1, *begin_free2, *end_free1, *end_free2);
 
  free(str_time);

}


FILE *initialize_output_file(char *output_path){

  FILE *fp; 
  fp = fopen(output_path, "w");

  if (fp == NULL){
   fprintf(stderr, "Can not open file %s, exiting.\n", output_path);
   exit(EXIT_FAILURE);
  }

  fprintf(fp, "read_id");

  fprintf(fp, ",M_start,M_end,M_match,M_mismatch,M_deletion,M_insertion\n");
  
  return fp;
}


void add_alignment_result_to_file(FILE* fp, int *cigar_stats, kseq_t *seq){

  int i;
  
  fputs(seq->name.s, fp);
  fputs(",", fp);
  
  fprintf(fp, "%d,", cigar_stats[0]);
  fprintf(fp, "%d,", (int)seq->seq.l - cigar_stats[1]);

  for(i=2; i<6; i++){
    if (cigar_stats != NULL){
      fprintf(fp, "%d",cigar_stats[i]);
    }
    else fprintf(fp, "0");
    if (i<5)fprintf(fp, ",");
  }
  fprintf(fp, "\n");
}     


void fputseq(FILE *fp, kseq_t seq_record, struct Seq merged_seq, int start,
             int end){

 assertm(start<=end, __FILE__, __FUNCTION__, __LINE__,
        "start (%d) can not be bigger than end (%d)",
        start,end);

 fputs("@", fp);
 fputs(seq_record.name.s, fp);
 if (seq_record.comment.l) {fputs(" ", fp); fputs(seq_record.comment.s, fp);}
 fputs("\n", fp);

 merged_seq.s[end]='\0';
 merged_seq.q[end]='\0';

 fputs(merged_seq.s+start, fp);
 fputs("\n+\n", fp);
 fputs(merged_seq.q+start, fp);
 fputs("\n", fp);

}


void initialize_files(gzFile *fp1, gzFile *fp2, FILE **fp3, FILE **fp4,
                      char *seq_file1_path, char *seq_file2_path,
                      char * output_merged_path, char *output_path){

  *fp1 = gzopen(seq_file1_path, "r");
  
  if (seq_file2_path != NULL){
    *fp2 = gzopen(seq_file2_path, "r");
  }
  
  if (output_merged_path != NULL){
    *fp3 = fopen(output_merged_path, "w");

    if (*fp3 == NULL){
      fprintf(stderr, "Can not open file %s, exiting.\n", output_merged_path);
      exit(EXIT_FAILURE);
    }

  }
  else{
    *fp3 = stdout;
  }

  if (output_path != NULL){
    *fp4 = initialize_output_file(output_path);
  }
  

}


int main(int argc, char **argv)
{
	/*
	Summary: 
	Merges paired reads from two files. 
	If trim=1 only overlapping parts are reported back. 
	If output_merged_path=NULL result is reported to stdout.
	If output_path is not NULL statistics about the merge is reported to here
	If merged_sequence quality is less than mask_threshold, the merged nt is
		returned as a lower-case letter (which is a guaranteed mismatch for most
		common alignment tools in downstream analysis)
	*/
	
  clock_t begin = clock();
  double time_spent = 0;

  char *seq_file1_path = NULL;
  char *seq_file2_path = NULL;
  char *output_merged_path = NULL;
  char *output_path = NULL;
  char *qtype = "Phred+33";
  
  int cigar_size = 1001;
  char *cigar = (char*)safe_malloc(sizeof(char)*cigar_size, "cigar", 
                                   NULL, NULL, -1);
  cigar[cigar_size-1]='\0';
  
  int cigar_stats[6] = {0, 0, 0, 0, 0, 0};
  struct Seq rc_seq_buffer;
  struct Seq extended_seqs[2];
  struct Seq merged_seq;
  initialize_seq(&rc_seq_buffer, 1000);
  initialize_seq(&extended_seqs[0], 1000);
  initialize_seq(&extended_seqs[1], 1000);
  initialize_seq(&merged_seq, 1000);

	int N0 = 0;
  int N1 = -1;
  int counter = 0;
  int l1, l2;
  int begin_free1=0;
  int begin_free2=0;
  int end_free1=0;
  int end_free2=0;
  int trim=1;
  int mask_threshold=30;
	int merge_debug=0;
	int min_len=0;
	int min_ends_free=0;
	
  gzFile fp1, fp2;
  FILE *fp3 = NULL;
  FILE *fp4 = NULL;
  
  kseq_t *seq_record1;
  kseq_t *seq_record2;
  wavefront_aligner_t* wf_aligner;
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine_2p;
  attributes.affine2p_penalties.mismatch = 1;
  attributes.affine2p_penalties.match = 0;
  attributes.affine2p_penalties.gap_opening1 = 6;
  attributes.affine2p_penalties.gap_extension1 = 4;

  parse_args(argc, argv, &seq_file1_path, &seq_file2_path, &output_merged_path,
           &output_path, &attributes.affine2p_penalties.mismatch, &attributes.affine2p_penalties.match,
           &attributes.affine2p_penalties.gap_opening1, &attributes.affine2p_penalties.gap_extension1, 
           &begin_free1, &begin_free2, &end_free1, &end_free2, &N0, &N1, &qtype, &trim,
           &mask_threshold, &merge_debug);

  initialize_files(&fp1, &fp2, &fp3, &fp4, seq_file1_path, seq_file2_path,
                   output_merged_path, output_path); 
 
  seq_record1 = kseq_init(fp1);
  seq_record2 = kseq_init(fp2);

  attributes.alignment_form.span = alignment_endsfree;
  attributes.affine2p_penalties.gap_opening2 = attributes.affine2p_penalties.gap_opening1;
  attributes.affine2p_penalties.gap_extension2 = attributes.affine2p_penalties.gap_extension1;
  wf_aligner = wavefront_aligner_new(&attributes);

  while ((l1 = kseq_read(seq_record1)) >= 0 && (counter<N1 || N1 == -1))
  {
  	
		l2 = kseq_read(seq_record2);
		min_len = min(seq_record1->seq.l, seq_record2->seq.l);
		
		if (counter<N0){
			counter += 1;
  		continue;
  	}
  
		check_read(seq_record1->name.s, seq_record1->seq.s, seq_record1->qual.s, 
							 seq_record1->seq.l, seq_record1->qual.l);
		check_read(seq_record2->name.s, seq_record2->seq.s, seq_record2->qual.s, 
							 seq_record2->seq.l, seq_record2->qual.l);
		
		if (min_len < begin_free1 || min_len < begin_free2 || min_len < end_free1
				|| min_len < end_free2){
				
			min_ends_free = (int)((float)min_len/2);
			
			wavefront_aligner_set_alignment_free_ends(wf_aligner, 
																								min(begin_free1, min_ends_free), 
																								min(begin_free2, min_ends_free), 
																								min(end_free1, min_ends_free), 
																								min(end_free2, min_ends_free));
				
		}
		else{
			wavefront_aligner_set_alignment_free_ends(wf_aligner, begin_free1, begin_free2, 
                                            end_free1, end_free2);
		}
		
		
    
    if (l2<0) break;
    
    assertm(seq_record2->seq.l==seq_record2->qual.l, __FILE__, __FUNCTION__, __LINE__,
         "%s R2: seq length is %d but qual length is %d\n", seq_record2->name.s, 
         seq_record2->seq.l, seq_record2->qual.l);
    
    assertm(strcmp(seq_record1->name.s,seq_record2->name.s)==0, __FILE__, __FUNCTION__, __LINE__,
           "seq R1 name is %s but seq R2 name is %s\n", seq_record1->name.s, 
           seq_record2->name.s);

    align_pairs(seq_record1, seq_record2, &rc_seq_buffer, &cigar, 
                &cigar_size, wf_aligner, extended_seqs);
    
    get_cigar_stats(cigar, 1, 0, begin_free1, begin_free2, end_free1, end_free2,
                    cigar_stats);
                    
    merge_seqs(extended_seqs[0].s, extended_seqs[0].q, extended_seqs[1].s, 
               extended_seqs[1].q, qtype, &merged_seq, mask_threshold);
               
    if (merge_debug==1){
		  fprintf(stderr, "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n+\n", seq_record1->name.s,
		  				extended_seqs[0].s, cigar, extended_seqs[1].s, merged_seq.s, 
		  				extended_seqs[0].q, extended_seqs[1].q, merged_seq.q);
	 	}
	 
	 
    fputseq(fp3, *seq_record1, merged_seq, (trim==1)?cigar_stats[0]:0, 
        (trim==1)?(merged_seq.l-cigar_stats[1]):merged_seq.l);

    if (fp4 != NULL){
        add_alignment_result_to_file(fp4, cigar_stats, seq_record1);
    }
                                                          
    counter +=1;
    
  }
  
  kseq_destroy(seq_record1);
  kseq_destroy(seq_record2);
  gzclose(fp1);
  
  if (seq_file2_path != NULL){
    gzclose(fp2);
  }
  
  wavefront_aligner_delete(wf_aligner);
  
  if (fp3 != NULL){
    fclose(fp3);
  }
  
  if (fp4 != NULL){
  	fclose(fp4);
  }

  if (cigar != NULL){
    free(cigar);
  }
  
  free_seq(&merged_seq);
  free_seq(&extended_seqs[0]);
  free_seq(&extended_seqs[1]);
  free_seq(&rc_seq_buffer);
 
  clock_t end = clock();
  time_spent += (double)(end - begin) / CLOCKS_PER_SEC;

  fprintf(stderr, "%d sequences in T = %f seconds.\n", counter-N0, time_spent);
  
  return 0;
}
