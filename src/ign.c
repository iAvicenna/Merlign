#include "../include/alignment.h"
#include "../include/io.h"
#include "../include/structures.h"
#include "../include/utils.h"
#include "../include/error.h"
#include <stdlib.h>
#include <getopt.h>
#define BUFFERSIZE 10000

void parse_args(int argc, char **argv, char **reads_path, char **output_path, 
                char **references_path, char **mids_file_path, int *mid_start,
                int *mid_buffer, int *match, int *mismatch, int *gap_opening, int *gap_extension,
                int *N0, int *N1, int *mid_stop_when, int *mid_max_distance,
                float *reference_stop_when, int *read_begin_free, int *read_end_free,
                int *ref_begin_free, int *ref_end_free, char **qtype, int *dif_threshold,
                int *realign_qual_threshold, int *realign_error_threshold,
                int *alignment_debug)
{
  int c;
  
  while (1)
    {
      static struct option long_options[] =
        {
          {"reads",     required_argument, 0, 'F'},
          {"output", required_argument, 0, 'G'},
          {"references", required_argument, 0, 'H'},
          {"mids",    required_argument, 0, 'A'},
          {"mid_start", required_argument, 0, 'S'},
          {"mid_buffer", required_argument, 0, 'E'},
          {"match", required_argument, 0, 'q'},
          {"mismatch", required_argument, 0, 'w'},
          {"gap_open", required_argument, 0, 'e'},
          {"gap_extend", required_argument, 0, 'r'},
          {"read_begin_free", required_argument, 0, 't'}, 
          {"read_end_free", required_argument, 0, 'y'},
          {"ref_begin_free", required_argument, 0, 'u'},
          {"ref_end_free", required_argument, 0, 'i'},
          {"N0", required_argument, 0, 'B'},
          {"N1", required_argument, 0, 'N'},
          {"mid_stop_when", required_argument, 0, 'W'},
          {"mid_max_distance", required_argument, 0, 'M'},
          {"reference_stop_when", required_argument, 0, 'R'},
          {"qtype", required_argument, 0, 'Q'},
          {"dif_threshold", required_argument, 0, 'T'},
          {"realign_qual_threshold", required_argument, 0, 'V'},
          {"realign_error_threshold", required_argument, 0, 'C'},
          {"alignment_debug", required_argument, 0, 'D'},

          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "F:G:H:A:S:E:q:w:e:r:t:y:u:i:B:N:W:M:R:Q:T:V:C:D:",
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
          fprintf (stderr, "option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;
        case 'F':
          *reads_path = optarg;
          break;
        case 'G':
          *output_path = optarg;
          break;
        case 'H':
          *references_path = optarg;
          break;
        case 'B':
          *N0 = atoi(optarg);
          break;
        case 'N':
          *N1 = atoi(optarg);
          break;
        case 'S':
          *mid_start = atoi(optarg);
          break;
        case 'E':
          *mid_buffer = atoi(optarg);
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
          *read_begin_free = atoi(optarg);
          break;
        case 'y':
          *read_end_free = atoi(optarg);
          break;
        case 'u':
          *ref_begin_free = atoi(optarg);
          break;
        case 'i':
          *ref_end_free = atoi(optarg);
          break;
        case 'A':
          *mids_file_path = optarg;
          break;
        case 'W':
          *mid_stop_when = atoi(optarg);
          break;
        case 'M':
          *mid_max_distance = atoi(optarg);
          break;
        case 'R':
          *reference_stop_when = atof(optarg);
          break;
        case 'Q':
          *qtype = optarg;
          break;
        case 'T':
          *dif_threshold = atoi(optarg);
          break;
        case 'V':
        	*realign_qual_threshold = atoi(optarg);
        	break;
        case 'C':
      	  *realign_error_threshold = atoi(optarg);
      	  break;
      	case 'D':
      		*alignment_debug = atoi(optarg);
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

  char *str_time;
  time_t rawtime;
  struct tm * timeinfo;
  
  check_parameters_ign(*match, *mismatch, *gap_opening, *gap_extension, *read_end_free,
                       *ref_begin_free, *read_begin_free, *ref_begin_free, *mid_start, 
                       *mid_buffer, *alignment_debug); 
                       
  assertm(*references_path !=NULL, NULL, NULL, -1,  
  			 "Minimal usage is ign --references arg1"
         " (if --reads is not supplied stdin is used if output is not "
         "supplied stout is used)");
  if (*reads_path != NULL) check_file(*reads_path, "--reads");
  check_file(*references_path, "--references");
  if (*mids_file_path != NULL)check_file(*mids_file_path, "--mids");

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  str_time = (char*)calloc(strlen(asctime(timeinfo)), sizeof(char));
  strncpy(str_time, asctime(timeinfo), strlen(asctime(timeinfo))-1);

  fprintf(stderr, "\nStart time: %s (UT: %u).\n", str_time, (unsigned)time(NULL));

  fprintf(stderr, "reads path: %s \nreferences path: %s \noutput path: %s\n"
                  "N0 %d, N1 %d, dif threshold %d, qtype %s, alignment_debug %d\n", 
                  *reads_path, *references_path, *output_path, *N0, *N1, 
                  *dif_threshold, *qtype, *alignment_debug);

  if (*mids_file_path != NULL){
    fprintf(stderr, "mids path: %s\n", *mids_file_path);
  }

  fprintf(stderr, "\npenalties: \n");
  fprintf(stderr, "match=%d, mismatch=%d, ", *match, *mismatch);
  fprintf(stderr, "gap open=%d, gap extension=%d, ", *gap_opening, *gap_extension);
  fprintf(stderr, "read begin/end free=%d/%d ref begin/end free=%d/%d\n\n", 
          *read_begin_free, *read_end_free, *ref_begin_free, *ref_end_free);
  
  fprintf(stderr, "mid_start: %d, mid_end: %d, mid_max_distance: %d, "
                  "mid_stop_when: %d, reference_stop_when %f, \n", *mid_start, 
                  *mid_buffer, *mid_max_distance, *mid_stop_when, *reference_stop_when);  
  
  fprintf(stderr, "dif_threshold %d, realign_qual_threshold %d, realign_error_threshold %d\n\n",
          *dif_threshold, *realign_error_threshold, *realign_qual_threshold);

  free(str_time);

}


FILE *initialize_output_file(char *output_path, int include_mids){

    FILE * fp; 
    
    if (output_path != NULL){
    	fp = fopen(output_path, "w");
			assertm(fp != NULL, __FILE__, __FUNCTION__, __LINE__, 
					 "Can not open file %s supplied as --output", output_path);
		}
		else{
			fp = stdout;
		}

    fprintf(fp, "read_id");

    if (include_mids==1){
      fprintf(fp, ",M1_index,M1_start,M1_end,M1_match,M1_mismatch,M1_deletion,"
                   "M1_insertion,");
      fprintf(fp, "M2_index,M2_start,M2_end,M2_match,M2_mismatch,M2_deletion,"
                   "M2_insertion");

    }
                            
    fprintf(fp, ",align_index,is_rc,align_cigar,align_mods,align_mod_quals,align_av_err\n");
                 
    return fp;

}


void add_insert_alignment_result_to_file(struct InsertAlignmentResult iar, 
                                         FILE* fp)
{

    fprintf(fp, "%d,%d,%s,%s,%s,%f", iar.index, iar.is_rc, iar.condensed_cigar, iar.seq_mods,
                                   iar.qual_mods, iar.aligned_err_prob);
  
}


void add_alignment_result_to_file(struct AlignmentResult ar, FILE* fp,
                                  int *cigar_stats, int put_soft){

  int i,start_index;
  start_index = 0;
  
  if (put_soft==0) start_index=2;
  
  fprintf(fp, "%d,", ar.index);
  fprintf(fp, "%d,", ar.start);
  fprintf(fp, "%d,", ar.end);


  for (i=start_index; i<6; i++){
    if (cigar_stats != NULL){
      fprintf(fp, "%d", cigar_stats[i]);
      }
    else fprintf(fp, "0");
    if (i<5)fprintf(fp, ",");
  }

}


int main(int argc, char **argv){

  char *seq_file_path = NULL;
  char *mids_file_path = NULL;
  char *references_file_path = NULL;
  char *output_path = NULL;
  char buffer[BUFFERSIZE];
  int cigar_stats1[6] = {0, 0, 0, 0, 0, 0};
  int cigar_stats2[6] = {0, 0, 0, 0, 0, 0};
  char *qtype = "Phred+33";
  int buffer_size=101;
  char *cigar_buffer = (char*)safe_malloc(buffer_size*sizeof(char), "cigar_buffer", 
                                          __FILE__, __FUNCTION__, __LINE__);
 
 	cigar_buffer[100]='\0';
 
  int mid_start = 0;
  int mid_buffer = 3;
  int read_begin_free=15;
  int read_end_free=15;
  int ref_begin_free=5;
  int ref_end_free=5;
  int mid_stop_when = 1;
  int mid_max_distance = 2;
  int dif_threshold=50;
  int N0 = 0;
  int N1 = -1;
  int counter, insert_index;
  int realign_dist_threshold = 4;
  int realign_error_threshold = 2;
  int alignment_debug = 0;
  int min_end_free=0;
	float reference_stop_when = 0.9;
  
  double time_spent=0;

  clock_t begin = clock();

  counter = 0; 
  
  struct AlignmentResult mid_alignment_result1;
  struct AlignmentResult mid_alignment_result2;
  struct InsertAlignmentResult iar;
  struct SeqList mids;
  struct SeqList rc_mids;
  struct SeqList references;
  struct Seq extended_seqs[2];
  struct Seq seq;
  
  initialize_seq(&seq, 1000);
  initialize_seq(&extended_seqs[0], 1000);
  initialize_seq(&extended_seqs[1], 1000);  
  initialize_ins_alig_res(&iar, 1000, 1000, 1000);
  initialize_alignment_result(&mid_alignment_result1, 1000);
  initialize_alignment_result(&mid_alignment_result2, 1000);

  gzFile fp1 = NULL;
  gzFile fp2 = NULL;
	FILE *fp3=NULL;

  wavefront_aligner_t* wf_aligner;
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine_2p;
  attributes.affine2p_penalties.mismatch = 2;
  attributes.affine2p_penalties.match = 0;
  attributes.affine2p_penalties.gap_opening1 = 4;
  attributes.affine2p_penalties.gap_extension1 = 2;
 
  wavefront_aligner_t* wf_aligner_precise;
  wavefront_aligner_attr_t attributes_precise = wavefront_aligner_attr_default;
  attributes_precise.distance_metric = gap_affine_2p;
  attributes_precise.affine2p_penalties.mismatch = 2;
  attributes_precise.affine2p_penalties.match = -1;
  attributes_precise.affine2p_penalties.gap_opening1 = 4;
  attributes_precise.affine2p_penalties.gap_extension1 = 2;
  
  parse_args(argc, argv, &seq_file_path, &output_path, &references_file_path, &mids_file_path, &mid_start, &mid_buffer,
             &attributes.affine2p_penalties.match, &attributes.affine2p_penalties.mismatch,
             &attributes.affine2p_penalties.gap_opening1, &attributes.affine2p_penalties.gap_extension1,
             &N0, &N1, &mid_stop_when, &mid_max_distance, &reference_stop_when,
             &read_begin_free, &read_end_free, &ref_begin_free, &ref_end_free,
             &qtype, &dif_threshold, &realign_dist_threshold, &realign_error_threshold,
             &alignment_debug);

  if (seq_file_path != NULL){
    fp1 = gzopen(seq_file_path, "r");
		assertm(fp1 != NULL, __FILE__, __FUNCTION__, __LINE__, 
					 "Can not open file %s supplied as --reads", seq_file_path);
  }
  kseq_t *seq_record;
  seq_record = kseq_init(fp1);
  references = load_sequences(references_file_path, 1);
	
	if (output_path != NULL){
		fp3 = initialize_output_file(output_path, (int)(mids_file_path!=NULL));
	}
	
  if (mids_file_path != NULL){
     mids =  load_sequences(mids_file_path, 0);
     reverse_complement_list(&rc_mids, mids);
  }
  
  attributes.alignment_form.span = alignment_endsfree;
  attributes.affine2p_penalties.gap_opening2 = attributes.affine2p_penalties.gap_opening1;
  attributes.affine2p_penalties.gap_extension2 = attributes.affine2p_penalties.gap_extension1;
  
  attributes_precise.alignment_form.span = alignment_endsfree;
  attributes_precise.affine2p_penalties.gap_opening2 = attributes_precise.affine2p_penalties.gap_opening1;
  attributes_precise.affine2p_penalties.gap_extension2 = attributes_precise.affine2p_penalties.gap_extension1;
 

  wf_aligner = wavefront_aligner_new(&attributes);
  
  wf_aligner_precise = wavefront_aligner_new(&attributes_precise);
  wavefront_aligner_set_alignment_free_ends(wf_aligner_precise, ref_begin_free, 
  																		      ref_end_free, read_begin_free, 
  																		      read_end_free);    


  while (get_seq(&fp1, buffer,  &seq, seq_record, alignment_debug)==1){
  
  		if (counter<N0){
  			counter += 1;
  			continue;
  		}
  
      if (fp3 != NULL)fprintf(fp3, "%s,", seq.id);
      
      check_read(seq.id, seq.s, seq.q, seq.l, 
      					 (seq.q == NULL)?-1:(int)strlen(seq.q));
      
      

      if (seq.l<read_begin_free || seq.l<read_end_free){
      		
		    min_end_free = (int)((float)seq.l/2);
				wavefront_aligner_set_alignment_free_ends(wf_aligner, ref_begin_free,
																									ref_end_free, min_end_free, 
																									min_end_free);
				
				wavefront_aligner_set_alignment_free_ends(wf_aligner_precise, 
																									ref_begin_free,
																									ref_end_free, min_end_free, 
																									min_end_free);      		
      }
      else{
      	wavefront_aligner_set_alignment_free_ends(wf_aligner, ref_begin_free, 
      																						ref_end_free, read_begin_free, 
      																						read_end_free);   
      																						 
      	wavefront_aligner_set_alignment_free_ends(wf_aligner_precise, 
      																						ref_begin_free, ref_end_free, 
      																						read_begin_free, 
      																						read_end_free);   
      }
      
		  
      align_long(seq.s, seq.q, references,  wf_aligner, 
                 reference_stop_when, qtype, dif_threshold,
                 extended_seqs, &iar, &cigar_buffer, &buffer_size);
      

      if (realign_dist_threshold>0 && realign_error_threshold>0 && 
          attributes.affine2p_penalties.match == 0 &&
          iar.distance > realign_dist_threshold && 
          iar.aligned_err_prob < realign_error_threshold)        //realign more precisely but more slowly
      {
         align_long(seq.s, seq.q, references,  wf_aligner_precise, 
                    reference_stop_when, qtype, dif_threshold,
                    extended_seqs, &iar, &cigar_buffer, &buffer_size);
                 
      }
      
      if (alignment_debug==1){
		    fprintf(stderr, "%s\n%s\n%s\n%s\n%d\n+\n", extended_seqs[0].s,
		    				extended_seqs[1].s, extended_seqs[0].q,
		    				iar.condensed_cigar, iar.distance);
			}
      
      insert_index = iar.index;
      
      if (mids_file_path != NULL){

        align_short(seq.s, mids, mid_max_distance, mid_stop_when, 
                    "start", mid_start, mid_buffer, &mid_alignment_result1);
                                            
        align_short(seq.s, rc_mids, mid_max_distance, mid_stop_when, 
                    "end", seq.l - mid_start - mid_buffer, mid_buffer,
                    &mid_alignment_result2);
        
        get_cigar_stats(mid_alignment_result1.cigar, 0, 1, 0, 0, 'S', 'S', cigar_stats1);
        get_cigar_stats(mid_alignment_result2.cigar, 0, 1, 0, 0, 'S', 'S', cigar_stats2);  
        
        if (fp3 != NULL && insert_index%2==0){ // insert is 5 prime to 3 prime
          add_alignment_result_to_file(mid_alignment_result1, fp3, cigar_stats1, 0);     
          fprintf(fp3, ","); 
          add_alignment_result_to_file(mid_alignment_result2, fp3, cigar_stats2, 0);  
          fprintf(fp3, ",");    
        }
        else if (fp3 != NULL){ //insert is reverse_complement (insert is 3 prime to 5 prime)
          add_alignment_result_to_file(mid_alignment_result2, fp3, cigar_stats2, 0);     
          fprintf(fp3, ","); 
          add_alignment_result_to_file(mid_alignment_result1, fp3, cigar_stats1, 0);  
          fprintf(fp3, ",");    
        }

        
      }
      
      if (fp3 != NULL){

        if (insert_index%2==1){
          iar.is_rc = 1;  
        }
        else {iar.is_rc=0;}

			  iar.index = (int)((double)(iar.index - iar.index%2)/2); /*since reverse complements are also considered there are
                                                                  2x the original number of references with 2n+1th 
                                                                  being the rc of 2nth sequence*/
                                                                  
        if (strcmp(iar.seq_mods,"OTHER")==0)
        {
        	iar.index = -1;
        	iar.is_rc = -1;
        }                                                        
                                                                  
        add_insert_alignment_result_to_file(iar, fp3);
				fprintf(fp3, "\n");                                            
			}
			
      counter += 1;
      if (N1>0 && counter>=N1){
        break;
      }
      iar.distance = -1;
      
      		 
      
  
  }

  wavefront_aligner_delete(wf_aligner);                                          
  wavefront_aligner_delete(wf_aligner_precise);                                          
   
  kseq_destroy(seq_record);
  if (fp1 != NULL) gzclose(fp1);
  gzclose(fp2);
	
	if (fp3 != NULL) fclose(fp3);  
  
  if (mids_file_path != NULL){
    free_seq_list(&mids);
    free_seq_list(&rc_mids);
  }
  
  free_seq_list(&references);
  free_seq(&extended_seqs[0]);
  free_seq(&extended_seqs[1]);
  free_ins_alig_res(&iar);
  free(cigar_buffer);
  free_seq(&seq);
	free_alignment_result(&mid_alignment_result1);
	free_alignment_result(&mid_alignment_result2);

  clock_t end = clock();
  time_spent += (double)(end - begin) / CLOCKS_PER_SEC;

  fprintf(stderr, "%d sequences in T = %f seconds.\n", counter-N0, time_spent);
  
}
