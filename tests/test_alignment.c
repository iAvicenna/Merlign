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

START_TEST(test_get_start_end){
  
  int start,end,i;
	char *extended_refs[4] = {"    ", "  A  ", "ABCDEF ", " ABCDE"};
	int lengths[4] = {4,5,7,6};
	int starts[4] = {4,2,0,1};
	int ends[4] = {-1,2,5,5};
	
	
	for (i=0; i<4; i++){
		start = get_start(extended_refs[i], lengths[i]);
		end = get_end(extended_refs[i], lengths[i]);
		
		ck_assert_int_eq(start, starts[i]);
		ck_assert_int_eq(end, ends[i]);
	}
	
	
}END_TEST


START_TEST(test_get_insert_alignment_result){
  
  char *qtype = "Phred+33";
  struct InsertAlignmentResult iar;
	initialize_ins_alig_res(&iar, 0, 0, 0);
  int thresholds[14] = {10,10,10,10,10,10,10,10,10,4, 8, 6, 50, 50};
  int min_lens[14] =   {17,14,14,17,14,14,20,8, 8, 20,20,16,233,233};
  int i;
  
  char *extended_reads[14] = {
    "AGTCGTAGCTAGCATCA",
    "   CGTAGCTAGCATCA",
    "   CGTAGCTAGCATCA",
    "AGTCGTAGCTAGCATCA",
    "   CGTAGCTAGGATCA",
    "   CGTAGCTAGGATCA",
    "AGC   CGTAGCTAGGATCAAGC",
    "AGCTGAAG",
    "AGCTGAAG",
    "AGC   CGTAGCTAGGATCAAGC", 
    "AGC   CGTAGCTAGGATCAAGC",
    "   CGTAGCTAGGATCAA A",
    "CGACTACAGGACTCTCCCT ATGCTCGAACATATAGACTGGTTTGGTCACTGTCCGTGCTCGGGTGGTGAACCCCCCAAATGTACAATTTGTCAAATTTGCCATTGTTTGGCATAGTCGCGTTCAGCGCTGGATATTTGTATTCTAATTTGTGCAACCAATTCAATCTACTAAAGAAACTAACAACAGATCCCCTTTTGCAAGCATAGCTTTTCCCATCCTGAGCCTGTAGTCG",
    "CGACTACAGGACTCTCCCT ATGCTCGAACATATAGACTGGTTTGGTCACTGTCCGTGCTCGGGTGGTGAACCCCCCAAATGTACAATTTGTCAAATTTGCCATTGTTTGGCATAGTCGCGTTCAGCGCTGGATATTTGTATTCTAATTTGTGCAACCAATTCAATCTACTAAAGAAACTAACAACAGATCCCCTTTTGCAAGCATAGCTTTTCCCATCCTGAGCCTGTAGTCG"

  };
  
  
  
  char *extended_quals[14] = {
    "ABCDEFGHIJKLMNOPR",
    "   DEFGHIJKLMNOPR",
    "   \"EFGHIJKLMNOPR",
    "ABCDEFGHIJKLMNOPR",
    "   DEFGHIJKLMNOPR",
    "   DEFGHIJKLMNOPR",
    "123   DEFGHIJKLMNOPRPTY",
    "ABCDEFGH",
    "ABCDEFGH",
    "123   DEFGHIJKLMNOPRPTY",
    "123   DEFGHIJKLMNOPRPTY",
    "   DEFGHIJKLMNOPQR S",
    "PjjjjjjjjjjPjjjPjjj jjj\\jjjjPjjjjjjjjjjPjjPjjjjjjj\\jjjjjjjjj\\\\jjjjjjjjjjjjj\\jjjj\\jjjjjjjjjjjjj\\jjjjjjjjjjjjjPjjjjjjjjjjjjjjjB:jjPjjj\\jP\\jjj\\jjjPjjjjjjjjjPPj\\jj\\:jjPjjPj\\\\\\N\\jjjjj\\jjN\\\\\\:\\P\\\\PN\\\\\\\\\\\\\\P\\\\P\\\\\\\\\\PN\\\\j\\\\\\\\j\\\\\\P\\\\\\\\\\\\j\\\\\\jN",
    "PjjjjjjjjjjPjjjPjjj jjj\\jjjjPjjjjjjjjjjPjjPjjjjjjj\\jjjjjjjjj\\\\jjjjjjjjjjjjj\\jjjj\\jjjjjjjjjjjjj\\jjjjjjjjjjjjjPjjjjjjjjjjjjjjjB:jjPjjj\\jP\\jjj\\jjjPjjjjjjjjjPPj\\jj\\:jjPjjPj\\\\\\N\\jjjjj\\jjN\\\\\\:\\P\\\\PN\\\\\\\\\\\\\\P\\\\P\\\\\\\\\\PN\\\\j\\\\\\\\j\\\\\\P\\\\\\\\\\\\j\\\\\\jN"

  };

  char *extended_references[14] = {
    "AGTCGTAGCTAGCATCA",
    "AGTCGTAGCTAGCATCA",
    "AGT GTAGCTAGCATCA",
    "AGTCGTAGTTAGCATCA",
    "AGTCGTAGCTAGTATCA",
    "AGTCG AGCTAGTATCA",
    "   AGTCG AGCTAGTATCA   ",
    "     G  ",
    "     A  ",
    "   AGTCG AGCTAGTATCA   ", 
    "   AGTCG AGCTAGTATCA   ",
    "AGTCGTAGCTAGTATCA A ",
    "         GACTCTCCCTGATGCTCGAACATATATGCTGGTTTGGTCACTGTCCGTGCTCGGGTGGTGAACCCCCCAAATGTACAATTTGTCAAATTTGCCATTGTTTGGCATAGTCACGTTCAGCGCTGGATATTTGTATTCTAATTTGTGCAACCAATTCAATCTACTAAAGAAACTGTTAACAGATCCCCTTTTGCAAGCATAGCTTTTCCCATCCTGAGCXXXXXXXXX",
    "         GACTCTCCCTGATGCTCGAACATATATGCTGGTTTGGTCACTGTCCGTGCTCGGGTGGTGAACCCCCCAAATGTACAATTTGTCAAATTTGCCATTGTTTGGCATAGTCACGTTCAGCGCTGGATATTTGTATTCTAATTTGTGCAACCAATTCAATCTACTAAAGAAACTGTTAACAGATCCCCTTTTGCAAGCATAGCTTTTCCCATCCTGAGC         "

    
};
  
  char *cigar_arrays[14] = {
    "MMMMMMMMMMMMMMMMM",
    "IIIMMMMMMMMMMMMMM",
    "IIIDMMMMMMMMMMMMM",
    "MMMMMMMMXMMMMMMMM",
    "IIIMMMMMMMMMXMMMM",
    "IIIMMDMMMMMMXMMMM",
    "DDDIIIMMDMMMMMMXMMMMDDD",
    "DDDDDXDD",
    "DDDDDMDD",
    "DDDIIIMMDMMMMMMXMMMMDDD", 
    "DDDIIIMMDMMMMMMXMMMMDDD",
    "IIIMMMMMMMMMXMMMMDID",
    "DDDDDDDDDMMMMMMMMMMDMMMMMMMMMMMMMMMXXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMXXXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMXXXXXXXXX",  
    "DDDDDDDDDMMMMMMMMMMDMMMMMMMMMMMMMMMXXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMXXXMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDDDDDDDDD"    
  };
  
  char *iars0[14] = {
    "",
    "   ",
    "   C",
    "C",
    "   G",
    "   TG",
    "AGC   TGAGC",
    "AGCTGAAG",
    "AGCTGAG",
    "OTHER",
    "AGC   TGAGC",
    "   GA A",
    "CGACTACAG GAGAACCTGTAGTCG",
    "CGACTACAG GAGAACCTGTAGTCG"

  };
  
  char *iars2[14] = {
    "17M",
    "3I14M",
    "3I1D13M",
    "8M1X8M",
    "3I9M1X4M",
    "3I2M1D6M1X4M",
    "3D3I2M1D6M1X4M3D",
    "5D1X2D",
    "5D1M2D",
    "NA",
    "3D3I2M1D6M1X4M3D",
    "3I9M1X4M1D1I1D",
    "9D10M1D15M2X81M1X61M3X42M9X",
    "9D10M1D15M2X81M1X61M3X42M9D"  
  };
  
  //do cases with insertions at the beginning and endto make sure
  //start and end are handled correctly. check that threshold is 
  //is handled correctly. 
  
  char *iars1[14] = {
    "\"\"",
    "\"   \"",
    "\"   \"\"\"",
    "\"I\"",
    "\"   M\"",
    "\"   FM\"",
    "\"123   FMPTY\"",
    "\"ABCDEFGH\"",
    "\"ABCDEGH\"",
    "NA",
    "\"123   FMPTY\"",
    "\"   MR S\"",
    "\"Pjjjjjjjj jjjjN\\\\\\\\j\\\\\\jN\"",
		"\"Pjjjjjjjj jjjjN\\\\\\\\j\\\\\\jN\""

  };
  

  for (i=0; i<14; i++){

    get_insert_alignment_result(extended_reads[i], extended_references[i], 
												 extended_quals[i], cigar_arrays[i], qtype, 
												 thresholds[i], &iar, min_lens[i]);

    ck_assert_str_eq(iar.seq_mods, iars0[i]);
    ck_assert_str_eq(iar.qual_mods, iars1[i]);
    ck_assert_str_eq(iar.condensed_cigar, iars2[i]);

    
  }
  
  free_ins_alig_res(&iar);
}
END_TEST



START_TEST(test_align_pairs){

	int l1,l2,counter=0;

  struct Seq rc_seq_buffer;
  struct Seq extended_seqs[2];
  
  initialize_seq(&rc_seq_buffer, 1000);
  initialize_seq(&extended_seqs[0], 1000);
  initialize_seq(&extended_seqs[1], 1000);

	char *cigar_buffer = (char*)safe_malloc(11*sizeof(char), "cigar_buffer", 
																					__FILE__, __FUNCTION__, __LINE__);
	int buffer_size=11;
 
  kseq_t *seq_record1;
  kseq_t *seq_record2;
  
  char *path1 = "./test_files/simulated_reads_R1.fastq";
  char *path2 = "./test_files/simulated_reads_R2.fastq";
  gzFile fp1, fp2;

  fp1 = gzopen(path1, "r");
  fp2 = gzopen(path2, "r");
	
  seq_record1 = kseq_init(fp1);
  seq_record2 = kseq_init(fp2);


	wavefront_aligner_t* wf_aligner;
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine_2p;
  attributes.affine2p_penalties.mismatch = 1;
  attributes.affine2p_penalties.match = 0;
  attributes.affine2p_penalties.gap_opening1 = 4;
  attributes.affine2p_penalties.gap_extension1 = 2;
	wf_aligner = wavefront_aligner_new(&attributes);
  wavefront_aligner_set_alignment_free_ends(wf_aligner, 0, 0, 0, 0);
  
  while ((l1 = kseq_read(seq_record1)) >= 0 && 
         (l2 = kseq_read(seq_record2)) >= 0 )
  {
      counter += 1;
      
      align_pairs(seq_record1, seq_record2, &rc_seq_buffer, &cigar_buffer, 
                  &buffer_size, wf_aligner, extended_seqs); 

      ck_assert(strcmp(seq_record1->seq.s, extended_seqs[0].s)
                || strcmp(seq_record2->seq.s, extended_seqs[1].s)); 

  }

  ck_assert_int_le(0, counter);

	wavefront_aligner_delete(wf_aligner);
	free(cigar_buffer);
	free_seq(&extended_seqs[0]);
	free_seq(&extended_seqs[1]);
	free_seq(&rc_seq_buffer);
  kseq_destroy(seq_record1);
  kseq_destroy(seq_record2);
  gzclose(fp1);
  gzclose(fp2);
               	
}


START_TEST(test_align_long_simple){

	int i,j;
	int dif_threshold = 20;
	int stop_when = 1;
	char *cigar_buffer = (char*)safe_malloc(11*sizeof(char), "cigar_buffer", 
																					__FILE__, __FUNCTION__, __LINE__);
	int buffer_size=11;
	
	int index[3] = {0, 2, 4};
	char *cigars[3] = {
		"46M",
		"71M",
		"37M"
	};
	int lengths[3] = {0, 2, 3};

	
	char *qtype = "Phred+33";
	
	char *quals[3]=
	{
		"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
		"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
	};
	
	char *seqs[3] = {
		"AGATCGGAAGAGGTCGACTACGATCGATCGATCGATCGCATCGCAT",
		"TGAGCATCGACGATCGCGGCATCGATCGTACGTACGTACGTCGATGTGTCGACTATCGCATGCACGCATCA",
		"AGGTATCATCGAGCGTACGATGCATCGATCTATCGAT"
	};
	
	struct SeqList references;
	struct Seq extended_seqs[2];
	struct InsertAlignmentResult iar;
	initialize_ins_alig_res(&iar, 1, 1, 1);
	
	initialize_seq(&extended_seqs[0], 10);
	initialize_seq(&extended_seqs[1], 10);

	references = load_sequences("./test_files/align_long_test_references.fasta", 1);

	wavefront_aligner_t* wf_aligner;
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine_2p;
  attributes.affine2p_penalties.mismatch = 1;
  attributes.affine2p_penalties.match = 0;
  attributes.affine2p_penalties.gap_opening1 = 4;
  attributes.affine2p_penalties.gap_extension1 = 2;
	wf_aligner = wavefront_aligner_new(&attributes);
  wavefront_aligner_set_alignment_free_ends(wf_aligner, 0, 0, 0, 0);
  
	
	for (i=0; i<3; i++){
		
			align_long(seqs[i], quals[i], references, wf_aligner, stop_when, qtype, 
								 dif_threshold, extended_seqs, &iar, &cigar_buffer,
								 &buffer_size);
			
			ck_assert_int_eq(iar.index, index[i]);
			ck_assert_str_eq(iar.condensed_cigar, cigars[i]);
			ck_assert_str_eq(iar.seq_mods, "");
			ck_assert_str_eq(iar.qual_mods, "\"\"");
			
			for (j=0; j<3; j++){
				ck_assert_int_eq(iar.l[j], lengths[j]);
			}

      iar.distance = -1;
			
			
	}
	
	wavefront_aligner_delete(wf_aligner);
	free_ins_alig_res(&iar);
	free(cigar_buffer);
	free_seq_list(&references);
	free_seq(&extended_seqs[0]);
	free_seq(&extended_seqs[1]);
	
}


START_TEST(test_align_long_complex){

	int i,j;
	int dif_threshold = 20;
	int stop_when = 1;
	char *cigar_buffer = (char*)safe_malloc(101*sizeof(char), "cigar_buffer", 
																					__FILE__, __FUNCTION__, __LINE__);
	int buffer_size=101;
	
	int index[3] = {0, 2, 4}; /*0 2 4 because for each reference its reverse is also added to the list*/
	char *cigars[3] = {
		"4M1X40M1X",
		"5M1I65M1D",
		"9M3D21M3I4M"
	};
	char *seq_mods[3] = {
		"AC",
		"C ",
		"   AGC"
	};
	char *qual_mods[3] = {
		"\"ED\"",
		"\"F \"",
		"\"   JKL\""
	};
	int lengths[3][3] = {{2, 4, 9}, {2, 4, 9}, {6, 8, 11}};

	
	char *qtype = "Phred+33";
	
	char *quals[3]=
	{
		"ABCDEFGHIJKLMNOPRSQTUVYZ12345567890JDHALJDAKLD",
		"ABCDEFGHIJKLMNOPRSQTUVYZ12345567890JDHALJDAKLDLCNIOQEUIOQEUEIUNCNZMCNAS",
		"ABCDEFGHIMNOPRSQTUVYZ123455678JKL90JD"
	};
	
	char *seqs[3] = {
		"AGATAGGAAGAGGTCGACTACGATCGATCGATCGATCGCATCGCAC",
		"TGAGCCATCGACGATCGCGGCATCGATCGTACGTACGTACGTCGATGTGTCGACTATCGCATGCACGCATC",
		"AGGTATCATGCGTACGATGCATCGATCTATAGCCGAT"
	};
	
	struct SeqList references;
	struct Seq extended_seqs[2];
	struct InsertAlignmentResult iar;
	initialize_ins_alig_res(&iar, 1, 1, 1);
	
	initialize_seq(&extended_seqs[0], 1);
	initialize_seq(&extended_seqs[1], 1);

	references = load_sequences("./test_files/align_long_test_references.fasta", 1);

	wavefront_aligner_t* wf_aligner;
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine_2p;
  attributes.affine2p_penalties.mismatch = 1;
  attributes.affine2p_penalties.match = 0;
  attributes.affine2p_penalties.gap_opening1 = 4;
  attributes.affine2p_penalties.gap_extension1 = 2;
	wf_aligner = wavefront_aligner_new(&attributes);
  wavefront_aligner_set_alignment_free_ends(wf_aligner, 0, 0, 0, 0);
  

	for (i=0; i<3; i++){
		
			align_long(seqs[i], quals[i], references, wf_aligner, stop_when, qtype, 
								 dif_threshold, extended_seqs, &iar, &cigar_buffer,
								 &buffer_size);
			ck_assert_int_eq(iar.index, index[i]);
			ck_assert_str_eq(iar.condensed_cigar, cigars[i]);
			ck_assert_str_eq(iar.seq_mods, seq_mods[i]);
			ck_assert_str_eq(iar.qual_mods, qual_mods[i]);
			
			for (j=0; j<3; j++){
				ck_assert_int_eq(iar.l[j], lengths[i][j]);
			}

      iar.distance = -1;
			
	}
	
	wavefront_aligner_delete(wf_aligner);
	free_ins_alig_res(&iar);
	free(cigar_buffer);
	free_seq_list(&references);
	free_seq(&extended_seqs[0]);
	free_seq(&extended_seqs[1]);
	
}


START_TEST(test_align_short_simple){
	
	int i;
	char *adapter_seqs[] = {"AGCTGTGCTC", "GCTAGTGTCGTC", "TGAGCTGCTGTA"};
	struct AlignmentResult ar;
	

	initialize_alignment_result(&ar, 1000);
	struct SeqList adapters;
	initialize_seq_list(&adapters, 3);
	for (i=0; i<3; i++){
		add_seq_to_seqlist(&adapters, i, adapter_seqs[i], NULL, NULL);
	}

	char *reads[] = {"TCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCATAGCTGTGCTC",
									 "GCTAGTGTCGTCTCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCAT",
									 "TCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCATTGAGCTGCTGTA"};
	int starts[] = {43, 0, 43};
	int buffers[]   = {2, 2, 2};
	char *cigars[] = {"10=", "12=", "12="};
	for (i=0; i<3; i++){
  
		align_short(reads[i], adapters, 1, 0, "start", starts[i], buffers[i], &ar);
		ck_assert_int_eq(starts[i], ar.start);
		ck_assert_int_eq(i, ar.index);
		ck_assert_int_eq(starts[i]+strlen(adapter_seqs[ar.index])-1, ar.end);
    ck_assert_int_eq(0, ar.distance);
		ck_assert_str_eq(cigars[i], ar.cigar);
		
	}
	free_alignment_result(&ar);
	free_seq_list(&adapters);
									 

}
END_TEST


START_TEST(test_align_short_complex){

	/*
	Tests of aligning adapters where the max edit distance allowed is 1
	All reads except fourth have the adapter with either one deletion, insertion
	or mismatch at various locations including extremities so that is checked for.
	It is also checked that the fourth read which has two errors returns no results.
	*/
	
	struct SeqList adapters;
	initialize_seq_list(&adapters, 3);
	
	int i;
	char *adapter_seqs[] = {"AGCTGTGCTC", "GCTAGTGTCGTC", "TGAGCTGCTGTA"};
	struct AlignmentResult ar;
	initialize_alignment_result(&ar, 1000);
	for (i=0; i<3; i++){
		add_seq_to_seqlist(&adapters, i, adapter_seqs[i], NULL, NULL);
	}

	char *reads[] = {"TCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCATAGCTGTCTC",
									 "GTAGTGTCGTCTCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCAT",
									 "TCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCATTGAGCTGCGTA",
									 "TCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCATAGCTGGTC",
									 "TCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCATAGCTGGGCTC",
									 "TCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCATTGAAGCTGCTGTA",
									 "CTAGTGTCGTCTCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCAT",
									 "TCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCATAGCTGTGCT"};
									 
									 //do more cases where sequence has two adapters with 1 distance
									 
									 
	int starts[] = {40, 0, 40, 0, 39, 35, 0, 39};
  int buffers[] = {2, 2, 2, 2, 5, 9, 2, 5};

	int distances[] = {1, 1, 1, strlen(reads[3]), 1, 1, 1, 1};
	int indices[] = {0, 1, 2, -1, 0, 2, 1, 0};
	char *cigars[] = {"6=1D3=","1=1D10=", "8=1D3=", "", "5=1X4=", "3=1I9=",
										"1D11=", "9=1D"};

  int start_answers[] = {43, 0, 43, -1, 43, 43, 0, 43};
  int end_answers[] = {51, 10, 53, -1, 52, 55, 10, 51};
	
	for (i=0; i<8; i++){
		align_short(reads[i], adapters, 1, 0, "start", starts[i], buffers[i], &ar);
		
		ck_assert_int_eq(start_answers[i], ar.start);
		ck_assert_int_eq(end_answers[i], ar.end);
		ck_assert_int_eq(indices[i], ar.index);
		ck_assert_int_eq(distances[i], ar.distance);
		
		if (ar.cigar == NULL){
			ck_assert(cigars[i] == NULL);
		}
		else{
			ck_assert_str_eq(cigars[i], ar.cigar);
		}
		
	}
  free_alignment_result(&ar);
	free_seq_list(&adapters);
									 
}
END_TEST


START_TEST(test_align_short_more_complex){
	/*
	Test of aligning adapters where an edit distance of 2 is allowed
	or there are double adapters in the sequence with varying distances.

	*/
	
	int i;
	char *adapter_seqs[] = {"AGCTGTGCTC", "GCTAGTGTCGTC", "TGAGCTGCTGTA"};
	struct AlignmentResult ar;
	initialize_alignment_result(&ar, 1);
	struct SeqList adapters;
	initialize_seq_list(&adapters, 3);
	
	for (i=0; i<3; i++){
		add_seq_to_seqlist(&adapters, i, adapter_seqs[i], NULL, NULL);
	}
	
  char *reads[] = {"TAGTGTCGTCTCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCAT", //distance 2
  						 		 "GCTGTGTGGTCTCGACTCAGCACTAGCATCGCGATCGATCGATCGATCGATCAT", //distance2
   						 		 "AGCTGTGCTCGCTAGTGTCGTCGCTCGATCAGCATCTACGTCGATCATCGACTA", //2 adapters
  						 		 "AGCCGTGCTCGCTAGTGTCGTCGCTCGATCAGCATCTACGTCGATCATCGACTA", //2 adapters
  						 		 "GCTGTGCTCTGAGCGCAGTAGTCGCTGCTCTCTCGCTCGCTCGCTCGCTCG", //2 adapters
  						 		 "GCTGCTCGATGCTACGATCGATCGATCGATCGATCGACACTGACATCGATCGA",
  								};
  								
  								//do a case when stop distance = 1 and there is 
  								//one at the beginning with 1 distance and enxt with 0mak
  								
	int starts[] = {0, 0, 0, 8, 0, 0};
	int buffers[] = {2, 2, 2, 4, 2, 2};
	int start_answers[] = {0, 0, 0, 10, 0, -1};
	int end_answers[] = {9, 10, 9, 21, 8, -1};
	int max_distances[] = {2, 2, 1, 0, 2, 2, 2};
	
	int distances[] = {2, 2, 0, 0, 1, strlen(reads[5])};
	int indices[] = {1, 1, 0, 1, 0, -1};
	char *cigars[] = {"2D10=", "3=1D4=1X3=", "10=", "12=", "1D9=", NULL};
	 
	for (i=0; i<6; i++){
		align_short(reads[i], adapters, max_distances[i], 0, "start", starts[i], 
								buffers[i], &ar);

		ck_assert_int_eq(start_answers[i], ar.start);
		ck_assert_int_eq(end_answers[i], ar.end);
		ck_assert_int_eq(indices[i], ar.index);
		ck_assert_int_eq(distances[i], ar.distance);
		if (i !=5) ck_assert_str_eq(cigars[i], ar.cigar);
		else ck_assert_str_eq(ar.cigar, "");
	}
  free_alignment_result(&ar);
	free_seq_list(&adapters);
}
END_TEST


START_TEST(test_short_align_multiple_alignments){
	
	int i;
	char *adapter_seqs[] = {"AGCTGTGCTC", "GCTAGTGTCGTC", "TGAGCTGCTGTA"};
	struct AlignmentResult ar;
	
	initialize_alignment_result(&ar, 1000);
	struct SeqList adapters;
	initialize_seq_list(&adapters, 3);
	for (i=0; i<3; i++){
		add_seq_to_seqlist(&adapters, i, adapter_seqs[i], NULL, NULL);
	}

	char *reads[6] = {"AGCTGTGCTCYYYYYYYYAGCTGTGCTC",
                    "AGCTGTGCTCYYYYYYYYAGCTGTGCTC",
                    "AGCTGTGCTCYYYYYYYYAGCTGTGCTC",
                    "AGCTGTGCTCGCTAGTGTCGTC",
                    "ACCTGTGCTCGCTAGTGTCGTC",
                    "ACCTGTGCTCGCTAGCGTCGTCTGAGCTGCTGTA",
                    };

  int starts[] = {0, 18, 0, 0, 0, 0};
	int buffers[]   = {27, 9, 27, 22, 22, 34};
  
  char *best_locs[6] = {"start", "start", "end", "start", 
                        "end", "end"};
  int start_answers[] = {0, 18, 18, 0, 10, 22}; 
  int index_answers[] = {0, 0, 0, 0, 1, 2};
  int distance_answers[] = {0, 0, 0, 0, 0, 0};
  
	for (i=0; i<6; i++){
  
		align_short(reads[i], adapters, 1, 0, best_locs[i], starts[i],
                buffers[i], &ar);
		ck_assert_int_eq(start_answers[i], ar.start);
    ck_assert_int_eq(distance_answers[i], ar.distance);
    ck_assert_int_eq(index_answers[i], ar.index);

  }
	free_alignment_result(&ar);
	free_seq_list(&adapters);
									 

}
END_TEST

START_TEST(test_misc){
  
  int nums1[] = {0, 2, 5};
  int nums2[] = {2, 5, 0};
  
  ck_assert_int_eq(argextremum(nums1, 3, "min"), 0);
  ck_assert_int_eq(argextremum(nums1, 3, "max"), 2);
  ck_assert_int_eq(argextremum(nums2, 3, "min"), 2);
  ck_assert_int_eq(argextremum(nums2, 3, "max"), 1);



}

Suite *alignment_suite(void)
{
    Suite *s;
    TCase *tc_core;
    
    s = suite_create("ALIGNMENT");
    tc_core = tcase_create("Core");
    
		tcase_add_test(tc_core, test_get_insert_alignment_result);
    tcase_add_test(tc_core, test_get_start_end);
		tcase_add_test(tc_core, test_align_long_simple);
		tcase_add_test(tc_core, test_align_long_complex);
    tcase_add_test(tc_core, test_align_short_simple);
    tcase_add_test(tc_core, test_align_short_complex);
    tcase_add_test(tc_core, test_align_short_more_complex);
    tcase_add_test(tc_core, test_short_align_multiple_alignments);
    tcase_add_test(tc_core, test_align_pairs);
    tcase_add_test(tc_core, test_misc);    

    suite_add_tcase(s, tc_core);
    
    return s;
}


int main(void)
{

  int number_failed=0;
  Suite *s;
  SRunner *sr;

  s = alignment_suite();
  sr = srunner_create(s);
  srunner_set_fork_status(sr, CK_NOFORK);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
