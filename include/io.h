#ifndef IO_H
#define IO_H

#include <zlib.h>
#include <stdio.h>
#include <string.h>

#include "../ext/kseq.h"
#include "../include/structures.h"
#include "../include/utils.h"
KSEQ_INIT(gzFile, gzread);

UTILITY_STATIC_HEADER(int count_sequences(char *file_name));

struct SeqList load_sequences(char *file_name, int add_rc);

int get_seq(gzFile* fp, char* buffer, struct Seq *seq, kseq_t *seq_record,
						int debug);

#endif

