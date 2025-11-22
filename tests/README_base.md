# Mer|ign

## Introduction
Merlign (or Mer|ign) is a combination of two command line tools (mer and ign)
written in C for MERging, demultiplexing, alIGNing and trimming of ngs reads.
The main external dependancies are the blazingly fast alignment libraries
[Edlib](https://github.com/Martinsos/edlib#building) (used for aligning MIDs)
and [WFA2](https://github.com/smarco/WFA2-lib)
(used for aligning paired reads to each other and merged reads to references).
Merging is achieved via the Bayesian framework described in
[Error filtering, pair assembly and error correction for next-generation sequencing reads](https://academic.oup.com/bioinformatics/article/31/21/3476/194979?login=false). The optimal
framework is achieved by piping the results of mer into ign (hence the name
mer|ign) however each tool is able to operate on its own and able to produce
output into or read input from files. There is a bash script in examples
showing how to do this in a parallelized fashion.

The components of this pipeline in detail are:

### mer

Bayesian merge paired reads and optionally trim them.

This is achieved first by aligning paired reads to each other via WFA2 library
and then merging them using functionalities within this library which employs
the Bayesian algorithm referenced above. If trim=True then any part of the paired
reads that is out of the alignment is trimmed. For instance if you have two
paired reads with adapters of the form (where FW=forward and BW=backward)


    FW MID - FW INSERT - FW MID - ADAPTER
    BW MID - BW INSERT - BW MID - ADAPTER


then the alignment will be


    ADAPTER - FW MID - FW INSERT - FW MID
              FW MID - FW INSERT - FW MID - ADAPTER


and the adapters will be automatically trimmed.

mer then can produce two outputs, one a fastq.gz file of merged_sequences
(can also direct this to stdout so it can be piped to ign) and another
output csv file which describes for each sequence properties of the merge
(such as number of mismatches, indels etc). The output format is
described in [Outputs](#outputs) section.

When merging, for any position whose merged quality score is less than the mask_threshold
is recorded with lower letters which can be analysed seperately downstream.

### ign

Align references to reads and optionally identify MIDs.

Given references, they and their reverse complements are aligned
to the merged reads to identify which references align to a given read the best.
This part is achieved mainly through WFA2 and some boiler plate code
and is the first step of ign. If for a particular reference
percent of match to a read exceeds stop_when (default 0.9), then it terminates and does not check
other references against the given read.

Each read's ends are compared against a list of MIDs and their
reverse complements to demultiplex the reads. mid_start and mid_buffer determine
which part of the read is checked for MID alignment. See section Parameters.
These extended ends are aligned to the MIDs using edlib which is a blazingly fast alignment
library that uses Levenshtein distance. mid_stop_when works in a similar way to
stop_when (default 1) above to stop looking at MIDs when the number of mismatches is lower or equal
to this value.

ign produces one output file csv which describes propoerties of MID and
reference alignments for each sequence. The output format is described in
[Outputs](#outputs) section.

No quality filtering is carried out during this analysis however relevant
quality information for each read is outputted and can be integrated into
downstream analysis for filtering low quality reads. See [Outputs](#outputs) section
for reported quality propoerties.

This package was put together for the paper [PAPER](link) therefore its main target
are libraries made from relatively short inserts (~200bp) so that reference sequences
are also short and map completely to the reads.


## Building

You can build merlign by calling

    make purge
    make all

in the root directory, which will also build external dependencies. If you want to rebuild Merlign without
rebuilding external dependencies (faster), then you can call

    make merlign


## Testing

You can test merlign by calling

    make all

in the /tests/ directory. You can also call

    make tests

from the root directory in which case tests will be run and README will be
updated based on the results of the tests. Merlign has only been tested on
Ubuntu 20.04


## Build Flags

DEBUG and VERBOSE turn on and off certain printing functionalities such as error stack.
In particular if both VERBOSE and DEBUG are defined, error stacks and error messages are printed.
If only VERBOSE is defined, error messages are printed but not stacks.

If UNIT_TESTING is defined then this exposes certain internal functions in various files
to the testing unit.

When testing merlign is built with DEBUG=1, VERBOSE=0, UNIT_TESTING=1. Regularly it is
built with DEBUG=0, VERBOSE=1, UNIT_TESTING=0.

If a test is silently failing or some error message seems to be missing compile tests with VERBOSE=1.
If you set DEBUG=0 some coverage tests may fail as some functions only run when DEBUG=1.


## Examples

The examples folder contain a stand alone bash script for calling
merlign (with mer piping results to ign) from bash with parallelization
over multiple files.


## Outputs

### mer

The .csv file produced has the same number of rows as the number of reads
processed. It has 7 columns:

<b>read_id</b>: common id of the paired read aligned and merged

<b>M_start/M_end</b>: starting and ending positions of the merge (see example in
introduction)

<b>M_match/M_mismatch/M_deletion/M_insertion</b>: Number of matches, mismatches, deletions and insertions.

The other output which is either directed to a file (if supplied as input) or to stdout
is a list of merged sequences and their merged quality scores. See [Introduction](#introduction)
for details.

### ign

The .csv file produced has the same number of rows as the number of reads processed.
It has 21 columns.

<b>read_id</b>: id of the read processed

<b>align_index</b>: index of the best matching reference (-1 if nothing matches).

<b>is_rc</b>: Whether or not the reference or its reverse complement matches the read.

<b>align_cigar</b>: cigar string for the alignment of the read to the best matching reference.

<b>align_mods</b>: The substitutions and indels that are required to transform the best matchong r
eference to the read.

<b>align_mod_quals</b>: quality scores at position where there is a mismatch or an indel. If it is
a read indel then quality is indicated as an empty space.

<b>align_av_err</b>: The expected number of sequencing errors on the part of the read that aligns to
the reference, which is a sum of per position quality based errror rates.

<b>M1_index/M2_index</b>: Indices of the MIDs identified at either ends of the read. -1 if nothing is identified.

<b>M1_start/M1_end/M2_start/M2_end</b>: positions on the read where the MID alignments start and end.

<b>M1_match/M1_mismatch/M1_deletion/M1_insertion/M2_match/...</b>: Number of matches, mismatches, deletions and insertions for each identified MID.


## Parameters

### Mer

<b>--file1 (-F), --file2 (-G)</b>: paths of paired reads to merge

<b>--merged_path (-M)</b>: If supplied, writes merged sequences and qualities to here. Otherwise it sends it to stdout.

<b>--output_path (-0)</b>: Path of the .csv outputs file (see Outputs section).

<b>--match (-q), --mismatch (-w), --gap_open/extend (-e/-r), --end_free1/free2 (-t/-y),
--begin_free1/free2 (-u/-i)</b>: alignment penalty parameters supplied to WFA2 where 1 and 2
refers to paired reads 1 and 2  gap_open and gap_extend apply to both
reads (default: 0, 1, 6/4, 0/0, 0/0).

<b>--N0/1 (-B/-N)</b>: A debugging option to only consider subset of sequences between integers N0 and N1
(default: 0 and -1).

<b>--qtype (-Q)</b>: quality type of the file could be one of Solexa+64, Phred+33, Phred+64 (default: Phred+33).

<b>--trim (-T)</b>: If 1, trims merged paired reads or not (see Introduction for details, default: 1).

<b>--mask_threshold (-K)</b>: If merged quality of the paired read is lower than this value, the letter
at this position is recorded as a lower case letter for identification in downstream analysis (default: 30)

<b>--merge_debug (-D)</b>: If 1, prints some extra information about the alignment (to merged_path
if supplied. otherwise to stdout, default: 0)

### Ign

--reads (-F): path to the reads file

--output (-G): path of the .csv file for outputs (see section Outputs).

--mids (-A): path for the MID files. If not supplied MID identification is skipped.

--mid_start (-S) / mid_buffer (-E): Given a read and a MID, it defines the two flanking regions corresponding
to

    [mid_start, mid_start + MID length + mid_buffer]

and

   [read length - mid_start - mid_buffer, read length - mid_start + max MID length].

If reads are not trimmed it is advised to use a larger mid_buffer (default: 0, 3).

<b>--match (-q), --mismatch (-w), --gap_open/extend (-e/-r), --read_begin/end_free (-t/-y), --ref_begin/end_free (-u/-i)</b>: alignment penalty parameters supplied to WFA2. gap_open/gap_extend applies to the reference (default: 0, -1, 4/2, 15/15, 5/5)

<b>--N0/1 (-B/-N)</b>: A debugging option to only consider subset of sequences between integers N0 and N1 (default: 0 and -1).

<b>mid_stop_when (-W)</b>: When aligning multiple MIDs, if any of the MIDs has number of mismatches lower or equal to
this value, ign will stop and not attempt to align remaining MIDs (default: 1).

<b>mid_max_distance (-M)</b>: A MID alignment is only accepted if number of mismatches is lower than this. Otherwise
matched MID index and other values in the output are set to -1 (default: 2).

<b>--reference_stop_when (-R)</b>: When aligning multiple references, if any reference (or its reverse) aligns with a
%match (compared to minimum of read and reference) is greater or equal to this value, then ign stops and
does not attempt to align remaining references (default: 0.9).

<b>--qtype (-Q)</b>: quality type of the file could be one of Solexa+64, Phred+33, Phred+64 (default: Phred+33).

<b>--dif_threshold (-T) </b>: Similary to mid_max_distance, if number of mismatches is higher than this value
read is recorded as OTHER in the output file rather than with the index of the reference it best matches (default: 50)

<b>--realign_qual_threshold (-V), --realign_error_threshold (-C)</b>: By default
ign creates two aligners:

A primary aligner with penalties mismatch: 2, gap_open1: 4, gap_extension1: 2 and match: 0.

A more precise aligner with penalties mismatch: 2, gap_open1: 4, gap_extension1: 2 and match: -1.

The match, mismatch, reference gap_open/gap_extend scores of the first aligner is then overwritten with the
parameters above, if supplied. The second aligner is only used when supplied match parameter is 0.

Note that the main difference between two aligners is the match score and the former case is much faster than the latter.
However for a sequence and read which has lots of mismatches and relatively high values begin/end_free parameters,
this might result in short alignments (where begin/end_free are completely exhausted). If match != 0 then
it only uses the primary aligner. If match=0 and if the alignment distance is higher than realign_threshold (default: 4)
and also the aligned error probability is less than realign_error_threshold (default: 2) then the precise aligner
is used. The latter is not to realign for reads that have low quality and more expected errors. For the
fastest option set either of the thresholds to 0 and then precise aligner is skipped and only an aligner with match=0
is used. If you want a non-fast but single aligner, then just set for instance match=1.

<b>--alignment_debug (-D)</b>: If 1, prints some extra information about the alignment to the output path (default: 0)


## Testing Status

Unit tests done via [Check](https://libcheck.github.io/check/),
which uses [gcov](https://gcc.gnu.org/onlinedocs/gcc/Gcov.html) and
[Valgrind](https://valgrind.org/) to carry out coverage and memory tests.
