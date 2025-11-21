CC = gcc
CXX = c++
MEM = valgrind
COV = gcov

TEST=1
DEBUG=1
VERBOSE=0

LDFLAGS = -lcheck -lpthread -lcheck -lgcov -lm -lz
CFLAGS = -Wall -Wextra -rdynamic -g -Wno-unused-function
WFAFLAGS = -fopenmp
COVFLAGS = --coverage
SRC_MAKE_ARGS = TEST=$(TEST) DEBUG=$(DEBUG) VERBOSE=$(VERBOSE)
MEM_ARGS = -s --leak-check=full --show-leak-kinds=all --log-fd=9 9>

EXT_FOLDER = ../ext
SRC = ../src
TEST_FROM_SRC = ../tests
LOGS_FOLDER = ./logs
WFA_FOLDER=$(EXT_FOLDER)/WFA2

EDLIB_BUILD_FOLDER = $(EXT_FOLDER)/edlib/edlib/build
EDLIB_INCLUDE_FOLDER = $(EXT_FOLDER)/edlib/edlib/include
EDLIB_OBJ = $(EDLIB_BUILD_FOLDER)/edlib.o
LEDLIB = -L$(EDLIB_BUILD_FOLDER) -ledlib
IEDLIB = -I$(EDLIB_INCLUDE_FOLDER)

IWFA = -I$(WFA_FOLDER)
LWFA =  -L$(WFA_FOLDER)/lib -lwfa

SRC_NAMES = utils structures alignment io error
TEST_NAMES = utils structures alignment io badinputs merlign cov mem
SRC_OBJS = $(foreach wrd, $(SRC_NAMES), $(SRC)/$(wrd).o)
TESTS = $(foreach wrd, $(TEST_NAMES), test_$(wrd))
TEST_OBJS = $(foreach wrd, $(TEST_NAMES), test_$(wrd).o)

VERBOSE_FLAG = 
ifeq ($(VERBOSE),1)
	VERBOSE_FLAG = -DVERBOSE
endif

N=100

TARGET=merlign

.PHONY: tests clean test_alignment test_io test_utils test_structures test_badinputs test_merlign test_mem test_cov

all: readme clean clean_logs clean_cov test_build_ext test_build tests rebuild\
		   clean_all

readme:
	@cat README_base.md > README.md
	@echo "Tests performed at $(shell date) (UT:$(shell date +%s))" $ '\n' >> README.md

tests: test_io test_utils clean_structure_gcda test_structures test_badinputs\
			   test_alignment test_cov test_mem
	
test_build_ext:
	$(MAKE) -s clean all -C $(WFA_FOLDER) MARCH_FLAG="" BUILD_WFA_PARALLEL="0"

test_build:
	$(MAKE) clean $(TARGET) -C $(SRC) $(SRC_MAKE_ARGS)

rebuild:
	$(MAKE) -s clean all -C $(WFA_FOLDER)
	$(MAKE) -s clean $(TARGET) -C $(SRC)

test_io.o:
	$(CC) -c test_io.c -o test_io.o $(CFLAGS) -DUNIT_TESTING $(VERBOSE_FLAG)
	$(CC) $(foreach wrd,error utils io structures,$(SRC)/$(wrd).o) test_io.o\
	 	$(CFLAGS)  $(COVFLAGS) $(LDFLAGS) -o test_io

test_io: test_io.o
	@echo "TESTING IO.C"
	@./test_io | sed -z '$ s/\n$$/\n<br \/>\n/' >> README.md
	@tail -n 2 README.md| sed -z '$ s/\n<br \/>\n/\n/'

test_utils.o:
	$(CC) -c test_utils.c -o test_utils.o  $(CFLAGS) $(COVFLAGS) -DUNIT_TESTING\
		$(VERBOSE_FLAG)
	$(CC)  $(foreach wrd,error utils structures, $(SRC)/$(wrd).o)\
	 	test_utils.o $(LDFLAGS) $(CFLAGS) $(COVFLAGS) -o test_utils

test_utils: test_utils.o
	@echo "TESTING UTILS.C"
	@./test_utils | sed -z '$ s/\n$$/\n<br \/>\n/' >> README.md
	@tail -n 2 README.md| sed -z '$ s/\n<br \/>\n/\n/'

test_alignment.o:
	$(CC) -c -g test_alignment.c -o test_alignment.o  $(IEDLIB) $(IWFA)\
		$(CFLAGS) $(COVFLAGS) -DUNIT_TESTING $(VERBOSE_FLAG)
	$(CXX) -g $(foreach wrd,error io utils alignment structures,$(SRC)/$(wrd).o)\
	 	test_alignment.o $(CFLAGS) $(WFAFLAGS) $(COVFLAGS) $(LDFLAGS)\
	  $(LWFA) $(LEDLIB) -o test_alignment

test_alignment: test_alignment.o
	@echo "TESTING ALIGNMENT.C"
	@./test_alignment | sed -z '$ s/\n$$/\n<br \/>\n/' >> README.md
	@tail -n 2 README.md| sed -z '$ s/\n<br \/>\n/\n/'

test_structures.o:
	$(CC) -c test_structures.c -o test_structures.o $(CFLAGS) $(VERBOSE_FLAG)
	$(CC) $(foreach wrd,error utils structures,$(SRC)/$(wrd).o)\
	 	test_structures.o $(CFLAGS)  $(COVFLAGS) $(LDFLAGS) -o test_structures

test_structures: test_structures.o
	@echo "TESTING STRUCTURES.C"
	@./test_structures | sed -z '$ s/\n$$/\n<br \/>\n/' >> README.md
	@tail -n 2 README.md| sed -z '$ s/\n<br \/>\n/\n/'

test_badinputs.o:
	$(CC) -c -g test_badinputs.c -o test_badinputs.o $(IEDLIB) $(CFLAGS) $(IWFA)\
	 	$(WFAFLAGS) -DUNIT_TESTING $(VERBOSE_FLAG)
	$(CXX) $(foreach wrd,error utils io structures alignment ,$(SRC)/$(wrd).o)\
	 	test_badinputs.o $(EDLIB_OBJECT) $(CLAGS) $(COVFLAGS) $(LDFLAGS) $(LEDLIB)\
	 	$(LWFA) $(WFAFLAGS) -o test_badinputs

test_badinputs: test_badinputs.o
	@echo "TESTING BAD INPUTS"
	@./test_badinputs | sed -z '$ s/\n$$/\n<br \/>\n/' >> README.md
	@tail -n 2 README.md| sed -z '$ s/\n<br \/>\n/\n$$/'

test_merlign.o:
	$(CC) -c test_merlign.c -o test_merlign.o $(IEDLIB) $(IWFA) $(CFLAGS)\
	 	$(COVFLAGS) -DUNIT_TESTING $(VERBOSE_FLAG)
	$(CXX) -g $(foreach wrd,error io utils alignment structures,$(SRC)/$(wrd).o)\
	 	test_merlign.o $(CFLAGS) $(WFAFLAGS) $(COVFLAGS) $(LDFLAGS) $(LWFA)\
	 	$(LEDLIB) -o test_merlign

test_merlign:
	@echo "TESTING MERLIGN"
	@$(SRC)/mer --match 0 --file1 ./test_files/simulated_reads_R1.fastq --file2\
	 	./test_files/simulated_reads_R2.fastq --adapters\
	 	./test_files/adapters_for_simulated_reads.fasta --adapter_start=200\
	 	--adapter_buffer=120 --output_path ./test_files/simulated_mer_output_test.csv\
	 	--merged_path ./test_files/simulated_merged_seqs_test.fastq --N1 1000\
	 	> /dev/null 2>&1	
	@$(SRC)/ign  --match -1 --reads ./test_files/simulated_merged_seqs_test.fastq\
	 	--references ./test_files/refs_for_simulated_reads.fasta --mids\
	 	./test_files/MIDs.fasta --mid_start 0 --mid_buffer 15 --output\
	 	./test_files/simulated_ign_output_test.csv --N1 1000 > /dev/null 2>&1
	@$(SRC)/mer --match 0 --file1 ./test_files/simulated_reads_R1.fastq --file2\
	 	./test_files/simulated_reads_R2.fastq --adapters\
	 	./test_files/adapters_for_simulated_reads.fasta --adapter_start=0\
	 	--adapter_buffer=120 --output_path ./test_files/simulated_mer_output_test.csv\
	 	--merged_path ./test_files/simulated_merged_seqs_test.fastq --N1 1000\
	@$(SRC)/mer --match 0 --file1 ./test_files/simulated_reads_R1.fastq --file2\
	 	./test_files/simulated_reads_R2.fastq --adapters\
	 	./test_files/adapters_for_simulated_reads.fasta --adapter_start=1\
	 	--adapter_buffer=120 --output_path ./test_files/simulated_mer_output_test.csv\
	 	--merged_path ./test_files/simulated_merged_seqs_test.fastq --N1 1000\
	


	@./test_merlign | sed -z '$ s/\n$$/\n<br \/>\n/' >> README.md
	@tail -n 2 README.md| sed -z '$ s/\n<br \/>\n/\n/'

test_cov: $(TEST_OBJS) $(SRC_OBJS)
	$(CC) -c test_cov.c -o test_cov.o $(CFLAGS)
	$(CC) test_cov.o $(CFLAGS)  $(COVFLAGS) $(LDFLAGS) -o test_cov
	cd $(SRC); gcov $(SRC)/utils.c | head -n 2 > \
		${TEST_FROM_SRC}/${LOGS_FOLDER}/utils.cov.log
	cd $(SRC); gcov $(SRC)/structures.c | head -n 2 >\
	 	${TEST_FROM_SRC}/${LOGS_FOLDER}/structures.cov.log
	cd $(SRC); gcov $(SRC)/io.c | head -n 2 >\
	 	${TEST_FROM_SRC}/${LOGS_FOLDER}/io.cov.log
	cd $(SRC); gcov $(SRC)/alignment.c | head -n 2 >\
	 	${TEST_FROM_SRC}/${LOGS_FOLDER}/alignment.cov.log     
	
	@echo "TESTING COVERAGE"
	@./test_cov | sed -z '$ s/\n$$/\n<br \/>\n/' >> README.md
	@tail -n 2 README.md| sed -z '$ s/\n<br \/>\n/\n/'

test_mem:
	$(CC) -c test_mem.c -o test_mem.o $(CFLAGS)
	$(CC) test_mem.o $(CFLAGS) $(LDFLAGS) -o test_mem
	@echo "TESTING MEMORY (AND YOUR PATIENCE)"
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/test_structures.mem.log ./test_structures\
		> /dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/test_utils.mem.log ./test_utils > /dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/test_io.mem.log ./test_io > /dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/test_alignment.mem.log ./test_alignment >\
	 	/dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/mer1.mem.log $(SRC)/mer\
		--file1 ./test_files/example_file1_R1.fastq.gz\
		--file2 ./test_files/example_file1_R2.fastq.gz\
		--N1 $(N) > /dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/mer2.mem.log $(SRC)/mer\
		--file1 ./test_files/example_file1_R1.fastq.gz\
		--file2 ./test_files/example_file1_R2.fastq.gz --N1 $(N)\
	 	--end_free1 100 --end_free2 100 --begin_free1 100 --begin_free2 100\
	 	> /dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/mer3.mem.log $(SRC)/mer --file1\
		 	./test_files/example_file1_R1.fastq.gz --file2\
		 	./test_files/example_file1_R2.fastq.gz --N1 $(N)\
		 	--end_free1 100 --end_free2 100 --begin_free1 100 --begin_free2 100\
		 	--trim 0 > /dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/mer4.mem.log $(SRC)/mer --file1\
		 	./test_files/example_file1_R1.fastq.gz --file2\
		 	./test_files/example_file1_R2.fastq.gz --N1 $(N)\
		 	--trim 0 > /dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/mer5.mem.log $(SRC)/mer --match 0\
		--file1 ./test_files/simulated_reads_R1.fastq\
	 	--file2 ./test_files/simulated_reads_R2.fastq\
	 	--output_path ./test_files/simulated_mer_output_test.csv\
	 	--merged_path ./test_files/simulated_merged_seqs_test.fastq --N1 $(N) >\
	 	/dev/null 2>&1	
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/ign1.mem.log $(SRC)/ign\
		--reads ./test_files/ign_example_seqs1.fastq.gz\
	 	--mids ./test_files/MIDs.fasta\
	 	--references ./test_files/BE92.fasta\
	 	--mid_start 0 --mid_buffer 4 --N1 $(N) > /dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/ign2.mem.log $(SRC)/ign\
		--reads ./test_files/ign_example_seqs1.fastq.gz\
		--references ./test_files/BE92.fasta\
	 	--N1 $(N) > /dev/null 2>&1
	@$(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/ign3.mem.log $(SRC)/ign  --match -1\
		--reads ./test_files/simulated_merged_seqs_test.fastq\
		--references ./test_files/BE92.fasta --mids ./test_files/MIDs.fasta\
	 	--mid_start 0 --mid_buffer 15 --output\
	 	./test_files/simulated_ign_output_test.csv --N1 $(N) > /dev/null 2>&1
	@$(SRC)/mer --file1 ./test_files/example_file1_R1.fastq.gz --file2\
	 	./test_files/example_file1_R2.fastq.gz --adapters ./test_files/adapters.fasta\
	 	--adapter_start=200 --adapter_buffer=100 --N1 $(N) 2> /dev/null\
	 	| $(MEM) $(MEM_ARGS) $(LOGS_FOLDER)/merlign.mem.log $(SRC)/ign  --references\
	 	./test_files/BE92.fasta --mids ./test_files/MIDs.fasta --mid_start 0\
		--mid_buffer 15 --output\
	 	./test_files/simulated_ign_output_test.csv --N1 $(N) > /dev/null 2>&1
	@./test_mem | sed -z '$ s/\n$$/\n<br \/>\n/' >> README.md
	@tail -n 2 README.md| sed -z '$ s/\n<br \/>\n/\n/'


clean:
	@echo "Cleaning tests, test objects and source objects."
	@rm -rf $(TEST_OBJS) $(TESTS) $(SRC_OBJS)

clean_cov:
	@echo "Cleaning gcov, gcno, gcda files."
	@rm -rf *.gcov *.gcno *.gcda $(SRC)/*.gcov $(SRC)/*.gcno $(SRC)/*.gcda

clean_logs:
	@echo "Cleaning logs."
	@rm -rf $(LOGS_FOLDER)/*.log

clean_structure_gcda:
	@rm -rf $(SRC)/structures.gcda

clean_all:
	$(MAKE) clean
	$(MAKE) clean_cov
	$(MAKE) clean_structure_gcda
