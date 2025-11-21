CC = gcc
CXX = g++

EXT_FOLDER = ../ext
WFA_FOLDER = $(EXT_FOLDER)/WFA2
EDLIB_INCLUDE_FOLDER = $(EXT_FOLDER)/edlib/edlib/include
EDLIB_BUILD_FOLDER = $(EXT_FOLDER)/edlib/edlib/build

CFLAGS = -Wall -Wextra -Wno-unused-function
CXXFLAGS = -fopenmp
LDFLAGS = -lwfa -lm -lrt -lz -ledlib  -L$(WFA_FOLDER)/lib -L$(EDLIB_BUILD_FOLDER)
INC = -I$(EDLIB_INCLUDE_FOLDER) -I$(WFA_FOLDER)

DEBUG_FLAGS =
DEBUG = 0

ifeq ($(DEBUG),1)
	DEBUG_FLAGS = -DDEBUG -rdynamic -g
	CFLAGS := $(CFLAGS) -O0
	CXXFLAGS := $(CXXFLAGS) -O0
else
	CFLAGS := $(CFLAGS) -O3
	CXXFLAGS := $(CXXFLAGS) -O3
endif

VERBOSE = 1
VERBOSE_FLAG = 
ifeq ($(VERBOSE),1)
	VERBOSE_FLAG = -DVERBOSE
endif

TEST_FLAGS = 
TEST = 0
ifeq ($(TEST),1)
	TEST_FLAGS =  -DUNIT_TESTING  
	CFLAGS :=  $(CFLAGS) -Wno-unused-parameter --coverage
	CXXFLAGS := $(CXXFLAGS) --coverage
	LDFLAGS := $(LDFLAGS) -lgcov
endif

FILES = structures io utils error alignment mer ign
OBJFILES = structures.o io.o utils.o error.o alignment.o
TARGET = merlign
EXE = mer ign
.PHONY: clean tidy

$(TARGET): $(OBJFILES)
	$(CC) -c mer.c -o mer.o $(CFLAGS) $(INC) $(DEBUG_FLAGS) $(VERBOSE_FLAG)          
	
	$(CXX) $(OBJFILES) mer.o -o mer $(CXXFLAGS) $(INC) $(LDFLAGS)\
		$(LINKS) $(DEBUG_FLAGS) 
	
	$(CC) -c ign.c -o ign.o $(CFLAGS) $(INC) $(DEBUG_FLAGS) $(VERBOSE_FLAG)
	
	$(CXX) $(OBJFILES) ign.o -o ign $(CXXFLAGS) $(INC) $(LDFLAGS)\
	 	$(LINKS) $(DEBUG_FLAGS)  
	

error.o: error.c
	$(CC) -c error.c -o error.o $(CFLAGS) $(DEBUG_FLAGS) $(TEST_FLAGS) $(VERBOSE_FLAG)

structures.o: structures.c
	$(CC) -c structures.c -o structures.o $(CFLAGS) $(DEBUG_FLAGS) $(VERBOSE_FLAG)

io.o: io.c
	$(CC) -c io.c -o io.o $(CFLAGS) $(DEBUG_FLAGS) $(TEST_FLAGS) $(VERBOSE_FLAG)
	
utils.o: utils.c
	$(CC) -c utils.c -o utils.o $(CFLAGS) $(DEBUG_FLAGS) $(TEST_FLAGS)\
		$(VERBOSE_FLAG)
	
alignment.o: alignment.c
	$(CC) -c alignment.c -o alignment.o $(CFLAGS) $(INC) $(DEBUG_FLAGS)\
		$(TEST_FLAGS) $(VERBOSE_FLAG)
	
clean:
	@echo "Cleaning source object and executable files"
	@rm -f $(OBJFILES) $(EXE) $(foreach wrd, $(EXE), $(wrd).o) 
	@echo "Cleaning source object files, test files and debugging files"
	@rm -f $(OBJFILES) $(foreach wrd, $(EXE), $(wrd).o) $(foreach wrd, $(FILES), $(wrd).gcda) $(foreach wrd, $(FILES), $(wrd).gcno) $(foreach wrd, $(FILES), $(wrd).*.gcov) 
