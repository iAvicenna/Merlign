SRC_FOLDER = src
EXT_FOLDER = ext
WFA2_FOLDER = $(EXT_FOLDER)/WFA2
EDLIB_FOLDER = $(EXT_FOLDER)/edlib
TESTS_FOLDER = tests
TARGET = merlign

DEBUG=0
VERBOSE=1

.PHONY: ext all merlign purge clean tidy tests

ext: 
	git submodule update --init --recursive
	cp $(EXT_FOLDER)/Makefile_WFA2 $(EXT_FOLDER)/WFA2/Makefile
	$(MAKE) all -C $(EXT_FOLDER)
	
all:
	$(MAKE) ext
	$(MAKE) merlign -C $(SRC_FOLDER) DEBUG=$(DEBUG) VERBOSE=$(VERBOSE)
	
$(TARGET):
	$(MAKE) $(TARGET) -C $(SRC_FOLDER) 

purge:
	$(MAKE) clean -C $(EXT_FOLDER)
	$(MAKE) clean -C $(SRC_FOLDER) 
	$(MAKE) clean -C $(TESTS_FOLDER)

clean:
	$(MAKE) clean -C $(SRC_FOLDER) 
	$(MAKE) clean -C $(TESTS_FOLDER)


tests:
	@echo "make tests needs to rebuild libraries with test specific flags and then\
		rebuild them with default flags."
	@$(MAKE) -s all -C $(TESTS_FOLDER) 
	@cp $(TESTS_FOLDER)/README.md README.md
