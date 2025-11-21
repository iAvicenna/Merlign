C++ = c++
MAKE = make
WFA_FOLDER = WFA2
EDLIB_FOLDER = edlib
EDLIB_INCLUDE_FOLDER = $(EDLIB_FOLDER)/edlib/include
EDLIB_SRC_FOLDER = $(EDLIB_FOLDER)/edlib/src
EDLIB_BUILD_FOLDER = $(EDLIB_FOLDER)/build
EDLIB_OBJ = $(EDLIB_BUILD_FOLDER)/edlib.o
EDLIB_A = $(EDLIB_BUILD_FOLDER)/edlib.o
MARCH_FLAG=-march=native

.PHONY: all

all: edlib.o wfalib.o

edlib.o: 
	$(C++) -c $(EDLIB_SRC_FOLDER)/edlib.cpp -o $(EDLIB_OBJ) -I $(EDLIB_INCLUDE_FOLDER) 
	$(MAKE) -C $(EDLIB_FOLDER) BUILD_DIR=$(EDLIB_BUILD_FOLDER) clean all
	
wfalib.o:
	$(MAKE) -C $(WFA_FOLDER) clean all MARCH_FLAG=$(MARCH_FLAG)
	
clean:
	rm -rf $(WFA_FOLDER) $(EDLIB_FOLDER)
