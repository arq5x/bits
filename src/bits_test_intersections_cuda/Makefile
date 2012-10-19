include ../../defs.cuda
include ../../defs.local
UTILITIES_DIR = ../utils/
OBJ_DIR = ../../obj/
BIN_DIR = ../../bin/

LFLAGS_A=-fPIC
LFLAGS_B=-L$(SDK_PATH)/C/lib \
		 -L$(SDK_PATH)/C/common/lib/$(PLATFORM) \
		 -L$(SDK_PATH)/shared/lib \
		 -L$(CUDPP_PATH)/lib \
		 -L$(CUDA_LIB)/ \
		 -lcudpp \
		 -lcudart \
		 -lcuda

# -------------------
# define our includes
# -------------------
INCLUDES = -I$(UTILITIES_DIR)/bedFile/ \
           -I$(UTILITIES_DIR)/version/ \
           -I$(UTILITIES_DIR)/gzstream/ \
           -I$(UTILITIES_DIR)/genomeFile/ \
           -I$(UTILITIES_DIR)/lineFileUtilities/ \
           -I$(UTILITIES_DIR)/fileType/ \
           -I$(UTILITIES_DIR)/bits/interval_intersection/lib/seq/ \
           -I$(UTILITIES_DIR)/bits/interval_intersection/lib/cuda/ \
		   -I$(UTILITIES_DIR)/bits/utils/c/mt/ \
		   -I$(UTILITIES_DIR)/bits/utils/c/timer/ \
		   -I$(UTILITIES_DIR)/bits/utils/c/file_read/ \
		   -I$(CUDPP_INC)
		   #-I$(UTILITIES_DIR)/bits/utils/cuda/radixsort/ \
		   #-I$(UTILITIES_DIR)/bits/utils/cuda/reduce/ 
# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= 
OBJECTS= $(SOURCES:.cpp=.o)
CUDA_SOURCES=bits_test_cuda_main.cu \
			 bits_test_cuda.cu 
CUDA_OBJECTS=bits_test_cuda_main.cu.o \
			 bits_test_cuda.cu.o
_EXT_OBJECTS=bedFile.o \
			 lineFileUtilities.o \
			 gzstream.o \
			 fileType.o \
			 interval.o \
			 bsearch.o \
			 mt.o \
			 timer.o \
			 file_read.o \
			 genomeFile.o \
			 interval_cuda.cu.o 
			 #radixsort.cu.o \
			 #radixsort.o \
			 #reduce.cu.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= bits_test_cuda
BUILT_CUDA_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(CUDA_OBJECTS))

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS) $(BUILT_CUDA_OBJECTS) $(EXT_OBJECTS) 
	@echo " * linking $(PROGRAM)"
	g++ $(LFLAGS_A) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ $(LIBS) \
		$(LOCAL_LIBS) -lz $(LFLAGS_B) -lgsl -lgslcblas -lm

$(BUILT_OBJECTS): $(SOURCES)
	@echo " * compiling" $(*F).cpp
	$(NVCC) -m64 -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) \
		$(LOCAL_INCLUDES)

$(BUILT_CUDA_OBJECTS): $(CUDA_SOURCES)
	@echo "  * compiling" $(*F)
	$(NVCC) -m64 -c -o $@ $(*F) $(INCLUDES) $(LFLAGS_B) 


$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/lineFileUtilities/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bedFile/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/gzstream/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/genomeFile/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/fileType/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bits/interval_intersection/lib/seq/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bits/interval_intersection/lib/cuda/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bits/binary_search/lib/seq/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bits/binary_search/lib/cuda/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bits/utils/c/mt/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bits/utils/c/timer/
	@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bits/utils/c/file_read/

#@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bits/utils/cuda/reduce/
#@$(MAKE) --no-print-directory -C $(UTILITIES_DIR)/bits/utils/cuda/radixsort/

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
