include ../../../../../../defs.cuda

OBJ_DIR = ../../../../../../obj/
UTILITIES_DIR = ../../../../../utils/

# SDK_PATH=$(HOME)/NVIDIA_GPU_Computing_SDK
# NVCC=/usr/local/cuda/bin/nvcc
# SDK_COMMON_INC=$(SDK_PATH)/C/common/inc
# SDK_SHARED_INC=$(SDK_PATH)/C/shared/inc
# CUDA_INCLUDE=/usr/local/cuda/include

LFLAGS_A=-fPIC
LFLAGS_B=-L$(SDK_PATH)/C/lib \
		 -L$(SDK_PATH)/C/common/lib/$(PLATFORM) \
		 -L$(SDK_PATH)/shared/lib \
		 -L$(CUDPP_PATH)/shared/lib \
		 -lcudpp \
		 -lcudart

# -------------------
# define our includes
# -------------------
INCLUDES=-I$(CUDA_INCLUDE) \
		 -I$(SDK_SHARED_INC) \
		 -I$(SDK_COMMON_INC) \
		 -I$(UTILITIES_DIR)/bits/utils/c/mt \
		 -I$(CUDPP_PATH)/include \
		 -I$(UTILITIES_DIR)/bits/binary_search/lib/cuda \
		 -I$(UTILITIES_DIR)/bits/binary_search/lib/seq \
		 -I$(UTILITIES_DIR)/bits/interval_intersection/lib/seq \
		 -I$(UTILITIES_DIR)/bits/utils/cuda/random \
		 -I$(UTILITIES_DIR)/bits/utils/c/timer
		 #-I$(UTILITIES_DIR)/bits/utils/cuda/radixsort \
		 #-I$(UTILITIES_DIR)/bits/utils/cuda/reduce \

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= 
OBJECTS= $(SOURCES:.c=.o)
CUDA_SOURCES=interval_cuda.cu
CUDA_OBJECTS=interval_cuda.cu.o
_EXT_OBJECTS=bsearch.o mt.o timer.o bsearch.cu.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
BUILT_CUDA_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(CUDA_OBJECTS))

all: $(BUILT_OBJECTS) $(BUILT_CUDA_OBJECTS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).c
	@$(CXX) -c -o $@ $(*F).c $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -fopenmp -O3

$(BUILT_CUDA_OBJECTS): $(CUDA_SOURCES)
	@echo "  * compiling" $(*F)
	$(NVCC) -m64 -c -o $@ $(*F) $(INCLUDES) $(LFLAGS_B) -O3

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C -W $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(BUILT_OBJECTS) $(BUILT_CUDA_OBJECTS)

.PHONY: clean


