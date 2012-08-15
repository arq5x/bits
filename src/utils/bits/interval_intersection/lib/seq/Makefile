OBJ_DIR = ../../../../../../obj/
UTILITIES_DIR = ../../../../../utils/

#SDK_PATH=$(HOME)/NVIDIA_GPU_Computing_SDK
#CUDPP_PATH=$(HOME)/Downloads/cudpp_src_2.0
#NVCC=/usr/local/cuda/bin/nvcc
#SDK_COMMON_INC=$(SDK_PATH)/C/common/inc
#SDK_SHARED_INC=$(SDK_PATH)/C/shared/inc
#CUDA_INCLUDE=/usr/local/cuda/include
#CUDPP_PATH=$(HOME)/Downloads/cudpp_src_2.0

#LFLAGS_A=-fPIC
#LFLAGS_B=-L/usr/local/cuda/lib64 \
		 #-L$(SDK_PATH)/C/lib \
		 #-L$(SDK_PATH)/C/common/lib/linux \
		 #-L$(SDK_PATH)/shared/lib \
		 #-L$(CUDPP_PATH)/shared/lib \
		 #-lcudpp \
		 #-lcudart \
		 #-lcudpp_x86_64 

# -------------------
# define our includes
# -------------------
INCLUDES=-I$(UTILITIES_DIR)/bits/utils/c/mt \
		 -I$(UTILITIES_DIR)/bits/binary_search/lib/seq \
		 -I$(UTILITIES_DIR)/bits/utils/c/timer

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= interval.c
OBJECTS= $(SOURCES:.c=.o)
_EXT_OBJECTS=bsearch.o mt.o timer.o
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
BUILT_CUDA_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(CUDA_OBJECTS))

all: $(BUILT_OBJECTS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).c
	@$(CXX) -O3 -c -o $@ $(*F).c $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C -W $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(BUILT_OBJECTS) $(BUILT_CUDA_OBJECTS)

.PHONY: clean
