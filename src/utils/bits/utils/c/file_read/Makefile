include ../../../../../../defs.local
OBJ_DIR = ../../../../../../obj
UTILITIES_DIR = ../../../../../utils/
INCLUDES = -I$(UTILITIES_DIR)/bedFile/ \
			-I$(UTILITIES_DIR)/version/ \
			-I$(UTILITIES_DIR)/gzstream/ \
			-I$(UTILITIES_DIR)/genomeFile/ \
			-I$(UTILITIES_DIR)/lineFileUtilities/ \
			-I$(UTILITIES_DIR)/fileType/ \
			-I$(UTILITIES_DIR)/bits/utils/c/timer/ \
			-I$(UTILITIES_DIR)/bits/utils/c/mt/ \
			-I$(UTILITIES_DIR)/bits/interval_intersection/lib/seq 

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= file_read.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS=
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -O3 -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -lm \
		$(LOCAL_INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

PHONY: clean
