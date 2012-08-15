# ==========================
# BEDTools Makefile
# (c) 2009 Aaron Quinlan
# ==========================

# define our object and binary directories
export OBJ_DIR	= obj
export BIN_DIR	= bin
export SRC_DIR	= src
export CXX		= g++
export CXXFLAGS = -Wall -O2 -fPIC
export LIBS		= -lz
export BT_ROOT  = src/utils/BamTools/

SUBDIRS = $(SRC_DIR)/bits_count_intersections \
		  $(SRC_DIR)/bits_count_intersections_cuda \
		  $(SRC_DIR)/bits_count_intersections_per_interval \
		  $(SRC_DIR)/bits_count_intersections_per_interval_cuda \
		  $(SRC_DIR)/bits_count_significance_cuda \
		  $(SRC_DIR)/bits_enumerate_intersections \
		  $(SRC_DIR)/bits_enumerate_intersections_cuda \
		  $(SRC_DIR)/bits_test_intersections \
		  $(SRC_DIR)/bits_test_intersections_cuda 

UTIL_SUBDIRS =	$(SRC_DIR)/utils/lineFileUtilities \
				$(SRC_DIR)/utils/bedFile \
				$(SRC_DIR)/utils/genomeFile \
				$(SRC_DIR)/utils/gzstream \
				$(SRC_DIR)/utils/fileType \
				$(SRC_DIR)/utils/bits/utils/c/mt \
				$(SRC_DIR)/utils/bits/utils/c/timer \
				$(SRC_DIR)/utils/bits/binary_search/lib/cuda \
				$(SRC_DIR)/utils/bits/binary_search/lib/seq \
				$(SRC_DIR)/utils/bits/interval_intersection/lib/cuda \
				$(SRC_DIR)/utils/bits/interval_intersection/lib/seq

all:
	[ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR)
	[ -d $(BIN_DIR) ] || mkdir -p $(BIN_DIR)
	
	@echo "Building BEDTools:"
	@echo "========================================================="
	
	@for dir in $(UTIL_SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done


.PHONY: all

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
	@rm -Rf $(BT_ROOT)/lib
	@rm -f $(BT_ROOT)/src/api/*.o
	@rm -f $(BT_ROOT)/src/api/internal/*.o
	@rm -Rf $(BT_ROOT)/include

.PHONY: clean
