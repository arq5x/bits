/*****************************************************************************
file_read.h
(c) 2012 - Ryan M. Layer
Hall Laboratory
Quinlan Laboratory
Department of Computer Science
Department of Biochemistry and Molecular Genetics
Department of Public Health Sciences and Center for Public Health Genomics,
University of Virginia
rl6sf@virginia.edu

Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#ifndef __FILE_READ___
#include "bedFile.h"
#include "genomeFile.h"
#include "interval.h"
#include <vector>
#include <map>
#include <string>
void read_and_map_files(GenomeFile *_genome,
				   map<string,CHRPOS> *_offsets,
				   BedFile *_bedA,
				   BedFile *_bedB,
				   vector<struct interval> *_A,
				   vector<struct interval> *_B);
void read_and_map_files_to_array(GenomeFile *_genome,
								 map<string,CHRPOS> *_offsets,
								 BedFile *_bedA,
								 BedFile *_bedB,
								 struct interval **_A,
								 unsigned int *A_size,
								 unsigned int **_B_starts,
								 unsigned int **_B_ends,
								 unsigned int *B_size);

void read_and_map_files_to_array_skip_vector(GenomeFile *_genome,
								 map<string,CHRPOS> *_offsets,
								 BedFile *_bedA,
								 BedFile *_bedB,
								 struct interval **_A,
								 unsigned int *A_size,
								 unsigned int **_B_starts,
								 unsigned int **_B_ends,
								 unsigned int *B_size);

void read_and_map_files_to_interval_arrays_skip_vector(GenomeFile *_genome,
								 map<string,CHRPOS> *_offsets,
								 BedFile *_bedA,
								 BedFile *_bedB,
								 struct interval **_A,
								 unsigned int *A_size,
								 struct interval **_B,
								 unsigned int *B_size);

#if 1
void rand_human_chr(char *c, unsigned int *c_id);
void rand_human_interval(unsigned int len,
                         char *chr,
                         unsigned int *start,
                         unsigned int *end);
void read_and_map_rand(GenomeFile *_genome,
				   map<string,CHRPOS> *_offsets,
				   unsigned int num_intervals,
				   unsigned int interval_len,
				   vector<struct interval> *_A,
				   vector<struct interval> *_B);
#endif
#endif
