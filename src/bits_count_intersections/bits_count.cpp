/*****************************************************************************
bits_count.cpp
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
#include "lineFileUtilities.h"
#include "bits_count.h"
#include "timer.h"
#include "file_read.h"


/*
Constructor
*/
BitsCount::BitsCount(string bedAFile, string bedBFile, string genomeFile) {

    _bedAFile = bedAFile;
    _bedBFile = bedBFile;
    _genomeFile = genomeFile;
    
    // create new BED file objects for A and B
    _bedA = new BedFile(bedAFile);
    _bedB = new BedFile(bedBFile);
    //_genome = new BedFile(genomeFile);
    _genome = new GenomeFile(genomeFile);
    
    CountOverlaps();
}


/*
Destructor
*/
BitsCount::~BitsCount(void) {
}

void BitsCount::CountOverlaps() {
    struct interval *A;
	unsigned int *B_starts, *B_ends;
	unsigned int A_size, B_size;

	read_and_map_files_to_array(_genome,
			           &_offsets,
					   _bedA,
					   _bedB,
					   &A,
					   &A_size,
					   &B_starts,
					   &B_ends,
					   &B_size);


    uint32_t tot_overlaps = 
		count_intersections_bsearch_seq_mem( A,
											 A_size,
											 B_starts,
											 B_ends,
											 B_size);

	cout << tot_overlaps << endl;
}
