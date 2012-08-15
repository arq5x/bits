/*****************************************************************************
bits_count_cuda.cu
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
#include "bits_count_cuda.h"
#include "interval_cuda.h"
#include "timer.h"
#include "file_read.h"
#include "cutil.h"


/*
Constructor
*/
BitsCountCUDA::BitsCountCUDA(string bedAFile,
							 string bedBFile,
							 string genomeFile) {

    _bedAFile = bedAFile;
    _bedBFile = bedBFile;
    _genomeFile = genomeFile;
    
    // create new BED file objects for A and B
    _bedA = new BedFile(bedAFile);
    _bedB = new BedFile(bedBFile);
    //_genome = new BedFile(genomeFile);
    _genome = new GenomeFile(genomeFile);
    
    CountOverlapsCUDA();
}


/*
Destructor
*/
BitsCountCUDA::~BitsCountCUDA(void) {
}

void BitsCountCUDA::CountOverlapsCUDA() {
	int *prt;
	cudaMalloc(&prt, 0);
	vector<struct interval> A, B;
	read_and_map_files(_genome, &_offsets, _bedA, _bedB, &A, &B);

    uint32_t tot_overlaps = count_intersections_bsearch_cuda(&A[0],
															 A.size(),
															 &B[0],
															 B.size());


	cout << tot_overlaps << endl;
}
