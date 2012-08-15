/*****************************************************************************
bits_enumerate_cuda.cu
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
#include "bits_enumerate_cuda.h"
#include "interval_cuda.h"
#include "timer.h"
#include "file_read.h"


/*
Constructor
*/
BitsEnumerateCUDA::BitsEnumerateCUDA(string bedAFile,
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
    
    EnumerateOverlapsCUDA();
}


/*
Destructor
*/
BitsEnumerateCUDA::~BitsEnumerateCUDA(void) {
}





void BitsEnumerateCUDA::EnumerateOverlapsCUDA()
{
	int *prt;
	cudaMalloc(&prt, 0);

	vector<struct interval> A, B;
	read_and_map_files(_genome, &_offsets, _bedA, _bedB, &A, &B);
    										     
	unsigned int *R = (unsigned int *) malloc (sizeof(unsigned int) * A.size());
	unsigned int *E;

    uint32_t tot_overlaps = enumerate_intersections_bsearch_cuda(&A[0],
															 A.size(),
															 &B[0],
															 B.size(),
															 R,
															 &E);
	unsigned int start = 0;
	
    for (size_t i = 0; i < _bedA->bedList.size(); ++i) {
		while (start < R[i]) {
			printf("%s\t%d\t%d\t%s\t%d\t%d\n",
					_bedA->bedList[i].chrom.c_str(),
					_bedA->bedList[i].start,
					_bedA->bedList[i].end,
					_bedB->bedList[ E[start] ].chrom.c_str(),
					_bedB->bedList[ E[start] ].start,
					_bedB->bedList[ E[start] ].end);	
			++start;
		}
	}
}
