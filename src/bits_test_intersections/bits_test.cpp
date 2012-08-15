/*****************************************************************************
bits_test_intersections.cpp
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
#include "bits_test.h"
#include "timer.h"
#include "file_read.h"
#include <math.h>


/*
Constructor
*/
BitsTest::BitsTest(string bedAFile,
				   string bedBFile,
				   string genomeFile,
				   unsigned int N) {

    _bedAFile = bedAFile;
    _bedBFile = bedBFile;
    _genomeFile = genomeFile;
	_N = N;
    
    // create new BED file objects for A and B
    _bedA = new BedFile(bedAFile);
    _bedB = new BedFile(bedBFile);
    //_genome = new BedFile(genomeFile);
    _genome = new GenomeFile(genomeFile);
    
    TestOverlaps();
}


/*
Destructor
*/
BitsTest::~BitsTest(void) {
}





void BitsTest::TestOverlaps() 
{
    vector<struct interval> A, B;
    read_and_map_files(_genome,
                       &_offsets,
                       _bedA,
                       _bedB,
                       &A,
                       &B);

	CHRPOS max_offset = 0;
	map<string,CHRPOS>::const_iterator itr;
	for (itr = _offsets.begin(); itr != _offsets.end(); ++itr)
		max_offset += itr->second;

	unsigned int O;
	double mean,sd,p;

    test_intersections_bsearch_seq(&A[0],
								   A.size(),
								   &B[0],
								   B.size(),
								   _N,
								   max_offset,
								   &O,
								   &mean,
								   &sd,
								   &p);
	printf("O:%u\tE:%f\tsd:%f\tp:%f\n", O, mean, sd, p);
}
