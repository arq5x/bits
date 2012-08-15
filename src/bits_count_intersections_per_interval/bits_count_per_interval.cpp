/*****************************************************************************
bits_count_per_interval.ccp
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
#include "bits_count_per_interval.h"
#include "timer.h"
#include "file_read.h"


/*
Constructor
*/
BitsCountPerInterval::
BitsCountPerInterval(string bedAFile, string bedBFile, string genomeFile) {

    _bedAFile = bedAFile;
    _bedBFile = bedBFile;
    _genomeFile = genomeFile;
    
    // create new BED file objects for A and B
    _bedA = new BedFile(bedAFile);
    _bedB = new BedFile(bedBFile);
    //_genome = new BedFile(genomeFile);
    _genome = new GenomeFile(genomeFile);
    
    CountOverlapsPerInterval();
}


/*
Destructor
*/
BitsCountPerInterval::~BitsCountPerInterval(void) {
}





void BitsCountPerInterval::CountOverlapsPerInterval() 
{
	vector<struct interval> A, B;
	read_and_map_files(_genome,
			           &_offsets,
					   _bedA,
					   _bedB,
					   &A,
					   &B);


	uint32_t *R = (uint32_t *) malloc (A.size() * sizeof(uint32_t));
	uint32_t tot_overlaps =
			per_interval_count_intersections_bsearch_seq(&A[0],
														 A.size(),
														 &B[0],
														 B.size(),
														 R); 
    for (size_t i = 0; i < _bedA->bedList.size(); ++i) {
		cout << _bedA->bedList[i].chrom <<  "\t" <<
				_bedA->bedList[i].start << "\t" <<
				_bedA->bedList[i].end   << "\t" <<
				R[i] << endl;
	}
}
