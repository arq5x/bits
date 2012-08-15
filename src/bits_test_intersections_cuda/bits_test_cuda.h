/*****************************************************************************
bits_test_cuda.h
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
#ifndef BITS_COUNT_CUDA_RAND_H
#define BITS_COUNT_CUDA_RAND_H

#include "bedFile.h"
#include "genomeFile.h"
#include "interval.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

class BitsTestCUDA {

public:
    // constructor
    // destructor
	BitsTestCUDA(string bedAFile,
				 string bedBFile,
				 string genomeFile,
				 unsigned int N);
    ~BitsTestCUDA(void);

private:

    string _bedAFile;
    string _bedBFile;
    string _genomeFile;
	unsigned int _N;
    BedFile *_bedA, *_bedB;
    GenomeFile *_genome;
    
    map<string,CHRPOS> _offsets;

    void TestOverlapsCUDA(void);
};

#endif /* BITS_COUNT_H */

