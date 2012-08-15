/*****************************************************************************
bits_test_intersections.h
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
#ifndef BITS_COUNT_H
#define BITS_COUNT_H

#include "bedFile.h"
#include "genomeFile.h"
#include "interval.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;



class BitsTest {

public:
    // constructor
    BitsTest(string bedAFile,
			 string bedBFile,
			 string genomeFile,
			 unsigned int N);
    // destructor
    ~BitsTest(void);

private:

    //------------------------------------------------
    // private attributes
    //------------------------------------------------
    string _bedAFile;
    string _bedBFile;
    string _genomeFile;

	unsigned int _N;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;
    GenomeFile *_genome;
    
    map<string,CHRPOS> _offsets;

    //------------------------------------------------
    // private methods
    //------------------------------------------------
    void TestOverlaps(void);
};

#endif /* BITS_COUNT_H */

