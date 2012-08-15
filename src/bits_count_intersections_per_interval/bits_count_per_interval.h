/*****************************************************************************
bits_count_per_interval.h
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
#ifndef BITS_COUNT_PER_INTERVAL_H
#define BITS_COUNT_PER_INTERVAL_H

#include "bedFile.h"
#include "genomeFile.h"
#include "interval.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;



class BitsCountPerInterval {

public:
    // constructor
    BitsCountPerInterval(string bedAFile, string bedBFile, string genomeFile);
    // destructor
    ~BitsCountPerInterval(void);

private:

    //------------------------------------------------
    // private attributes
    //------------------------------------------------
    string _bedAFile;
    string _bedBFile;
    string _genomeFile;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;
    GenomeFile *_genome;
    
    map<string,CHRPOS> _offsets;

    //------------------------------------------------
    // private methods
    //------------------------------------------------
    void CountOverlapsPerInterval(void);
};

#endif /* BITS_COUNT_H */

