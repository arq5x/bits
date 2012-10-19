/*****************************************************************************
file_read.cpp
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
#include "file_read.h"
#include "bedFile.h"
#include "genomeFile.h"
#include "interval.h"
#include "timer.h"
#include "mt.h"
#include "math.h"
#include <vector>
#include <string>


//{{{ void read_and_map_files(GenomeFile *_genome,
void read_and_map_files(GenomeFile *_genome,
				   map<string,CHRPOS> *_offsets,
				   BedFile *_bedA,
				   BedFile *_bedB,
				   vector<struct interval> *_A,
				   vector<struct interval> *_B)
{
    vector<string> chromList =  _genome->getChromList();
    CHRPOS curr_offset = 0;
    for (size_t c = 0; c < chromList.size(); ++c) {
        string currChrom = chromList[c];
        CHRPOS currChromSize = _genome->getChromSize(currChrom);
		(*_offsets)[currChrom] = curr_offset;
		curr_offset += currChromSize;
	}
    // pload A and B into vectors
    _bedA->loadBedFileIntoVector();

    _bedB->loadBedFileIntoVector();

	string last_chr;
	CHRPOS last_proj = 0;
    // project A into U
	struct interval ivl;
	CHRPOS projected_start;
	CHRPOS projected_end;
    for (size_t i = 0; i < _bedA->bedList.size(); ++i) {
			if (_bedA->bedList[i].chrom.compare(last_chr) != 0)
				last_proj = (*_offsets)[_bedA->bedList[i].chrom];

			projected_start = last_proj + _bedA->bedList[i].start;
			projected_end = last_proj +  _bedA->bedList[i].end;
            ivl.start = projected_start + 1;
            ivl.end   = projected_end;
			_A->push_back(ivl);
    }

    // project B into U
	last_chr = "";
    for (size_t i = 0; i < _bedB->bedList.size(); ++i) {
		//struct interval ivl;
		if (_bedB->bedList[i].chrom.compare(last_chr) != 0)
			last_proj = (*_offsets)[_bedB->bedList[i].chrom];

		projected_start = last_proj + _bedB->bedList[i].start;
		projected_end = last_proj +  _bedB->bedList[i].end;

        ivl.start = projected_start + 1;
        ivl.end   = projected_end;
        _B->push_back(ivl);
    }
}
//}}}

//{{{ void read_and_map_files_to_array(GenomeFile *_genome,
void read_and_map_files_to_array(GenomeFile *_genome,
								 map<string,CHRPOS> *_offsets,
								 BedFile *_bedA,
								 BedFile *_bedB,
								 struct interval **_A,
								 unsigned int *A_size,
								 unsigned int **_B_starts,
								 unsigned int **_B_ends,
								 unsigned int *B_size)
{
    vector<string> chromList =  _genome->getChromList();
    CHRPOS curr_offset = 0;
    for (size_t c = 0; c < chromList.size(); ++c) {
        string currChrom = chromList[c];
        CHRPOS currChromSize = _genome->getChromSize(currChrom);
		(*_offsets)[currChrom] = curr_offset;
		curr_offset += currChromSize;
	}
    // pload A and B into vectors
    _bedA->loadBedFileIntoVector();
	*A_size = _bedA->bedList.size();

    _bedB->loadBedFileIntoVector();
	*B_size = _bedB->bedList.size();

	*_A = (struct interval *) malloc( 
			sizeof(struct interval) * (*A_size));
	string last_chr;
	CHRPOS last_proj = 0;
    for (size_t i = 0; i < (*A_size); ++i) {
			if (_bedA->bedList[i].chrom.compare(last_chr) != 0)
				last_proj = (*_offsets)[_bedA->bedList[i].chrom];

			(*_A)[i].start = last_proj + _bedA->bedList[i].start + 1;
				(*_A)[i].end = last_proj +  _bedA->bedList[i].end;

			last_chr = _bedA->bedList[i].chrom;
    }
    
	*_B_starts = (unsigned int*) malloc( 
				sizeof(struct interval) * (*B_size));
	*_B_ends = (unsigned int*) malloc( 
				sizeof(struct interval) * (*B_size));
    // project B into U
	last_chr = "";
    for (size_t i = 0; i < (*B_size); ++i) {
		if (_bedB->bedList[i].chrom.compare(last_chr) != 0)
			last_proj = (*_offsets)[_bedB->bedList[i].chrom];

		(*_B_starts)[i] = last_proj + _bedB->bedList[i].start + 1;
		(*_B_ends)[i] = last_proj +  _bedB->bedList[i].end;

		last_chr = _bedB->bedList[i].chrom;
    }
}
//}}}

//{{{ void read_and_map_files_to_array(GenomeFile *_genome,
void read_and_map_files_to_array_skip_vector(GenomeFile *_genome,
								 map<string,CHRPOS> *_offsets,
								 BedFile *_bedA,
								 BedFile *_bedB,
								 struct interval **_A,
								 unsigned int *A_size,
								 unsigned int **_B_starts,
								 unsigned int **_B_ends,
								 unsigned int *B_size)
{
    vector<string> chromList =  _genome->getChromList();
    CHRPOS curr_offset = 0;
    for (size_t c = 0; c < chromList.size(); ++c) {
        string currChrom = chromList[c];
        CHRPOS currChromSize = _genome->getChromSize(currChrom);
		(*_offsets)[currChrom] = curr_offset;
		curr_offset += currChromSize;
	}

	_bedA->loadBedFileIntoIntervalArray(_A, A_size, _offsets);
	_bedB->loadBedFileIntoStartEndArrays( _B_starts,
										  _B_ends,
										  B_size ,
										  _offsets );
}
//}}}

//{{{ Random generation only works for hg19
#define min_v(a,b) \
	    ({ __typeof__ (a) _a = (a); \
		        __typeof__ (b) _b = (b); \
		         _a < _b ? _a : _b; })

//{{{unsigned int HG19[] = {
unsigned int HG19[] = {
	249250621, // chr1
	243199373, // chr1
	198022430, // chr1
	191154276, // chr1
	180915260, // chr1
	171115067, // chr1
	159138663, // chr1
	146364022, // chr1
	141213431, // chr1
	135534747, // chr1
	135006516, // chr1
	133851895, // chr1
	115169878, // chr1
	107349540, // chr1
	102531392, // chr1
	90354753 , // chr1
	81195210 , // chr1
	78077248 , // chr1
	59128983 , // chr1
	63025520 , // chr1
	48129895 , // chr1
	51304566 , // chr1
	155270560, // chr1
	59373566  // chr1
};
//}}}

//{{{void rand_human_chr(char *c, unsigned int *c_id) {
void rand_human_chr(char *c, unsigned int *c_id) {
    int max = 23;
    //unsigned int bits = (int)( ceil(log(max_offset)/log(2) ) );
    //unsigned int mask = (2 << (bits-1)) - 1;
    //unsigned int bits = 5;
    unsigned int mask = 31;
    
    *c_id = get_rand(max,mask);

    if ( (*c_id) == 23 ) 
        sprintf(c, "chrY");
    else if ( (*c_id) == 22) 
        sprintf(c, "chrX"); 
    else 
        sprintf(c, "chr%d", *c_id + 1);
}   
//}}}

//{{{void rand_human_interval(unsigned int len,
void rand_human_interval(unsigned int len,
                         char *chr,
                         unsigned int *start,
                         unsigned int *end)
{
    unsigned int c_id;
    rand_human_chr(chr, &c_id);
   
    unsigned int max_offset = HG19[c_id];
   
    unsigned int bits = (int)( ceil(log(max_offset)/log(2) ) );
    unsigned int mask = (2 << (bits-1)) - 1;
    
    *start = get_rand(max_offset,mask);
    *end = min_v(max_offset, *start+len+1);
}   
//}}}

//{{{ void read_and_map_files(GenomeFile *_genome,
void read_and_map_rand(GenomeFile *_genome,
				   map<string,CHRPOS> *_offsets,
				   unsigned int num_intervals,
				   unsigned int interval_len,
				   vector<struct interval> *_A,
				   vector<struct interval> *_B)
{
	init_genrand((unsigned) time(NULL));

    vector<string> chromList =  _genome->getChromList();
    CHRPOS curr_offset = 0;
    for (size_t c = 0; c < chromList.size(); ++c) {
        string currChrom = chromList[c];
        CHRPOS currChromSize = _genome->getChromSize(currChrom);
		(*_offsets)[currChrom] = curr_offset;
		curr_offset += currChromSize;
	}

	string last_chr;
	struct interval ivl;
	//CHRPOS projected_start;
	//CHRPOS projected_end;
	char chrom[6];
	unsigned int start, end;
	for (unsigned int i = 0; i < num_intervals; ++i) {
		rand_human_interval(interval_len, chrom, &start, &end);
		string c(chrom);
		CHRPOS projected_start = (*_offsets)[c] + start;
		CHRPOS projected_end   = (*_offsets)[c] + end;
		ivl.start = projected_start + 1;
		ivl.end   = projected_end;
		_A->push_back(ivl);
    }

	for (unsigned int i = 0; i < num_intervals; ++i) {
		rand_human_interval(interval_len, chrom, &start, &end);
		string c(chrom);
		CHRPOS projected_start = (*_offsets)[c] + start;
		CHRPOS projected_end   = (*_offsets)[c] + end;
   
		ivl.start = projected_start + 1;
		ivl.end   = projected_end;
		_B->push_back(ivl);
    }
}
//}}}
//}}}
