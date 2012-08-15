/*****************************************************************************
interval_cuda.h
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
#ifndef __INTERVAL_CUDA_H__
#define __INTERVAL_CUDA_H__
#include <cudpp.h>


void per_interval_count_intersections_bsearch_cuda(struct interval *A,
												  unsigned int size_A,
												  struct interval *B,
												  unsigned int size_B,
												  unsigned int *R);

unsigned int count_intersections_bsearch_cuda(struct interval *A,
										      unsigned int size_A,
											  struct interval *B,
										      unsigned int size_B);

unsigned int count_intersections_sort_bsearch_cuda(struct interval *A,
										      unsigned int size_A,
											  struct interval *B,
										      unsigned int size_B);

unsigned int count_intersections_i_gm_bsearch_cuda(struct interval *A,
										      unsigned int size_A,
											  struct interval *B,
										      unsigned int size_B,
										      unsigned int size_I);

void allocate_and_move( struct interval *A,
						unsigned int **A_starts_h,
						unsigned int **A_starts_d,
					   	unsigned int **A_lens_h ,
						unsigned int **A_lens_d,
						unsigned int size_A,
						struct interval *B,
						unsigned int **B_starts_h ,
						unsigned int **B_starts_d,
						unsigned int **B_ends_h ,
						unsigned int **B_ends_d,
						unsigned int size_B,
						unsigned int **R_d);

__global__
void count_bsearch_cuda (	unsigned int *A_start,
							unsigned int *A_len,
							int A_size,
							unsigned int *B_start,
							unsigned int *B_end,
							int B_size,
							unsigned int *R,
							int n);
__global__
void count_i_gm_bsearch_cuda (	unsigned int *A_start,
							unsigned int *A_len,
							int A_size,
							unsigned int *B_start,
							unsigned int *B_end,
							int B_size,
							unsigned int *I_start,
							unsigned int *I_end,
							int I_size,
							unsigned int *R,
							int n);

__global__
void enumerate_bsearch_cuda (unsigned int *A_start,
							unsigned int *A_len,
							int A_size,
							unsigned int *B_start,
							unsigned int *B_end,
							int B_size,
							unsigned int *B_starts_id_d,
							unsigned int *R,
							unsigned int *E,
							int n);

unsigned int enumerate_intersections_bsearch_cuda(struct interval *A,
										      unsigned int size_A,
											  struct interval *B,
										      unsigned int size_B,
											  unsigned int *R,
											  unsigned int **E);

unsigned int test_intersections_bsearch_cuda(struct interval *A,
										     unsigned int size_A,
											 struct interval *B,
										     unsigned int size_B,
											 unsigned int n,
											 unsigned int max_offset,
											 unsigned int *O,
											 double *mean,
											 double *sd,
											 double *p);

__global__
void set_len_cuda (	unsigned int *start,
					unsigned int *end,
					unsigned int *len,
					int size,
					int n);

__global__
void set_end_cuda (	unsigned int *start,
					unsigned int *end,
					unsigned int *len,
					int size,
					int n);

void cudpp_sort_by_key( unsigned int *keys_d,
						unsigned int size);

void cudpp_sort_by_key_value( unsigned int *keys_d,
							  unsigned int *values_d,
							  unsigned int size);

void cudpp_sum(unsigned int *list_d,
			   unsigned int *sum_d,
			   unsigned int size);

void cudpp_rand(unsigned int *out_d,
				unsigned int seed,
				unsigned int size);
	
void cudpp_planned_rand(CUDPPHandle *rand_plan,
						unsigned int *out_d,
						unsigned int size);

void cudpp_rand_init(CUDPPHandle *rand_cudpp,
					 CUDPPHandle *rand_plan,
					 unsigned int seed,
					 unsigned int size);

__global__
void map_list_cuda (unsigned int *list,
					unsigned int size,
					unsigned int old_max,
					unsigned int new_max,
					int n);

void bits_cuda(int block_size,
			   unsigned int per_thread,
			   unsigned int *A_starts_d,
			   unsigned int *A_lens_d,
			   unsigned int size_A,
			   unsigned int *B_starts_d,
			   unsigned int *B_ends_d,
			   unsigned int size_B,
			   unsigned int *R_d);

void cudpp_scan(unsigned int *R_d,
			    unsigned int *Ro_d,
			    unsigned int size);
__global__
void set_id_cuda (	unsigned int *list,
					int size,
					int n);

void cuda_free();
#endif
