/*****************************************************************************
interval_cuda.cu
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
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cutil.h>
#include "cudpp.h"
#include "interval.h"
#include "timer.h"
#include "interval_cuda.h"
#include "bsearch_cuda.h"
#include "bsearch_cuda.cu"
#include <gsl/gsl_statistics_int.h>

//{{{ void per_interval_count_intersections_bsearch_cuda(struct interval *A,
void per_interval_count_intersections_bsearch_cuda(struct interval *A,
												  unsigned int size_A,
												  struct interval *B,
												  unsigned int size_B,
												  unsigned int *R)
{
	cudaFree(NULL);

	unsigned int *A_starts_h, *A_lens_h, *B_starts_h, *B_ends_h;
	unsigned int *A_starts_d, *A_lens_d, *B_starts_d, *B_ends_d;
	unsigned int *R_d;
	allocate_and_move(A,
					 &A_starts_h,
					 &A_starts_d,
					 &A_lens_h ,
					 &A_lens_d,
					 size_A,
					 B,
					 &B_starts_h ,
					 &B_starts_d,
					 &B_ends_h ,
					 &B_ends_d,
					 size_B,
					 &R_d);

	cudpp_sort_by_key(B_starts_d, size_B);

	cudpp_sort_by_key(B_ends_d, size_B);

	bits_cuda(256,
			  1,
			  A_starts_d,
			  A_lens_d,
			  size_A,
			  B_starts_d,
			  B_ends_d,
			  size_B,
			  R_d);

	cudaMemcpy(R, R_d, size_A* sizeof(unsigned int), cudaMemcpyDeviceToHost);

	cudaThreadSynchronize();
	cudaError_t err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Result move: %s.\n", cudaGetErrorString( err) );


	cudaFree(A_starts_d);
	cudaFree(A_lens_d);
	cudaFree(B_starts_d);
	cudaFree(B_ends_d);
	cudaFree(R_d);
	free(A_starts_h);
	free(A_lens_h);
	free(B_starts_h);
	free(B_ends_h);
}
//}}}

//{{{ unsigned int count_intersections_bsearch_cuda(struct interval *A,
unsigned int count_intersections_bsearch_cuda(struct interval *A,
										      unsigned int size_A,
											  struct interval *B,
										      unsigned int size_B)
{
	//cudaFree(NULL);

	unsigned int *A_starts_h, *A_lens_h, *B_starts_h, *B_ends_h;
	unsigned int *A_starts_d, *A_lens_d, *B_starts_d, *B_ends_d;
	unsigned int *R_d;
	allocate_and_move(A,
					 &A_starts_h,
					 &A_starts_d,
					 &A_lens_h ,
					 &A_lens_d,
					 size_A,
					 B,
					 &B_starts_h ,
					 &B_starts_d,
					 &B_ends_h ,
					 &B_ends_d,
					 size_B,
					 &R_d);

	cudpp_sort_by_key(B_starts_d, size_B);

	cudpp_sort_by_key(B_ends_d, size_B);

	bits_cuda(256,
			  1,
			  A_starts_d,
			  A_lens_d,
			  size_A,
			  B_starts_d,
			  B_ends_d,
			  size_B,
			  R_d);

	unsigned int *Ro_d;
	cudaMalloc((void **)&Ro_d, sizeof(unsigned int));
	cudaMemset(Ro_d, 0, sizeof(unsigned int));

	cudpp_sum(R_d, Ro_d, size_A);

	unsigned int R;
	cudaMemcpy(&R, Ro_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);

	cudaThreadSynchronize();
	cudaError_t err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Result move: %s.\n", cudaGetErrorString( err) );

	cudaFree(A_starts_d);
	cudaFree(A_lens_d);
	cudaFree(B_starts_d);
	cudaFree(B_ends_d);
	cudaFree(R_d);
	cudaFree(Ro_d);
	free(A_starts_h);
	free(A_lens_h);
	free(B_starts_h);
	free(B_ends_h);

	return R;
}
//}}}

//{{{ unsigned int test_intersections_bsearch_cuda(struct interval *A,
unsigned int test_intersections_bsearch_cuda(struct interval *A,
										     unsigned int size_A,
											 struct interval *B,
										     unsigned int size_B,
											 unsigned int n,
											 unsigned int max_offset,
											 unsigned int *O,
											 double *mean,
											 double *sd,
											 double *p)
{
	cudaError_t err;
	int block_size = 256;
	dim3 dimBlock(block_size);
	int grid_size;

	//{{{ allocate and move A_starts, A_lens, B_starts, B_ends, B_lens
	struct timeval t_start = in();
	unsigned int *A_starts_h, *A_lens_h, *B_starts_h, *B_ends_h;
	unsigned int *A_starts_d, *A_lens_d, *B_starts_d, *B_ends_d;
	unsigned int *R_d;
	struct timeval start = in();
	allocate_and_move(A,
					&A_starts_h,
					&A_starts_d,
					&A_lens_h ,
					&A_lens_d,
					size_A,
					B,
					&B_starts_h ,
					&B_starts_d,
					&B_ends_h ,
					&B_ends_d,
					size_B,
					&R_d);
	unsigned int *Ro_d;
	cudaMalloc((void **)&Ro_d, sizeof(unsigned int));
	cudaMemset(Ro_d, 0, sizeof(unsigned int));
	unsigned int *B_lens_d;
	cudaMalloc((void **)&B_lens_d, size_B*sizeof(unsigned int));

	unsigned int i_mem_time = out(t_start);
	//fprintf(stderr, "i_mem:%u\t", i_mem_time);
	//}}}
	
	cudpp_sort_by_key_value( A_starts_d,
							 A_lens_d,
							 size_A);

	cudpp_sort_by_key_value( B_starts_d,
							 B_ends_d,
							 size_B);

	//{{{ set B_len, ordered by B_start
	grid_size = ( size_B + block_size - 1) / (block_size * 1);
	dim3 dimGridSet( grid_size );
	set_len_cuda <<<dimGridSet, dimBlock >>> (B_starts_d,
											  B_ends_d,
											  B_lens_d,
											  size_B,
											  1);
	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "set_len_cuda: %s.\n", cudaGetErrorString( err) );

	//}}}

	cudpp_sort_by_key(B_ends_d,
					  size_B);

	bits_cuda(256,
			  1,
			  A_starts_d,
			  A_lens_d,
			  size_A,
			  B_starts_d,
			  B_ends_d,
			  size_B,
			  R_d);
	
	unsigned int i_do_time = out(t_start);
	cudpp_sum(R_d, Ro_d, size_A);


	//{{{ Move R_d to R
	t_start = in();
	unsigned int R;
	cudaMemcpy(&R, Ro_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Result move: %s.\n", cudaGetErrorString( err) );

	i_mem_time += out(t_start);
	//fprintf(stderr, "i_mem:%u\t", i_mem_time);
	//}}}
	
	*O = R;

	// Generate random A_start and B_start
	srand(1);
	unsigned int *sims = (unsigned int *) malloc(n * sizeof(unsigned int));
	int i, N = 0;
	for(i = 0; i < n; i++) {

		//{{{  Rands for A_start and B_start
		cudpp_rand(A_starts_d, rand(), size_A);
		cudpp_rand(B_starts_d, rand(), size_B);

		unsigned int shift = 4294967295;

		grid_size = ( size_B + block_size - 1) / (block_size * 1);
		dim3 dimGridMapB( grid_size );
		map_list_cuda <<<dimGridMapB, dimBlock >>> (B_starts_d,
													size_B,
													shift,
													max_offset,
													1);
		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "map_list_cuda B: %s.\n", cudaGetErrorString( err) );

		grid_size = ( size_A + block_size - 1) / (block_size * 1);
		dim3 dimGridMapA( grid_size );
		map_list_cuda <<<dimGridMapA, dimBlock >>> (A_starts_d,
													size_A,
													shift,
													max_offset,
													1);
		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "map_list_cuda A: %s.\n", cudaGetErrorString( err) );
		//}}}
		
		cudpp_sort_by_key(A_starts_d, size_A);
		cudpp_sort_by_key(B_starts_d, size_B);
			
		//{{{ set B_end, ordered by B_start
		grid_size = ( size_B + block_size - 1) / (block_size * 1);
		//dimGridSet( grid_size );
		set_end_cuda <<<dimGridSet, dimBlock >>> (B_starts_d,
												  B_ends_d,
												  B_lens_d,
												  size_B,
												  1);
		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "map_list_cuda A: %s.\n", cudaGetErrorString( err) );

		//}}}

		bits_cuda(256,
				  1,
				  A_starts_d,
				  A_lens_d,
				  size_A,
				  B_starts_d,
				  B_ends_d,
				  size_B,
				  R_d);

		unsigned int i_do_time = out(t_start);
		cudpp_sum(R_d, Ro_d, size_A);

		//{{{ Move R_d to R
		t_start = in();
		cudaMemcpy(&R, Ro_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Result move: %s.\n", cudaGetErrorString( err) );

		i_mem_time += out(t_start);
		//fprintf(stderr, "i_mem:%u\t", i_mem_time);
		//}}}
		
		sims[i] = R;	
	
		if (R >= *O)
			N = N + 1;
	}

	*mean = gsl_stats_int_mean((const int*)sims, 1, n);
	*sd = gsl_stats_int_sd_m((const int*)sims, 1, n, *mean);
	*p =  ( (double) N + 1) / ( (double) n + 1);

	cudaFree(A_starts_d);
	cudaFree(A_lens_d);
	cudaFree(B_starts_d);
	cudaFree(B_ends_d);
	cudaFree(B_lens_d);
	cudaFree(R_d);
	cudaFree(Ro_d);
	free(A_starts_h);
	free(A_lens_h);
	free(B_starts_h);
	free(B_ends_h);
	free(sims);
	
	return 0;
}
//}}}

//{{{ unsigned int enumerate_intersections_bsearch_cuda(struct interval *A,
unsigned int enumerate_intersections_bsearch_cuda(struct interval *A,
										      unsigned int size_A,
											  struct interval *B,
										      unsigned int size_B,
											  unsigned int *R,
											  unsigned int **E)
{
	cudaFree(NULL);

	unsigned int *A_starts_h, *A_lens_h, *B_starts_h, *B_ends_h;
	unsigned int *A_starts_d, *A_lens_d, *B_starts_d, *B_ends_d;
	unsigned int *R_d;

	allocate_and_move(A,
					 &A_starts_h,
					 &A_starts_d,
					 &A_lens_h ,
					 &A_lens_d,
					 size_A,
					 B,
					 &B_starts_h ,
					 &B_starts_d,
					 &B_ends_h ,
					 &B_ends_d,
					 size_B,
					 &R_d);

	cudpp_sort_by_key(B_ends_d, size_B);

	// B_starts_id is the index of the interval in the orginial bed file, when
	// the intersecting intervals are enumerated, the B_starts_id value will be
	// stored so that the correct intervals can be displayed 
	unsigned int *B_starts_id_d;
	cudaMalloc((void **)&B_starts_id_d, size_B*sizeof(unsigned int));

	int block_size = 256;
	dim3 dimBlock(block_size);
	int grid_size = ( size_B + block_size - 1) / (block_size * 1);
	dim3 dimGridIdSet( grid_size );

	set_id_cuda <<<dimGridIdSet, dimBlock>>> (B_starts_id_d, size_B, 1);

	cudaThreadSynchronize();
	cudaError_t err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "set_id_cuda: %s.\n", cudaGetErrorString( err) );

	cudpp_sort_by_key_value(B_starts_d, B_starts_id_d, size_B);

	bits_cuda(256,
			  1,
			  A_starts_d,
			  A_lens_d,
			  size_A,
			  B_starts_d,
			  B_ends_d,
			  size_B,
			  R_d);

	unsigned int *Ro_d;
	cudaMalloc((void **)&Ro_d, size_A*sizeof(unsigned int));

	// Do a prefix scan of the results to be used to identify the offset of the
	// intervals within B
	cudpp_scan(R_d, Ro_d, size_A);

	cudaMemcpy(R, Ro_d , size_A*sizeof(unsigned int), cudaMemcpyDeviceToHost);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Result move: %s.\n", cudaGetErrorString( err) );
	
	unsigned int N = (R)[size_A - 1];

	cudaFree(R_d);

	cudaMemcpy(B_starts_d, B_starts_h, size_B * sizeof(unsigned int),
			cudaMemcpyHostToDevice);
	cudaMemcpy(B_ends_d, B_ends_h, size_B * sizeof(unsigned int),
			cudaMemcpyHostToDevice);

	cudpp_sort_by_key_value( B_starts_d,
							 B_ends_d,
							 size_B);

	unsigned int *E_d;
	cudaMalloc((void **)&E_d, N*sizeof(unsigned int));

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Malloc E: %s.\n", cudaGetErrorString( err) );

	grid_size = ( size_A + block_size - 1) / (block_size * 1);
	dim3 dimGridSearch( grid_size );

	enumerate_bsearch_cuda <<< dimGridSearch, dimBlock >>> (
			A_starts_d, A_lens_d, size_A,
			B_starts_d, B_ends_d, size_B,
			B_starts_id_d,
			Ro_d,
			E_d,
			1);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "enumerate_bsearch_cuda: %s.\n",
					cudaGetErrorString( err) );

	*E = (unsigned int *) malloc( N*sizeof(unsigned int) );
	cudaMemcpy(*E, E_d, N*sizeof(unsigned int), cudaMemcpyDeviceToHost);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Result move: %s.\n", cudaGetErrorString( err) );
	return N;
}
//}}}

//{{{void cudpp_rand_init(CUDPPHandle *rand_cudpp,
void cudpp_rand_init(CUDPPHandle *rand_cudpp,
					 CUDPPHandle *rand_plan,
					 unsigned int seed,
					 unsigned int size)
{
	cudppCreate(rand_cudpp);

	CUDPPConfiguration rand_config;
	rand_config.op = CUDPP_ADD;
	rand_config.datatype = CUDPP_UINT;
	rand_config.algorithm = CUDPP_RAND_MD5;
	rand_config.options = 0;

	CUDPPResult res = cudppPlan(*rand_cudpp, 
								rand_plan,
								rand_config,
								size,
								1,
								0);;
	if (CUDPP_SUCCESS != res) {
		printf("Error creating rand CUDPPPlan\n");
        exit(-1);
	}

	res = cudppRandSeed(*rand_plan, seed);
	if (CUDPP_SUCCESS != res) {
		printf("Error in cudppRandSeed\n");
        exit(-1);
	}

}
//}}}	

//{{{ void cudpp_rand(unsigned int *keys_d,
void cudpp_planned_rand(CUDPPHandle *rand_plan,
						unsigned int *out_d,
						unsigned int size)
{

	CUDPPResult res = cudppRand(*rand_plan, out_d, size);
    if (CUDPP_SUCCESS != res) {
        printf("Error in cudppRand()\n");
        exit(-1);
    }
}
//}}}

//{{{ void cudpp_rand(unsigned int *keys_d,
void cudpp_rand(unsigned int *out_d,
				unsigned int seed,
				unsigned int size)
{

	CUDPPHandle rand_cudpp;
	cudppCreate(&rand_cudpp);

	CUDPPConfiguration rand_config;
	rand_config.op = CUDPP_ADD;
	rand_config.datatype = CUDPP_UINT;
	rand_config.algorithm = CUDPP_RAND_MD5;
	rand_config.options = 0;

	CUDPPHandle rand_plan= 0;
	CUDPPResult res = cudppPlan(rand_cudpp, 
								&rand_plan,
								rand_config,
								size,
								1,
								0);;

	if (CUDPP_SUCCESS != res) {
		printf("Error creating rand CUDPPPlan\n");
        exit(-1);
	}

	res = cudppRandSeed(rand_plan, seed);
	if (CUDPP_SUCCESS != res) {
		printf("Error in cudppRandSeed\n");
        exit(-1);
	}

	res = cudppRand(rand_plan, out_d, size);
    if (CUDPP_SUCCESS != res) {
        printf("Error in cudppRand()\n");
        exit(-1);
    }

	res = cudppDestroy(rand_cudpp);
    if (CUDPP_SUCCESS != res) {
		printf("Error shutting down CUDPP Library.\n");
        exit(-1);
    }
}
//}}}

//{{{void cudpp_sort_by_key( unsigned int *keys_d,
void cudpp_sort_by_key( unsigned int *keys_d,
						unsigned int size)
{
	CUDPPHandle sort_cudpp;
	cudppCreate(&sort_cudpp);

	CUDPPConfiguration sort_config;
	sort_config.datatype = CUDPP_UINT;
	sort_config.algorithm = CUDPP_SORT_RADIX;
	sort_config.options = CUDPP_OPTION_KEYS_ONLY;

	CUDPPHandle sort_plan= 0;
	CUDPPResult res = cudppPlan(sort_cudpp, 
								&sort_plan,
								sort_config,
								size,
								1,
								0);;

	if (CUDPP_SUCCESS != res) {
		printf("Error creating sort CUDPPPlan\n");
        exit(-1);
	}

    res = cudppSort(sort_plan, keys_d, NULL, size);
    if (CUDPP_SUCCESS != res) {
        printf("Error in cudppSort()\n");
        exit(-1);
    }

	res = cudppDestroy(sort_cudpp);
    if (CUDPP_SUCCESS != res) {
		printf("Error shutting down CUDPP Library.\n");
        exit(-1);
    }
}
//}}}

//{{{void cudpp_sort_by_key_value( unsigned int *keys_d,
void cudpp_sort_by_key_value(unsigned int *keys_d,
							 unsigned int *values_d,
							 unsigned int size)
{
	CUDPPHandle sort_cudpp;
	cudppCreate(&sort_cudpp);

	CUDPPConfiguration sort_config;
	sort_config.datatype = CUDPP_UINT;
	sort_config.algorithm = CUDPP_SORT_RADIX;
	sort_config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;

	CUDPPHandle sort_plan= 0;
	CUDPPResult res = cudppPlan(sort_cudpp, 
								&sort_plan,
								sort_config,
								size,
								1,
								0);;

	if (CUDPP_SUCCESS != res) {
		printf("Error creating sort CUDPPPlan\n");
        exit(-1);
	}

    res = cudppSort(sort_plan, keys_d, values_d, size);
    if (CUDPP_SUCCESS != res) {
        printf("Error in cudppSort()\n");
        exit(-1);
    }

	res = cudppDestroy(sort_cudpp);
    if (CUDPP_SUCCESS != res) {
		printf("Error shutting down CUDPP Library.\n");
        exit(-1);
    }
}
//}}}

//{{{ void cudpp_sum(unsigned int *list_d,
void cudpp_sum(unsigned int *list_d,
			   unsigned int *sum_d,
			   unsigned int size)
{
	CUDPPHandle sum_cudpp;
	cudppCreate(&sum_cudpp);

	CUDPPConfiguration sum_config;
	sum_config.datatype = CUDPP_UINT;
	sum_config.algorithm = CUDPP_REDUCE;
	sum_config.options = 0;

	CUDPPHandle sum_plan = 0;
	CUDPPResult res = cudppPlan(sum_cudpp, 
								&sum_plan,
								sum_config,
								size,
								1,
								0);

	if (CUDPP_SUCCESS != res) {
		printf("Error creating sort CUDPPPlan\n");
        exit(-1);
	}

    res = cudppReduce(sum_plan, sum_d, list_d, size);
    if (CUDPP_SUCCESS != res) {
        printf("Error in cudppReduce() R_d\n");
        exit(-1);
    }

	res = cudppDestroy(sum_plan);
    if (CUDPP_SUCCESS != res) {
		printf("Error shutting down CUDPP Library.\n");
        exit(-1);
    }
}
//}}}

//{{{void cudpp_scan(unsigned int *R_d,
void cudpp_scan(unsigned int *R_d,
			    unsigned int *Ro_d,
			    unsigned int size)
{
	CUDPPHandle scan_cudpp;
	cudppCreate(&scan_cudpp);

	CUDPPConfiguration scan_config;
	scan_config.op = CUDPP_ADD;
	scan_config.datatype = CUDPP_UINT;
	scan_config.algorithm = CUDPP_SCAN;
	scan_config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_INCLUSIVE;

	//unsigned int *Ro_d;
	//cudaMalloc((void **)&Ro_d, size_A*sizeof(unsigned int));

	CUDPPHandle scan_plan = 0;
	CUDPPResult res = cudppPlan(scan_cudpp, 
								&scan_plan,
								scan_config,
								size,
								1,
								0);

	if (CUDPP_SUCCESS != res) {
		printf("Error creating scan CUDPPPlan\n");
        exit(-1);
	}

    res = cudppScan(scan_plan, Ro_d, R_d, size);
    if (CUDPP_SUCCESS != res) {
        printf("Error in cudppScan() R_d\n");
        exit(-1);
    }

	res = cudppDestroy(scan_plan);
    if (CUDPP_SUCCESS != res) {
		printf("Error shutting down CUDPP Library.\n");
        exit(-1);
    }
}
//}}}	

//{{{ unsigned int count_intersections_sort_bsearch_cuda(struct interval *A,
unsigned int count_intersections_sort_bsearch_cuda(struct interval *A,
										      unsigned int size_A,
											  struct interval *B,
										      unsigned int size_B)
{
	int block_size = 256;
	dim3 dimBlock(block_size);
	int grid_size = ( size_A + block_size - 1) / (block_size * 1);
	dim3 dimGridSearch( grid_size );
	cudaError_t err;

	start(); //data_prep_time
	//{{{ Allocate and move 
	unsigned int *A_starts_h, *A_lens_h, *B_starts_h, *B_ends_h;
	unsigned int *A_starts_d, *A_lens_d, *B_starts_d, *B_ends_d;
	unsigned int *R_d;
	allocate_and_move(A,
					&A_starts_h,
					&A_starts_d,
					&A_lens_h ,
					&A_lens_d,
					size_A,
					B,
					&B_starts_h ,
					&B_starts_d,
					&B_ends_h ,
					&B_ends_d,
					size_B,
					&R_d);
	//}}}
	stop(); //data_prep_time
	unsigned long data_prep_time = report();

	start(); //sort_time
	//{{{ Sort B_starts and B_ends
	// Sort B by start
	//nvRadixSort::RadixSort radixsortB_starts(size_B, true);
	//radixsortB_starts.sort((unsigned int*)B_starts_d, 0, size_B, 32);

	//cudaThreadSynchronize();
	//err = cudaGetLastError();
	//if(err != cudaSuccess)
		//fprintf(stderr, "Sort B_starts: %s.\n", cudaGetErrorString( err) );

	// Sort B by end
	//nvRadixSort::RadixSort radixsortB_ends(size_B, true);
	//radixsortB_ends.sort((unsigned int*)B_ends_d, 0, size_B, 32);

	//cudaThreadSynchronize();
	//err = cudaGetLastError();
	//if(err != cudaSuccess)
		//fprintf(stderr, "Sort B_ends: %s.\n", cudaGetErrorString( err) );
	//}}}
	stop();	//sort_time
	unsigned long sort_time = report();
	
	start();
	//{{{ Compute and count intersections
	count_bsearch_cuda <<<dimGridSearch, dimBlock >>> (
			A_starts_d, A_lens_d, size_A,
			B_starts_d, B_ends_d, size_B,
			R_d,
			1);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Sort B_ends: %s.\n", cudaGetErrorString( err) );

	//parallel_sum(R_d, block_size, size_A, 1024);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Parallel sum: %s.\n", cudaGetErrorString( err) );


	unsigned int R;
	cudaMemcpy(&R, R_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Result move: %s.\n", cudaGetErrorString( err) );

	//}}}
	stop(); //intersect_time
	/*
	unsigned long intersect_time = report();

	unsigned long total_time = data_prep_time + 
							   sort_time +
							   pre_sort_time +
							   intersect_time;
	printf("sort\t"
		   "total:%lu\t"
		   "prep:%lu,%f\t"
		   "sort:%lu,%f\t"
		   "presort:%lu,%f\t"
		   "intersect:%lu,%f\n",
		   total_time,
		   data_prep_time,  (double)data_prep_time / (double)total_time,
		   sort_time, (double)sort_time / (double)total_time,
		   pre_sort_time, (double)pre_sort_time / (double)total_time,
		   intersect_time, (double)intersect_time / (double)total_time);
	*/
	cudaFree(A_starts_d);
	cudaFree(A_lens_d);
	cudaFree(B_starts_d);
	cudaFree(B_ends_d);
	cudaFree(R_d);

	return R;
}
//}}}

//{{{ unsigned int count_intersections_i_bsearch_cuda(struct interval *A,
unsigned int count_intersections_i_gm_bsearch_cuda(struct interval *A,
										      unsigned int size_A,
											  struct interval *B,
										      unsigned int size_B,
										      unsigned int size_I)
{
	int block_size = 256;
	dim3 dimBlock(block_size);
	int grid_size = ( size_A + block_size - 1) / (block_size);
	dim3 dimGridSearch( grid_size );
	cudaError_t err;

	start(); //data_prep_time

	//{{{ Allocate and move 
	unsigned int *A_starts_h, *A_lens_h, *B_starts_h, *B_ends_h;
	unsigned int *A_starts_d, *A_lens_d, *B_starts_d, *B_ends_d;
	unsigned int *R_d;
	allocate_and_move(A,
					&A_starts_h,
					&A_starts_d,
					&A_lens_h ,
					&A_lens_d,
					size_A,
					B,
					&B_starts_h ,
					&B_starts_d,
					&B_ends_h ,
					&B_ends_d,
					size_B,
					&R_d);
	//}}}
	//
	stop(); //data_prep_time
	unsigned long data_prep_time = report();

	start();//sort_time
	//{{{ Sort B_starts and B_ends
	// Sort B by start
	//nvRadixSort::RadixSort radixsortB_starts(size_B, true);
	//radixsortB_starts.sort((unsigned int*)B_starts_d, 0, size_B, 32);

	//cudaThreadSynchronize();
	//err = cudaGetLastError();
	//if(err != cudaSuccess)
		//fprintf(stderr, "Sort B_starts: %s.\n", cudaGetErrorString( err) );

	// Sort B by end
	//nvRadixSort::RadixSort radixsortB_ends(size_B, true);
	//radixsortB_ends.sort((unsigned int*)B_ends_d, 0, size_B, 32);

	//cudaThreadSynchronize();
	//err = cudaGetLastError();
	//if(err != cudaSuccess)
		//fprintf(stderr, "Sort B_ends: %s.\n", cudaGetErrorString( err) );
	//}}}
	stop();	//sort_time
	unsigned long sort_time = report();

	start();//index_time
	//{{{ Generate index
	unsigned int *I_starts_d, *I_ends_d;
	cudaMalloc((void **)&I_starts_d, (size_I)*sizeof(unsigned int));
	cudaMalloc((void **)&I_ends_d, (size_I)*sizeof(unsigned int));

	int index_grid_size = ( size_I + block_size - 1) / (block_size);
	dim3 index_dimGrid( index_grid_size );

	gen_index <<<index_dimGrid, dimBlock>>> ( B_starts_d, size_B, I_starts_d, size_I);
	gen_index <<<index_dimGrid, dimBlock>>> ( B_ends_d, size_B, I_ends_d, size_I);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Count i bsearch: %s.\n", cudaGetErrorString( err) );
	//}}}
	stop();	//index_time
	unsigned long index_time = report();

	//{{{ Compute and count intersections
	count_i_gm_bsearch_cuda <<<dimGridSearch, dimBlock >>> (
			A_starts_d, A_lens_d, size_A,
			B_starts_d, B_ends_d, size_B,
			I_starts_d, I_ends_d, size_I,
			R_d,
			1);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Count i bsearch: %s.\n", cudaGetErrorString( err) );

	//parallel_sum(R_d, block_size, size_A, 1024);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Parallel sum: %s.\n", cudaGetErrorString( err) );


	unsigned int R;
	cudaMemcpy(&R, R_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Result move: %s.\n", cudaGetErrorString( err) );

	//}}}
	stop(); //intersect_time
	unsigned long intersect_time = report();

	unsigned long total_time = data_prep_time + 
							   sort_time +
							   index_time +
							   intersect_time;
	printf("index gm\t"
		   "total:%lu\t"
		   "prep:%lu,%f\t"
		   "sort:%lu,%f\t"
		   "index:%lu,%f\t"
		   "intersect:%lu,%f\n",
		   total_time,
		   data_prep_time,  (double)data_prep_time / (double)total_time,
		   sort_time, (double)sort_time / (double)total_time,
		   index_time, (double)index_time / (double)total_time,
		   intersect_time, (double)intersect_time / (double)total_time);

	cudaFree(A_starts_d);
	cudaFree(A_lens_d);
	cudaFree(B_starts_d);
	cudaFree(B_ends_d);
	cudaFree(I_starts_d);
	cudaFree(I_ends_d);
	cudaFree(R_d);

	return R;
}
//}}}

//{{{ __global__ void count_bsearch_cuda (	unsigned int *A_start,
/*
 * @param A_start list of start positions to query, does not need to be sorted
 * @param A_len list of lengths that correspond to A_start
 * @param A_size size of A_start and A_len
 * @param B_start list of sorted start positions to be queried
 * @param B_end list of sorted end positions to be queired 
 * @param B_size size of B_start and B_end
 * @param R number of intersections for each interval in A
 * @param n number of intervals per thread
 */
__global__
void count_bsearch_cuda (	unsigned int *A_start,
							unsigned int *A_len,
							int A_size,
							unsigned int *B_start,
							unsigned int *B_end,
							int B_size,
							unsigned int *R,
							int n)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	while ( i < (n * grid_size) ) {

		if (i < A_size) {
			unsigned int start = A_start[i];
			unsigned int end = start + A_len[i];

			int cant_before = bound_binary_search(B_end,
														   B_size,
														   start,
														   -1,
														   B_size);

			int cant_after = bound_binary_search(B_start,
														  B_size,
														  end,
														  -1,
														  B_size);

			while ( end == B_start[cant_after] )
				++cant_after;

			cant_after = A_size - cant_after;	

			R[i] = A_size - cant_before - cant_after;
		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void count_i_gm_bsearch_cuda (	unsigned int *A_start,
/*
 * @param A_start list of start positions to query, does not need to be sorted
 * @param A_len list of lengths that correspond to A_start
 * @param A_size size of A_start and A_len
 * @param B_start list of sorted start positions to be queried
 * @param B_end list of sorted end positions to be queired 
 * @param B_size size of B_start and B_end
 * @param R number of intersections for each interval in A
 * @param n number of intervals per thread
 */
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
							int n)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	while ( i < (n * grid_size) ) {
		if (i < A_size) {
			unsigned int start = A_start[i];
			unsigned int end = start + A_len[i];

			int cant_before = i_binary_search(B_end,
											  B_size,
											  start,
											  I_end,
											  I_size);
	
			int cant_after = i_binary_search(B_start,
											 B_size,
											 end,
											 I_start,
											 I_size);

			while ( end == B_start[cant_after] )
				++cant_after;

			cant_after = A_size - cant_after;	

			R[i] = A_size - cant_before - cant_after;
		}
		i += grid_size;
	}
}
//}}}

//{{{void allocate_and_move( struct interval *A,
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

						unsigned int **R_d)
{

	cudaError_t err;
	//{{{ Move intervals to unsigned int arrays
	*A_starts_h = (unsigned int *) malloc( (size_A) * sizeof(unsigned int));
	*A_lens_h = (unsigned int *) malloc( (size_A) * sizeof(unsigned int));

	*B_starts_h = (unsigned int *) malloc( (size_B) * sizeof(unsigned int));
	*B_ends_h = (unsigned int *) malloc( (size_B) * sizeof(unsigned int));

	int i;
	for (i = 0; i < size_B; i++) {
		(*B_starts_h)[i] = B[i].start;
		(*B_ends_h)[i] = B[i].end;
	}

	for (i = 0; i < size_A; i++) {
		(*A_starts_h)[i] = A[i].start;
		(*A_lens_h)[i] = A[i].end - A[i].start;
	}
	//}}}

	//{{{ Move inteval arrays to device
	cudaMalloc((void **)A_starts_d, (size_A)*sizeof(unsigned int));
	cudaMalloc((void **)A_lens_d, (size_A)*sizeof(unsigned int));
	cudaMalloc((void **)B_starts_d, (size_B)*sizeof(unsigned int));
	cudaMalloc((void **)B_ends_d, (size_B)*sizeof(unsigned int));

	cudaMemcpy(*A_starts_d, *A_starts_h, (size_A) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(*A_lens_d, *A_lens_h, (size_A) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);
	cudaMemcpy(*B_starts_d, *B_starts_h, (size_B) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(*B_ends_d, *B_ends_h, (size_B) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Interval move: %s.\n", cudaGetErrorString( err) );
	//}}}
	
	//{{{ Alocate space for result on device
	cudaMalloc((void **)R_d, (size_A)*sizeof(unsigned int));
	unsigned long memup_time = report();

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "R_d malloc: %s.\n", cudaGetErrorString( err) );
	//}}}
}
//}}}

//{{{ __global__ void enumerate_bsearch_cuda (	unsigned int *A_start,
/*
 * @param A_start list of start positions to query, does not need to be sorted
 * @param A_len list of lengths that correspond to A_start
 * @param A_size size of A_start and A_len
 * @param B_start list of sorted start positions to be queried
 * @param B_end list of sorted end positions to be queired 
 * @param B_size size of B_start and B_end
 * @param R number of intersections for each interval in A
 * @param n number of intervals per thread
 */
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
							 int n)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	while ( i < (n * grid_size) ) {
		if (i < A_size) {
			unsigned int A_i_start = A_start[i];
			unsigned int A_i_end = A_i_start + A_len[i];

			unsigned int start= 0, end;

			if (i != 0)
				start = R[i - 1];

			end = R[i];

			if (end - start > 0) {
				unsigned int from = bound_binary_search(B_start,
													    B_size,
													    A_i_end,
													    -1,
													    B_size);

				while ( ( B_start[from] == A_i_end) && (from < B_size))
					++from;

				while (  (end - start) > 0 ) {
					if ( (A_i_start <= B_end[from]) && 
							(A_i_end >= B_start[from]) ) {
						 E[start] = B_starts_id_d[from];
						 ++start;
					}

					--from;
				}
			}
		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void set_len_cuda (	unsigned int *start,
__global__
void set_len_cuda (	unsigned int *start,
					unsigned int *end,
					unsigned int *len,
					int size,
					int n)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	while ( i < (n * grid_size) ) {
		if (i < size) {
			len[i] = end[i] - start[i];
		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void set_end_cuda (	unsigned int *start,
__global__
void set_end_cuda (	unsigned int *start,
					unsigned int *end,
					unsigned int *len,
					int size,
					int n)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	while ( i < (n * grid_size) ) {
		if (i < size) {
			end[i] = start[i] + len[i];
		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void map_list_cuda (unsigned int *list,
__global__
void map_list_cuda (unsigned int *list,
					unsigned int size,
					unsigned int old_max,
					unsigned int new_max,
					int n)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;
	float norm = (float) new_max/(float) old_max;

	while ( i < (n * grid_size) ) {
		if (i < size) {
			list[i] = (unsigned int)(list[i] * norm);
		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void set_len_cuda (	unsigned int *start,
__global__
void set_id_cuda (unsigned int *list,
				  int size,
				  int n)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	while ( i < (n * grid_size) ) {
		if (i < size) {
			list[i] = i;
		}
		i += grid_size;
	}
}
//}}}

//{{{void bits_cuda(int block_size,
void bits_cuda(int block_size,
			   unsigned int per_thread,
			   unsigned int *A_starts_d,
			   unsigned int *A_lens_d,
			   unsigned int size_A,
			   unsigned int *B_starts_d,
			   unsigned int *B_ends_d,
			   unsigned int size_B,
			   unsigned int *R_d) 
{
	dim3 dimBlock(block_size);
	int grid_size = ( size_A + block_size - 1) / (block_size * 1);
	dim3 dimGridSearch( grid_size );

	count_bsearch_cuda <<< dimGridSearch, dimBlock >>> (
			A_starts_d, A_lens_d, size_A,
			B_starts_d, B_ends_d, size_B,
			R_d,
			1);

	cudaThreadSynchronize();
	cudaError_t err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "count_bsearch_cuda: %s.\n", cudaGetErrorString( err) );
}
//}}}

//{{{ void cuda_free() {
void cuda_free() {
	cudaFree(NULL);
}
//}}}
