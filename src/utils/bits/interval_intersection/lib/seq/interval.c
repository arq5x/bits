/*****************************************************************************
interval.c
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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "interval.h"
#include "bsearch.h"
#include "mt.h"
#include "timer.h"
#include <gsl/gsl_statistics_int.h>

//{{{ int compare_interval_by_start (const void *a, const void *b)
int compare_interval_by_start (const void *a, const void *b)
{  
	struct interval *a_i = (struct interval *)a;
	struct interval *b_i = (struct interval *)b;
	if (a_i->start < b_i->start)
		return -1;
	else if (a_i->start > b_i->start)
		return 1;
	else
		return 0;
}
//}}}

//{{{ void enumerate_intersections_bsearch_seq(struct interval *A,

/**
  * @param A intervals in set A
  * @param size_A size of set A
  * @param B intervals in set B
  * @param size_B size of set B
  * @param R prefix sum of the intersection between A and B, A[0] is the number
  * of intervals in B that intersect A[0], A[1] is A[0] + the number of
  * intervals in B that intersect A[1], and so on
  * @param E array that will hold the enumberated interval intersections
  * @param size_E 
  */
unsigned int enumerate_intersections_bsearch_seq(struct interval *A,
												 unsigned int size_A,
												 struct interval *B,
												 unsigned int size_B,
												 unsigned int **R,
												 unsigned int **E)
{

	//{{{ i_do
	*R = (unsigned int *) malloc(size_A * sizeof(unsigned int));

	unsigned int O = per_interval_count_intersections_bsearch_seq(A,
																  size_A,
																  B,
																  size_B,
																  *R);
	//}}} i_do
	
	//{{{ e_scan
	int i;
	// Scan R, the resulting array will define the location within E of the
	// enumerated intersecting intervals for each interval in A 
	// A[i] is allocated E[ R[i-1] ] ... E[ R[i] ]
	for (i = 1; i < size_A; i++)
		(*R)[i] = (*R)[i] + (*R)[i-1];
	//}}}

	//{{{ e_sort
	for (i = 0; i < size_B; i++)
		B[i].order = i;

	qsort(B, size_B, sizeof(struct interval), compare_interval_by_start); 
	//}}} e_sort

	//{{{ e_mem
	unsigned int *B_starts =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));
	for (i = 0; i < size_B; i++)
		B_starts[i] = B[i].start;

	*E = (unsigned int *) malloc(O * sizeof(unsigned int));
	//}}} e_mem

	//{{{ e_do
	unsigned int start = 0, end;
	for (i = 0; i < size_A; i++) {
		if (i != 0)
			start = (*R)[i - 1];

		end = (*R)[i];
		if (end - start > 0) {
			unsigned int from = bsearch_seq(A[i].end,
											B_starts,
											size_B,
											-1,
											size_B);

			while ( ( B_starts[from] == A[i].end) && from < size_B)
				++from;

			while (  (end - start) > 0 ) {
				if ( (A[i].start <= B[from].end) &&
					 (A[i].end >= B[from].start) ) {
					(*E)[start] = B[from].order;
					start++;
				} 

				--from;
			}
		}
	}
	//}}} e_do

	//{{{ e_mem
	free(B_starts);
	//}}} e_mem
	
	return O;
}
//}}}

//{{{ unsigned int per_interval_count_intersections_bsearch_seq(struct interval *A,
unsigned int per_interval_count_intersections_bsearch_seq(struct interval *A,
														   unsigned int size_A,
														   struct interval *B,
														   unsigned int size_B,
														   unsigned int *R)
{

	unsigned int i, O = 0;


	//{{{ i_mem
	unsigned int *B_starts =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));
	unsigned int *B_ends =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));

	for (i = 0; i < size_B; i++) {
		B_starts[i] = B[i].start;
		B_ends[i] = B[i].end;
	}
	//}}} i_mem

	//{{{ i_sort
	qsort(B_starts, size_B, sizeof(unsigned int), compare_unsigned_int); 
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 
	//}}} i_sort

	//{{{ i_do
	for (i = 0; i < size_A; i++) {
		unsigned int num_cant_before = bsearch_seq(A[i].start,
												   B_ends,
												   size_B,
												   -1,
												   size_B);
		unsigned int b = bsearch_seq(A[i].end,
								     B_starts,
								     size_B,
								     -1,
								     size_B);

		while ( ( B_starts[b] == A[i].end) && b < size_B)
			++b;

		unsigned int num_cant_after = size_B - b;

		unsigned int num_left = size_B - num_cant_before - num_cant_after;
		O += num_left;
		R[i] = num_left;

	}
	//}}} i_do

	//{{{ i_mem
	free(B_starts);
	free(B_ends);
	//}}}


	return O;
}
//}}}

//{{{ unsigned int per_interval_count_intersections_bsearch_seq_fast(
unsigned int per_interval_count_intersections_bsearch_seq_sorted_int(
														  struct interval *A,
														  unsigned int size_A,
														  unsigned int *B_start,
														  unsigned int *B_end,
														  unsigned int size_B,
														  unsigned int *R)
{

	unsigned int i, O = 0;

	for (i = 0; i < size_A; i++) {
		unsigned int num_cant_before = bsearch_seq(A[i].start,
												   B_end,
												   size_B,
												   -1,
												   size_B);
		unsigned int b = bsearch_seq(A[i].end,
								     B_start,
								     size_B,
								     -1,
								     size_B);

		while ( ( B_start[b] == A[i].end) && b < size_B)
			++b;

		unsigned int num_cant_after = size_B - b;

		unsigned int num_left = size_B - num_cant_before - num_cant_after;
		O += num_left;
		R[i] = num_left;
	}
	return O;
}
//}}}

//{{{ unsigned int count_intersections_bsearch_seq(struct interval *A,
unsigned int count_intersections_bsearch_seq(struct interval *A,
										     unsigned int size_A,
											 struct interval *B,
										     unsigned int size_B)
{

	unsigned int i, O = 0;

	//{{{ i_mem
	unsigned int *B_starts =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));
	unsigned int *B_ends =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));

	for (i = 0; i < size_B; i++) {
		B_starts[i] = B[i].start;
		B_ends[i] = B[i].end;
	}
	//}}} i_mem

	//{{{ i_sort
	qsort(B_starts, size_B, sizeof(unsigned int), compare_unsigned_int); 
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 
	//}}} i_sort

	//{{{ i_do
	O = count_seq(A,
				  size_A,
				  B_starts,
				  B_ends,
				  size_B);
	//}}} i_do

	//{{{ i_mem
	free(B_starts);
	free(B_ends);
	//}}} i_mem
	return O;
}
//}}}

//{{{ unsigned int count_intersections_i_bsearch_seq(struct interval *A,
unsigned int count_intersections_i_bsearch_seq(struct interval *A,
										       unsigned int size_A,
											   struct interval *B,
										       unsigned int size_B,
											   unsigned int size_I)
{

	unsigned int i, O = 0;

	unsigned int *B_starts =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));
	unsigned int *B_ends =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));

	for (i = 0; i < size_B; i++) {
		B_starts[i] = B[i].start;
		B_ends[i] = B[i].end;
	}

	qsort(B_starts, size_B, sizeof(unsigned int), compare_unsigned_int); 
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 

	unsigned int *I_starts = (unsigned int *)
			malloc(size_I * sizeof(unsigned int));
	unsigned int *I_ends = (unsigned int *)
			malloc(size_I * sizeof(unsigned int));

	create_index(B_starts, size_B, I_starts, size_I);
	create_index(B_ends, size_B, I_ends, size_I);

	for (i = 0; i < size_A; i++) {
		unsigned int num_cant_before = i_bsearch_seq(A[i].start,
												B_ends,
												size_B,
												I_ends,
												size_I);
		unsigned int b = i_bsearch_seq(A[i].end,
									   B_starts,
									   size_B,
									   I_starts,
									   size_I);

		while ( ( B_starts[b] == A[i].end) && b < size_B)
			++b;

		unsigned int num_cant_after = size_B - b;

		unsigned int num_left = size_B - num_cant_before - num_cant_after;

		O += num_left;
	}

	free(B_starts);
	free(B_ends);
	return O;
}
//}}}

//{{{ unsigned int count_intersections_t_bsearch_seq(struct interval *A,
unsigned int count_intersections_t_bsearch_seq(struct interval *A,
										       unsigned int size_A,
											   struct interval *B,
										       unsigned int size_B,
											   unsigned int size_T)
{

	unsigned int i, O = 0;

	unsigned int *B_starts =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));
	unsigned int *B_ends =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));

	for (i = 0; i < size_B; i++) {
		B_starts[i] = B[i].start;
		B_ends[i]   = B[i].end;
	}

	qsort(B_starts, size_B, sizeof(unsigned int), compare_unsigned_int); 
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 

	unsigned int *T_starts = (unsigned int *)
			malloc(size_T * sizeof(unsigned int));
	unsigned int *T_ends = (unsigned int *)
			malloc(size_T * sizeof(unsigned int));

	create_tree(B_starts, size_B, T_starts, size_T);
	create_tree(B_ends,   size_B, T_ends,   size_T);

	for (i = 0; i < size_A; i++) {

		unsigned int num_cant_before = t_bsearch_seq(A[i].start,
												B_ends,
												size_B,
												T_ends,
												size_T);



		unsigned int b = t_bsearch_seq(A[i].end,
									   B_starts,
									   size_B,
									   T_starts,
									   size_T);

		//unsigned int x = b;

		while ( ( B_starts[b] == A[i].end) && b < size_B)
			++b;

		unsigned int num_cant_after = size_B - b;

		unsigned int num_left = size_B - num_cant_before - num_cant_after;

		O += num_left;

	}

	free(B_starts);
	free(B_ends);
	return O;
}
//}}}

//{{{ unsigned int count_intersections_bsearch_seq(struct interval *A,
unsigned int count_intersections_sort_bsearch_seq(struct interval *A,
										          unsigned int size_A,
											      struct interval *B,
										          unsigned int size_B)
{

	unsigned int i, O = 0;

	unsigned int *B_starts =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));
	unsigned int *B_ends =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));

	for (i = 0; i < size_B; i++) {
		B_starts[i] = B[i].start;
		B_ends[i] = B[i].end;
	}

	qsort(B_starts, size_B, sizeof(unsigned int), compare_unsigned_int); 
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 

	qsort(A, size_A, sizeof(struct interval), compare_interval_by_start); 

	for (i = 0; i < size_A; i++) {
		unsigned int num_cant_before = bsearch_seq(A[i].start,
												   B_ends,
												   size_B,
												   -1,
												   size_B);
		unsigned int b = bsearch_seq(A[i].end,
								     B_starts,
								     size_B,
								     -1,
								     size_B);

		while ( ( B_starts[b] == A[i].end) && b < size_B)
			++b;

		unsigned int num_cant_after = size_B - b;

		unsigned int num_left = size_B - num_cant_before - num_cant_after;
		O += num_left;
	}

	free(B_starts);
	free(B_ends);
	return O;
}
//}}}

//{{{unsigned int count_intersections_brute_force_seq(struct interval *A,
unsigned int count_intersections_brute_force_seq(struct interval *A,
												 unsigned int size_A,
											     struct interval *B,
										         unsigned int size_B)
{
	unsigned int i, j, O = 0;

	for (i = 0; i < size_A; i++) 
		for (j = 0; j < size_B; j++) 
			if ( ( A[i].start <= B[j].end ) &&
				 ( A[i].end >= B[j].start ) )
				++O;
	return O;
}
//}}}

//{{{void generate_interval_sets(struct interval *A,
// must run init_genrand(seed) first
void generate_interval_sets(struct interval *A,
							unsigned int size_A,
							unsigned int len_A,
							struct interval *B,
							unsigned int size_B,
							unsigned int len_B,
							unsigned int P)
{

	int i;

    for (i = 0; i < size_A; i++) {
		A[i].start = genrand_int32();
		A[i].end = A[i].start + len_A;
	}
	
	qsort(A, size_A, sizeof(struct interval), compare_interval_by_start); 
	
	/*
	 * Draw a number to see if the next interval will intersect or not.
	 * Draw a number to get the next interval, make new interval intersect or
	 * not with the drawn interval based on the first draw.
	 */
	int p_max = 100;
	unsigned int p_mask = get_mask(p_max);
	unsigned int i_mask = get_mask(size_A);
	unsigned int l_mask = get_mask(len_A);

    for (i = 0; i < size_B; i++) {
		unsigned int next_i = get_rand(size_A, i_mask);
		unsigned int next_p = get_rand(p_max, p_mask);

		if (P >= next_p) // intersect
			// Pick an rand between start and end to start from
			B[i].start = A[next_i].start + get_rand(len_A, l_mask);
		else  // do not intersect
			B[i].start = A[next_i].end + get_rand(len_A, l_mask);

		B[i].end = B[i].start + len_B;
	}
}
//}}}

//{{{void generate_interval_sets(struct interval *A,
// must run init_genrand(seed) first
void generate_ind_interval_sets(struct interval *A,
							unsigned int size_A,
							unsigned int len_A,
							struct interval *B,
							unsigned int size_B,
							unsigned int len_B)
{

	int i;

    for (i = 0; i < size_A; i++) {
		A[i].start = genrand_int32();
		A[i].end = A[i].start + len_A;
	}
	
    for (i = 0; i < size_B; i++) {
		B[i].start = genrand_int32();
		B[i].end = B[i].start + len_B;
	}
}
//}}}

//{{{ unsigned int test_intersections_bsearch_seq(struct interval *A,
void test_intersections_bsearch_seq(struct interval *A,
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
	unsigned int i,j;
	int R = 0;
	struct timeval t_start;
	unsigned int bits = (int)( ceil(log(max_offset)/log(2) ) );
	unsigned int mask = (2 << (bits-1)) - 1;

	//{{{ i_sort sort intervals by start
	qsort(B, size_B, sizeof(struct interval), compare_interval_by_start); 
	//}}}

	//{{{ i_mem split intervals into start, end, lenght, intervals
	unsigned int *A_starts =
			(unsigned int *) malloc(size_A * sizeof(unsigned int));
	unsigned int *B_starts =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));
	unsigned int *B_ends =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));
	unsigned int *B_lens =
			(unsigned int *) malloc(size_B * sizeof(unsigned int));

	for (i = 0; i < size_B; i++) {
		B_starts[i] = B[i].start;
		B_ends[i] = B[i].end;
		B_lens[i] = B[i].end - B[i].start;
	}
	//}}} i_mem

	//{{{ i_sort ends
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 
	//}}} i_sort

	//{{{ i_do find observed intersections
	*O  = count_seq(A,
					size_A,
					B_starts,
					B_ends,
					size_B);

	//}}} i_do
	
	init_genrand((unsigned) time(NULL));

	unsigned int *sims = (unsigned int *) malloc(n * sizeof(unsigned int));
	for (j = 0; j < n; j++) {

            /*
		two_side_simple_permute(mask,
							 max_offset,
							 A,
							 A_starts,
							 size_A,
							 B_starts,
							 B_ends,
							 B_lens,
							 size_B);
            */

		one_side_simple_permute(mask,
							 max_offset,
							 A,
							 A_starts,
							 size_A,
							 B_starts,
							 B_ends,
							 B_lens,
							 size_B);



		unsigned int T  = count_seq(A,
									size_A,
									B_starts,
									B_ends,
									size_B);

		sims[j] = T;
		if (T >= *O) 
			R = R + 1;
	}

	*mean = gsl_stats_int_mean((const int*)sims, 1, n);
	*sd = gsl_stats_int_sd_m((const int*)sims, 1, n, *mean);
	*p =  ( (double) R + 1) / ( (double) n + 1);


	//{{{ i_mem free B_starts B_ends B_lens
	free(B_starts);
	free(B_ends);
	free(B_lens);
	free(sims);
	//}}}
}
//}}}

//{{{ unsigned int count_seq(struct interval *A,
unsigned int count_seq(struct interval *A,
					   unsigned int size_A,
					   unsigned int *B_starts,
					   unsigned int *B_ends,
					   unsigned int size_B)
{
	unsigned int i, O = 0;
	unsigned int num_cant_before, b, num_cant_after, num_left;

	for (i = 0; i < size_A; i++) {

		num_cant_before = bsearch_seq(A[i].start,
												   B_ends,
												   size_B,
												   -1,
												   size_B);
		b = bsearch_seq(A[i].end,
								     B_starts,
								     size_B,
								     -1,
								     size_B);

		while ( ( B_starts[b] == A[i].end) && b < size_B)
			++b;

		num_cant_after = size_B - b;

		num_left = size_B - num_cant_before - num_cant_after;
		O += num_left;
	}

	return O;
}
//}}}

//{{{ unsigned int count_intersections_bsearch_seq_mem(struct interval *A,
unsigned int count_intersections_bsearch_seq_mem(struct interval *A,
												 unsigned int size_A,
												 unsigned int *B_starts,
												 unsigned int *B_ends,
												 unsigned int size_B)
{

	unsigned int i, O = 0;

	struct timeval t_start = in();
	qsort(B_starts, size_B, sizeof(unsigned int), compare_unsigned_int); 
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 
	unsigned int sort_time = out(t_start);


	O = count_seq(A,
				  size_A,
				  B_starts,
				  B_ends,
				  size_B);

	return O;
}
//}}}


//{{{ void two_side_simple_permute(unsigned int mask,
void two_side_simple_permute(unsigned int mask,
							 unsigned int max_offset,
							 struct interval *A,
							 unsigned int *A_starts,
							 unsigned int size_A,
							 unsigned int *B_starts,
							 unsigned int *B_ends,
							 unsigned int *B_lens,
							 unsigned int size_B)
{
	int i;
	for (i = 0; i < size_A; i++) 
		A_starts[i] = get_rand(max_offset, mask);

	qsort(A_starts, size_A, sizeof(unsigned int), compare_unsigned_int); 

	int len;
	for (i = 0; i < size_A; i++) {
		len = A[i].end - A[i].start;
		A[i].start = A_starts[i];
		A[i].end = A[i].start + len;
	}


	for (i = 0; i < size_B; i++) 
		B_starts[i] = get_rand(max_offset, mask);

	qsort(B_starts, size_B, sizeof(unsigned int), compare_unsigned_int); 

	for (i = 0; i < size_B; i++) 
		B_ends[i] = B_starts[i] + B_lens[i];
	
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 
}
//}}}

//{{{ void one_side_simple_permute(unsigned int mask,
void one_side_simple_permute(unsigned int mask,
							 unsigned int max_offset,
							 struct interval *A,
							 unsigned int *A_starts,
							 unsigned int size_A,
							 unsigned int *B_starts,
							 unsigned int *B_ends,
							 unsigned int *B_lens,
							 unsigned int size_B)
{
	int i;

	for (i = 0; i < size_B; i++) 
		B_starts[i] = get_rand(max_offset, mask);

	qsort(B_starts, size_B, sizeof(unsigned int), compare_unsigned_int); 

	for (i = 0; i < size_B; i++) 
		B_ends[i] = B_starts[i] + B_lens[i];
	
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 
}
//}}}

//{{{ void one_side_gap_permute(unsigned int mask,
void one_side_gap_permute(unsigned int mask,
							 unsigned int max_offset,
							 struct interval *A,
							 unsigned int *A_starts,
							 unsigned int size_A,
							 unsigned int *B_starts,
							 unsigned int *B_ends,
							 unsigned int *B_lens,
							 unsigned int size_B)
{
	int i;

	for (i = 0; i < size_B; i++) 
		B_starts[i] = get_rand(max_offset, mask);

	qsort(B_starts, size_B, sizeof(unsigned int), compare_unsigned_int); 

	for (i = 0; i < size_B; i++) 
		B_ends[i] = B_starts[i] + B_lens[i];
	
	qsort(B_ends, size_B, sizeof(unsigned int), compare_unsigned_int); 
}
//}}}

