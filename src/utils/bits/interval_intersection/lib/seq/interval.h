/*****************************************************************************
interval.h
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
#ifndef __INTERVAL_H__
#define __INTERVAL_H__

struct interval {
	// order is hack to maintin the position of the interval in the orginal
	// data file
	unsigned int start, end, order;
};

unsigned int per_interval_count_intersections_bsearch_seq_sorted_int(
														  struct interval *A,
														  unsigned int size_A,
														  unsigned int *B_start,
														  unsigned int *B_end,
														  unsigned int size_B,
														  unsigned int *R);

int compare_interval_by_start (const void *a, const void *b);

unsigned int per_interval_count_intersections_bsearch_seq(struct interval *A,
														   unsigned int size_A,
														   struct interval *B,
														   unsigned int size_B,
														   unsigned int *R);
unsigned int enumerate_intersections_bsearch_seq(struct interval *A,
										 unsigned int size_A,
										 struct interval *B,
										 unsigned int size_B,
										 unsigned int **R,
										 unsigned int **E);

unsigned int count_intersections_bsearch_seq(struct interval *A,
										     unsigned int size_A,
											 struct interval *B,
										     unsigned int size_B);

unsigned int count_intersections_i_bsearch_seq(struct interval *A,
										       unsigned int size_A,
											   struct interval *B,
										       unsigned int size_B,
											   unsigned int size_I);

unsigned int count_intersections_t_bsearch_seq(struct interval *A,
										       unsigned int size_A,
											   struct interval *B,
										       unsigned int size_B,
											   unsigned int size_T);

unsigned int count_intersections_sort_bsearch_seq(struct interval *A,
										          unsigned int size_A,
											      struct interval *B,
										          unsigned int size_B);

unsigned int count_intersections_brute_force_seq(struct interval *A,
												 unsigned int size_A,
											     struct interval *B,
										         unsigned int size_B);
void generate_interval_sets(struct interval *A,
							unsigned int size_A,
							unsigned int len_A,
							struct interval *B,
							unsigned int size_B,
							unsigned int len_B,
							unsigned int P);

void generate_ind_interval_sets(struct interval *A,
							unsigned int size_A,
							unsigned int len_A,
							struct interval *B,
							unsigned int size_B,
							unsigned int len_B);

unsigned int count_seq(struct interval *A,
					   unsigned int size_A,
					   unsigned int *B_starts,
					   unsigned int *B_ends,
					   unsigned int size_B);

void test_intersections_bsearch_seq(struct interval *A,
									unsigned int size_A,
									struct interval *B,
									unsigned int size_B,
									unsigned int n,
									unsigned int max_offset,
								    unsigned int *O,
									double *mean,
									double *sd,
									double *p);

unsigned int count_intersections_bsearch_seq_mem(struct interval *A,
												 unsigned int size_A,
												 unsigned int *B_starts,
												 unsigned int *B_ends,
												 unsigned int size_B);

void two_side_simple_permute(unsigned int mask,
							 unsigned int max_offset,
							 struct interval *A,
							 unsigned int *A_starts,
							 unsigned int size_A,
							 unsigned int *B_starts,
							 unsigned int *B_ends,
							 unsigned int *B_lens,
							 unsigned int size_B);

void one_side_simple_permute(unsigned int mask,
							 unsigned int max_offset,
							 struct interval *A,
							 unsigned int *A_starts,
							 unsigned int size_A,
							 unsigned int *B_starts,
							 unsigned int *B_ends,
							 unsigned int *B_lens,
							 unsigned int size_B);

#endif
