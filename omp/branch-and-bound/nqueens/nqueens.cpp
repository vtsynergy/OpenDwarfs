/* N-Queens (branch-and-bound dwarf) -- OpenMP target-offload translation of
 * opencl/branch-and-bound/nqueens (Ping-Che Chen's OpenCL N-queen solver).
 *
 * Algorithm (unchanged): the host enumerates partial boards of the first (board_size - level) rows
 * (using the symmetry-reduced first-row loop j = 0 .. board_size/2-1) into batches; each batch is
 * uploaded and a kernel (`nqueen`, or `nqueen1` for j == 0) completes the remaining `level` rows with an
 * iterative bit-mask search, counting solutions and symmetry-unique solutions; the host sums the
 * per-task results.
 *
 * Translation of the *default NVIDIA configuration* of the OpenCL app (scalar kernels, global
 * atomics = work-stealing persistent threads, WORK_ITEMS=256, nThreads = 256 x 80 SMs x 3 = 61440),
 * i.e. the code selected by -D WORK_ITEMS=256 -D USE_ATOMICS -D ENABLE_CHAR. The other OpenCL build
 * variants (vectorised, FORCE_LOCAL, CPU without atomics) are not translated. The kernel bodies below
 * are the OpenCL kernel text; only get_global_id / atom_inc / address-space qualifiers changed.
 * The OpenCL host issues clFinish after every launch, so its device/host overlap is already
 * serialised; the translation reads the results back right after each launch, which gives
 * identical sums.
 *
 * Usage: nqueens [options] N      (1 <= N <= 32)   options: -cpu  -threads #  -blocksize #
 * Timed region: the whole Compute() call (host enumeration + uploads + kernels + result downloads).
 * Validation (added): compared with the original single-threaded CPU solver (nqueen_cpu.cpp, run
 * untimed) and, where known, the published solution counts.
 */
#include <iostream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <omp.h>
#include "nqueen_cpu.h"
#include "dwarf_timer.h"

typedef unsigned int uint;
#define PREDICATED
#define USE_ATOMICS
#define WORK_ITEMS 256
/* default macros of the scalar, non-FORCE_LOCAL kernel configuration */
#define ARRAY_DECL uint ns[12];
#define NS(x) ns[x]
#define FORBIDDEN_DECL
#define FORBIDDEN(x) forbidden_data[x]
#define BOARD_ARRAY_DECL int board_array[32];
#define BOARD_ARRAY(x) board_array[x]
#define BOARD(x) ((x) < board_size - level ? params[idx + pitch * (4 + (x))] : NS((x) - board_size + level))

#pragma omp declare target
static inline int bit_scan(unsigned int x) {
	int res = 0;
	res |= (x & 0xaaaaaaaa) ? 1 : 0;
	res |= (x & 0xcccccccc) ? 2 : 0;
	res |= (x & 0xf0f0f0f0) ? 4 : 0;
	res |= (x & 0xff00ff00) ? 8 : 0;
	res |= (x & 0xffff0000) ? 16 : 0;
	return res;
}
static inline int atom_inc_d(int *p) {
	int v;
#pragma omp atomic capture
	{ v = *p; (*p)++; }
	return v;
}

/* ---- kernel `nqueen` (one work-item = persistent worker; atomics = work stealing) ---- */
static void nqueen_item(int idx, int board_size, int level, int threads, int pitch, const uint *params, uint *results,
                        const uint *forbidden_data, int *global_index) {
	int tid = 0; (void)tid;
	
	

	BOARD_ARRAY_DECL

		uint ms;
	ARRAY_DECL
		FORBIDDEN_DECL

		uint solutions = 0;
	uint unique_solutions = 0;
	int i;
	int j;

	uint mask = params[idx];
	uint left_mask = params[idx + pitch];
	uint right_mask = params[idx + pitch * 2];
	uint index = params[idx + pitch * 3];
	uint board_mask = (1 << board_size) - 1;
	uint left_mask_big = 0;
	uint right_mask_big = 0;

	results[idx + pitch * 2] = params[idx];	// for checking kernel exeuction completion

	for(j = 0; j < board_size - level; j++) {
		BOARD_ARRAY(j) = bit_scan(params[idx + pitch * (4 + j)]);
	}	

#ifdef FORCE_LOCAL
	if(tid < 32) {
		for(i = 0; i <= level; i++) {
			FORBIDDEN(i) = forbidden_data[i];
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);
#endif	

	i = 0;

	ms = mask | left_mask | right_mask | FORBIDDEN(0);
	NS(0) = ((ms + 1) & ~ms);

	for(;;) {
		if(i >= 0) {
			if((NS(i) & board_mask) != 0) {
				mask |= NS(i);
				left_mask_big = (left_mask_big << 1) | (left_mask >> 31);
				left_mask = (left_mask | NS(i)) << 1;
				right_mask_big = (right_mask_big >> 1) | (right_mask << 31);
				right_mask = (right_mask | NS(i)) >> 1;
				ms = mask | left_mask | right_mask | FORBIDDEN(i + 1);
				NS(i + 1) = ((ms + 1) & ~ms);
				i++;
			}
			else {
				if(i == level) {
					int repeat_times = 8;
					bool equal = true;

					bool rotate1 = (BOARD(index) == (1 << (board_size - 1)));
					bool rotate2 = (BOARD(board_size - index - 1) == 1);
					bool rotate3 = (NS(level - 1) == (1 << (board_size - index - 1)));

					if(rotate1 || rotate2 || rotate3) {
						for(j = 0; j < level; j++) {
							BOARD_ARRAY(j + board_size - level) = bit_scan(NS(j));
						}

						int min_pos = board_size;
						int relation = 0;

						// rotate cw
						if(rotate1) {
							equal = true;
							//						for(j = 0; j < board_size; j++) {
							//							t_array[board_size - board_array[j] - 1] = j;
							//						}
							//							
							//						for(j = 0; j < board_size; j++) {
							//							if(t_array[j] < board_array[j]) {
							//								repeat = true;
							//								equal = false;
							//								break;
							//							}
							//							else if(t_array[j] > board_array[j]) {
							//								equal = false;
							//								break;
							//							}
							//						}

							relation = 0;
							for(j = 0; j < board_size; j++) {
#ifndef PREDICATED							
								if(BOARD_ARRAY(board_size - BOARD_ARRAY(j) - 1) != j) {
									equal = false;
									if(min_pos > board_size - BOARD_ARRAY(j) - 1) {
										relation = BOARD_ARRAY(board_size - BOARD_ARRAY(j) - 1) - j;
										min_pos = board_size - BOARD_ARRAY(j) - 1;
									}							
								}
#else
								int p = board_size - BOARD_ARRAY(j) - 1;
								bool cond = (BOARD_ARRAY(p) != j);
								equal = equal && !cond;
								cond = cond && (min_pos > p);
								relation = cond ? BOARD_ARRAY(p) - j : relation;
								min_pos = cond ? p : min_pos;
#endif								
							}

							repeat_times = equal ? 2 : repeat_times;
							//							if(equal) {
							//								repeat_times = 2;
							//							}
						}

						if(relation >= 0 && rotate2) {
							// rotate ccw
							equal = true;
							//						for(j = 0; j < board_size; j++) {
							//							t_array[board_array[j]] = (board_size - j - 1);
							//						}
							//						for(j = 0; j < board_size; j++) {
							//							if(t_array[j] < board_array[j]) {
							//								repeat = true;
							//								equal = false;
							//								break;
							//							}
							//							else if(t_array[j] > board_array[j]) {
							//								equal = false;
							//								break;
							//							}
							//						}

							min_pos = board_size;
							relation = 0;
							for(j = 0; j < board_size; j++) {
#ifndef PREDICATED							
								if(BOARD_ARRAY(BOARD_ARRAY(j)) != board_size - j - 1) {
									equal = false;
									if(min_pos > BOARD_ARRAY(j)) {
										relation = BOARD_ARRAY(BOARD_ARRAY(j)) - (board_size - j - 1);
										min_pos = BOARD_ARRAY(j);
									}
								}
#else
								int p = BOARD_ARRAY(j);
								bool cond = (BOARD_ARRAY(p) != board_size - j - 1);
								equal = equal && !cond;
								cond = cond && (min_pos > p);
								relation = cond ? BOARD_ARRAY(p) - (board_size - j - 1) : relation;
								min_pos = cond ? p : min_pos;
#endif

							}

							repeat_times = equal ? 2 : repeat_times;
							//							if(equal) {
							//								repeat_times = 2;
							//							}
						}

						if(relation >= 0 && repeat_times == 8 && rotate3) {
							// rotate 180
							equal = true;
							//						for(j = 0; j < board_size; j++) {
							//							t_array[board_size - j - 1] = (board_size - board_array[j] - 1);
							//						}
							//							
							//						for(j = 0; j < board_size; j++) {
							//							if(t_array[j] < board_array[j]) {
							//								repeat = true;
							//								equal = false;
							//								break;
							//							}
							//							else if(t_array[j] > board_array[j]) {
							//								equal = false;
							//								break;
							//							}
							//						}

							min_pos = board_size;
							relation = 0;
							for(j = board_size - 1; j >= board_size / 2; j--) {
								if(BOARD_ARRAY(board_size - j - 1) != board_size - BOARD_ARRAY(j) - 1) {
									equal = false;
									relation = BOARD_ARRAY(board_size - j - 1) - (board_size - BOARD_ARRAY(j) - 1);
									break;
								}
							}
							//						for(j = 0; j < board_size; j++) {
							//							if(BOARD_ARRAY(board_size - j - 1) != board_size - BOARD_ARRAY(j) - 1) {
							//								equal = false;
							//								if(min_pos > board_size - j - 1) {
							//									relation = BOARD_ARRAY(board_size - j - 1) - (board_size - BOARD_ARRAY(j) - 1);
							//									min_pos = board_size - j - 1;
							//								}
							//							}
							//						}

							repeat_times = equal ? 4 : repeat_times;
							//							if(equal) {
							//								repeat_times = 4;
							//							}
						}

						solutions += (relation >= 0) ? repeat_times : 0;
						unique_solutions += (relation >= 0) ? 1 : 0;
						//						if(relation >= 0) {
						//							solutions += repeat_times;
						//							unique_solutions++;
						//						}
					}
					else {
						solutions += 8;
						unique_solutions ++;
					}
				}

				i--;
				if(i >= 0) {
					mask &= ~NS(i);
					left_mask = ((left_mask >> 1) | (left_mask_big << 31)) & ~NS(i);
					left_mask_big >>= 1;
					right_mask = ((right_mask << 1) | (right_mask_big >> 31)) & ~NS(i);
					right_mask_big <<= 1;
					ms = mask | left_mask | right_mask | NS(i) | FORBIDDEN(i);
					NS(i) = ((ms + NS(i)) & ~ms);
				}
			}
		}
		else {
			results[idx] = solutions;
			results[idx + pitch] = unique_solutions;
			results[idx + pitch * 3] = params[idx];	// for checking kernel exeuction completion


#ifdef USE_ATOMICS
			idx = atom_inc_d(global_index);
			if(idx < threads) {
				solutions = 0;
				unique_solutions = 0;
				i = 0;

				mask = params[idx];
				left_mask = params[idx + pitch];
				right_mask = params[idx + pitch * 2];
				index = params[idx + pitch * 3];
				left_mask_big = 0;
				right_mask_big = 0;

				results[idx + pitch * 2] = params[idx];	// for checking kernel exeuction completion			

				for(j = 0; j < board_size - level; j++) {
					BOARD_ARRAY(j) = bit_scan(params[idx + pitch * (4 + j)]);
				}

				ms = mask | left_mask | right_mask | FORBIDDEN(0);
				NS(0) = ((ms + 1) & ~ms);
			}
			else {
				break;
			}
#else
			break;
#endif
		}
	}
}

/* ---- kernel `nqueen1` (first row = 0 case) ---- */
static void nqueen1_item(int idx, int board_size, int level, int threads, int pitch, const uint *params, uint *results,
                         const uint *forbidden, int *global_index) {
	int tid = 0; (void)tid;
	
	

	uint ms;
	ARRAY_DECL

		uint mask = params[idx];
	uint left_mask = params[idx + pitch];
	uint right_mask = params[idx + pitch * 2];
	int second_row = params[idx + pitch * 3];
	uint board_mask = (1 << board_size) - 1;
	uint left_mask_big = 0;
	uint right_mask_big = 0;
	uint solutions = 0;
	int i = 0;

	ms = mask | left_mask | right_mask | (i < second_row ? 2 : 0);
	NS(0) = ((ms + 1) & ~ms);

	results[idx + pitch * 2] = params[idx];	// for checking kernel exeuction completion

	for(;;) {
		if(i >= 0) {
			if((NS(i) & board_mask) != 0) {
				mask |= NS(i);
				left_mask_big = (left_mask_big << 1) | (left_mask >> 31);
				left_mask = (left_mask | NS(i)) << 1;
				right_mask_big = (right_mask_big >> 1) | (right_mask << 31);
				right_mask = (right_mask | NS(i)) >> 1;
				ms = mask | left_mask | right_mask | (i + 1 < second_row ? 2 : 0);
				NS(i + 1) = ((ms + 1) & ~ms);
				i++;
			}
			else {
				if(i == level) {
					solutions++;
				}

				i--;
				if(i >= 0) {
					mask &= ~NS(i);
					left_mask = ((left_mask >> 1) | (left_mask_big << 31)) & ~NS(i);
					left_mask_big >>= 1;
					right_mask = ((right_mask << 1) | (right_mask_big >> 31)) & ~NS(i);
					right_mask_big <<= 1;
					ms = mask | left_mask | right_mask | NS(i) | (i < second_row ? 2 : 0);
					NS(i) = ((ms + NS(i)) & ~ms);
				}
			}

			/*
			   bool n = ((NS(i) & board_mask) != 0);
			   uint nss = NS(i);
			   mask = n ? (mask | nss) : mask;
			   left_mask_big = n ? ((left_mask_big << 1) | (left_mask >> 31)) : left_mask_big;
			   left_mask = n ? ((left_mask | nss) << 1) : left_mask;
			   right_mask_big = n ? ((right_mask_big >> 1) | (right_mask << 31)) : right_mask_big;
			   right_mask = n ? ((right_mask | nss) >> 1) : right_mask;
			   ms = mask | left_mask | right_mask | (i + 1 < second_row ? 2 : 0);
			   nss = NS(i + 1);
			   NS(i + 1) = n ? ((ms + 1) & ~ms) : nss;
			//			i = n ? i + 1 : i;

			solutions = (!n && i == level) ? solutions + 1 : solutions;

			i = n ? i + 1 : i - 1;
			//			i = !n ? i - 1 : i;
			n = (!n && i >= 0);
			nss = ~NS(max(i, 0));
			mask = n ? (mask & nss) : mask;
			left_mask = n ? (((left_mask >> 1) | (left_mask_big << 31)) & nss) : left_mask;
			left_mask_big = n ? (left_mask_big >> 1) : left_mask_big;
			right_mask = n ? (((right_mask << 1) | (right_mask_big >> 31)) & nss) : right_mask;
			right_mask_big = n ? (right_mask_big << 1) : right_mask_big;
			ms = mask | left_mask | right_mask | ~nss | (i < second_row ? 2 : 0);
			if(n) NS(i) = (ms + ~nss) & ~ms;
			 */
		}
		else {
			results[idx] = solutions * 8;
			results[idx + pitch] = solutions;
			results[idx + pitch * 3] = params[idx];	// for checking kernel exeuction completion

#ifdef USE_ATOMICS
			idx = atom_inc_d(global_index);
			if(idx < threads) {
				mask = params[idx];
				left_mask = params[idx + pitch];
				right_mask = params[idx + pitch * 2];
				second_row = params[idx + pitch * 3];
				left_mask_big = 0;
				right_mask_big = 0;
				solutions = 0;
				i = 0;

				results[idx + pitch * 2] = params[idx];	// for checking kernel exeuction completion

				ms = mask | left_mask | right_mask | (i < second_row ? 2 : 0);
				NS(0) = ((ms + 1) & ~ms);
			}
			else {
				break;
			}
#else
			break;
#endif
		}
	}
}
#pragma omp end declare target

static inline int bit_scan_h(unsigned int x) { return bit_scan(x); }

/* One launch of nqueen / nqueen1 over `work_dim` work-items (work-groups of 256 -> flat loop: the
 * kernels use no barriers or local memory in this configuration). */
__attribute__((noinline)) static void launch(bool first_row_zero, int work_dim, int board_size, int level, int threads, int pitch, uint *params,
                   uint *results, uint *forbidden, int *gindex, size_t nparams, size_t nres) {
	if (first_row_zero) {
#pragma omp target teams distribute parallel for thread_limit(256) map(present, alloc : params[0 : nparams], results[0 : nres], forbidden[0 : 32], gindex[0 : 1])
		for (int idx = 0; idx < work_dim; idx++) nqueen1_item(idx, board_size, level, threads, pitch, params, results, forbidden, gindex);
	} else {
#pragma omp target teams distribute parallel for thread_limit(256) map(present, alloc : params[0 : nparams], results[0 : nres], forbidden[0 : 32], gindex[0 : 1])
		for (int idx = 0; idx < work_dim; idx++) nqueen_item(idx, board_size, level, threads, pitch, params, results, forbidden, gindex);
	}
}

struct SolverState {
	int nThreads = 61440;       // m_nThreads (OpenCL: 256 work-items x 80 CUs x (4-1) on V100)
	int nMaxWorkItems = 256;    // m_nMaxWorkItems
};

struct BatchCtx {
	uint *mask_vector, *results, *forb_buf;
	int *gindex_buf;
	size_t nparams, nres;
	int board_size, level, max_pitch, nMaxWorkItems;
	long long solutions, unique_solutions;
};

/* One kernel launch + immediate download/accumulate of the results. (A plain function, not a lambda:
 * nvc++ mis-computes OpenMP array sections of variables captured by a lambda.) */
static void run_batch(BatchCtx &c, int j, int t_size, bool tail, const unsigned int *forbidden, int work_dim) {
	uint *mv = c.mask_vector, *rs = c.results, *fb = c.forb_buf;
	int *gi = c.gindex_buf;
	size_t np = c.nparams, nr = c.nres;
	int board_size = c.board_size, level = c.level, max_pitch = c.max_pitch;
	int n = c.nMaxWorkItems;
	if (tail && work_dim % n != 0) work_dim += n - work_dim % n;
	// forbidden data (written once per j in the OpenCL host code; rewriting is equivalent)
	for (int k = 0; k <= level; k++) fb[k] = forbidden[board_size - level - 1 + k];
	int nfb = level + 1;
#pragma omp target update to(fb[0 : nfb])
#pragma omp target update to(mv[0 : np])
	gi[0] = work_dim;
#pragma omp target update to(gi[0 : 1])
	launch(j == 0, work_dim, board_size, level, t_size, max_pitch, mv, rs, fb, gi, np, nr);
#pragma omp target update from(rs[0 : nr])
	for (int k = 0; k < t_size; k++) {
		if (rs[k + max_pitch * 2] != rs[k + max_pitch * 3]) { std::cerr << "kernel execution check failed\n"; exit(2); }
		c.solutions += rs[k];
		c.unique_solutions += rs[k + max_pitch];
	}
}

static long long Compute(SolverState &S, int board_size, long long *unique) {
	// estimate amount of levels need to be computed on the device
	long long total = 10000000000LL;
	total /= 10;  // for atomics
	int level = 0;
	int i = board_size;
	while (total > 0 && i > 0) { total /= ((i + 1) / 2); i--; level++; }
	if (level > board_size - 2) level = board_size - 2;
	if (level > 11) level = 11;

	int threads = S.nThreads * 16;  // atomics: many more tasks than workers
	int max_threads = threads;
	int max_pitch = (max_threads + 15) & ~0xf;

	// device buffers (ParamBuffer, ResultBuffer, ForbiddenBuffer, GlobalIndex); host arrays carry the mapping
	size_t nparams = (size_t)max_pitch * (4 + 32), nres = (size_t)max_pitch * 4;
	std::vector<unsigned int> mask_vector_v(nparams), results_v(nres);
	uint *mask_vector = mask_vector_v.data(), *results = results_v.data();
	uint forb_buf[32] = {0};
	int gindex_buf[1] = {0};
#pragma omp target enter data map(alloc : mask_vector[0 : nparams], results[0 : nres], forb_buf[0 : 32], gindex_buf[0 : 1])

	long long solutions = 0, unique_solutions = 0;
	unsigned int board_mask = (1 << board_size) - 1;
	int total_size = 0;
	int last_total_size = 0; (void)last_total_size;

	BatchCtx ctx = {mask_vector, results, forb_buf, gindex_buf, nparams, nres, board_size, level, max_pitch, S.nMaxWorkItems, 0, 0};
	for (int j = 0; j < board_size / 2; j++) {
		unsigned int masks[32], left_masks[32], right_masks[32], ms[32], ns[32], forbidden[32];
		unsigned int border_mask = 0;
		int idx = 0;
		int ii = 0;
		masks[0] = (1 << j);
		left_masks[0] = 1 << (j + 1);
		right_masks[0] = (1 << j) >> 1;
		ms[0] = masks[0] | left_masks[0] | right_masks[0];
		ns[0] = (1 << j);
		for (int k = 0; k < j; k++) {
			border_mask |= (1 << k);
			border_mask |= (1 << (board_size - k - 1));
		}
		for (int k = 0; k < board_size; k++) {
			if (k == board_size - 2) forbidden[k] = border_mask;
			else if ((k + 1) < j || (k + 1) > board_size - j - 1) forbidden[k] = 1 | (1 << (board_size - 1));
			else forbidden[k] = 0;
		}
		forbidden[board_size - 1] = 0xffffffff;

		while (ii >= 0) {
			if (j == 0) {
				if (ii >= 1) {
					unsigned int m = ms[ii] | (ii + 1 < idx ? 2 : 0);
					ns[ii + 1] = (m + 1) & ~m;
				} else {
					ns[ii + 1] = ((ms[ii] + 1) & ~ms[ii]);
					if (ii == 0) idx = bit_scan_h(ns[ii + 1]);
				}
			} else {
				unsigned int m = ms[ii] | forbidden[ii];
				ns[ii + 1] = (m + 1) & ~m;
			}

			if (ii == board_size - level - 1) {
				mask_vector[total_size] = masks[ii];
				mask_vector[total_size + max_pitch] = left_masks[ii];
				mask_vector[total_size + max_pitch * 2] = right_masks[ii];
				if (j == 0) mask_vector[total_size + max_pitch * 3] = idx - ii < 0 ? 0 : idx - ii;
				else mask_vector[total_size + max_pitch * 3] = j;  // check rotation
				for (int k = 0; k <= ii; k++) mask_vector[total_size + max_pitch * (k + 4)] = ns[k];
				total_size++;

				if (total_size == max_threads) {
					run_batch(ctx, j, threads, false, forbidden, S.nThreads);
					if (total_size > threads)
						for (int k = 0; k < 4 + board_size - level; k++)
							memcpy(&mask_vector[max_pitch * k], &mask_vector[threads + max_pitch * k], (total_size - threads) * sizeof(int));
					total_size -= threads;
				}
				ii--;
			} else if ((ns[ii + 1] & board_mask) != 0) {
				ms[ii] |= ns[ii + 1];
				masks[ii + 1] = masks[ii] | ns[ii + 1];
				left_masks[ii + 1] = (left_masks[ii] | ns[ii + 1]) << 1;
				right_masks[ii + 1] = (right_masks[ii] | ns[ii + 1]) >> 1;
				ms[ii + 1] = masks[ii + 1] | left_masks[ii + 1] | right_masks[ii + 1];
				ii++;
			} else {
				ii--;
			}
		}

		while (total_size > 0) {
			for (int k = total_size; k < max_threads; k++) {
				mask_vector[k] = 0xffffffff;
				mask_vector[k + max_pitch] = 0xffffffff;
				mask_vector[k + max_pitch * 2] = 0xffffffff;
				mask_vector[k + max_pitch * 3] = 0;
				mask_vector[k + max_pitch * 4] = 0;
			}
			int t_size = total_size > threads ? threads : total_size;
			int work_dim = (t_size < S.nThreads) ? t_size : S.nThreads;
			run_batch(ctx, j, t_size, true, forbidden, work_dim);
			if (total_size > t_size)
				for (int k = 0; k < 4 + board_size - level; k++)
					memcpy(&mask_vector[max_pitch * k], &mask_vector[t_size + max_pitch * k], (total_size - t_size) * sizeof(int));
			total_size -= t_size;
		}
	}

#pragma omp target exit data map(delete : mask_vector[0 : nparams], results[0 : nres], forb_buf[0 : 32], gindex_buf[0 : 1])
	solutions = ctx.solutions; unique_solutions = ctx.unique_solutions;
	if (unique != 0) *unique = unique_solutions;
	return solutions;
}

int main(int argc, char **argv) {
	std::cerr << "N-Queen solver for OpenCL (OpenMP target-offload translation)\n";
	std::cerr << "Ping-Che Chen\n\n";
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " [options] N\n";
		std::cerr << "\tN: board size (1 ~ 32)\n";
		std::cerr << "\t-cpu: use the original single-threaded CPU solver\n";
		std::cerr << "\t-threads #: set number of persistent threads to #\n";
		std::cerr << "\t-blocksize #: set size of thread blocks to #\n";
		return 0;
	}
	bool force_cpu = false;
	int threads = 0, block_size = 0, start = 1;
	while (start < argc - 1) {
		if (std::strcmp(argv[start], "-cpu") == 0) force_cpu = true;
		else if (std::strcmp(argv[start], "-threads") == 0 && start < argc - 2) { threads = std::atoi(argv[start + 1]); start++; }
		else if (std::strcmp(argv[start], "-blocksize") == 0 && start < argc - 2) { block_size = std::atoi(argv[start + 1]); start++; }
		else std::cerr << "Unknown option " << argv[start] << "\n";
		start++;
	}
	int board_size = std::atoi(argv[start]);
	if (board_size < 1 || board_size > 32) { std::cerr << "Inalid board size (only 1 ~ 32 allowed)\n"; return 0; }

	long long solutions = 0, unique_solutions = 0;
	SolverState S;
	if (threads) S.nThreads = threads;
	if (block_size) S.nMaxWorkItems = block_size;
	dwarf_init("nqueens", "branch-and-bound", "omp");
	dwarf_set_size(board_size, S.nThreads);
	if (force_cpu) {
		for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
			dwarf_region_begin();
			solutions = nqueen_cpu(board_size, &unique_solutions);
			dwarf_region_end(rep);
		}
	} else {
		std::cerr << "\tUsing " << S.nThreads << " threads\n\tBlock size = " << S.nMaxWorkItems << " threads\n\tUsing global atomics\n";
		for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
			dwarf_region_begin();
			solutions = Compute(S, board_size, &unique_solutions);
			dwarf_region_end(rep);
		}
	}
	std::cout << board_size << "-queen has " << solutions << " solutions (" << unique_solutions << " unique)\n";

	/* validation (added): original CPU solver, run untimed */
	long long ref_unique = 0, ref_sol = nqueen_cpu(board_size, &ref_unique);
	bool ok = (ref_sol == solutions && ref_unique == unique_solutions);
	std::cout << (ok ? "SUCCESS" : "FAILURE") << ": " << solutions << "/" << unique_solutions << " vs CPU reference " << ref_sol << "/" << ref_unique << "\n";
	return dwarf_report(ok);
}
