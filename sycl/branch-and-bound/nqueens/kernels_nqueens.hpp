// Generated from opencl/branch-and-bound/nqueens/kernels_nqueens.cl: kernel bodies verbatim, preprocessor branches
// resolved as the OpenCL host builds them (CPU device: -D FORCE_CPU; NVIDIA GPU: scalar kernels with -D USE_ATOMICS,
// PREDICATED rotation check, private arrays). -local (FORCE_LOCAL) and the vectorized kernels are not translated.
#ifndef KERNELS_NQUEENS_SYCL_HPP
#define KERNELS_NQUEENS_SYCL_HPP
#include "dwarf_sycl.hpp"
namespace nqk {
typedef unsigned int uint;
inline int bit_scan(unsigned int x)
{
	int res = 0;
	res |= (x & 0xaaaaaaaa) ? 1 : 0;
	res |= (x & 0xcccccccc) ? 2 : 0;
	res |= (x & 0xf0f0f0f0) ? 4 : 0;
	res |= (x & 0xff00ff00) ? 8 : 0;
	res |= (x & 0xffff0000) ? 16 : 0;
	return res;
}
inline int nq_atom_inc(int *p) {
	return sycl::atomic_ref<int, sycl::memory_order::relaxed, sycl::memory_scope::device, sycl::access::address_space::global_space>(*p).fetch_add(1);
}

/* ---- CPU (FORCE_CPU) kernels ---- */
#define BOARD(x) ((x) < board_size - level ? params[idx + pitch * (4 + (x))] : ns[(x) - board_size + level])
inline void cpu_nqueen(sycl::queue &q, sycl::nd_range<1> ndr, int board_size, int level, int threads, int pitch, uint* params, uint* results, const uint* forbidden)
{
	q.parallel_for(ndr, [=](sycl::nd_item<1> it) {
		int idx = it.get_global_id(0);
		int tid = it.get_local_id(0); (void)tid;
		uint mask = params[idx];
		uint left_mask = params[idx + pitch];
		uint right_mask = params[idx + pitch * 2];
		uint index = params[idx + pitch * 3];
		uint board_mask = (1 << board_size) - 1;
		uint left_mask_big = 0;
		uint right_mask_big = 0;

		uint left_masks[32];
		uint right_masks[32];
		uint ms[32];
		uint ns[33];
		uint solutions = 0;
		uint unique_solutions = 0;
		int i = 0;
		int j;

		int t_array[32];
		int board_array[32];

		ms[0] = mask | left_mask | right_mask | forbidden[0];
		ns[0] = ((ms[0] + 1) & ~ms[0]);
		left_masks[0] = left_mask;
		right_masks[0] = right_mask;

		results[idx + pitch * 2] = params[idx];	// for checking kernel exeuction completion

		while(i >= 0) {
			if((ns[i] & board_mask) != 0) {
				mask |= ns[i];
				left_masks[i+1] = (left_masks[i] | ns[i]) << 1;
				right_masks[i+1] = (right_masks[i] | ns[i]) >> 1;
				ms[i+1] = mask | left_masks[i+1] | right_masks[i+1] | forbidden[i + 1];
				ns[i+1] = ((ms[i+1] + 1) & ~ms[i+1]);
				i++;
			}
			else {
				if(i == level) {
					int repeat_times = 8;
					bool repeat = false;
					bool equal = true;

					bool rotate1 = (BOARD(index) == (1 << (board_size - 1)));
					bool rotate2 = (BOARD(board_size - index - 1) == 1);
					bool rotate3 = (ns[level - 1] == (1 << (board_size - index - 1)));

					if(rotate1 || rotate2 || rotate3) {
						for(j = 0; j < board_size - level; j++) {
							board_array[j] = bit_scan(params[idx + pitch * (4 + j)]);
						}	
						for(j = 0; j < level; j++) {
							board_array[j + board_size - level] = bit_scan(ns[j]);
						}

						if(rotate3) {
							// rotate 180
							equal = true;
							for(j = 0; j < board_size; j++) {
								t_array[board_size - j - 1] = (board_size - board_array[j] - 1);
							}

							for(j = 0; j < board_size; j++) {
								if(t_array[j] < board_array[j]) {
									repeat = true;
									equal = false;
									break;
								}
								else if(t_array[j] > board_array[j]) {
									equal = false;
									break;
								}
							}

							if(equal) {
								repeat_times = 4;
							}
						}

						// rotate cw
						if(!repeat && rotate1) {
							equal = true;
							for(j = 0; j < board_size; j++) {
								t_array[board_size - board_array[j] - 1] = j;
							}

							for(j = 0; j < board_size; j++) {
								if(t_array[j] < board_array[j]) {
									repeat = true;
									equal = false;
									break;
								}
								else if(t_array[j] > board_array[j]) {
									equal = false;
									break;
								}
							}

							if(equal) {
								repeat_times = 2;
							}
						}

						if(!repeat && rotate2) {
							// rotate ccw
							equal = true;
							for(j = 0; j < board_size; j++) {
								t_array[board_array[j]] = (board_size - j - 1);
							}
							for(j = 0; j < board_size; j++) {
								if(t_array[j] < board_array[j]) {
									repeat = true;
									equal = false;
									break;
								}
								else if(t_array[j] > board_array[j]) {
									equal = false;
									break;
								}
							}

							if(equal) {
								repeat_times = 2;
							}
						}

						if(!repeat) {
							solutions += repeat_times;
							unique_solutions++;
						}
					}
					else {
						solutions += 8;
						unique_solutions ++;
					}
				}

				i--;
				if(i >= 0) {
					mask &= ~ns[i];
					ms[i] |= ns[i];
					ns[i] = ((ms[i] + 1) & ~ms[i]);
				}
			}
		}

		results[idx] = solutions;
		results[idx + pitch] = unique_solutions;
		results[idx + pitch * 3] = params[idx];	// for checking kernel exeuction completion
}).wait_and_throw();
}

inline void cpu_nqueen1(sycl::queue &q, sycl::nd_range<1> ndr, int board_size, int level, int threads, int pitch, uint* params, uint* results, const uint* forbidden)
{
	q.parallel_for(ndr, [=](sycl::nd_item<1> it) {
		int idx = it.get_global_id(0);
		int tid = it.get_local_id(0); (void)tid;
		uint mask = params[idx];
		uint left_mask = params[idx + pitch];
		uint right_mask = params[idx + pitch * 2];
		int second_row = params[idx + pitch * 3];
		uint board_mask = (1 << board_size) - 1;
		uint left_mask_big = 0;
		uint right_mask_big = 0;

		uint left_masks[32];
		uint right_masks[32];
		uint ms[32];
		uint ns[33];
		uint solutions = 0;
		int i = 0;

		ms[0] = mask | left_mask | right_mask | (i < second_row ? 2 : 0);
		ns[0] = ((ms[0] + 1) & ~ms[0]);
		left_masks[0] = left_mask;
		right_masks[0] = right_mask;

		results[idx + pitch * 2] = params[idx];	// for checking kernel exeuction completion

		while(i >= 0) {
			if((ns[i] & board_mask) != 0) {
				mask |= ns[i];
				left_masks[i+1] = (left_masks[i] | ns[i]) << 1;
				right_masks[i+1] = (right_masks[i] | ns[i]) >> 1;
				ms[i+1] = mask | left_masks[i+1] | right_masks[i+1] | (i + 1 < second_row ? 2 : 0);
				ns[i + 1] = ((ms[i+1] + 1) & ~ms[i+1]);
				i++;
			}
			else {
				if(i == level) {
					solutions++;
				}

				i--;
				if(i >= 0) {
					mask &= ~ns[i];
					ms[i] |= ns[i];
					ns[i] = ((ms[i] + 1) & ~ms[i]);
				}
			}
		}

		results[idx] = solutions * 8;
		results[idx + pitch] = solutions;
		results[idx + pitch * 3] = params[idx];	// for checking kernel exeuction completion
}).wait_and_throw();
}

#undef BOARD

/* ---- GPU (persistent-thread) kernels ---- */
#define NS(x) ns[x]
#define FORBIDDEN(x) forbidden_data[x]
#define BOARD_ARRAY(x) board_array[x]
#define BOARD(x) ((x) < board_size - level ? params[idx + pitch * (4 + (x))] : NS((x) - board_size + level))
template <bool ATOMICS> inline void gpu_nqueen(sycl::queue &q, sycl::nd_range<1> ndr, int board_size, int level, int threads, int pitch, uint* params, uint* results, const uint* forbidden_data, int* global_index)
{
	q.parallel_for(ndr, [=](sycl::nd_item<1> it) {
		int idx = it.get_global_id(0);
		int tid = it.get_local_id(0); (void)tid;

		int board_array[32];

			uint ms;
		uint ns[12];
		

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
									int p = board_size - BOARD_ARRAY(j) - 1;
									bool cond = (BOARD_ARRAY(p) != j);
									equal = equal && !cond;
									cond = cond && (min_pos > p);
									relation = cond ? BOARD_ARRAY(p) - j : relation;
									min_pos = cond ? p : min_pos;
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
									int p = BOARD_ARRAY(j);
									bool cond = (BOARD_ARRAY(p) != board_size - j - 1);
									equal = equal && !cond;
									cond = cond && (min_pos > p);
									relation = cond ? BOARD_ARRAY(p) - (board_size - j - 1) : relation;
									min_pos = cond ? p : min_pos;

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


	if constexpr (ATOMICS) {
				idx = nq_atom_inc(global_index);
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
	}
	else break;
			}
		}
}).wait_and_throw();
}

template <bool ATOMICS> inline void gpu_nqueen1(sycl::queue &q, sycl::nd_range<1> ndr, int board_size, int level, int threads, int pitch, uint* params, uint* results, const uint* forbidden, int* global_index)
{
	q.parallel_for(ndr, [=](sycl::nd_item<1> it) {
		int idx = it.get_global_id(0);
		int tid = it.get_local_id(0); (void)tid;

		uint ms;
		uint ns[12];

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

	if constexpr (ATOMICS) {
				idx = nq_atom_inc(global_index);
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
	}
	else break;
			}
		}
}).wait_and_throw();
}

#undef NS
#undef FORBIDDEN
#undef BOARD_ARRAY
#undef BOARD
} // namespace nqk
#endif
