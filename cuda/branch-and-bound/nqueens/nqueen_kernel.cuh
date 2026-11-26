/* CUDA translation of kernels_nqueens.cl for the variant the OpenCL host builds on the NVIDIA GPU:
 * -D WORK_ITEMS=256 -D USE_ATOMICS -D ENABLE_CHAR, not vectorised, FORCE_LOCAL/FORCE_CPU undefined. The CPU-only and vector
 * kernels are not translated. atom_inc -> atomicAdd(p, 1); __constant buffer -> const device pointer. */
#define USE_ATOMICS
#define ENABLE_CHAR
#define PREDICATED
typedef unsigned int uint;

__device__ static inline int bit_scan_dev(unsigned int x)
{
	int res = 0;
	res |= (x & 0xaaaaaaaa) ? 1 : 0;
	res |= (x & 0xcccccccc) ? 2 : 0;
	res |= (x & 0xf0f0f0f0) ? 4 : 0;
	res |= (x & 0xff00ff00) ? 8 : 0;
	res |= (x & 0xffff0000) ? 16 : 0;
	return res;
}


#ifdef FORCE_LOCAL
#define ARRAY_DECL __shared__ uint ns[12][WORK_ITEMS];
#define NS(x) ns[x][tid]
#define FORBIDDEN_DECL __shared__ uint forbidden[12][32];
#define FORBIDDEN(x) forbidden[x][tid % 32]
//#define BOARD_ARRAY_DECL __shared__ char board_array[32][WORK_ITEMS];
//#define BOARD_ARRAY(x) board_array[x][tid]
#ifdef ENABLE_CHAR
#define BOARD_ARRAY_DECL char board_array[32];
#else
#define BOARD_ARRAY_DECL int board_array[32];
#endif
#define BOARD_ARRAY(x) board_array[x]
#else
#define ARRAY_DECL uint ns[12];
#define NS(x) ns[x]
#define FORBIDDEN_DECL
#define FORBIDDEN(x) forbidden_data[x]
#define BOARD_ARRAY_DECL int board_array[32];
#define BOARD_ARRAY(x) board_array[x]
#endif

#define BOARD(x) ((x) < board_size - level ? params[idx + pitch * (4 + (x))] : NS((x) - board_size + level))

__global__ void nqueen(int board_size, int level, int threads, int pitch, uint* params, uint* results, const uint* forbidden_data, int* global_index)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;

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
		BOARD_ARRAY(j) = bit_scan_dev(params[idx + pitch * (4 + j)]);
	}	

#ifdef FORCE_LOCAL
	if(tid < 32) {
		for(i = 0; i <= level; i++) {
			FORBIDDEN(i) = forbidden_data[i];
		}
	}

	__syncthreads();
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
							BOARD_ARRAY(j + board_size - level) = bit_scan_dev(NS(j));
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
			idx = atomicAdd(global_index, 1);
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
					BOARD_ARRAY(j) = bit_scan_dev(params[idx + pitch * (4 + j)]);
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

__global__ void nqueen1(int board_size, int level, int threads, int pitch, uint* params, uint* results, const uint* forbidden, int* global_index)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;

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
			idx = atomicAdd(global_index, 1);
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


