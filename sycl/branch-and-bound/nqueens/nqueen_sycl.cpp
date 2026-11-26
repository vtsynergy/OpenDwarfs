// nqueen_sycl.cpp -- N-queen solver, translated from nqueen_cl.cpp (Ping-Che Chen; adapted for OpenDwarfs)
#include <cstring>
#include <vector>
#include "nqueen_sycl.h"
#include "kernels_nqueens.hpp"

namespace {
template <class T> struct DevBuf {
	sycl::queue &q; T *p;
	DevBuf(sycl::queue &q_, size_t n) : q(q_), p(dsycl::alloc<T>(q_, n)) {}
	~DevBuf() { dsycl::release(q, p); }
};
}

NQueenSolver::NQueenSolver(sycl::queue &q_, int threads, int block_size, bool force_no_atomics) : q(q_)
{
	m_bCPU = q.get_device().is_cpu();
	/* the OpenCL host enables global atomics on non-CPU devices that advertise them (NVIDIA does) unless -noatomics */
	m_bAtomics = !m_bCPU && !force_no_atomics;
	/* [FIX-OCL-34] algorithmic work-group 256 and persistent-thread count 61440 */
	m_nWorkItems = block_size != 0 ? block_size : DWARF_WG_ALG;
	m_nThreads = threads != 0 ? threads : 61440;
}

void NQueenSolver::Launch(int j, size_t work_dim, int board_size, int level, int threads, int pitch, unsigned int *params,
		unsigned int *results, unsigned int *forbidden, int *global_index)
{
	sycl::nd_range<1> ndr(work_dim, m_nWorkItems);
	if(m_bCPU) {
		if(j == 0) nqk::cpu_nqueen1(q, ndr, board_size, level, threads, pitch, params, results, forbidden);
		else       nqk::cpu_nqueen(q, ndr, board_size, level, threads, pitch, params, results, forbidden);
	} else if(m_bAtomics) {
		if(j == 0) nqk::gpu_nqueen1<true>(q, ndr, board_size, level, threads, pitch, params, results, forbidden, global_index);
		else       nqk::gpu_nqueen<true>(q, ndr, board_size, level, threads, pitch, params, results, forbidden, global_index);
	} else {
		if(j == 0) nqk::gpu_nqueen1<false>(q, ndr, board_size, level, threads, pitch, params, results, forbidden, global_index);
		else       nqk::gpu_nqueen<false>(q, ndr, board_size, level, threads, pitch, params, results, forbidden, global_index);
	}
}

void NQueenSolver::Collect(const unsigned int *result_buffer, std::vector<unsigned int> &results, int max_pitch, int count,
		long long &solutions, long long &unique_solutions, int errcode)
{
	dsycl::read(q, &results[0], result_buffer, (size_t)max_pitch * sizeof(int) * 4);
	for(int k = 0; k < count; k++) {
		if(results[k + max_pitch * 2] != results[k + max_pitch * 3]) {
			throw NQError(errcode);
		}
		solutions += results[k];
		unique_solutions += results[k + max_pitch];
	}
}

inline int bit_scan_host(unsigned int x) { return nqk::bit_scan(x); }

long long NQueenSolver::Compute(int board_size, long long* unique)
{
	long long total = 10000000000LL;
	total /= 10;	// for atomics
	int level = 0;
	int i = board_size;
	while(total > 0 && i > 0) {
		total /= ((i + 1) / 2);
		i--;
		level++;
	}


	if(level > board_size - 2) {
		level = board_size - 2;
	}

	if(level > 11) {
		level = 11;
	}

	const int threads = m_bAtomics ? m_nThreads * 16 : m_nThreads;
	const int max_threads = threads;
	const int max_pitch = (max_threads + 15) & ~0xf;

	/* device buffers (inside the timed region, as in the OpenCL tree) */
	DevBuf<unsigned int> param_buffer(q, (size_t)max_pitch * (32 + 32));
	DevBuf<unsigned int> result_buffer(q, (size_t)max_pitch * 64);
	DevBuf<unsigned int> forbidden_buffer(q, 32);
	DevBuf<int> global_index_buffer(q, 1);

	std::vector<unsigned int> mask_vector(max_pitch * (4 + 32));
	std::vector<unsigned int> results(max_pitch * 4);
	bool forbidden_written = false;
	bool pending = false;     /* a launched batch whose results have not been collected (ocdTempEvent != 0) */
	int last_total = 0;

	long long solutions = 0;
	long long unique_solutions = 0;
	bool has_data = false;

	unsigned int board_mask = (1 << board_size) - 1;
	int total_size = 0;
	for(int j = 0; j < board_size / 2; j++) {
		unsigned int masks[32];
		unsigned int left_masks[32];
		unsigned int right_masks[32];
		unsigned int ms[32];
		unsigned int ns[32];
		unsigned int forbidden[32];
		unsigned int border_mask = 0;
		int idx = 0;
		int i = 0;

		masks[0] = (1 << j);
		left_masks[0] = 1 << (j + 1);
		right_masks[0] = (1 << j) >> 1;
		ms[0] = masks[0] | left_masks[0] | right_masks[0];
		ns[0] = (1 << j);

		for(int k = 0; k < j; k++) {
			border_mask |= (1 << k);
			border_mask |= (1 << (board_size - k - 1));
		}

		for(int k = 0; k < board_size; k++) {
			if(k == board_size - 2) {
				forbidden[k] = border_mask;
			}
			else if((k + 1) < j || (k + 1) > board_size - j - 1) {
				forbidden[k] = 1 | (1 << (board_size - 1));
			}
			else {
				forbidden[k] = 0;
			}
		}
		forbidden[board_size - 1] = 0xffffffff;

		forbidden_written = false;
		while(i >= 0) {
			if(j == 0) {
				if(i >= 1) {
					unsigned int m = ms[i] | (i + 1 < idx ? 2 : 0);
					ns[i + 1] = (m + 1) & ~m;
				}
				else {
					ns[i + 1] = ((ms[i] + 1) & ~ms[i]);
					if(i == 0) {
						idx = bit_scan_host(ns[i + 1]);
					}
				}
			}
			else {
				unsigned int m = ms[i] | forbidden[i];
				ns[i + 1] = (m + 1) & ~m;
			}

			if(i == board_size - level - 1) {
				mask_vector[total_size] = masks[i];
				mask_vector[total_size + max_pitch] = left_masks[i];
				mask_vector[total_size + max_pitch * 2] = right_masks[i];
				if(j == 0) {
					mask_vector[total_size + max_pitch * 3] = idx - i < 0 ? 0 : idx - i;
				}
				else {
					mask_vector[total_size + max_pitch * 3] = j;
				}
				for(int k = 0; k <= i; k++) {
					mask_vector[total_size + max_pitch * (k + 4)] = ns[k];
				}
				total_size++;
				if(total_size == max_threads) {
					if(pending) { Collect(result_buffer.p, results, max_pitch, last_total, solutions, unique_solutions, 1); pending = false; }
					if(!forbidden_written) {
						dsycl::write(q, forbidden_buffer.p, forbidden + board_size - level - 1, (level + 1) * sizeof(int));
						forbidden_written = true;
					}
					dsycl::write(q, param_buffer.p, &mask_vector[0], (size_t)max_pitch * sizeof(int) * (4 + 32));
					size_t work_dim = m_nThreads;
					if(m_bAtomics) { int num_threads = (int)work_dim; dsycl::write(q, global_index_buffer.p, &num_threads, sizeof(int)); }
					Launch(j, work_dim, board_size, level, threads, max_pitch, param_buffer.p, result_buffer.p, forbidden_buffer.p, global_index_buffer.p);
					pending = true;
					last_total = threads;

					if(total_size > threads) {
						for(int k = 0; k < 4 + board_size - level; k++) {
							memcpy(&mask_vector[max_pitch * k], &mask_vector[threads + max_pitch * k], (total_size - threads) * sizeof(int));
						}
					}
					total_size -= threads;
				}

				i--;
			}
			else if((ns[i + 1] & board_mask) != 0) {
				ms[i] |= ns[i + 1];
				masks[i+1] = masks[i] | ns[i + 1];
				left_masks[i+1] = (left_masks[i] | ns[i + 1]) << 1;
				right_masks[i+1] = (right_masks[i] | ns[i + 1]) >> 1;
				ms[i+1] = masks[i+1] | left_masks[i+1] | right_masks[i + 1];
				i++;
			}
			else {
				i--;
			}
		}

		while(total_size > 0) {
			for(int k = total_size; k < max_threads; k++) {
				mask_vector[k] = 0xffffffff;
				mask_vector[k + max_pitch] = 0xffffffff;
				mask_vector[k + max_pitch * 2] = 0xffffffff;
				mask_vector[k + max_pitch * 3] = 0;
				mask_vector[k + max_pitch * 4] = 0;
			}
			if(pending) { Collect(result_buffer.p, results, max_pitch, last_total, solutions, unique_solutions, 2); pending = false; }

			int t_size = total_size > threads ? threads : total_size;
			if(!forbidden_written) {
				dsycl::write(q, forbidden_buffer.p, forbidden + board_size - level - 1, (level + 1) * sizeof(int));
				forbidden_written = true;
			}
			dsycl::write(q, param_buffer.p, &mask_vector[0], (size_t)max_pitch * sizeof(int) * (4 + 32));

			size_t work_dim = t_size < m_nThreads ? t_size : m_nThreads;
			size_t n = m_nWorkItems;
			if(work_dim % n != 0) work_dim += n - work_dim % n;
			if(m_bAtomics) { int num_thread = (int)work_dim; dsycl::write(q, global_index_buffer.p, &num_thread, sizeof(int)); }
			Launch(j, work_dim, board_size, level, t_size, max_pitch, param_buffer.p, result_buffer.p, forbidden_buffer.p, global_index_buffer.p);
			pending = true;
			last_total = t_size;

			if(total_size > t_size) {
				for(int k = 0; k < 4 + board_size - level; k++) {
					memcpy(&mask_vector[max_pitch * k], &mask_vector[t_size + max_pitch * k], (total_size - t_size) * sizeof(int));
				}
			}
			total_size -= t_size;
		}
	}

	if(pending) { Collect(result_buffer.p, results, max_pitch, last_total, solutions, unique_solutions, 3); pending = false; }

	if(unique != 0) {
		*unique = unique_solutions;
	}
	return solutions;
}
