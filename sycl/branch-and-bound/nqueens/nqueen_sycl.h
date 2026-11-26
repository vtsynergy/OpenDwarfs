// N-queen solver for SYCL (interface of opencl/branch-and-bound/nqueens/nqueen_cl.h)
#ifndef NQUEEN_SYCL_H
#define NQUEEN_SYCL_H
#include <vector>
#include "dwarf_sycl.hpp"

struct NQError { int code; NQError(int c) : code(c) {} };

class NQueenSolver
{
public:
	NQueenSolver(sycl::queue &q, int threads, int block_size, bool force_no_atomics);
	long long Compute(int board_size, long long* unique);
	int GetThreads() const { return m_nThreads; }
	int GetBlockSize() const { return m_nWorkItems; }
	bool AtomicsEnabled() const { return m_bAtomics; }
private:
	void Launch(int j, size_t work_dim, int board_size, int level, int threads, int pitch, unsigned int *params,
			unsigned int *results, unsigned int *forbidden, int *global_index);
	void Collect(const unsigned int *result_buffer, std::vector<unsigned int> &results, int max_pitch, int count,
			long long &solutions, long long &unique_solutions, int errcode);
	sycl::queue &q;
	bool m_bCPU, m_bAtomics;
	int m_nThreads, m_nWorkItems;
};
#endif
