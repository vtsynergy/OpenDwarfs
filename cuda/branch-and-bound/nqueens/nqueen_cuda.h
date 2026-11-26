// N-queen solver for OpenCL
// Ping-Che Chen


#ifndef NQUEEN_CL_H
#define NQUEEN_CL_H


#include <iostream>
#include <vector>


class CLError
{
public:
	CLError(int err, int line = 0) : m_ErrNo(err), m_Line(line) {}
	int GetErrorNo() const { return m_ErrNo; }
	int GetErrorLine() const { return m_Line; }
private:
	int m_ErrNo;
	int m_Line;
};

inline std::ostream& operator<<(std::ostream& stream, const CLError& x)
{
	stream << "CUDA error: " << x.GetErrorNo();
	return stream;
}


class NQueenSolver
{
	struct solver_info
	{
		solver_info() : m_ParamBuffer(0), m_ResultBuffer(0), m_ForbiddenBuffer(0), m_GlobalIndex(0), m_nMaxWorkItems(0), m_nThreads(0), m_nLastTotalSize(0) {}
		int *m_ParamBuffer, *m_ResultBuffer, *m_ForbiddenBuffer, *m_GlobalIndex;
		size_t m_nMaxWorkItems;
		int m_nThreads;
		int m_nLastTotalSize;
	};


public:

	NQueenSolver(int threads = 0, int block_size = 0);
	~NQueenSolver();

	long long Compute(int board_size, long long* unique);
	int GetThreads(int idx) { return m_SolverInfo[idx].m_nThreads; }
	int GetBlockSize(int idx) { return m_SolverInfo[idx].m_nMaxWorkItems; }
	bool AtomicsEnabled(int idx) { return true; }
	bool VectorizationEnabled(int idx) { return false; }

private:

	std::vector<solver_info> m_SolverInfo;
};


#endif
