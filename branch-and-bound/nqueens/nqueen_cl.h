// N-queen solver for OpenCL
// Ping-Che Chen


#ifndef NQUEEN_CL_H
#define NQUEEN_CL_H


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include <iostream>
#include <vector>


class CLError
{
public:

	CLError(cl_int err, int line = 0) : m_ErrNo(err), m_Line(line) {}

	cl_int GetErrorNo() const { return m_ErrNo; }
	int GetErrorLine() const { return m_Line; }

private:

	cl_int m_ErrNo;
	int m_Line;
};


inline std::ostream& operator<<(std::ostream& stream, const CLError& x)
{
	stream << "OpenCL error: " << x.GetErrorNo();
	if(x.GetErrorLine() != 0) {
		stream << " (line: " << x.GetErrorLine() << ")";
	}

	return stream;
}


class NQueenSolver
{
	struct solver_info
	{
		solver_info() : m_Queue(0), m_Program(0), m_NQueen(0), m_NQueen1(0), m_ParamBuffer(0), m_ResultBuffer(0), m_ForbiddenBuffer(0), m_GlobalIndex(0), m_Event(0),
			m_bCPU(false), m_bEnableAtomics(false), m_bEnableVectorize(false), m_bEnableChar(false), m_bEnableLocal(false), m_nMaxWorkItems(0), m_nThreads(0)
		{
		}

		~solver_info()
		{
			if(m_Event != 0) {
				clReleaseEvent(m_Event);
			}

			if(m_ParamBuffer != 0) {
				clReleaseMemObject(m_ParamBuffer);
			}

			if(m_ResultBuffer != 0) {
				clReleaseMemObject(m_ResultBuffer);
			}

			if(m_ForbiddenBuffer != 0) {
				clReleaseMemObject(m_ForbiddenBuffer);
			}

			if(m_GlobalIndex != 0) {
				clReleaseMemObject(m_GlobalIndex);
			}

			if(m_NQueen != 0) {
				clReleaseKernel(m_NQueen);
			}

			if(m_NQueen1 != 0) {
				clReleaseKernel(m_NQueen1);
			}

			if(m_Program != 0) {
				clReleaseProgram(m_Program);
			}

			if(m_Queue != 0) {
				clReleaseCommandQueue(m_Queue);
			}
		}

		cl_command_queue m_Queue;

		cl_program m_Program;
		cl_kernel m_NQueen;
		cl_kernel m_NQueen1;

		cl_mem m_ParamBuffer;
		cl_mem m_ResultBuffer;
		cl_mem m_ForbiddenBuffer;
		cl_mem m_GlobalIndex;

		cl_event m_Event;

		bool m_bCPU;
		bool m_bEnableAtomics;
		bool m_bEnableVectorize;
		bool m_bEnableChar;
		bool m_bEnableLocal;

		size_t m_nMaxWorkItems;
		int m_nThreads;

		cl_ulong m_TotalTime;
		int m_nLastTotalSize;
	};


public:

	NQueenSolver(cl_context context, std::vector<cl_device_id> devices, bool profiling = false, int threads = 0, int block_size = 0, bool force_local = false, bool force_no_atomics = false, bool force_no_vec = false, bool force_vec4 = false);
	~NQueenSolver();

	long long Compute(int board_size, long long* unique);
	cl_ulong GetProfilingTime(int idx) { return m_SolverInfo[idx].m_TotalTime; }
	int GetThreads(int idx) { return m_SolverInfo[idx].m_nThreads; }
	int GetBlockSize(int idx) { return m_SolverInfo[idx].m_nMaxWorkItems / (m_SolverInfo[idx].m_bEnableVectorize ? (m_bForceVec4 ? 4 : 2) : 1); }
	bool AtomicsEnabled(int idx) { return m_SolverInfo[idx].m_bEnableAtomics; }
	bool VectorizationEnabled(int idx) { return m_SolverInfo[idx].m_bEnableVectorize; }

private:

	void InitKernels(int i, int block_size);
	void BuildProgram(int i, const std::string& program, int vector_width, int work_items);

	cl_context m_Context;
	std::vector<cl_device_id> m_Devices;

	bool m_bProfiling;
	bool m_bForceLocal;
	bool m_bForceNoAtomics;
	bool m_bForceNoVectorization;
	bool m_bForceVec4;
	bool m_bDoubleQueue;

	std::vector<solver_info> m_SolverInfo;
};


#endif
