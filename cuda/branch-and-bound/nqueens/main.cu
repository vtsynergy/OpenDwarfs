// N-queen solver for OpenCL
// Ping-Che Chen


#include "dwarf_cuda.h"
#include "dwarf_timer.h"
#include "nqueen_cpu.h"
#include "nqueen_cuda.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <cstdlib>

int main(int argc, char** argv)
{

	dwarf_cuda_init(&argc, &argv);

	std::cerr << "N-Queen solver for CUDA\n";
	std::cerr << "Ping-Che Chen\n\n";
	if(argc < 2) {
		std::cerr << "Usage: " << argv[0] << " [options] N\n";
		std::cerr << "\tN: board size (1 ~ 32)\n";
		std::cerr << "\t-cpu: use CPU (multi-threaded on Windows)\n";
		//std::cerr << "\t-clcpu: use OpenCL CPU device\n";
		std::cerr << "\t-prof: enable profiler\n";
		std::cerr << "\t-threads #: set number of threads to #\n";
		std::cerr << "\t-blocksize #: set size of thread blocks to #\n";
		//std::cerr << "\t-platform #: select platform #\n";
		//std::cerr << "\t-device #: select device # (default: use all devices)\n";
		std::cerr << "\t-local: use local memory for arrays (default: off)\n";
		std::cerr << "\t-noatomics: do not use global atomics\n";
		std::cerr << "\t-novec: do not use vectorization\n";
		std::cerr << "\t-vec4: use 4D vectors instead of 2D (only when vectorized- default: off)\n";
		return 0;
	}

	// handle options
	bool force_cpu = false;
	//bool use_clcpu = false;
	//cl_uint platform_index = 0;
	//int device_index = -1;
	bool profiling = false;
	int threads = 0;
	int block_size = 0;
	bool local = false;//default OFF (was true)
	bool noatomics = false;
	bool novec = false;
	bool use_vec4 = false;

	int start = 1;
	while(start < argc - 1) {
		if(std::strcmp(argv[start], "-cpu") == 0) {
			force_cpu = true;
		}
		//else if(std::strcmp(argv[start], "-clcpu") == 0) {
		//	use_clcpu = true;
		//}
		//else if(std::strcmp(argv[start], "-prof") == 0) {
		//	profiling = true;
		//}
		else if(std::strcmp(argv[start], "-threads") == 0 && start < argc - 2) {
			threads = std::atoi(argv[start + 1]);
			start++;
		}
		else if(std::strcmp(argv[start], "-blocksize") == 0 && start < argc - 2) {
			block_size = std::atoi(argv[start + 1]);
			start++;
		}
		//else if(std::strcmp(argv[start], "-platform") == 0 && start < argc - 2) {
		//	platform_index = std::atoi(argv[start + 1]);
		//	start++;
		//}
		//else if(std::strcmp(argv[start], "-device") == 0 && start < argc - 2) {
		//	device_index = std::atoi(argv[start + 1]);
		//	start++;
		//}
		else if(std::strcmp(argv[start], "-local") == 0) {
			local = true;
		}
		//else if(std::strcmp(argv[start], "-nolocal") == 0) {
		//	local = false;
		//}
		else if(std::strcmp(argv[start], "-noatomics") == 0) {
			noatomics = true;
		}
		else if(std::strcmp(argv[start], "-novec") == 0) {
			novec = true;
		}
		else if(std::strcmp(argv[start], "-vec4") == 0) {
			use_vec4 = true;
		}
		else {
			std::cerr << "Unknown option " << argv[start] << "\n";
		}

		start ++;
	}

	int board_size = std::atoi(argv[start]);
	if(board_size < 1 || board_size > 32) {
		std::cerr << "Inalid board size (only 1 ~ 32 allowed)\n";
		return 0;
	}

	clock_t start_time, end_time;
	dwarf_init("nqueens", "branch-and-bound", "cuda");
	dwarf_set_size(board_size, 0);
	long long solutions = 0;
	long long unique_solutions = 0;
	if(force_cpu) {
		start_time = std::clock();
		for(int rep = dwarf_first_rep(); rep < dwarf_reps(); rep++) { dwarf_region_begin(); solutions = nqueen_cpu(board_size, &unique_solutions); dwarf_region_end(rep); }
		end_time = std::clock();
	}
	else {
		try {
			NQueenSolver nqueen(threads, block_size);
			std::cerr << "\tUsing " << nqueen.GetThreads(0) << " threads\n";
			std::cerr << "\tBlock size = " << nqueen.GetBlockSize(0) << " threads\n";
			std::cerr << "\tUsing global atomics\n";
			(void)local; (void)noatomics; (void)novec; (void)use_vec4; (void)profiling;

			/* common timer: Compute() = host enumeration + device buffer creation/uploads + kernels + result reads */
			for(int rep = dwarf_first_rep(); rep < dwarf_reps(); rep++) { dwarf_region_begin(); solutions = nqueen.Compute(board_size, &unique_solutions); dwarf_region_end(rep); }
		}
		catch(CLError x)
		{
			if(x.GetErrorNo() >= 1 && x.GetErrorNo() <= 3) {
				std::cerr << x.GetErrorNo() << " CUDA kernel execution failed\n";
			}
			else {
				std::cerr << x << "\n";
			}
		}
	}

	std::cerr << board_size << "-queen has " << solutions << " solutions (" << unique_solutions << " unique)\n";
	//std::cerr << "Time used: " << std::setprecision(3) << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC << "s\n";

	/* [FIX-OCL-14] the app had no correctness check (and the error-reporting chain `if,if,if,else` printed the generic message
	 * after errors 1 and 2 as well). Added: comparison with the original single-threaded CPU solver, run untimed. */
	long long ref_unique = 0, ref_sol = nqueen_cpu(board_size, &ref_unique);
	bool ok = (ref_sol == solutions && ref_unique == unique_solutions);
	std::cout << (ok ? "SUCCESS" : "FAILURE") << ": " << solutions << "/" << unique_solutions << " vs CPU reference " << ref_sol << "/" << ref_unique << "\n";
	return dwarf_report(ok);
}
