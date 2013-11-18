// N-queen solver for OpenCL
// Ping-Che Chen


#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include "../../include/rdtsc.h"
#include "../../include/common_args.h"
#include "nqueen_cpu.h"
#include "nqueen_cl.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <cstdlib>

int main(int argc, char** argv)
{

	ocd_init(&argc, &argv, NULL);
	ocd_initCL();

	std::cerr << "N-Queen solver for OpenCL\n";
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
	long long solutions = 0;
	long long unique_solutions = 0;
	if(force_cpu) {
		start_time = std::clock();
		solutions = nqueen_cpu(board_size, &unique_solutions);
		end_time = std::clock();
	}
	else {
		cl_int err;
		//cl_uint num;
		//err = clGetPlatformIDs(0, 0, &num);
		//if(err != CL_SUCCESS) {
		//	std::cerr << "Unable to get platforms\n";
		//	return 0;
		//}

		//std::vector<cl_platform_id> platforms(num);
		//err = clGetPlatformIDs(num, &platforms[0], &num);
		//if(err != CL_SUCCESS) {
		//	std::cerr << "Unable to get platforms\n";
		//	return 0;
		//}

		// display all platform data
		//for(cl_uint i = 0; i < num; i++) {
		//	size_t size;
		//	clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, 0, &size);
		//	std::string name;
		//	name.reserve(size + 1);
		//	clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, size, &name[0], 0);
		//	std::cerr << "Platform [" << i << "]: " << name.c_str() << "\n";
		//}

		//if(platform_index >= num) {
		//	platform_index = num - 1;
		//}

		//std::cerr << "Select platform " << platform_index << "\n";

		//cl_context context;
		//cl_context_properties prop[] = { CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>(platforms[platform_index]), 0 };
		//context = clCreateContextFromType(prop, use_clcpu ? CL_DEVICE_TYPE_CPU : CL_DEVICE_TYPE_GPU, 0, 0, &err);
		//if(err != CL_SUCCESS) {
		//	std::cerr << "Unable to create context\n";
		//	return 0;
		//}

		//if(use_clcpu) {
		//	std::cerr << "Using CPU device\n";
		//}
		//else {
		//	std::cerr << "Using GPU device\n";
		//}

		// show device list
		size_t num_devices;
		//err = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, 0, &num_devices);
		//if(err != CL_SUCCESS) {
		//	std::cerr << "Unable to get device information\n";
		//	return 0;
		//}

		num_devices=1;//In OpenDwarfs we only work with one device at a time.
		std::vector<cl_device_id> devices(num_devices / sizeof(cl_device_id));
		//err = clGetContextInfo(context, CL_CONTEXT_DEVICES, num_devices, &devices[0], &num_devices);
		//if(err != CL_SUCCESS) {
		//	std::cerr << "Unable to get device information\n";
		//	return 0;
		//}

		//if(device_index >= 0 && device_index < devices.size()) {
		//	cl_device_id device_id = devices[device_index];
		devices.clear();
		devices.resize(1);
		devices[0] = device_id;
		//}

		//context and devices[0]=device_id are now as gotten from ocd_initCL();
		//ocd_initCL(); has also return the commands command queue.


		try {
			NQueenSolver nqueen(context, devices, profiling, threads, block_size, local, noatomics, novec, use_vec4);
			//std::cerr << "Number of devices:" << devices.size() << "\n";
			for(int i = 0; i < devices.size(); i++) {
				size_t name_length;
				err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 0, 0, &name_length);
				if(err == CL_SUCCESS) {
					std::string name;
					name.resize(name_length + 1);
					clGetDeviceInfo(devices[i], CL_DEVICE_NAME, name_length, &name[0], &name_length);
					name[name_length] = 0;
					std::cerr << "Device " << i << ": " << name.c_str() << "\n";
					std::cerr << "\tUsing " << nqueen.GetThreads(i) << " threads\n";
					std::cerr << "\tBlock size = " << nqueen.GetBlockSize(i) << " threads\n";
					if(nqueen.AtomicsEnabled(i)) {
						std::cerr << "\tUsing global atomics\n";
					}

					if(nqueen.VectorizationEnabled(i)) {
						std::cerr << "\tUsing vectorization\n";

						if(use_vec4) {
							std::cerr << "\tUse 4D vectors\n";
						}
						else {
							std::cerr << "\tUse 2D vectors\n";
						}
					}
				}
			}

			//start_time = std::clock();
			solutions = nqueen.Compute(board_size, &unique_solutions);
			//end_time = std::clock();

			//if(profiling) {
			//	for(int i = 0; i < devices.size(); i++) {
			//		std::cerr << "Profile time for device " << i << ": " << nqueen.GetProfilingTime(i) << "ns\n";
			//	}
			//}
		}
		catch(CLError x)
		{
			if(x.GetErrorNo() == 1) {
				std::cerr << "1 OpenCL kernel execution failed\n";
			}
			if(x.GetErrorNo() == 2) {
				std::cerr << "2 OpenCL kernel execution failed\n";
			}
			if(x.GetErrorNo() == 3) {
				std::cerr << "3 OpenCL kernel execution failed\n";
			}
			else {
				std::cerr << x << "\n";
			}
		}

		clReleaseContext(context);
	}

	std::cerr << board_size << "-queen has " << solutions << " solutions (" << unique_solutions << " unique)\n";
	//std::cerr << "Time used: " << std::setprecision(3) << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC << "s\n";

	ocd_finalize();
	return 0;
}
