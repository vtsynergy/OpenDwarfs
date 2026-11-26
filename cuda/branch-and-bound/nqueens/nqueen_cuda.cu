// N-queen solver for OpenCL
// Ping-Che Chen
// adapted for OpenDwarfs


#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include "nqueen_cuda.h"
#include "dwarf_cuda.h"
#include "nqueen_kernel.cuh"

//#define CHECK_ERROR(err) { if(err != CL_SUCCESS) throw CLError(err, __LINE__); }

/* Only the configuration the OpenCL host selects on the NVIDIA GPU exists here: global atomics, no vectorisation, 256-thread
 * blocks [FIX-OCL-34] and 61440 persistent threads (240 blocks). */
NQueenSolver::NQueenSolver(int threads, int block_size)
{
	m_SolverInfo.resize(1);
	m_SolverInfo[0].m_nMaxWorkItems = block_size != 0 ? block_size : DWARF_WG_ALG;
	m_SolverInfo[0].m_nThreads = threads != 0 ? threads : 61440;
	if(m_SolverInfo[0].m_nMaxWorkItems != DWARF_WG_ALG) {
		fprintf(stderr, "only block size %d is supported (the kernels require it, as reqd_work_group_size does in OpenCL)\n", DWARF_WG_ALG);
		exit(1);
	}
}


NQueenSolver::~NQueenSolver()
{
}


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


long long NQueenSolver::Compute(int board_size, long long* unique)
{
	// estimate amount of levels need to be computed on the device
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

//	std::cerr << "start: " << clock() << "\n";

//	if(m_bCPU) {
//		level = board_size - 2;
//	}

	std::vector<int> threads(m_SolverInfo.size()); //Will be always 1 in OpenDwarfs
	int max_threads = 0;
	int max_pitch;
	bool dev_pending = false; /* OpenCL: ocdTempEvent != 0 (a launch whose results have not been read) */
	for(int i = 0; i < threads.size(); i++) {
		if(true) { /* atomics */
			threads[i] = m_SolverInfo[i].m_nThreads * 16;
		}
		else {
			threads[i] = m_SolverInfo[i].m_nThreads;
		}

		if(max_threads < threads[i]) {
			max_threads = threads[i];
		}
	}

	max_pitch = (max_threads + 15) & ~0xf;

	for(int i = 0; i < threads.size(); i++) { //Will be always 1 in OpenDwarfs
		// create data buffer
		if(m_SolverInfo[i].m_ParamBuffer != 0) cudaFree(m_SolverInfo[i].m_ParamBuffer);
		CUDA_CHECK(cudaMalloc((void**)&m_SolverInfo[i].m_ParamBuffer, max_pitch * sizeof(int) * (32 + 32)));

		if(m_SolverInfo[i].m_ResultBuffer != 0) cudaFree(m_SolverInfo[i].m_ResultBuffer);
		CUDA_CHECK(cudaMalloc((void**)&m_SolverInfo[i].m_ResultBuffer, max_pitch * sizeof(int) * 64));

		if(m_SolverInfo[i].m_ForbiddenBuffer != 0) cudaFree(m_SolverInfo[i].m_ForbiddenBuffer);
		CUDA_CHECK(cudaMalloc((void**)&m_SolverInfo[i].m_ForbiddenBuffer, 32 * sizeof(int)));

		if(m_SolverInfo[i].m_GlobalIndex != 0) cudaFree(m_SolverInfo[i].m_GlobalIndex);
		CUDA_CHECK(cudaMalloc((void**)&m_SolverInfo[i].m_GlobalIndex, sizeof(int)));
	}

	std::vector<unsigned int> mask_vector(max_pitch * (4 + 32));
	std::vector<unsigned int> results(max_pitch * 4);
	std::vector<bool> forbidden_written(threads.size());

	long long solutions = 0;
	long long unique_solutions = 0;
	bool has_data = false;

	unsigned int board_mask = (1 << board_size) - 1;
	int total_size = 0;
	int last_total_size = 0;
	int device_idx = 0;
	for(int j = 0; j < board_size / 2; j++) {
	// only do nqueen1
//	int j = 1;
//	{
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

		//cl_kernel queen = (j == 0 ? m_NQueen1 : m_NQueen);

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

		for(int k = 0; k < threads.size(); k++) {
			forbidden_written[k] = false;
		}
		while(i >= 0) {
			if(j == 0) {
				if(i >= 1) {
					unsigned int m = ms[i] | (i + 1 < idx ? 2 : 0);
					ns[i + 1] = (m + 1) & ~m;
				}
				else {
					ns[i + 1] = ((ms[i] + 1) & ~ms[i]);
					if(i == 0) {
						idx = bit_scan(ns[i + 1]);
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
					// check rotation
					mask_vector[total_size + max_pitch * 3] = j;
				}
				for(int k = 0; k <= i; k++) {
					mask_vector[total_size + max_pitch * (k + 4)] = ns[k];
				}
				total_size++;
				if(total_size == max_threads) {
//					std::cerr << "device A[" << device_idx << "]: " << clock() << "\n";
//					if(has_data) {
//						err = clEnqueueReadBuffer(m_Queue, result_buffer, CL_TRUE, 0, threads * sizeof(int) * 2, &results[0], 0, 0, 0);
//						CHECK_ERROR(err);
//
//						if(m_bProfiling) {
//							cl_ulong start, end;
//							clGetEventProfilingInfo(profile_event, CL_PROFILING_COMMAND_START, sizeof(start), &start, 0);
//							clGetEventProfilingInfo(profile_event, CL_PROFILING_COMMAND_END, sizeof(end), &end, 0);
//
//							m_TotalTime += end - start;
//
//							clReleaseEvent(profile_event);
//							profile_event = 0;
//						}
//
//						for(int k = 0; k < last_total_size; k++) {
//							solutions += results[k];
//							unique_solutions += results[k + threads];
//						}
//					}

					int last_device_idx = device_idx; (void)last_device_idx;
					device_idx = 0; /* one device */

					if(dev_pending) {
//						std::cerr << "get data from device[" << device_idx << "]: " << clock() << "\n";

						// get data from the device
						CUDA_CHECK(cudaMemcpy(&results[0], m_SolverInfo[device_idx].m_ResultBuffer, max_pitch * sizeof(int) * 4, cudaMemcpyDeviceToHost));
						CUDA_CHECK(cudaDeviceSynchronize());
						//CHECK_ERROR(err);

						//if(m_bProfiling) {
						//	cl_ulong start, end;
						//	clGetEventProfilingInfo(m_SolverInfo[device_idx].m_Event, CL_PROFILING_COMMAND_START, sizeof(start), &start, 0);
						//	clGetEventProfilingInfo(m_SolverInfo[device_idx].m_Event, CL_PROFILING_COMMAND_END, sizeof(end), &end, 0);

						//	m_SolverInfo[device_idx].m_TotalTime += end - start;
//							std::cerr << "device time: " << end - start << "\n";
						//}

						dev_pending = false;

						for(int k = 0; k < m_SolverInfo[device_idx].m_nLastTotalSize; k++) {
							if(results[k + max_pitch * 2] != results[k + max_pitch * 3]) {
								throw CLError(1);
							}

							solutions += results[k];
							unique_solutions += results[k + max_pitch];
						}
					}

//					std::cerr << "device [" << device_idx << "]: " << " prepare launch: " << clock() << "\n";

					const bool use_nqueen1 = (j == 0);
					int arg_board_size = board_size;
					int arg_level = level;
					int arg_threads = threads[device_idx];
					int arg_pitch = max_pitch;

					if(!forbidden_written[device_idx]) {
						CUDA_CHECK(cudaMemcpy(m_SolverInfo[device_idx].m_ForbiddenBuffer, forbidden + board_size - level - 1, (level + 1) * sizeof(int), cudaMemcpyHostToDevice));
						CUDA_CHECK(cudaDeviceSynchronize());
						//CHECK_ERROR(err);
						forbidden_written[device_idx] = true;
					}

					CUDA_CHECK(cudaMemcpy(m_SolverInfo[device_idx].m_ParamBuffer, &mask_vector[0], max_pitch * sizeof(int) * (4 + 32), cudaMemcpyHostToDevice));
					CUDA_CHECK(cudaDeviceSynchronize());
					//CHECK_ERROR(err);

					size_t work_dim = m_SolverInfo[device_idx].m_nThreads;
					int num_threads_ = (int)work_dim;
					CUDA_CHECK(cudaMemcpy(m_SolverInfo[device_idx].m_GlobalIndex, &num_threads_, sizeof(int), cudaMemcpyHostToDevice));
					CUDA_CHECK(cudaDeviceSynchronize());

					{ /* 1-D launch of work_dim threads in blocks of 256 */
						dim3 gridDim_((unsigned int)(work_dim / DWARF_WG_ALG)), blockDim_(DWARF_WG_ALG);
						if(use_nqueen1) nqueen1<<<gridDim_, blockDim_>>>(arg_board_size, arg_level, arg_threads, arg_pitch, (uint*)m_SolverInfo[device_idx].m_ParamBuffer, (uint*)m_SolverInfo[device_idx].m_ResultBuffer, (const uint*)m_SolverInfo[device_idx].m_ForbiddenBuffer, m_SolverInfo[device_idx].m_GlobalIndex);
						else nqueen<<<gridDim_, blockDim_>>>(arg_board_size, arg_level, arg_threads, arg_pitch, (uint*)m_SolverInfo[device_idx].m_ParamBuffer, (uint*)m_SolverInfo[device_idx].m_ResultBuffer, (const uint*)m_SolverInfo[device_idx].m_ForbiddenBuffer, m_SolverInfo[device_idx].m_GlobalIndex);
						CUDA_SYNC_LAUNCH();
						dev_pending = true;
					}

//					has_data = true;
//					last_total_size = total_size;
					m_SolverInfo[device_idx].m_nLastTotalSize = threads[device_idx];

//					std::cerr << "device [" << device_idx << "]: " << " launch: " << clock() << "\n";

					if(total_size > threads[device_idx]) {
						// adjust the data array
						for(int k = 0; k < 4 + board_size - level; k++) {
							memcpy(&mask_vector[max_pitch * k], &mask_vector[threads[device_idx] + max_pitch * k], (total_size - threads[device_idx]) * sizeof(int));
						}
					}

					total_size -= threads[device_idx];
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

//		if(has_data) {
//			err = clEnqueueReadBuffer(m_Queue, result_buffer, CL_TRUE, 0, threads * sizeof(int) * 2, &results[0], 0, 0, 0);
//			CHECK_ERROR(err);
//
//			if(m_bProfiling) {
//				cl_ulong start, end;
//				clGetEventProfilingInfo(profile_event, CL_PROFILING_COMMAND_START, sizeof(start), &start, 0);
//				clGetEventProfilingInfo(profile_event, CL_PROFILING_COMMAND_END, sizeof(end), &end, 0);
//
//				m_TotalTime += end - start;
//
//				clReleaseEvent(profile_event);
//				profile_event = 0;
//			}
//
//			for(int k = 0; k < last_total_size; k++) {
//				solutions += results[k];
//				unique_solutions += results[k + threads];
//			}
//
//			std::cout << solutions << "\n";
//			std::cout << unique_solutions << "\n";
//
//			has_data = false;
//		}

		while(total_size > 0) {
//			std::cerr << "device B[" << device_idx << "]: " << clock() << "\n";

			for(int k = total_size; k < max_threads; k++) {
				mask_vector[k] = 0xffffffff;
				mask_vector[k + max_pitch] = 0xffffffff;
				mask_vector[k + max_pitch * 2] = 0xffffffff;
				mask_vector[k + max_pitch * 3] = 0;
				mask_vector[k + max_pitch * 4] = 0;
			}

			int last_device_idx = device_idx; (void)last_device_idx;
			device_idx = 0; /* one device */

			if(dev_pending) {
//				std::cerr << "get data from device[" << device_idx << "]: " << clock() << "\n";

				// get data from the device
				CUDA_CHECK(cudaMemcpy(&results[0], m_SolverInfo[device_idx].m_ResultBuffer, max_pitch * sizeof(int) * 4, cudaMemcpyDeviceToHost));
				CUDA_CHECK(cudaDeviceSynchronize());
				//CHECK_ERROR(err);

				//if(m_bProfiling) {
				//	cl_ulong start, end;
				//	clGetEventProfilingInfo(m_SolverInfo[device_idx].m_Event, CL_PROFILING_COMMAND_START, sizeof(start), &start, 0);
				//	clGetEventProfilingInfo(m_SolverInfo[device_idx].m_Event, CL_PROFILING_COMMAND_END, sizeof(end), &end, 0);

				//	m_SolverInfo[device_idx].m_TotalTime += end - start;
//					std::cerr << "device time: " << end - start << "\n";
				//}

				dev_pending = false;

				for(int k = 0; k < m_SolverInfo[device_idx].m_nLastTotalSize; k++) {
					if(results[k + max_pitch * 2] != results[k + max_pitch * 3]) {
						throw CLError(2);
					}

					solutions += results[k];
					unique_solutions += results[k + max_pitch];
				}
			}

//			std::cerr << "device [" << device_idx << "]: " << " prepare launch: " << clock() << "\n";

			const bool use_nqueen1 = (j == 0);
			int arg_board_size = board_size;
			int arg_level = level;
			int t_size = total_size > threads[device_idx] ? threads[device_idx] : total_size;
			int arg_threads = t_size;
			int arg_pitch = max_pitch;

			if(!forbidden_written[device_idx]) {
				CUDA_CHECK(cudaMemcpy(m_SolverInfo[device_idx].m_ForbiddenBuffer, forbidden + board_size - level - 1, (level + 1) * sizeof(int), cudaMemcpyHostToDevice));
				CUDA_CHECK(cudaDeviceSynchronize());
				//CHECK_ERROR(err);
				forbidden_written[device_idx] = true;
			}

			CUDA_CHECK(cudaMemcpy(m_SolverInfo[device_idx].m_ParamBuffer, &mask_vector[0], max_pitch * sizeof(int) * (4 + 32), cudaMemcpyHostToDevice));
			CUDA_CHECK(cudaDeviceSynchronize());
			//CHECK_ERROR(err);

			size_t work_dim = t_size < m_SolverInfo[device_idx].m_nThreads ? t_size : m_SolverInfo[device_idx].m_nThreads;
			if(work_dim % DWARF_WG_ALG != 0) work_dim += DWARF_WG_ALG - work_dim % DWARF_WG_ALG;

			int num_threads_ = (int)work_dim;
			CUDA_CHECK(cudaMemcpy(m_SolverInfo[device_idx].m_GlobalIndex, &num_threads_, sizeof(int), cudaMemcpyHostToDevice));
			CUDA_CHECK(cudaDeviceSynchronize());

			{ /* 1-D launch of work_dim threads in blocks of 256 */
				dim3 gridDim_((unsigned int)(work_dim / DWARF_WG_ALG)), blockDim_(DWARF_WG_ALG);
				if(use_nqueen1) nqueen1<<<gridDim_, blockDim_>>>(arg_board_size, arg_level, arg_threads, arg_pitch, (uint*)m_SolverInfo[device_idx].m_ParamBuffer, (uint*)m_SolverInfo[device_idx].m_ResultBuffer, (const uint*)m_SolverInfo[device_idx].m_ForbiddenBuffer, m_SolverInfo[device_idx].m_GlobalIndex);
				else nqueen<<<gridDim_, blockDim_>>>(arg_board_size, arg_level, arg_threads, arg_pitch, (uint*)m_SolverInfo[device_idx].m_ParamBuffer, (uint*)m_SolverInfo[device_idx].m_ResultBuffer, (const uint*)m_SolverInfo[device_idx].m_ForbiddenBuffer, m_SolverInfo[device_idx].m_GlobalIndex);
				CUDA_SYNC_LAUNCH();
				dev_pending = true;
			}

//			std::cerr << "device [" << device_idx << "]: " << " launch: " << clock() << "\n";

//			has_data = true;
			m_SolverInfo[device_idx].m_nLastTotalSize = t_size;
			//total_size = 0;

			if(total_size > t_size) {
				// adjust the data array
				for(int k = 0; k < 4 + board_size - level; k++) {
					memcpy(&mask_vector[max_pitch * k], &mask_vector[t_size + max_pitch * k], (total_size - t_size) * sizeof(int));
				}
			}

			total_size -= t_size;
		}
	}

//	if(has_data) {
//		err = clEnqueueReadBuffer(m_Queue, result_buffer, CL_TRUE, 0, threads * sizeof(int) * 2, &results[0], 0, 0, 0);
//		CHECK_ERROR(err);
//
//		if(m_bProfiling) {
//			cl_ulong start, end;
//			clGetEventProfilingInfo(profile_event, CL_PROFILING_COMMAND_START, sizeof(start), &start, 0);
//			clGetEventProfilingInfo(profile_event, CL_PROFILING_COMMAND_END, sizeof(end), &end, 0);
//
//			m_TotalTime += end - start;
//
//			clReleaseEvent(profile_event);
//			profile_event = 0;
//		}
//
//		for(int k = 0; k < last_total_size; k++) {
//			solutions += results[k];
//			unique_solutions += results[k + threads];
//		}
//
//		std::cout << solutions << "\n";
//		std::cout << unique_solutions << "\n";
//	}

//	std::cerr << "final gathering: " << clock() << "\n";

	// gather data from all devices

//	for(int i = 0; i< m_Devices.size(); i++) {
//		err = clFlush(m_SolverInfo[i].m_Queue);
//		CHECK_ERROR(err);
//	}

 //Only one device a time for OpenDwarfs
	while(dev_pending) {
		{
			int idx=0;
			// get data from the device
			CUDA_CHECK(cudaMemcpy(&results[0], m_SolverInfo[idx].m_ResultBuffer, max_pitch * sizeof(int) * 4, cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaDeviceSynchronize());
			//CHECK_ERROR(err);

			//if(m_bProfiling) {
			//	cl_ulong start, end;
			//	clGetEventProfilingInfo(m_SolverInfo[idx].m_Event, CL_PROFILING_COMMAND_START, sizeof(start), &start, 0);
			//	clGetEventProfilingInfo(m_SolverInfo[idx].m_Event, CL_PROFILING_COMMAND_END, sizeof(end), &end, 0);

			//	m_SolverInfo[idx].m_TotalTime += end - start;
			//}

			dev_pending = false;

			for(int k = 0; k < m_SolverInfo[idx].m_nLastTotalSize; k++) {
				if(results[k + max_pitch * 2] != results[k + max_pitch * 3]) {
					throw CLError(3);
				}

				solutions += results[k];
				unique_solutions += results[k + max_pitch];
			}
		}
	}

	for(int i = 0; i < m_SolverInfo.size(); i++) {
		if(m_SolverInfo[i].m_ParamBuffer != 0) { cudaFree(m_SolverInfo[i].m_ParamBuffer); m_SolverInfo[i].m_ParamBuffer = 0; }
		if(m_SolverInfo[i].m_ResultBuffer != 0) { cudaFree(m_SolverInfo[i].m_ResultBuffer); m_SolverInfo[i].m_ResultBuffer = 0; }
		if(m_SolverInfo[i].m_ForbiddenBuffer != 0) { cudaFree(m_SolverInfo[i].m_ForbiddenBuffer); m_SolverInfo[i].m_ForbiddenBuffer = 0; }
		if(m_SolverInfo[i].m_GlobalIndex != 0) { cudaFree(m_SolverInfo[i].m_GlobalIndex); m_SolverInfo[i].m_GlobalIndex = 0; }
	}

	if(unique != 0) {
		*unique = unique_solutions;
	}

//	std::cerr << "complete: " << clock() << "\n";

	return solutions;
}
