/* Main File for SPMV (aka CSR) application which represents Sparse-Linear-Algebra dwarf
 *
 * This application multiplies a sparse matrix, stored in CSR format, by a vector, adds it with another vector, and returns the resulting vector.
 *
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include "common_args.h"
#include "dwarf_cuda.h"
#include "spmv_kernel.cuh"
#include "dwarf_timer.h"
//#include "../../../include/common_util.h"
#include "../inc/common.h"
#include "../inc/sparse_formats.h"

#define START_GTOD_TIMER { \
	gettimeofday(tv,NULL); \
	start_time = 1000 * (tv->tv_sec*1000000L + tv->tv_usec); }

#define END_GTOD_TIMER { \
	gettimeofday(tv,NULL); \
	end_time = 1000 * (tv->tv_sec*1000000L + tv->tv_usec); }

static struct option long_options[] = {
	/* name, has_arg, flag, val */
	{"cpu", 0, NULL, 'c'},
	{"device", 1, NULL, 'd'},
	{"verbose", 0, NULL, 'v'},
	{"input_file",1,NULL,'i'},
	{"print",0,NULL,'p'},
	{"affirm",0,NULL,'a'},
	{"repeat",1,NULL,'r'},
	{"kernel_file",1,NULL,'k'},
	{"wg_size",1,NULL,'w'},
	{"enqueue",1,NULL,'e'},
	{0,0,0,0}
};

//int platform_id=PLATFORM_ID, n_device=DEVICE_ID;

/**
 * Compares N float values and prints error msg if any corresponding entries differ by greater than .001
 */
/* [FIX-OCL-04] now returns the number of mismatching rows (it returned void, so a failed -a check still exited with status 0) */
unsigned int float_array_comp(const float* control, const float* experimental, const unsigned int N, const unsigned int exec_num)
{
	unsigned int j, bad = 0;
	float diff,perc;
	for (j = 0; j < N; j++)
	{
		diff = experimental[j] - control[j];
		if(fabsf(diff) > .001)
		{
			perc = fabsf(diff/control[j]) * 100;
			fprintf(stderr,"Possible error on exec #%u, difference of %.3f (%.1f%% error) [control=%.3f, experimental=%.3f] at row %d \n",exec_num,diff,perc,control[j],experimental[j],j);
			bad++;
		}
	}
	return bad;
}

/**
 * Sparse Matrix-Vector Multiply
 *
 * Multiplies csr matrix by vector x, adds vector y, and stores output in vector out
 */
void spmv_csr_cpu(const csr_matrix* csr,const float* x,const float* y,float* out)
{
	unsigned int row,row_start,row_end,jj;
	float sum = 0;
	for(row=0; row < csr->num_rows; row++)
	{
		sum = y[row];
		row_start = csr->Ap[row];
		row_end   = csr->Ap[row+1];

		for (jj = row_start; jj < row_end; jj++){
			sum += csr->Ax[jj] * x[csr->Aj[jj]];
		}
		out[row] = sum;
	}
}

/*
 * Returns an array of work group sizes with only 1 element. The value is the largest possible
 * work-group size (i.e., fewest number of work-groups possible will be used), whether thats
 * limited by the device or the global size of the application
 */
size_t* default_wg_sizes(unsigned int* num_wg_sizes,const size_t max_wg_size,const size_t global_size)
{
	unsigned int num_wg;
	size_t* wg_sizes;
	(*num_wg_sizes)=1;
	wg_sizes = (size_t*)malloc(sizeof(size_t)*(*num_wg_sizes));
	check(wg_sizes != NULL,"csr.main() - Heap Overflow! Cannot allocate space for wg_sizes");
	wg_sizes[0] = max_wg_size;
	num_wg = global_size / wg_sizes[0];
	if(global_size % wg_sizes[0] != 0) //if wg_size is not a factor of global_size
	{							//use min num_wg such that wg_size < global_size
		num_wg++;
		wg_sizes[0] = global_size / (num_wg);
	}
	return wg_sizes;
}

/*
 * Main Method
 *
 * Reads input file, builds kernel and then executes it a number of times.
 * The number of executions is equal to num_execs*num_wg_sizes*num_kernels. These
 * values can be controlled from the command line with the -r,-w, and -k options,
 * respectively, and are implemented via nested for-loops. As the input is scaled
 * to larger size, reading in the input file becomes very time consuming and this
 * allows one to test multiple variables without repeating that process.
 */
int main(int argc, char** argv)
{
	int num_wg,verbosity = 0,do_print=0,do_affirm=0,do_mem_align=0,opt, option_index=0;
	unsigned long density_ppm = 500000;
	int rep; unsigned int spmv_bad = 0; /* common timer repetition index / -a mismatch counter */
	unsigned int N = 512,num_execs=1,num_matrices,i,ii,iii,j,k,num_wg_sizes=0,num_kernels=0;
	unsigned long start_time, end_time;
	struct timeval *tv;
	char* file_path = NULL,*optptr;
	void* tmp;

	const char* usage = "Usage: %s -i <file_path> [-v] [-c] [-p] [-a] [-r <num_execs>] [-k <kernel_file-1>][-k <kernel_file-2>]...[-k <kernel_file-n>] [-w <wg_size-1>][-w <wg_size-2>]...[-w <wg_size-m>]\n\n \
			     -i: Read CSR Matrix from file <file_path>\n \
			     -k: Test SPMV 'n' times, once with each kernel_file-'1..n' - Default is 1 kernel named './spmv_csr_kernel.xxx' where xxx is 'aocx' if USE_AFPGA is defined, 'cl' otherwise.\n \
			     -v: Increase verbosity level by 1 - Default is 0 - Max is 2 \n \
			     -p: Print matrices to stdout in standard (2-D Array) format - Warning: lots of output\n \
			     -a: Affirm results with serial C code on CPU\n \
			     -r: Execute program with same data exactly <num_execs> times to increase sample size - Default is 1\n \
			     -w: Loop through each kernel execution 'm' times, once with each wg_size-'1..m' - Default is 1 iteration with wg_size set to the maximum possible (limited either by the device or the size of the input)\n\n";

	size_t global_size;
	size_t* wg_sizes = NULL;

	char **kernel_files=NULL;

	//Does NOT use ocd_init() because we use TIMER_INIT (the 3rd thing in ocd_init) MANY TIMES in LOOOP! (nz-ocl)
	dwarf_cuda_init(&argc, &argv);

	while ((opt = getopt_long(argc, argv, "::vmw:k:i:par:::", long_options, &option_index)) != -1 )
	{
		switch(opt)
		{
			case 'v':
				verbosity++;
				break;
				//case 'c':
				//	printf("using cpu\n");
				//	_deviceType = CL_DEVICE_TYPE_CPU;
				//	break;
			case 'i':
				if(optarg != NULL)
					file_path = optarg;
				else
					file_path = argv[optind];
				printf("Reading Input from '%s'\n",file_path);
				break;
			case 'p':
				do_print = 1;
				break;
			case 'a':
				do_affirm = 1;
				break;
			case 'r':
				if(optarg != NULL)
					num_execs = atoi(optarg);
				else
					num_execs = atoi(argv[optind]);
				printf("Executing %d times\n",num_execs);
				break;
			case 'k':
				if(optarg != NULL)
					optptr = optarg;
				else
					optptr = argv[optind];
				num_kernels++;
				tmp = realloc(kernel_files,sizeof(char*)*num_kernels);
				check(tmp != NULL,"csr.main() - Heap Overflow! Cannot allocate space for kernel_files");
				kernel_files = (char**)tmp;
				kernel_files[num_kernels-1] = optptr;
				printf("Testing with Kernel File: '%s'\n",kernel_files[num_kernels-1]);
				break;
			case 'w':
				if(optarg != NULL)
					optptr = optarg;
				else
					optptr = argv[optind];
				num_wg_sizes++;
				tmp = realloc(wg_sizes,sizeof(size_t)*num_wg_sizes);
				check(tmp != NULL,"csr.main() - Heap Overflow! Cannot allocate space for wg_sizes");
				wg_sizes = (size_t*)tmp;
				wg_sizes[num_wg_sizes-1] = atoi(optptr);
				break;
			default:
				fprintf(stderr, usage,argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	if(!file_path)
	{
		fprintf(stderr,"-i Option must be supplied\n\n");
		fprintf(stderr, usage,argv[0]);
		exit(EXIT_FAILURE);
	}

	csr_matrix* csr = read_csr(&num_matrices,file_path);

	if(do_print) print_csr_arr_std(csr,num_matrices,stdout);
	else if(verbosity) {printf("Number of input matrices: %d\nMatrix 0 Metadata:\n",num_matrices); print_csr_metadata(&csr[0],stdout);}

	unsigned int *csr_ap[num_matrices], *csr_aj[num_matrices];
	float *csr_ax[num_matrices], *x_loc[num_matrices], *y_loc[num_matrices];

	//The other arrays
	float *x_host = NULL, *y_host = NULL, *device_out[num_matrices], *host_out=NULL;
	unsigned int max_row_len=0,max_col_len=0;
	for(ii=0; ii<num_matrices; ii++)
	{
		device_out[ii] = (float*)float_new_array(csr[ii].num_rows,"csr.main() - Heap Overflow! Cannot Allocate Space for device_out");
		if(max_row_len < csr[ii].num_rows)
		{
			max_row_len = csr[ii].num_rows;
			y_host = (float*)float_array_realloc(y_host,csr[ii].num_rows,"csr.main() - Heap Overflow! Cannot Allocate Space for y_host");
			if(do_affirm)
			{
				host_out = (float*)realloc(host_out,sizeof(float)*max_row_len);
				check(host_out != NULL,"csr.main() - Heap Overflow! Cannot Allocate Space for 'host_out'");
			}
		}
		if(max_col_len < csr[ii].num_cols)
		{
			max_col_len = csr[ii].num_cols;
			x_host = (float*)float_array_realloc(x_host,csr[ii].num_cols,"csr.main() - Heap Overflow! Cannot Allocate Space for x_host");
		}
	}

	for(ii = 0; ii < max_col_len; ii++)
	{
		x_host[ii] = rand() / (RAND_MAX + 1.0);
		if(do_print) printf("x[%d] = %6.2f\n",ii,x_host[ii]);
	}
	for(ii = 0; ii < max_row_len; ii++)
	{
		y_host[ii] = rand() / (RAND_MAX + 2.0);
		if(do_print) printf("y[%d] = %6.2f\n",ii,y_host[ii]);
	}

	dwarf_init("spmv", "sparse-linear-algebra", "cuda");
	dwarf_set_size(csr[0].num_rows, csr[0].num_nonzeros);
	if(verbosity) printf("Input Generated.\n");

#ifdef ENABLE_TIMER
	tv = (struct timeval*)malloc(sizeof(struct timeval));
	check(tv != NULL,"csr.main() - Heap Overflow! Cannot allocate space for tv");
#endif

	for(k=0; k<num_matrices; k++)
	{
		if(verbosity >= 2) printf("Creating Data Buffers for Matrix #%d of %d...\n",k+1,num_matrices);
		CUDA_CHECK(cudaMalloc((void**)&csr_ap[k], sizeof(int)*(csr[k].num_rows+1)));
		CUDA_CHECK(cudaMalloc((void**)&x_loc[k], sizeof(float)*csr[k].num_cols));
		CUDA_CHECK(cudaMalloc((void**)&y_loc[k], sizeof(float)*csr[k].num_rows));
		CUDA_CHECK(cudaMalloc((void**)&csr_aj[k], sizeof(int)*csr[k].num_nonzeros));
		CUDA_CHECK(cudaMalloc((void**)&csr_ax[k], sizeof(float)*csr[k].num_nonzeros));
	}

	if(!kernel_files) //use default if no kernel files were given on commandline
	{
		num_kernels = 1;
		kernel_files = (char**)malloc(sizeof(char*)*num_kernels);
		kernel_files[0] = (char*)"spmv_kernel";

	}

	for(iii=0; iii<num_kernels; iii++) //loop through all kernels that need to be tested
	{
		printf("Kernel #%d: '%s'\n\n",iii+1,kernel_files[iii]);

		/* OpenCL queries the device's max work-group size here; it only labels the loop below, the launch itself is pinned */
		if(!wg_sizes)
		{
			global_size = csr[0].num_rows;
			wg_sizes = default_wg_sizes(&num_wg_sizes,1024,global_size);
		}

		for(ii=0; ii<num_wg_sizes; ii++) //loop through all wg_sizes that need to be tested
		{
			num_wg = global_size / wg_sizes[ii];
			printf("Executing with WG Size #%d of %d: %d...\n",ii+1,num_wg_sizes,wg_sizes[ii]);

			if(num_execs > 1) dwarf_set_repeats((int)num_execs);
			for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++) //repeat Host-Device transfer, kernel execution, and device-host transfer (warm-ups first)
			{						//to gather multiple samples of data
				i = rep < 0 ? 0 : (unsigned int)rep;
				if(verbosity) printf("Beginning execution #%d of %d\n",i+1,num_execs);

#ifdef ENABLE_TIMER
				TIMER_INIT
#endif
				dwarf_region_begin();

					for(k=0; k<num_matrices; k++)
					{
						if(verbosity >= 2) printf("Enqueuing Matrix #%d of %d into pipeline...\n",k+1,num_matrices);

						/* Write our data set into the input array in device memory (non-blocking, as in OpenCL) */
						CUDA_CHECK(cudaMemcpyAsync(csr_ap[k], csr[k].Ap, sizeof(unsigned int)*csr[k].num_rows+4, cudaMemcpyHostToDevice));
						CUDA_CHECK(cudaMemcpyAsync(csr_aj[k], csr[k].Aj, sizeof(unsigned int)*csr[k].num_nonzeros, cudaMemcpyHostToDevice));
						CUDA_CHECK(cudaMemcpyAsync(csr_ax[k], csr[k].Ax, sizeof(float)*csr[k].num_nonzeros, cudaMemcpyHostToDevice));
						CUDA_CHECK(cudaMemcpyAsync(x_loc[k], x_host, sizeof(float)*csr[k].num_cols, cudaMemcpyHostToDevice));
						CUDA_CHECK(cudaMemcpyAsync(y_loc[k], y_host, sizeof(float)*csr[k].num_rows, cudaMemcpyHostToDevice));

						global_size = csr[k].num_rows;
						/* one thread per row: DWARF_LS threads per block, ceil(rows/DWARF_LS) blocks (kernel guards row < num_rows) */
						csr_kernel<<<DWARF_NBLOCKS(global_size), DWARF_LS>>>((unsigned int)global_size, csr_ap[k], csr_aj[k], csr_ax[k], x_loc[k], y_loc[k]);
						CUDA_CHECK_LAUNCH();

						/* Read back the results from the device to verify the output */
						CUDA_CHECK(cudaMemcpyAsync(device_out[k], y_loc[k], sizeof(float)*csr[k].num_rows, cudaMemcpyDeviceToHost));
					}
				CUDA_CHECK(cudaDeviceSynchronize());
				dwarf_region_end(rep);

#ifdef ENABLE_TIMER
				TIMER_STOP
#endif

					for(k=0; k<num_matrices; k++)
					{
							if(do_print)
							{
								printf("\nMatrix #%d of %d:\n",k+1,num_matrices);
								for(j = 0; j < csr[k].num_rows; j++)
									printf("\trow: %d	output: %6.2f \n", j, device_out[k][j]);
							}
					}

#ifdef ENABLE_TIMER
				TIMER_PRINT
#endif

					if(do_affirm)
					{
						if(verbosity) printf("Validating results with serial C code on CPU...\n");
						for(k=0; k<num_matrices; k++)
						{
							spmv_csr_cpu(&csr[k],x_host,y_host,host_out);
							spmv_bad += float_array_comp(host_out,device_out[k],csr[k].num_rows,i+1);
						}
					}
			}
		}
	}
#ifdef ENABLE_TIMER
	TIMER_DEST
#endif

		/* Shutdown and cleanup */
		for(k=0; k<num_matrices; k++)
		{
			CUDA_CHECK(cudaFree(csr_ap[k])); CUDA_CHECK(cudaFree(csr_aj[k])); CUDA_CHECK(cudaFree(csr_ax[k]));
			CUDA_CHECK(cudaFree(x_loc[k])); CUDA_CHECK(cudaFree(y_loc[k]));
			free(device_out[k]);
		}

	free(kernel_files);
	free(wg_sizes);
	free(tv);
	free(x_host);
	free(y_host);
	if(do_affirm) free(host_out);
	free_csr(csr,num_matrices);
	if(do_affirm) printf("%s: %u rows differ by more than 0.001\n", spmv_bad ? "FAILURE" : "SUCCESS", spmv_bad);
	return dwarf_report(spmv_bad == 0);
}

