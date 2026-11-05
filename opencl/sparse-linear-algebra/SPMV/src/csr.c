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

#include "../../../include/rdtsc.h"
#include "../../../include/common_args.h"
//#include "../../../include/common_util.h"
#include "../inc/common.h"
#include "../inc/sparse_formats.h"
#ifdef __FPGA__
    #include "cl_ext.h"
#endif

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
void float_array_comp(const float* control, const float* experimental, const unsigned int N, const unsigned int exec_num)
{
	unsigned int j;
	float diff,perc;
	for (j = 0; j < N; j++)
	{
		diff = experimental[j] - control[j];
		if(fabsf(diff) > .001)
		{
			perc = fabsf(diff/control[j]) * 100;
			fprintf(stderr,"Possible error on exec #%u, difference of %.3f (%.1f%% error) [control=%.3f, experimental=%.3f] at row %d \n",exec_num,diff,perc,control[j],experimental[j],j);
		}
	}
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
	wg_sizes = malloc(sizeof(size_t)*(*num_wg_sizes));
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
 * stores a valid cl_mem buffer in address *ptr using the given flags and num_bytes.
 */
void csrCreateBuffer(const cl_context* p_context, cl_mem* ptr, const size_t num_bytes, const cl_mem_flags flags, const char* buff_name, int verbosity)
{
	cl_int err;
	char err_msg[128];
	if(verbosity >= 2) printf("Allocating %zu bytes for %s...\n",num_bytes,buff_name);
	*ptr = clCreateBuffer(*p_context, flags,num_bytes, NULL, &err);
	snprintf(err_msg,88,"Failed to allocate device memory for %s!",buff_name);
	CHKERR(err, err_msg);
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
	cl_int err;
	int num_wg,verbosity = 0,do_print=0,do_affirm=0,do_mem_align=0,opt, option_index=0;
	unsigned long density_ppm = 500000;
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
	size_t max_wg_size,kernelLength,items_read;

	//cl_device_id device_id;
	//cl_int dev_type;
	//cl_context context;
	cl_command_queue write_queue;
	//,kernel_queue, => commands
	cl_command_queue read_queue;
	cl_program program;
	cl_kernel kernel;

	FILE *kernel_fp;
	char *kernelSource,**kernel_files=NULL,*kernel_file_mode;

	//Does NOT use ocd_init() because we use TIMER_INIT (the 3rd thing in ocd_init) MANY TIMES in LOOOP! (nz-ocl)
	ocd_parse(&argc, &argv);
	ocd_check_requirements(NULL);
	ocd_initCL();

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
				kernel_files = tmp;
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
				wg_sizes = tmp;
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

	cl_mem csr_ap[num_matrices],csr_aj[num_matrices],csr_ax[num_matrices],x_loc[num_matrices],y_loc[num_matrices];
	cl_event kernel_exec[num_matrices],ap_write[num_matrices],aj_write[num_matrices],ax_write[num_matrices],x_loc_write[num_matrices],y_loc_write[num_matrices],y_read[num_matrices];

	//The other arrays
	float *x_host = NULL, *y_host = NULL, *device_out[num_matrices], *host_out=NULL;
	unsigned int max_row_len=0,max_col_len=0;
	for(ii=0; ii<num_matrices; ii++)
	{
		device_out[ii] = float_new_array(csr[ii].num_rows,"csr.main() - Heap Overflow! Cannot Allocate Space for device_out");
		if(max_row_len < csr[ii].num_rows)
		{
			max_row_len = csr[ii].num_rows;
			y_host = float_array_realloc(y_host,csr[ii].num_rows,"csr.main() - Heap Overflow! Cannot Allocate Space for y_host");
			if(do_affirm)
			{
				host_out = realloc(host_out,sizeof(float)*max_row_len);
				check(host_out != NULL,"csr.main() - Heap Overflow! Cannot Allocate Space for 'host_out'");
			}
		}
		if(max_col_len < csr[ii].num_cols)
		{
			max_col_len = csr[ii].num_cols;
			x_host = float_array_realloc(x_host,csr[ii].num_cols,"csr.main() - Heap Overflow! Cannot Allocate Space for x_host");
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

	if(verbosity) printf("Input Generated.\n");

#ifdef ENABLE_TIMER
	tv = malloc(sizeof(struct timeval));
	check(tv != NULL,"csr.main() - Heap Overflow! Cannot allocate space for tv");
#endif

	/* Retrieve an OpenCL platform */
	//device_id = GetDevice(platform_id, n_device,dev_type);

	//if(verbosity) ocd_print_device_info(device_id);

	/* Create a compute context */
	//context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	//CHKERR(err, "Failed to create a compute context!");

	for(k=0; k<num_matrices; k++)
	{
		if(verbosity >= 2) printf("Creating Data Buffers for Matrix #%d of %d...\n",k+1,num_matrices);
		if (_deviceType == 3){
#if defined(CL_MEM_BANK_1_ALTERA) && defined(CL_MEM_BANK_2_ALTERA)
			csrCreateBuffer(&context,&csr_ap[k],sizeof(int)*(csr[k].num_rows+1),CL_MEM_BANK_1_ALTERA | CL_MEM_READ_ONLY,"csr_ap",verbosity);
			csrCreateBuffer(&context,&x_loc[k],sizeof(float)*csr[k].num_cols,CL_MEM_BANK_1_ALTERA | CL_MEM_READ_ONLY,"x_loc",verbosity);
			csrCreateBuffer(&context,&y_loc[k],sizeof(float)*csr[k].num_rows,CL_MEM_BANK_2_ALTERA | CL_MEM_READ_WRITE,"y_loc",verbosity);
			csrCreateBuffer(&context,&csr_aj[k],sizeof(int)*csr[k].num_nonzeros,CL_MEM_BANK_1_ALTERA | CL_MEM_READ_ONLY,"csr_aj",verbosity);
			csrCreateBuffer(&context,&csr_ax[k],sizeof(float)*csr[k].num_nonzeros,CL_MEM_BANK_2_ALTERA | CL_MEM_READ_ONLY,"csr_ax",verbosity);
#else
			fprintf(stderr, "Must use Altera OpenCL SDK to be able to run with the FPGA option!\n");
			exit(-1);
#endif
		}
		else{
			csrCreateBuffer(&context,&csr_ap[k],sizeof(int)*(csr[k].num_rows+1), CL_MEM_READ_ONLY,"csr_ap",verbosity);
			csrCreateBuffer(&context,&x_loc[k],sizeof(float)*csr[k].num_cols, CL_MEM_READ_ONLY,"x_loc",verbosity);
			csrCreateBuffer(&context,&y_loc[k],sizeof(float)*csr[k].num_rows, CL_MEM_READ_WRITE,"y_loc",verbosity);
			csrCreateBuffer(&context,&csr_aj[k],sizeof(int)*csr[k].num_nonzeros, CL_MEM_READ_ONLY,"csr_aj",verbosity);
			csrCreateBuffer(&context,&csr_ax[k],sizeof(float)*csr[k].num_nonzeros, CL_MEM_READ_ONLY,"csr_ax",verbosity);
		}
	}

	if(!kernel_files) //use default if no kernel files were given on commandline
	{
		num_kernels = 1;
		kernel_files = malloc(sizeof(char*)*num_kernels);
		kernel_files[0] = "spmv_kernel";

	}

	for(iii=0; iii<num_kernels; iii++) //loop through all kernels that need to be tested
	{
		printf("Kernel #%d: '%s'\n\n",iii+1,kernel_files[iii]);
		program = ocdBuildProgramFromFile(context,device_id,kernel_files[iii], NULL);

		if(!wg_sizes) //use default work-group size if none was specified on command line
		{
			/* Get the maximum work group size for executing the kernel on the device */
			kernel = clCreateKernel(program, "csr", &err);
			CHKERR(err, "Failed to create a compute kernel!");
			err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &max_wg_size, NULL);
			if(verbosity) printf("Kernel Max Work Group Size: %d\n",max_wg_size);
			CHKERR(err, "Failed to retrieve kernel work group info!");
			global_size = csr[0].num_rows; //Preconditions: all matrices in input file are same size
			//				all kernels have same max workgroup size
			wg_sizes = default_wg_sizes(&num_wg_sizes,max_wg_size,global_size);
			clReleaseKernel(kernel);
		}

		for(ii=0; ii<num_wg_sizes; ii++) //loop through all wg_sizes that need to be tested
		{
			num_wg = global_size / wg_sizes[ii];
			printf("Executing with WG Size #%d of %d: %d...\n",ii+1,num_wg_sizes,wg_sizes[ii]);

			for(i=0; i<num_execs; i++) //repeat Host-Device transfer, kernel execution, and device-host transfer num_execs times
			{						//to gather multiple samples of data
				if(verbosity) printf("Beginning execution #%d of %d\n",i+1,num_execs);

				/* Create command queues, one for each stage in the write-execute-read pipeline */
				write_queue = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);
				CHKERR(err, "Failed to create a command queue!");
				commands = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);//re-assigned every time (ocd_initCL has already done this)
				CHKERR(err, "Failed to create a command queue!");
				read_queue = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);
				CHKERR(err, "Failed to create a command queue!");

				/* Get the maximum work group size for executing the kernel on the device */
				kernel = clCreateKernel(program, "csr", &err);
				CHKERR(err, "Failed to create a compute kernel!");

#ifdef ENABLE_TIMER
				TIMER_INIT
#endif

					for(k=0; k<num_matrices; k++)
					{
						if(verbosity >= 2) printf("Enqueuing Matrix #%d of %d into pipeline...\n",k+1,num_matrices);

						/* Write our data set into the input array in device memory */
						err = clEnqueueWriteBuffer(write_queue, csr_ap[k], CL_FALSE, 0, sizeof(unsigned int)*csr[k].num_rows+4, csr[k].Ap, 0, NULL, &ap_write[k]);
						CHKERR(err, "Failed to write to source array!");

						err = clEnqueueWriteBuffer(write_queue, csr_aj[k], CL_FALSE, 0, sizeof(unsigned int)*csr[k].num_nonzeros, csr[k].Aj, 0, NULL, &aj_write[k]);
						CHKERR(err, "Failed to write to source array!");

						err = clEnqueueWriteBuffer(write_queue, csr_ax[k], CL_FALSE, 0, sizeof(float)*csr[k].num_nonzeros, csr[k].Ax, 0, NULL, &ax_write[k]);
						CHKERR(err, "Failed to write to source array!");

						err = clEnqueueWriteBuffer(write_queue, x_loc[k], CL_FALSE, 0, sizeof(float)*csr[k].num_cols, x_host, 0, NULL, &x_loc_write[k]);
						CHKERR(err, "Failed to write to source array!");

						err = clEnqueueWriteBuffer(write_queue, y_loc[k], CL_FALSE, 0, sizeof(float)*csr[k].num_rows, y_host, 0, NULL, &y_loc_write[k]);
						CHKERR(err, "Failed to write to source array!");

						/* Set the arguments to our compute kernel */
						global_size = csr[k].num_rows;
						err = clSetKernelArg(kernel, 0, sizeof(int), &global_size);
						err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &csr_ap[k]);
						err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &csr_aj[k]);
						err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &csr_ax[k]);
						err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &x_loc[k]);
						err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &y_loc[k]);
						CHKERR(err, "Failed to set kernel arguments!");

						/* Enqueue Kernel */
						err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global_size, &wg_sizes[ii], 1, &y_loc_write[k], &kernel_exec[k]);
						CHKERR(err, "Failed to execute kernel!");

						/* Read back the results from the device to verify the output */
						err = clEnqueueReadBuffer(read_queue, y_loc[k], CL_FALSE, 0, sizeof(float)*csr[k].num_rows, device_out[k], 1, &kernel_exec[k], &y_read[k]);
						CHKERR(err, "Failed to read output array!");
					}
				clFinish(write_queue);
				clFinish(commands);
				clFinish(read_queue);

#ifdef ENABLE_TIMER
				TIMER_STOP
#endif

					for(k=0; k<num_matrices; k++)
					{
							START_TIMER(ap_write[k], OCD_TIMER_H2D, "CSR Data Copy", ocdTempTimer)
							END_TIMER(ocdTempTimer)

							START_TIMER(aj_write[k], OCD_TIMER_H2D, "CSR Data Copy", ocdTempTimer)
							END_TIMER(ocdTempTimer)

							START_TIMER(ax_write[k], OCD_TIMER_H2D, "CSR Data Copy", ocdTempTimer)
							END_TIMER(ocdTempTimer)

							START_TIMER(x_loc_write[k], OCD_TIMER_H2D, "CSR Data Copy", ocdTempTimer)
							END_TIMER(ocdTempTimer)

							START_TIMER(y_loc_write[k], OCD_TIMER_H2D, "CSR Data Copy", ocdTempTimer)
							END_TIMER(ocdTempTimer)

							START_TIMER(kernel_exec[k], OCD_TIMER_KERNEL, "CSR Kernel", ocdTempTimer)
							END_TIMER(ocdTempTimer)

							START_TIMER(y_read[k], OCD_TIMER_D2H, "CSR Data Copy", ocdTempTimer)
							END_TIMER(ocdTempTimer)

							if(do_print)
							{
								printf("\nMatrix #%d of %d:\n",k+1,num_matrices);
								for(j = 0; j < csr[k].num_rows; j++)
									printf("\trow: %d	output: %6.2f \n", j, device_out[k][j]);
							}
					}

				clReleaseCommandQueue(write_queue);
				CHKERR(err,"Failed to release write_queue!");
				clReleaseCommandQueue(commands);
				CHKERR(err,"Failed to release kernel_queue!");
				clReleaseCommandQueue(read_queue);
				CHKERR(err,"Failed to release read_queue!");
				clReleaseKernel(kernel);
				CHKERR(err,"Failed to release kernel!");

#ifdef ENABLE_TIMER
				TIMER_PRINT
#endif

					if(do_affirm)
					{
						if(verbosity) printf("Validating results with serial C code on CPU...\n");
						for(k=0; k<num_matrices; k++)
						{
							spmv_csr_cpu(&csr[k],x_host,y_host,host_out);
							float_array_comp(host_out,device_out[k],csr[k].num_rows,i+1);
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
			err = clReleaseMemObject(csr_ap[k]);
			CHKERR(err,"Failed to release csr_ap!");
			err = clReleaseMemObject(csr_aj[k]);
			CHKERR(err,"Failed to release csr_aj!");
			err = clReleaseMemObject(csr_ax[k]);
			CHKERR(err,"Failed to release csr_ax!");
			err = clReleaseMemObject(x_loc[k]);
			CHKERR(err,"Failed to release x_loc!");
			err = clReleaseMemObject(y_loc[k]);
			CHKERR(err,"Failed to release y_loc!");
			//		err = clReleaseEvent(aj_write[k]);	//releasing of any of these events is throwing an error with the altera sdk.
			//		if(verbosity) printf("k: %d\terr: %d\n",k,err); //Perhaps because the command-queue was already released?
			//		CHKERR(err,"Failed to release aj_write!");
			//		err = clReleaseEvent(ap_write[k]);
			//		if(verbosity) printf("k: %d\terr: %d\n",k,err);
			//		CHKERR(err,"Failed to release ap_write!");
			//		err = clReleaseEvent(ax_write[k]);
			//		CHKERR(err,"Failed to release ax_write!");
			//		err = clReleaseEvent(x_loc_write[k]);
			//		CHKERR(err,"Failed to release x_loc_write!");
			//		err = clReleaseEvent(y_loc_write[k]);
			//		if(verbosity) printf("k: %d\terr: %d\n",k,err);
			//		CHKERR(err,"Failed to release y_loc_write!");
			//		err = clReleaseEvent(kernel_exec[k]);
			//		CHKERR(err,"Failed to release kernel_exec!");
			free(device_out[k]);
		}

	clReleaseContext(context);
	CHKERR(err,"Failed to release context!");
	if(verbosity) printf("Released context\n");
	free(kernel_files);
	free(wg_sizes);
	free(tv);
	free(x_host);
	free(y_host);
	if(do_affirm) free(host_out);
	free_csr(csr,num_matrices);
	return 0;
}

