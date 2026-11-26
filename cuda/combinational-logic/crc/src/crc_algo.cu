/*Main file for CRC application which represents the combinational-logic dwarf.
 *
 * This application computes a 32-bit ethernet CRC on a number of input pages using
 * a "Slice-By-8" algorithm published by Intel.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "common_args.h"
#include "dwarf_cuda.h"
#include "dwarf_timer.h"
#include "../inc/crc_formats.h"
#include "../inc/eth_crc32_lut.h"
#include "crc_kernel.cuh"

#define DATA_SIZE 100000000

//The CRC algorithms used in this dwarf were copied and/or adapted from
//versions posted by Stephan Brumme on the website:
//http://create.stephan-brumme.com/crc32/
const uint32_t Polynomial = 0xEDB88320;

unsigned char verbosity=0;
//int platform_id=PLATFORM_ID, n_device=DEVICE_ID;

int64_t dataTransferTime = 0;
int64_t kernelExecutionTime = 0;

unsigned int* num_parallel_crcs;
unsigned int page_size=DATA_SIZE,num_wg_sizes=0,num_words,num_blocks,num_pages_last_block,num_block_sizes=0;
size_t* wg_sizes=NULL;

void printTimeDiff(struct timeval start, struct timeval end)
{
	printf("%ld microseconds\n", ((end.tv_sec * 1000000 + end.tv_usec)
				- (start.tv_sec * 1000000 + start.tv_usec)));
}

int64_t computeTimeDiff(struct timeval start, struct timeval end)
{
	int64_t diff = (end.tv_sec * 1000000 + end.tv_usec)
		- (start.tv_sec * 1000000 + start.tv_usec);
	return diff;
}

// /////Bitwise version of CRC///////////////////////
////////altered from the fastest version of crc32_bitwise() by Stephan Brumme////////
// Copyright (c) 2013 Stephan Brumme. All rights reserved.
// see http://create.stephan-brumme.com/disclaimer.html
//
unsigned int serialCRC(unsigned int* h_num, size_t size)
{
	unsigned int j;
	uint32_t crc = ~0x00;
	unsigned char* current = (unsigned char*) h_num;
	while (size--)
	{
		crc ^= *current++;
		for (j = 0; j < 8; j++)
			crc = (crc >> 1) ^ (-1*((int)(crc & 0x01)) & Polynomial);
	}
	return ~crc;
}

// /////Slice-by-8 version of CRC///////////////////////
////////altered from the version posted online by Stephan Brumme////////
// Copyright (c) 2013 Stephan Brumme. All rights reserved.
// see http://create.stephan-brumme.com/disclaimer.html
//
uint32_t crc32_8bytes(const void* data, size_t length)
{
	uint32_t* current = (uint32_t*) data;
	uint32_t crc = 0xFFFFFFFF;

	while (length >= 8) // process eight bytes at once
	{
		uint32_t one = *current++ ^ crc;
		uint32_t two = *current++;
		crc = crc32Lookup[7][ one      & 0xFF] ^
			crc32Lookup[6][(one>> 8) & 0xFF] ^
			crc32Lookup[5][(one>>16) & 0xFF] ^
			crc32Lookup[4][ one>>24        ] ^
			crc32Lookup[3][ two      & 0xFF] ^
			crc32Lookup[2][(two>> 8) & 0xFF] ^
			crc32Lookup[1][(two>>16) & 0xFF] ^
			crc32Lookup[0][ two>>24        ];
		length -= 8;
	}

	unsigned char* currentChar = (unsigned char*) current;
	while (length--) { // remaining 1 to 7 bytes
		crc = (crc >> 8) ^ crc32Lookup[0][(crc & 0xFF) ^ *currentChar++];
	}
	return ~crc;
}

void enqueueCRCDevice(unsigned int* h_num, unsigned int* h_answer, size_t global_size, size_t local_size, unsigned int* d_input, unsigned int* d_output)
{

	// Write our data set into the input array in device memory (non-blocking, as in OpenCL)
	CUDA_CHECK(cudaMemcpyAsync(d_input, h_num, sizeof(char)*page_size*global_size, cudaMemcpyHostToDevice));

	/* one thread per page: DWARF_LS threads per block, ceil(pages/DWARF_LS) blocks; the kernel guards gid >= num_items */
	unsigned int num_items = (unsigned int)global_size;
	(void)local_size;
	if(verbosity >=2) printf("enqueueCRCDevice(): items=%zd global=%zd local=%d\n",global_size,DWARF_ROUNDUP(global_size),DWARF_LS);
	crc32_slice8<<<DWARF_NBLOCKS(global_size), DWARF_LS>>>(d_input, page_size, num_words, d_output, num_items);
	CUDA_CHECK_LAUNCH();

	// Read back the results from the device to verify the output
	CUDA_CHECK(cudaMemcpyAsync(h_answer, d_output, sizeof(int)*global_size, cudaMemcpyDeviceToHost));
}

void setup_device(const char* kernel_file)
{
	(void)kernel_file; /* kernels are compiled into the binary */
	if(!wg_sizes)
	{
		num_wg_sizes = 1;
		wg_sizes = (size_t*)malloc(sizeof(size_t)*num_wg_sizes);
		wg_sizes[0] = 1;
	}
}

void usage()
{
	printf("crc -i <input_file> [hvp] [-r <num_execs>] [-w <wg_size-1>][-w <wg_size-2>]...[-w <wg_size-m>] [-k <kernel_file-1>][-k <kernel_file-2>]...[-k <kernel_file-n>]\n");
	printf("Common arguments:\n");
	printf("\t[-p N] [-d N] -- : device selection (-d)\n");
	printf("Program-specific arguments:\n");
	printf("\t-h | 'Print this help message'\n");
	printf("\t-v | 'Increase verbosity level by 1 - Default is 0 - Max is 2'\n");
	printf("\t-i | 'Input file name' [string]\n");
	printf("\t-a | 'Verify results on CPU'\n");
	printf("\t-p | 'Set the number of pages to CRC in parallel (i.e., the global size of each kernel) - Default is 16\n");
	printf("\t-r | 'Execute program with same data exactly <num_execs> times to increase sample size - Default is 1\n");
	printf("\t-w | 'Loop through each kernel execution 'm' times, once with each wg_size-'1..m' - Default is 1 iteration with wg_size set to the maximum possible (limited either by the device or the size of the input)\n");
	printf("\t-k | 'Test CRC 'n' times, once with each kernel_file-'1..n' - Default is 1 kernel named './crc_kernel.xxx' where xxx is 'aocx' if USE_AFPGA is defined, 'cl' otherwise.\n");

	printf("\nNOTE: Seperate common arguments and program specific arguments with the '--' delimeter\n");
	exit(0);
}

int main(int argc, char** argv)
{
	size_t maxSize=DATA_SIZE,global_size,local_size;
	FILE* fp=NULL;
	void* tmp;
	unsigned int *h_num,cpu_remainder;
	unsigned int run_serial=0,seed=time(NULL),h,ii,i,j,k,l,m,num_pages=1,num_execs=1,num_kernels=0;
	char* file=NULL,*optptr;
	char** kernel_files=NULL;
	int c;
	int rep; unsigned int crc_bad = 0; /* common timer: repetition index; serial-check mismatch counter */
	struct timeval start,end;

	dwarf_cuda_init(&argc, &argv);

	while((c = getopt (argc, argv, "avn:s:i:p:w:k:hr:")) != -1)
	{
		switch(c)
		{
			case 'h':
				usage();
				exit(0);
				break;
			case 'v':
				verbosity++;
				break;
			case 'a':
				run_serial = 1;
				break;
			case 'i':
				if(optarg != NULL)
					file = optarg;
				else
					file = argv[optind];
				printf("Reading Input from '%s'\n",file);
				break;
			case 'r':
				if(optarg != NULL)
					num_execs = atoi(optarg);
				else
					num_execs = atoi(argv[optind]);
				printf("Executing %d times\n",num_execs);
				break;
			case 'p':
				if(optarg != NULL)
					optptr = optarg;
				else
					optptr = argv[optind];
				num_block_sizes++;
				tmp = realloc(num_parallel_crcs,sizeof(size_t)*num_block_sizes);
				check(tmp != NULL,"csr.main() - Heap Overflow! Cannot allocate space for num_parallel_crcs");
				num_parallel_crcs = (unsigned int*)tmp;
				num_parallel_crcs[num_block_sizes-1] = atoi(optptr);
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
			default:
				fprintf(stderr, "Invalid argument: '%s'\n\n",optarg);
				usage();
		}	
	}

	check(file != NULL,"-i option must be supplied!");
	h_num = read_crc(&num_pages,&page_size,file);

	if(!num_block_sizes)
	{
		num_block_sizes=1;
		num_parallel_crcs = (unsigned int*)malloc(sizeof(int)*num_block_sizes);
		num_parallel_crcs[0] = 16;
	}

	num_words = page_size / 4;
	if(verbosity) printf("num_words = %u\n",num_words);

	//ocd_options opts = ocd_get_options();
	//platform_id = opts.platform_id;
	//n_device = opts.device_id;

	//#ifdef USEGPU
	//	 dev_type = CL_DEVICE_TYPE_GPU;
	//#elif defined(USE_AFPGA)
	//	 dev_type = CL_DEVICE_TYPE_ACCELERATOR;
	//#else
	//	dev_type = CL_DEVICE_TYPE_CPU;
	//#endif

	//if(verbosity) printf("Getting Device\n");
	//device_id = GetDevice(platform_id, n_device,dev_type);

	//context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	//CHKERR(err, "Failed to create a compute context!");

	if(!kernel_files) //use default if no kernel files were given on commandline
	{
		num_kernels = 1;
		kernel_files = (char**)malloc(sizeof(char*)*num_kernels);
		kernel_files[0] = (char*)"crc_kernel";
	}

	dwarf_init("crc", "combinational-logic", "cuda");
	dwarf_set_size(num_pages, page_size);
	if(num_execs > 1) dwarf_set_repeats((int)num_execs);
	for(h=0; h<num_block_sizes; h++)
	{
		if(verbosity) printf("Executing with block size #%u of %u: %u\n",h+1,num_block_sizes,num_parallel_crcs[h]);

		num_blocks = num_pages/num_parallel_crcs[h];
		if(num_pages % num_parallel_crcs[h] != 0)
		{
			num_blocks++;
			num_pages_last_block = num_pages % num_parallel_crcs[h];
		}
		else
		{
			num_pages_last_block = num_parallel_crcs[h];
		}

		if(verbosity) printf("Num Pages: %u - Num Parallel CRCs: %u - Num blocks = %u\n",num_pages,num_parallel_crcs[h],num_blocks);
		unsigned int *dev_input[num_blocks], *dev_output[num_blocks];
		unsigned int* ocl_remainders;
		ocl_remainders = (unsigned int*)int_new_array(num_pages,"crc_algo.main() - Heap Overflow! Cannot allocate space for ocl_remainders");

		for(i=0; i<num_blocks; i++)
		{
			CUDA_CHECK(cudaMalloc((void**)&dev_input[i], sizeof(char)*page_size*num_parallel_crcs[h]));
			CUDA_CHECK(cudaMalloc((void**)&dev_output[i], sizeof(int)*num_parallel_crcs[h]));
		}

		for(l=0; l<num_kernels; l++)
		{
			if(verbosity) printf("Executing with kernel #%u of %u: %s\n",l+1,num_kernels,kernel_files[l]);
			setup_device(kernel_files[l]);

			for(k=0; k<num_wg_sizes; k++)
			{
				if(verbosity) printf("Executing with Workgroup size #%u of %u: %zu\n",k+1,num_wg_sizes,wg_sizes[k]);

				/* common timer: repetitions (warm-ups first); the timed region is the H2D/kernel/D2H pipeline below */
				for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
				{
					ii = rep < 0 ? 0 : (unsigned int)rep;
					if(verbosity) printf("Beginning execution #%u of %u...\n",ii+1,num_execs);

#ifdef ENABLE_TIMER
					TIMER_INIT
#endif
					dwarf_region_begin();
						for(i=0; i<num_blocks; i++)
						{
							if(verbosity >= 2) printf("\tEnqueuing commmands for block #%d of %d...\n",i+1,num_blocks);
							if(i == num_blocks -1) //last iteration
							{
								global_size = num_pages_last_block;
								local_size = wg_sizes[k];
								if((global_size % local_size) != 0)
								{
									local_size = 1;
									while((global_size % local_size) == 0) local_size = local_size << 1;
									local_size = local_size >> 1;
								}
							}
							else
							{
								global_size = num_parallel_crcs[h];
								local_size = wg_sizes[k];
							}
							if(verbosity >= 2) printf("\tmain(): global_size=%zd - local_size=%zd\n",global_size,local_size);
							enqueueCRCDevice(&h_num[i*num_parallel_crcs[h]*num_words],&ocl_remainders[i*num_parallel_crcs[h]],global_size,local_size,dev_input[i],dev_output[i]);
						}
					CUDA_CHECK(cudaDeviceSynchronize());
					dwarf_region_end(rep);

#ifdef ENABLE_TIMER
					TIMER_STOP
#endif

						for(i=0; i<num_blocks; i++)
						{
							if(verbosity >= 2) printf("Parallel Computation: '%X'\n", ocl_remainders[i]);

						}

#ifdef ENABLE_TIMER
					TIMER_PRINT
#endif

						if(run_serial) // verify that we have the correct answer with regular C
						{
							printf("Validating results with serial CRC...\n");
							//						gettimeofday(&start,NULL);
							//						for(i=0; i<num_pages; i++)
							//						{
							//							cpu_remainder = serialCRC(&h_num[i*num_words], page_size);
							//							if(verbosity >= 2) printf("Bitwise Computation: '%X'\n", cpu_remainder);
							//							if(cpu_remainder != ocl_remainders[i])
							//								fprintf(stderr,"ERROR: OCL and bitwise remainders for page %u differ [OCL: '%X', Bitwise: '%X']\n",i+1,ocl_remainders[i],cpu_remainder);
							//						}
							//						gettimeofday(&end,NULL);
							//						printf("Bitwise CRC Time: ");
							//						printTimeDiff(start,end);

							gettimeofday(&start,NULL);
							for(i=0; i<num_pages; i++)
							{
								cpu_remainder = crc32_8bytes(&h_num[i*num_words], page_size);
								if(verbosity >= 3) printf("CPU - Slice-by-8 Computation: '%X'\n", cpu_remainder);
								if(cpu_remainder != ocl_remainders[i])
								{
									/* [FIX-OCL-01] mismatches were only printed to stderr and the exit status was always 0; they are now counted
									 * and reported as SUCCESS/FAILURE (and in the DWARF_RESULT valid= field). */
									crc_bad++;
									fprintf(stderr,"ERROR: OCL and CPU Slice-by-8 remainders for page %u differ [OCL: '%X', CPU: '%X']\n",i+1,ocl_remainders[i],cpu_remainder);
								}
							}
							gettimeofday(&end,NULL);
							printf("CPU Slice-by-8 CRC Time: ");
							printTimeDiff(start,end);
						}
				}
			}
		}


#ifdef ENABLE_TIMER
		TIMER_DEST
#endif
			for(i=0; i<num_blocks; i++)
			{
				CUDA_CHECK(cudaFree(dev_input[i]));
				CUDA_CHECK(cudaFree(dev_output[i]));
			}
	}
	free(h_num);

	if(run_serial) printf("%s: %u pages differ from serial slice-by-8\n", crc_bad ? "FAILURE" : "SUCCESS", crc_bad);
	return dwarf_report(crc_bad == 0);
}
