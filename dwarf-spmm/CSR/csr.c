#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <assert.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "common.h"
#include "sparse_formats.h"

#define CHECKERR(err) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "Error: %d\n", err);\
        exit(1); \
    }

#define BLOCK_SIZE 16

static int do_verify = 0;

static struct option long_options[] = {
      /* name, has_arg, flag, val */
      {"input", 1, NULL, 'i'},
      {"size", 1, NULL, 's'},
      {"verify", 0, NULL, 'v'},
      {0,0,0,0}
};

int
main ( int argc, char *argv[] )
{
  int matrix_dim = 32; /* default matrix_dim */
  int opt, option_index=0;
  func_ret_t ret;
  const char *input_file = NULL;
  float *m, *mm;
  stopwatch sw;

  cl_platform_id clPlatform;
  cl_device_id clDevice;
  cl_context clContext;
  cl_command_queue clCommands;
  cl_program clProgram;
  //cl_kernel clKernel_diagonal;
  cl_kernel clKernel_csr;

  cl_int errcode;

  FILE *kernelFile;
  char *kernelSource;
  size_t kernelLength;

  //cl_mem d_m;
  cl_mem x_loc;
  cl_mem y_loc;
  cl_mem x_csr;
  cl_mem y_csr;


  while ((opt = getopt_long(argc, argv, "::vs:i:", 
                            long_options, &option_index)) != -1 ) {
      switch(opt){
        case 'i':
          input_file = optarg;
	//fprintf(stderr, "shit\n");
          break;
        case 'v':
          do_verify = 1;
          break;
        case 's':
          matrix_dim = atoi(optarg);
          fprintf(stderr, "Currently not supported, use -i instead\n");
          fprintf(stderr, "Usage: %s [-v] [-s matrix_size|-i input_file]\n", argv[0]);
          exit(EXIT_FAILURE);
        case '?':
          fprintf(stderr, "invalid option\n");
          break;
        case ':':
          fprintf(stderr, "missing argument\n");
          break;
        default:
          fprintf(stderr, "Usage: %s [-v] [-s matrix_size|-i input_file]\n",
                  argv[0]);
          exit(EXIT_FAILURE);
      }
  }
  
  if ( (optind < argc) || (optind == 1)) {
      fprintf(stderr, "Usage: %s [-v] [-s matrix_size|-i input_file]\n", argv[0]);
      exit(EXIT_FAILURE);
  }

  if (input_file) {
      //needs to be completely redone for the random matrix
      //printf("Reading matrix from file %s\n", input_file);
      //ret = create_matrix_from_file(&m, input_file, &matrix_dim);
      //if (ret != RET_SUCCESS) {
      //    m = NULL;
      //    fprintf(stderr, "error create matrix from file %s\n", input_file);
      //    exit(EXIT_FAILURE);
      //}


  //modify to print my matrix
  if (do_verify){

    //printf("Before LUD\n");
    //print_matrix(m, matrix_dim);
    //matrix_duplicate(m, &mm, matrix_dim);
  }
printf("testing\n");
  errcode = clGetPlatformIDs(1, &clPlatform, NULL);
  CHECKERR(errcode);

  //errcode = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &clDevice, NULL);
  errcode = clGetDeviceIDs(clPlatform, CL_DEVICE_TYPE_GPU, 1, &clDevice, NULL);
  CHECKERR(errcode);

  clContext = clCreateContext(NULL, 1, &clDevice, NULL, NULL, &errcode);
  CHECKERR(errcode);

  clCommands = clCreateCommandQueue(clContext, clDevice, 0, &errcode);
  CHECKERR(errcode);

  kernelFile = fopen("lud_kernel.cl", "r");
  fseek(kernelFile, 0, SEEK_END);
  kernelLength = (size_t) ftell(kernelFile);
  kernelSource = (char *) malloc(sizeof(char)*kernelLength);
  rewind(kernelFile);
  fread((void *) kernelSource, kernelLength, 1, kernelFile);
  fclose(kernelFile);

  clProgram = clCreateProgramWithSource(clContext, 1, (const char **) &kernelSource, &kernelLength, &errcode);
  CHECKERR(errcode);

  free(kernelSource);

  errcode = clBuildProgram(clProgram, 1, &clDevice, NULL, NULL, NULL);

printf("testing2\n");
  if (errcode == CL_BUILD_PROGRAM_FAILURE)                                                                                                                                       
  {
printf("shouldnt execute\n");
    char *log;
    size_t logLength;
    errcode = clGetProgramBuildInfo(clProgram, clDevice, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLength);
    log = (char *) malloc(sizeof(char)*logLength);
    errcode = clGetProgramBuildInfo(clProgram, clDevice, CL_PROGRAM_BUILD_LOG, logLength, (void *) log, NULL);
    fprintf(stderr, "Kernel build error! Log:\n%s", log);
    free(log);
    return;
  }
printf("testing4\n");
  CHECKERR(errcode);



      
  } else {
    printf("No input file specified!\n");
    exit(EXIT_FAILURE);
  }
printf("testing3\n");
  //need to update for new kernels
  //clKernel_diagonal = clCreateKernel(clProgram, "lud_diagonal", &errcode);
  //CHECKERR(errcode);
  //clKernel_perimeter = clCreateKernel(clProgram, "lud_perimeter", &errcode);
  //CHECKERR(errcode);
  //clKernel_internal = clCreateKernel(clProgram, "lud_internal", &errcode);
  //CHECKERR(errcode);
  clKernel_csr = clCreateKernel(clProgram, "spmv_csr_scalar_kernel", &errcode);
  CHECKERR(errcode);

  //probably needs to be reworked
  //d_m = clCreateBuffer(clContext, CL_MEM_READ_WRITE, matrix_dim*matrix_dim*sizeof(float), NULL, &errcode);
  //CHECKERR(errcode);
  //d_m = clCreateBuffer(clContext, CL_MEM_READ_WRITE, matrix_dim*matrix_dim*sizeof(float), NULL, &errcode);
  //CHECKERR(errcode);

  /* beginning of timing point */
  stopwatch_start(&sw);
  //probably needs to be reworked
  //errcode = clEnqueueWriteBuffer(clCommands, d_m, CL_TRUE, 0, matrix_dim*matrix_dim*sizeof(float), (void *) m, 0, NULL, NULL);
  //CHECKERR(errcode);


  //create my random matrix
//printf("bleh");
      csr_matrix csr;
      csr = laplacian_5pt(512);
//printf("I need to give up on making a matrix");
      // fill matrix with random values: some matrices have extreme values, 
      // which makes correctness testing difficult, especially in single precision
      srand(13);
      unsigned int k = 0;
      for(k = 0; k < csr.num_nonzeros; k++){
         csr.Ax[k] = 1.0 - 2.0 * (rand() / (RAND_MAX + 1.0)); 
      }
printf("matrix made \n");

//ahhhh!
  //<IndexType,ValueType>
	//<unsigned int, float>
    float * x_host = float_new_array(csr.num_cols);
    float * y_host = float_new_array(csr.num_rows);
    unsigned int l = 0;
    for( l = 0; l < csr.num_cols; l++)
        x_host[l] = rand() / (RAND_MAX + 1.0);
    unsigned int o = 0;
    printf("god damnit \n");
    for( o = 0; o < csr.num_rows; o++)
        y_host[o] = 0; 


    //fill(y_host, y_host + csr.num_rows, 0);
   
    //initialize device arrays
    //float * x_loc = copy_array(x_host, csr.num_cols, HOST_MEMORY, loc);
    //float * y_loc = copy_array(y_host, csr.num_rows, HOST_MEMORY, loc);

   printf("construcing memory on device \n");
	x_loc = clCreateBuffer(clContext, CL_MEM_READ_WRITE, csr.num_cols*sizeof(float), NULL, &errcode);
  CHECKERR(errcode);
	y_loc = clCreateBuffer(clContext, CL_MEM_READ_WRITE, csr.num_rows*sizeof(float), NULL, &errcode);
  CHECKERR(errcode);
	x_csr = clCreateBuffer(clContext, CL_MEM_READ_WRITE, csr.num_cols*sizeof(float), NULL, &errcode);
  CHECKERR(errcode);
	y_csr = clCreateBuffer(clContext, CL_MEM_READ_WRITE, csr.num_rows*sizeof(float), NULL, &errcode);
  CHECKERR(errcode);

  printf("coppying to memory \n");
  errcode = clEnqueueWriteBuffer(clCommands, x_loc, CL_TRUE, 0, csr.num_cols*sizeof(float), (float *) x_host, 0, NULL, NULL);
  CHECKERR(errcode);
  errcode = clEnqueueWriteBuffer(clCommands, y_loc, CL_TRUE, 0, csr.num_rows*sizeof(float), (float *) y_host, 0, NULL, NULL);
  CHECKERR(errcode);
  errcode = clEnqueueWriteBuffer(clCommands, x_csr, CL_TRUE, 0, csr.num_cols*sizeof(float), (float *) csr.Ap, 0, NULL, NULL);
  CHECKERR(errcode);
  errcode = clEnqueueWriteBuffer(clCommands, y_csr, CL_TRUE, 0, csr.num_rows*sizeof(float), (float *) csr.Aj, 0, NULL, NULL);
  CHECKERR(errcode);


    //spmv_csr_scalar_kernel(csr.num_rows, x_csr, y_csr, x_loc, y_loc);

  printf("getting ready to benchmark \n");
  int i=0;
  size_t localWorkSize[2];
  size_t globalWorkSize[2];
  //float *m_debug = (float*)malloc(matrix_dim*matrix_dim*sizeof(float));

  //needs to be completely redone for controlling kernel

  //clKernel_csr(csr.num_rows, x_csr, y_csr, x_loc, y_loc);
      errcode = clSetKernelArg(clKernel_csr, 0, sizeof(unsigned int), (void *) &csr.num_rows);
      errcode = clSetKernelArg(clKernel_csr, 1, sizeof(cl_mem), (void *) &x_csr);
      errcode |= clSetKernelArg(clKernel_csr, 2, sizeof(cl_mem), (void *) &y_csr);
      errcode |= clSetKernelArg(clKernel_csr, 3, sizeof(cl_mem), (void *) &x_loc);
      errcode |= clSetKernelArg(clKernel_csr, 4, sizeof(cl_mem), (void *) &y_loc);
      CHECKERR(errcode);
      localWorkSize[0] = BLOCK_SIZE;
      globalWorkSize[0] = BLOCK_SIZE;
      errcode = clEnqueueNDRangeKernel(clCommands, clKernel_csr, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);

      clFinish(clCommands);



  /* end of timing point */
  stopwatch_stop(&sw);
  printf("Time consumed(ms): %lf\n", 1000*get_interval_by_sec(&sw));

  //clReleaseMemObject(d_m);

  if (do_verify){
    //printf("After LUD\n");
    //print_matrix(m, matrix_dim);
    //printf(">>>Verify<<<<\n");
    //lud_verify(mm, m, matrix_dim); 
    //free(mm);
  }

  //clReleaseKernel(clKernel_diagonal);
  //clReleaseKernel(clKernel_perimeter);
  //clReleaseKernel(clKernel_internal);
  //clReleaseProgram(clProgram);
  //clReleaseCommandQueue(clCommands);
  //clReleaseContext(clContext);

  //free(m);

  return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
