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
  cl_kernel clKernel_diagonal;
  cl_kernel clKernel_perimeter;
  cl_kernel clKernel_internal;

  cl_int errcode;

  FILE *kernelFile;
  char *kernelSource;
  size_t kernelLength;

  cl_mem d_m;

  while ((opt = getopt_long(argc, argv, "::vs:i:", 
                            long_options, &option_index)) != -1 ) {
      switch(opt){
        case 'i':
          input_file = optarg;
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
      printf("Reading matrix from file %s\n", input_file);
      ret = create_matrix_from_file(&m, input_file, &matrix_dim);
      if (ret != RET_SUCCESS) {
          m = NULL;
          fprintf(stderr, "error create matrix from file %s\n", input_file);
          exit(EXIT_FAILURE);
      }
  } else {
    printf("No input file specified!\n");
    exit(EXIT_FAILURE);
  }

  if (do_verify){
    printf("Before LUD\n");
    print_matrix(m, matrix_dim);
    matrix_duplicate(m, &mm, matrix_dim);
  }

  errcode = clGetPlatformIDs(1, &clPlatform, NULL);
  CHECKERR(errcode);

  errcode = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &clDevice, NULL);
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
  if (errcode == CL_BUILD_PROGRAM_FAILURE)                                                                                                                                       
  {
    char *log;
    size_t logLength;
    errcode = clGetProgramBuildInfo(clProgram, clDevice, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLength);
    log = (char *) malloc(sizeof(char)*logLength);
    errcode = clGetProgramBuildInfo(clProgram, clDevice, CL_PROGRAM_BUILD_LOG, logLength, (void *) log, NULL);
    fprintf(stderr, "Kernel build error! Log:\n%s", log);
    free(log);
    return;
  }
  CHECKERR(errcode);

  clKernel_diagonal = clCreateKernel(clProgram, "lud_diagonal", &errcode);
  CHECKERR(errcode);
  clKernel_perimeter = clCreateKernel(clProgram, "lud_perimeter", &errcode);
  CHECKERR(errcode);
  clKernel_internal = clCreateKernel(clProgram, "lud_internal", &errcode);
  CHECKERR(errcode);

  d_m = clCreateBuffer(clContext, CL_MEM_READ_WRITE, matrix_dim*matrix_dim*sizeof(float), NULL, &errcode);
  CHECKERR(errcode);

  /* beginning of timing point */
  stopwatch_start(&sw);
  errcode = clEnqueueWriteBuffer(clCommands, d_m, CL_TRUE, 0, matrix_dim*matrix_dim*sizeof(float), (void *) m, 0, NULL, NULL);
  CHECKERR(errcode);

  int i=0;
  size_t localWorkSize[2];
  size_t globalWorkSize[2];
  //float *m_debug = (float*)malloc(matrix_dim*matrix_dim*sizeof(float));

  for (i=0; i < matrix_dim-BLOCK_SIZE; i += BLOCK_SIZE) {
      errcode = clSetKernelArg(clKernel_diagonal, 0, sizeof(cl_mem), (void *) &d_m);
      errcode |= clSetKernelArg(clKernel_diagonal, 1, sizeof(int), (void *) &matrix_dim);
      errcode |= clSetKernelArg(clKernel_diagonal, 2, sizeof(int), (void *) &i);
      CHECKERR(errcode);
      localWorkSize[0] = BLOCK_SIZE;
      globalWorkSize[0] = BLOCK_SIZE;
      errcode = clEnqueueNDRangeKernel(clCommands, clKernel_diagonal, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
      CHECKERR(errcode);
      errcode = clSetKernelArg(clKernel_perimeter, 0, sizeof(cl_mem), (void *) &d_m);
      errcode |= clSetKernelArg(clKernel_perimeter, 1, sizeof(int), (void *) &matrix_dim);
      errcode |= clSetKernelArg(clKernel_perimeter, 2, sizeof(int), (void *) &i);
      CHECKERR(errcode);
      localWorkSize[0] = BLOCK_SIZE*2;
      globalWorkSize[0] = ((matrix_dim-i)/BLOCK_SIZE-1)*localWorkSize[0];
      errcode = clEnqueueNDRangeKernel(clCommands, clKernel_perimeter, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
      CHECKERR(errcode);
      errcode = clSetKernelArg(clKernel_internal, 0, sizeof(cl_mem), (void *) &d_m);
      errcode |= clSetKernelArg(clKernel_internal, 1, sizeof(int), (void *) &matrix_dim);
      errcode |= clSetKernelArg(clKernel_internal, 2, sizeof(int), (void *) &i);
      CHECKERR(errcode);
      localWorkSize[0] = BLOCK_SIZE;
      localWorkSize[1] = BLOCK_SIZE;
      globalWorkSize[0] = ((matrix_dim-i)/BLOCK_SIZE-1)*localWorkSize[0];
      globalWorkSize[1] = ((matrix_dim-i)/BLOCK_SIZE-1)*localWorkSize[1];
      errcode = clEnqueueNDRangeKernel(clCommands, clKernel_internal, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
      CHECKERR(errcode);
  }
  errcode = clSetKernelArg(clKernel_diagonal, 0, sizeof(cl_mem), (void *) &d_m);
  errcode |= clSetKernelArg(clKernel_diagonal, 1, sizeof(int), (void *) &matrix_dim);
  errcode |= clSetKernelArg(clKernel_diagonal, 2, sizeof(int), (void *) &i);
  CHECKERR(errcode);
  localWorkSize[0] = BLOCK_SIZE;
  globalWorkSize[0] = BLOCK_SIZE;
  errcode = clEnqueueNDRangeKernel(clCommands, clKernel_diagonal, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  CHECKERR(errcode);

  clFinish(clCommands);

  errcode = clEnqueueReadBuffer(clCommands, d_m, CL_TRUE, 0, matrix_dim*matrix_dim*sizeof(float), (void *) m, 0, NULL, NULL);

  /* end of timing point */
  stopwatch_stop(&sw);
  printf("Time consumed(ms): %lf\n", 1000*get_interval_by_sec(&sw));

  clReleaseMemObject(d_m);

  if (do_verify){
    printf("After LUD\n");
    print_matrix(m, matrix_dim);
    printf(">>>Verify<<<<\n");
    lud_verify(mm, m, matrix_dim); 
    free(mm);
  }

  clReleaseKernel(clKernel_diagonal);
  clReleaseKernel(clKernel_perimeter);
  clReleaseKernel(clKernel_internal);
  clReleaseProgram(clProgram);
  clReleaseCommandQueue(clCommands);
  clReleaseContext(clContext);

  free(m);

  return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
