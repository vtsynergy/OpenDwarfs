// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "srad.h"
#include <string.h>

int  BLOCK_SIZE = 16;
#include "../../include/rdtsc.h"
int platform_id = PLATFORM_ID, device_id = DEVICE_ID;
// includes, project


#define CHECKERR(err) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "Error: %d in line: %d\n", err, __LINE__);\
        exit(1); \
    }
//#define USEGPU 1
void random_matrix(float *I, int rows, int cols);
void runTest( int argc, char** argv);
void usage(int argc, char **argv)
{
	fprintf(stderr, "Usage: %s <rows> <cols> <y1> <y2> <x1> <x2> <lamda> <no. of iter> [platform & device]\n", argv[0]);
	fprintf(stderr, "\t<rows>   - number of rows\n");
	fprintf(stderr, "\t<cols>    - number of cols\n");
	fprintf(stderr, "\t<y1> 	 - y1 value of the speckle\n");
	fprintf(stderr, "\t<y2>      - y2 value of the speckle\n");
	fprintf(stderr, "\t<x1>       - x1 value of the speckle\n");
	fprintf(stderr, "\t<x2>       - x2 value of the speckle\n");
	fprintf(stderr, "\t<lamda>   - lambda (0,1)\n");
	fprintf(stderr, "\t<no. of iter>   - number of iterations\n");
	fprintf(stderr, "\t<platform>   - index of platform\n");
	fprintf(stderr, "\t<device>   - index of device\n");
	
	exit(1);
}
////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char** argv) 
{
    OCD_INIT
		runTest( argc, argv);
    OCD_FINISH
    return EXIT_SUCCESS;
}


void
runTest( int argc, char** argv) 
{
    int rows, cols, size_I, size_R, niter = 10, iter;
    float *I, *J, lambda, q0sqr, sum, sum2, tmp, meanROI,varROI ;
    int i, j, k;

#ifdef CPU
	float Jc, G2, L, num, den, qsqr;
	int *iN,*iS,*jE,*jW, k;
	float *dN,*dS,*dW,*dE;
	float cN,cS,cW,cE,D;
#endif

#ifdef GPU
	
    cl_device_id clDevice;
    cl_context clContext;
    cl_command_queue clCommands;
    cl_program clProgram;
    cl_kernel clKernel_srad1;
    cl_kernel clKernel_srad2;

    cl_int errcode;

	cl_mem J_cuda;
    cl_mem C_cuda;
	cl_mem E_C, W_C, N_C, S_C;

    FILE *kernelFile;
    char *kernelSource;
    size_t kernelLength;
	if(argc == 11){
	platform_id=atoi(argv[9]);
	device_id = atoi(argv[10]);
	}
	clDevice = GetDevice(platform_id, device_id);
    	size_t max_worksize[3];
errcode = clGetDeviceInfo(clDevice, CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t)*3, &max_worksize, NULL);
 CHECKERR(errcode);
        while(BLOCK_SIZE*BLOCK_SIZE>max_worksize[0])
                BLOCK_SIZE = BLOCK_SIZE/2;

    clContext = clCreateContext(NULL, 1, &clDevice, NULL, NULL, &errcode);
    CHECKERR(errcode);

    clCommands = clCreateCommandQueue(clContext, clDevice, CL_QUEUE_PROFILING_ENABLE, &errcode);
    CHECKERR(errcode);

    kernelFile = fopen("srad_kernel.cl", "r");
    fseek(kernelFile, 0, SEEK_END);
    kernelLength = (size_t) ftell(kernelFile);
    kernelSource = (char *) malloc(sizeof(char)*kernelLength);
    rewind(kernelFile);
    fread((void *) kernelSource, kernelLength, 1, kernelFile);
    fclose(kernelFile);

    clProgram = clCreateProgramWithSource(clContext, 1, (const char **) &kernelSource, &kernelLength, &errcode);
    CHECKERR(errcode);

    free(kernelSource);
	char arg[50];
	sprintf(arg,"-D BLOCK_SIZE=%d",BLOCK_SIZE);
    errcode = clBuildProgram(clProgram, 1, &clDevice, arg, NULL, NULL);
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

    clKernel_srad1 = clCreateKernel(clProgram, "srad_cuda_1", &errcode);
    CHECKERR(errcode);
    clKernel_srad2 = clCreateKernel(clProgram, "srad_cuda_2", &errcode);
    CHECKERR(errcode);

#endif

	unsigned int r1, r2, c1, c2;
	float *c;
    
	
 
	if (argc == 9 || argc == 11)
	{
		rows = atoi(argv[1]);  //number of rows in the domain
		cols = atoi(argv[2]);  //number of cols in the domain
		if ((rows%16!=0) || (cols%16!=0)){
		fprintf(stderr, "rows and cols must be multiples of 16\n");
		exit(1);
		}
		r1   = atoi(argv[3]);  //y1 position of the speckle
		r2   = atoi(argv[4]);  //y2 position of the speckle
		c1   = atoi(argv[5]);  //x1 position of the speckle
		c2   = atoi(argv[6]);  //x2 position of the speckle
		lambda = atof(argv[7]); //Lambda value
		niter = atoi(argv[8]); //number of iterations
	}
    else{
	usage(argc, argv);
    }



	size_I = cols * rows;
    size_R = (r2-r1+1)*(c2-c1+1);   

	I = (float *)malloc( size_I * sizeof(float) );
    J = (float *)malloc( size_I * sizeof(float) );
	c  = (float *)malloc(sizeof(float)* size_I) ;


#ifdef CPU

    iN = (int *)malloc(sizeof(unsigned int*) * rows) ;
    iS = (int *)malloc(sizeof(unsigned int*) * rows) ;
    jW = (int *)malloc(sizeof(unsigned int*) * cols) ;
    jE = (int *)malloc(sizeof(unsigned int*) * cols) ;    


	dN = (float *)malloc(sizeof(float)* size_I) ;
    dS = (float *)malloc(sizeof(float)* size_I) ;
    dW = (float *)malloc(sizeof(float)* size_I) ;
    dE = (float *)malloc(sizeof(float)* size_I) ;    
    

    for (i=0; i< rows; i++) {
        iN[i] = i-1;
        iS[i] = i+1;
    }    
    for (j=0; j< cols; j++) {
        jW[j] = j-1;
        jE[j] = j+1;
    }
    iN[0]    = 0;
    iS[rows-1] = rows-1;
    jW[0]    = 0;
    jE[cols-1] = cols-1;

#endif

#ifdef GPU

	//Allocate device memory
    J_cuda = clCreateBuffer(clContext, CL_MEM_READ_WRITE, sizeof(float)*size_I, NULL, &errcode);
    CHECKERR(errcode);
    C_cuda = clCreateBuffer(clContext, CL_MEM_READ_WRITE, sizeof(float)*size_I, NULL, &errcode);
    CHECKERR(errcode);
    E_C = clCreateBuffer(clContext, CL_MEM_READ_WRITE, sizeof(float)*size_I, NULL, &errcode);
    CHECKERR(errcode);
    W_C = clCreateBuffer(clContext, CL_MEM_READ_WRITE, sizeof(float)*size_I, NULL, &errcode);
    CHECKERR(errcode);
    S_C = clCreateBuffer(clContext, CL_MEM_READ_WRITE, sizeof(float)*size_I, NULL, &errcode);
    CHECKERR(errcode);
    N_C = clCreateBuffer(clContext, CL_MEM_READ_WRITE, sizeof(float)*size_I, NULL, &errcode);
    CHECKERR(errcode);

	
#endif 

	printf("Randomizing the input matrix\n");
	//Generate a random matrix
	random_matrix(I, rows, cols);

    for (k = 0;  k < size_I; k++ ) {
     	J[k] = (float)exp(I[k]) ;
    }
	printf("Start the SRAD main loop\n");
 for (iter=0; iter< niter; iter++){     
		sum=0; sum2=0;
        for (i=r1; i<=r2; i++) {
            for (j=c1; j<=c2; j++) {
                tmp   = J[i * cols + j];
                sum  += tmp ;
                sum2 += tmp*tmp;
            }
        }
        meanROI = sum / size_R;
        varROI  = (sum2 / size_R) - meanROI*meanROI;
        q0sqr   = varROI / (meanROI*meanROI);

#ifdef CPU
        
		for (i = 0 ; i < rows ; i++) {
            for (j = 0; j < cols; j++) { 
		
				k = i * cols + j;
				Jc = J[k];
 
				// directional derivates
                dN[k] = J[iN[i] * cols + j] - Jc;
                dS[k] = J[iS[i] * cols + j] - Jc;
                dW[k] = J[i * cols + jW[j]] - Jc;
                dE[k] = J[i * cols + jE[j]] - Jc;
			
                G2 = (dN[k]*dN[k] + dS[k]*dS[k] 
                    + dW[k]*dW[k] + dE[k]*dE[k]) / (Jc*Jc);

   		        L = (dN[k] + dS[k] + dW[k] + dE[k]) / Jc;

				num  = (0.5*G2) - ((1.0/16.0)*(L*L)) ;
                den  = 1 + (.25*L);
                qsqr = num/(den*den);
 
                // diffusion coefficent (equ 33)
                den = (qsqr-q0sqr) / (q0sqr * (1+q0sqr)) ;
                c[k] = 1.0 / (1.0+den) ;
                
                // saturate diffusion coefficent
                if (c[k] < 0) {c[k] = 0;}
                else if (c[k] > 1) {c[k] = 1;}
		}
	}
         for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {        

                // current index
                k = i * cols + j;
                
                // diffusion coefficent
					cN = c[k];
					cS = c[iS[i] * cols + j];
					cW = c[k];
					cE = c[i * cols + jE[j]];

                // divergence (equ 58)
                D = cN * dN[k] + cS * dS[k] + cW * dW[k] + cE * dE[k];
                
                // image update (equ 61)
                J[k] = J[k] + 0.25*lambda*D;
            }
	}

#endif // CPU


#ifdef GPU

	//Currently the input size must be divided by 16 - the block size
	int block_x = cols/BLOCK_SIZE ;
    int block_y = rows/BLOCK_SIZE ;

    size_t localWorkSize[2] = {BLOCK_SIZE, BLOCK_SIZE};
    size_t globalWorkSize[2] = {block_x*localWorkSize[0], block_y*localWorkSize[1]};


	//Copy data from main memory to device memory
	errcode = clEnqueueWriteBuffer(clCommands, J_cuda, CL_TRUE, 0, sizeof(float)*size_I, (void *) J, 0, NULL, &ocdTempEvent);

        clFinish(clCommands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "SRAD Data Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHECKERR(errcode);

	//Run kernels
    errcode = clSetKernelArg(clKernel_srad1, 0, sizeof(cl_mem), (void *) &E_C);
    errcode |= clSetKernelArg(clKernel_srad1, 1, sizeof(cl_mem), (void *) &W_C);
    errcode |= clSetKernelArg(clKernel_srad1, 2, sizeof(cl_mem), (void *) &N_C);
    errcode |= clSetKernelArg(clKernel_srad1, 3, sizeof(cl_mem), (void *) &S_C);
    errcode |= clSetKernelArg(clKernel_srad1, 4, sizeof(cl_mem), (void *) &J_cuda);
    errcode |= clSetKernelArg(clKernel_srad1, 5, sizeof(cl_mem), (void *) &C_cuda);
    errcode |= clSetKernelArg(clKernel_srad1, 6, sizeof(int), (void *) &cols);
    errcode |= clSetKernelArg(clKernel_srad1, 7, sizeof(int), (void *) &rows);
    errcode |= clSetKernelArg(clKernel_srad1, 8, sizeof(float), (void *) &q0sqr);
    CHECKERR(errcode);
    errcode = clEnqueueNDRangeKernel(clCommands, clKernel_srad1, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
    clFinish(clCommands);
        START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SRAD Kernels", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHECKERR(errcode);
    errcode = clSetKernelArg(clKernel_srad2, 0, sizeof(cl_mem), (void *) &E_C);
    errcode |= clSetKernelArg(clKernel_srad2, 1, sizeof(cl_mem), (void *) &W_C);
    errcode |= clSetKernelArg(clKernel_srad2, 2, sizeof(cl_mem), (void *) &N_C);
    errcode |= clSetKernelArg(clKernel_srad2, 3, sizeof(cl_mem), (void *) &S_C);
    errcode |= clSetKernelArg(clKernel_srad2, 4, sizeof(cl_mem), (void *) &J_cuda);
    errcode |= clSetKernelArg(clKernel_srad2, 5, sizeof(cl_mem), (void *) &C_cuda);
    errcode |= clSetKernelArg(clKernel_srad2, 6, sizeof(int), (void *) &cols);
    errcode |= clSetKernelArg(clKernel_srad2, 7, sizeof(int), (void *) &rows);
    errcode |= clSetKernelArg(clKernel_srad2, 8, sizeof(float), (void *) &lambda);
    errcode |= clSetKernelArg(clKernel_srad2, 9, sizeof(float), (void *) &q0sqr);
    CHECKERR(errcode);
    errcode = clEnqueueNDRangeKernel(clCommands, clKernel_srad2, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
    clFinish(clCommands);
	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SRAD Kernels", ocdTempTimer)
    	END_TIMER(ocdTempTimer)
		CHECKERR(errcode);

	//Copy data from device memory to main memory
	errcode = clEnqueueReadBuffer(clCommands, J_cuda, CL_TRUE, 0, sizeof(float)*size_I, (void *) J, 0, NULL, &ocdTempEvent);
        clFinish(clCommands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SRAD Data Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHECKERR(errcode);

#endif   
}

#ifdef GPU

    clFinish(clCommands);

#endif

#ifdef OUTPUT
    //Printing output	
		printf("Printing Output:\n"); 
    for( i = 0 ; i < rows ; i++){
		for ( j = 0 ; j < cols ; j++){
         printf("%.5f ", J[i * cols + j]); 
		}	
     printf("\n"); 
   }
#endif 

	printf("Computation Done\n");

	free(I);
	free(J);
#ifdef CPU
	free(iN); free(iS); free(jW); free(jE);
    free(dN); free(dS); free(dW); free(dE);
#endif
#ifdef GPU
    clReleaseMemObject(C_cuda);
    clReleaseMemObject(J_cuda);
    clReleaseMemObject(E_C);
    clReleaseMemObject(W_C);
    clReleaseMemObject(N_C);
    clReleaseMemObject(S_C);
    clReleaseKernel(clKernel_srad1);
    clReleaseKernel(clKernel_srad2);
    clReleaseProgram(clProgram);
    clReleaseCommandQueue(clCommands);
    clReleaseContext(clContext);
#endif 
	free(c);
  
}


void random_matrix(float *I, int rows, int cols){
    
    int i, j;
	srand(7);
	
	for( i = 0 ; i < rows ; i++){
		for ( j = 0 ; j < cols ; j++){
		 I[i * cols + j] = rand()/(float)RAND_MAX ;
		}
	}

}

