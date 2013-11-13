#include "global.h"
#include "functions.h"
#include "timeRec.h"
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#include "../../include/rdtsc.h"
#include "../../include/common_args.h"

/*************************************************************
 **************** Version 1 **********************************
32. This version is based on version 27.
	1) Constant memory will be used for scoring matrix
	2) Input string are stored in constant memory, it is only used for trace back, 
	   not be used for matrix filling.
db1.This version is based on version 32, modified for database search
db2.Add time record in the program.
verOpenCL4 Opencl implemenation, GPU lock-based sync
**************************************************************/


#define MIN(a, b) \
	(a < b ? a : b) \


char * loadSource(char *filePathName, size_t *fileSize)
{
	FILE *pfile;
	size_t tmpFileSize;
	char *fileBuffer;
	pfile = fopen(filePathName, "rb");

	if (pfile == NULL)
	{
		printf("Open file %s open error!\n", filePathName);
		return NULL;
	}

	fseek(pfile, 0, SEEK_END);
	tmpFileSize = ftell(pfile);

	fileBuffer = (char *)malloc(tmpFileSize);

	fseek(pfile, 0, SEEK_SET);
	fread(fileBuffer, sizeof(char), tmpFileSize, pfile);

	fclose(pfile);

	//debug================================
	//for (int i = 0; i < tmpFileSize; i++)
	//{
	//	printf("%c", fileBuffer[i]);
	//}
	//=====================================

	*fileSize = tmpFileSize;
	return fileBuffer;
}

int main(int argc, char ** argv)
{
	ocd_init(&argc, &argv, NULL);
	ocd_initCL();
		
	if (argc < 3)
	{
		printf("Calculate similarities between two strings.\n");
		printf("Maximum length of each string is: %d\n", MAX_LEN);
		printf("Usage: %s query database\n", argv[0]);
		printf("or: %s query database [openPenalty extensionPenalty block#]\n", argv[0]);
		printf("openPenalty (5.0), extensionPenalty (0.5)\n");
		return 1;
	}

	
	
	/////////////////////////////////////
	//      00 --> 01
	//		|	   |	
	//		10 --> 11
	////////////////////////////////////
	char queryFilePathName[255], dbDataFilePathName[255], dbLenFilePathName[255];
	int querySize, subSequenceNum, subSequenceSize;
	float openPenalty, extensionPenalty;
	int coalescedOffset = COALESCED_OFFSET;
	int nblosumWidth = 23;
	size_t blockSize = 64;
	size_t setZeroThreadNum, mfThreadNum;
	int blockNum = 14;

	cl_ulong maxLocalSize;

	int arraySize;

    struct timeval t1, t2;
	float tmpTime;
	FILE *pfile;

	//record time
	memset(&strTime, 0, sizeof(STRUCT_TIME));
	timerStart();

	openPenalty = 5.0f;
	extensionPenalty = 0.5;

	if (argc == 6)
	{
		openPenalty = atof(argv[3]);
		extensionPenalty = atof(argv[4]);
		blockNum = atoi(argv[5]);
	}
	
	//relocated to after MAX_COMPUTE_UNITS check
	//mfThreadNum = blockNum * blockSize;

	cl_program hProgram;
	cl_kernel hMatchStringKernel, hTraceBackKernel, hSetZeroKernel;
	size_t sourceFileSize;
	char *cSourceCL = NULL;

	//err = clGetPlatformIDs(1, &platformID, NULL);
	//CHKERR(err, "Get platform ID error!");

	cl_int err;

	//check to make sure the device supports this block count
	//then scale threads appropriately
	cl_uint devBlockNum = 0;
	CHKERR(clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS,\
		sizeof(cl_uint), &devBlockNum, 0), \
		"Error while querying CL_DEVICE_MAX_COMPUTE_UNITS.");
	if (devBlockNum == MIN(blockNum, devBlockNum)) {
		printf("Scaling blocks from %d to %d to fit on device\n",\
			blockNum, devBlockNum);
		blockNum = devBlockNum;
	}
	mfThreadNum = blockNum * blockSize;
	
	CHKERR(clGetDeviceInfo(device_id, CL_DEVICE_LOCAL_MEM_SIZE,\
		sizeof(cl_ulong), &maxLocalSize, 0), \
		"Error while querying CL_DEVICE_LOCAL_MEM_SIZE.");
	
	//load the source file
	char kernel_file[] = "kernels.cl";
	cSourceCL = loadSource(kernel_file, &sourceFileSize);

	hProgram = clCreateProgramWithSource(context, 1, (const char **)&cSourceCL, 
				&sourceFileSize, &err);
	CHKERR(err, "Create program with source error");

	err = clBuildProgram(hProgram, 0, 0, 0, 0, 0);
	//debug================================
	int logSize = 3000, i;
	size_t retSize;
	char logTxt[3000];
	err = clGetProgramBuildInfo(hProgram, device_id, CL_PROGRAM_BUILD_LOG, logSize, logTxt, &retSize);
	for (i = 0; i < retSize; i++)
	{
		printf("%c", logTxt[i]);
	}
	//===================================
	CHKERR(err, "Build program error");

	hMatchStringKernel = clCreateKernel(hProgram, "MatchStringGPUSync", &err);
	CHKERR(err, "Create MatchString kernel error");
	hTraceBackKernel = clCreateKernel(hProgram, "trace_back2", &err);
	CHKERR(err, "Create trace_back2 kernel error");
	hSetZeroKernel = clCreateKernel(hProgram, "setZero", &err);
	CHKERR(err, "Create setZero kernel error");

	sprintf(queryFilePathName, "%s", argv[1]);
	sprintf(dbDataFilePathName, "%s.data", argv[2]);
	sprintf(dbLenFilePathName, "%s.loc", argv[2]);

	char *allSequences, *querySequence, *subSequence;
	char *seq1, *seq2;
	cl_mem seq1D, seq2D;

	allSequences = new char[2 * (MAX_LEN)];
	if (allSequences == NULL)
	{
		printf("Allocate sequence buffer error!\n");
		return 1;
	}
	querySequence = allSequences;

	seq1D = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_char) * MAX_LEN, 0, &err);
	CHKERR(err, "Create seq1D memory");
	seq2D = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_char) * MAX_LEN, 0, &err);
	CHKERR(err, "Create seq2D memory");
	
	//read query sequence
	querySize = readQuerySequence(queryFilePathName, querySequence);
	if (querySize <= 0 || querySize > MAX_LEN)
	{
		printf("Query size %d is out of range (0, %d)\n",
				MAX_LEN,
				querySize);
		return 1;
	}
	encoding(querySequence, querySize);
	subSequence = allSequences + querySize;

	//allocate output sequence buffer
	char *outSeq1, *outSeq2;
	outSeq1 = new char[2 * MAX_LEN];
	outSeq2 = new char[2 * MAX_LEN];
	if (outSeq1 == NULL ||
		outSeq2 == NULL)
	{
		printf("Allocate output sequence buffer on host error!\n");
		return 1;
	}

	cl_mem outSeq1D, outSeq2D;
	outSeq1D = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_char) * MAX_LEN * 2, 0, &err);
	CHKERR(err, "Create outSeq1D memory");
	outSeq2D = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_char) * MAX_LEN * 2, 0, &err);
	CHKERR(err, "Create outSeq2D memory");

	//allocate thread number per launch and 
	//location difference information
	int *threadNum, *diffPos;
	threadNum = new int[2 * MAX_LEN];
	diffPos = new int[2 * MAX_LEN];
	if (threadNum == NULL ||
		diffPos == NULL)
	{
		printf("Allocate location buffer on host error!\n");
		return 1;
	}

	cl_mem threadNumD, diffPosD;
	threadNumD = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_int) * (2 * MAX_LEN), 0, &err);
	CHKERR(err, "Create threadNumD memory");
	diffPosD = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_int) * (2 * MAX_LEN), 0, &err);
	CHKERR(err, "Create diffPosD memory");

	//allocate matrix buffer
	char *pathFlag, *extFlag; 
	float *nGapDist, *hGapDist, *vGapDist;
	int maxElemNum = (MAX_LEN + 1) * (MAX_LEN + 1);
	pathFlag  = new char[maxElemNum];
	extFlag   = new char[maxElemNum];
	nGapDist = new float[maxElemNum];
	hGapDist = new float[maxElemNum];
	vGapDist = new float[maxElemNum];
	if (pathFlag  == NULL ||
		extFlag   == NULL ||
		nGapDist == NULL ||
		hGapDist == NULL ||
		vGapDist == NULL)
	{
		printf("Allocate DP matrices on host error!\n");
		return 1;
	}

	cl_mem pathFlagD, extFlagD,	nGapDistD, hGapDistD, vGapDistD;
	pathFlagD = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_char) * maxElemNum, 0, &err);
	CHKERR(err, "Create pathFlagD memory");
	extFlagD = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_char) * maxElemNum, 0, &err);
	CHKERR(err, "Create extFlagD memory");
	nGapDistD = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * maxElemNum, 0, &err);
	CHKERR(err, "Create nGapDistD memory");
	hGapDistD = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * maxElemNum, 0, &err);
	CHKERR(err, "Create hGapDistD memory");
	vGapDistD = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * maxElemNum, 0, &err);
	CHKERR(err, "Create vGapDistD memory");

	//Allocate the MAX INFO structure
	MAX_INFO *maxInfo;
	maxInfo = new MAX_INFO[1];
	if (maxInfo == NULL)
	{
		printf("Alloate maxInfo on host error!\n");
		return 1;
	}
	
	cl_mem maxInfoD;
	maxInfoD = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(MAX_INFO) * mfThreadNum, 0, &err);
	CHKERR(err, "Create maxInfoD memory");

	//allocate the distance table
	cl_mem blosum62D;
	int nblosumHeight = 23;
	blosum62D = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * nblosumWidth * nblosumHeight, 0, &err);
	err = clEnqueueWriteBuffer(commands, blosum62D, CL_TRUE, 0,
							   nblosumWidth * nblosumHeight * sizeof(cl_float), blosum62[0], 0, NULL, &ocdTempEvent);
        clFinish(commands);
        START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "SWAT Scoring Matrix Copy", ocdTempTimer)
        END_TIMER(ocdTempTimer)
	CHKERR(err, "copy blosum62 to device");
	cl_mem mutexMem;
	mutexMem = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_int), 0, &err);
	CHKERR(err, "create mutex mem error!");

	//copy the scoring matrix to the constant memory
	//copyScoringMatrixToConstant();

	//open the database
	pDBDataFile = fopen(dbDataFilePathName, "rb");
	if (pDBDataFile == NULL)
	{
		printf("DB data file %s open error!\n", dbDataFilePathName);
		return 1;
	}

	pDBLenFile = fopen(dbLenFilePathName, "rb");
	if (pDBLenFile == NULL)
	{
		printf("DB length file %s open error!\n", dbLenFilePathName);
		return 1;
	}

	//record time
	timerEnd();
	strTime.iniTime = elapsedTime();

	//read the total number of sequences
	fread(&subSequenceNum, sizeof(int), 1, pDBLenFile);

	//get the larger and smaller of the row and colum number
	int subSequenceNo, launchNum, launchNo;
	int rowNum, columnNum, matrixIniNum;
	int DPMatrixSize;
	int seq1Pos, seq2Pos, nOffset, startPos;

	for (subSequenceNo = 0; subSequenceNo < subSequenceNum; subSequenceNo++)
	{
		//record time
		timerStart();

		//read subject sequence
		fread(&subSequenceSize, sizeof(int), 1, pDBLenFile);
		if (subSequenceSize <= 0 || subSequenceSize > MAX_LEN)
		{
			printf("Size %d of bubject sequence %d is out of range!\n",
					subSequenceSize,
					subSequenceNo);
			break;
		}
		fread(subSequence, sizeof(char), subSequenceSize, pDBDataFile);

		gettimeofday(&t1, NULL);
		if (subSequenceSize > querySize)
		{
			seq1 = subSequence;
			seq2 = querySequence;
			rowNum = subSequenceSize + 1;
			columnNum = querySize + 1;
		}
		else
		{
			seq1 = querySequence;
			seq2 = subSequence;
			rowNum = querySize + 1;
			columnNum = subSequenceSize + 1;
		}

		launchNum = rowNum + columnNum - 1;

		//preprocessing for sequences
		DPMatrixSize = preProcessing(rowNum,
					  columnNum,
					  threadNum,
					  diffPos,
					  matrixIniNum);

		//record time
		timerEnd();
		strTime.preprocessingTime += elapsedTime();

		//record time
		timerStart();

		//use a kernel to initialize the matrix
		arraySize = DPMatrixSize * sizeof(char);
		setZeroThreadNum = ((arraySize - 1) / blockSize + 1) * blockSize;
		err  = clSetKernelArg(hSetZeroKernel, 0, sizeof(cl_mem), (void *)&pathFlagD);
		err |= clSetKernelArg(hSetZeroKernel, 1, sizeof(int), (void *)&arraySize);
		err |= clEnqueueNDRangeKernel(commands, hSetZeroKernel, 1, NULL, &setZeroThreadNum,
									 &blockSize, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SWAT DP Matrix Init", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		err |= clSetKernelArg(hSetZeroKernel, 0, sizeof(cl_mem), (void *)&extFlagD);
		err |= clEnqueueNDRangeKernel(commands, hSetZeroKernel, 1, NULL, &setZeroThreadNum,
									 &blockSize, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SWAT DP Matrix Init", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		CHKERR(err, "Initialize flag matrice");

		arraySize = matrixIniNum * sizeof(float);
		setZeroThreadNum = ((arraySize - 1) / blockSize + 1) * blockSize;
		err  = clSetKernelArg(hSetZeroKernel, 0, sizeof(cl_mem), (void *)&nGapDistD);
		err |= clSetKernelArg(hSetZeroKernel, 1, sizeof(int), (void *)&arraySize);
		err |= clEnqueueNDRangeKernel(commands, hSetZeroKernel, 1, NULL, &setZeroThreadNum,
									 &blockSize, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SWAT Distance Matrix Init", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		err |= clSetKernelArg(hSetZeroKernel, 0, sizeof(cl_mem), (void *)&hGapDistD);
		err |= clEnqueueNDRangeKernel(commands, hSetZeroKernel, 1, NULL, &setZeroThreadNum,
									 &blockSize, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SWAT Distance Matrix Init", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		err |= clSetKernelArg(hSetZeroKernel, 0, sizeof(cl_mem), (void *)&vGapDistD);
		err |= clEnqueueNDRangeKernel(commands, hSetZeroKernel, 1, NULL, &setZeroThreadNum,
									 &blockSize, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SWAT Distance Matrix Init", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		CHKERR(err, "Initialize dist matrice");

		arraySize = sizeof(MAX_INFO) * mfThreadNum;
		setZeroThreadNum = ((arraySize - 1) / blockSize + 1) * blockSize;
		err  = clSetKernelArg(hSetZeroKernel, 0, sizeof(cl_mem), (void *)&maxInfoD);
		err |= clSetKernelArg(hSetZeroKernel, 1, sizeof(int), (void *)&arraySize);
		err |= clEnqueueNDRangeKernel(commands, hSetZeroKernel, 1, NULL, &setZeroThreadNum,
									 &blockSize, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SWAT Max Info Matrix Init", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		CHKERR(err, "Initialize max info");

		arraySize = sizeof(int);
		setZeroThreadNum = ((arraySize - 1) / blockSize + 1) * blockSize;
		err  = clSetKernelArg(hSetZeroKernel, 0, sizeof(cl_mem), (void *)&mutexMem);
		err |= clSetKernelArg(hSetZeroKernel, 1, sizeof(int), (void *)&arraySize);
		err |= clEnqueueNDRangeKernel(commands, hSetZeroKernel, 1, NULL, &setZeroThreadNum,
									 &blockSize, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SWAT Mutex Init", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		CHKERR(err, "Initialize mutex variable");

		//copy input sequences to device
		err  = clEnqueueWriteBuffer(commands, seq1D, CL_FALSE, 0, (rowNum - 1) * sizeof(cl_char), seq1, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "SWAT Sequence Copy", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		err |= clEnqueueWriteBuffer(commands, seq2D, CL_FALSE, 0, (columnNum - 1) * sizeof(cl_char), seq2, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "SWAT Sequence Copy", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		CHKERR(err, "copy input sequence");

		err  = clEnqueueWriteBuffer(commands, diffPosD, CL_FALSE, 0, launchNum * sizeof(cl_int), diffPos, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "SWAT Mutex Info Copy", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		err |= clEnqueueWriteBuffer(commands, threadNumD, CL_FALSE, 0, launchNum * sizeof(cl_int), threadNum, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "SWAT Mutex Info Copy", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		CHKERR(err, "copy diffpos and/or threadNum mutexMem info error!");
		
		//record time
		timerEnd();
		strTime.copyTimeHostToDevice += elapsedTime();

		//record time
		timerStart();

		//set arguments
		err  = clSetKernelArg(hMatchStringKernel, 0, sizeof(cl_mem), (void *)&pathFlagD);
		err |= clSetKernelArg(hMatchStringKernel, 1, sizeof(cl_mem), (void *)&extFlagD);
		err |= clSetKernelArg(hMatchStringKernel, 2, sizeof(cl_mem), (void *)&nGapDistD);
		err |= clSetKernelArg(hMatchStringKernel, 3, sizeof(cl_mem), (void *)&hGapDistD);
		err |= clSetKernelArg(hMatchStringKernel, 4, sizeof(cl_mem), (void *)&vGapDistD);
		err |= clSetKernelArg(hMatchStringKernel, 5, sizeof(cl_mem), (void *)&diffPosD);
		err |= clSetKernelArg(hMatchStringKernel, 6, sizeof(cl_mem), (void *)&threadNumD);
		err |= clSetKernelArg(hMatchStringKernel, 7, sizeof(cl_int), (void *)&rowNum);
		err |= clSetKernelArg(hMatchStringKernel, 8, sizeof(cl_int), (void *)&columnNum);
		err |= clSetKernelArg(hMatchStringKernel, 9, sizeof(cl_mem), (void *)&seq1D);
		err |= clSetKernelArg(hMatchStringKernel, 10, sizeof(cl_mem), (void *)&seq2D);	
		err |= clSetKernelArg(hMatchStringKernel, 11, sizeof(cl_int), (void *)&nblosumWidth);
		err |= clSetKernelArg(hMatchStringKernel, 12, sizeof(cl_float), (void *)&openPenalty);
		err |= clSetKernelArg(hMatchStringKernel, 13, sizeof(cl_float), (void *)&extensionPenalty);
		err |= clSetKernelArg(hMatchStringKernel, 14, sizeof(cl_mem), (void *)&maxInfoD);
		err |= clSetKernelArg(hMatchStringKernel, 15, sizeof(cl_mem), (void *)&blosum62D);
		err |= clSetKernelArg(hMatchStringKernel, 16, sizeof(cl_mem), (void *)&mutexMem);
		//err |= clSetKernelArg(hMatchStringKernel, 17, maxLocalSize, NULL);
		CHKERR(err, "Set match string argument error!");

		err = clEnqueueNDRangeKernel(commands, hMatchStringKernel, 1, NULL, &mfThreadNum,
									 &blockSize, 0, NULL, &ocdTempEvent);
        clFinish(commands);
        START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SWAT Kernels", ocdTempTimer)
        END_TIMER(ocdTempTimer)
		CHKERR(err, "Launch kernel match string error");

		//record time
		timerEnd();
		strTime.matrixFillingTime += elapsedTime();

		//record time
		timerStart();
		err  = clSetKernelArg(hTraceBackKernel, 0, sizeof(cl_mem), (void *)&pathFlagD);
		err |= clSetKernelArg(hTraceBackKernel, 1, sizeof(cl_mem), (void *)&extFlagD);
		err |= clSetKernelArg(hTraceBackKernel, 2, sizeof(cl_mem), (void *)&diffPosD);
		err |= clSetKernelArg(hTraceBackKernel, 3, sizeof(cl_mem), (void *)&seq1D);
		err |= clSetKernelArg(hTraceBackKernel, 4, sizeof(cl_mem), (void *)&seq2D);	
		err |= clSetKernelArg(hTraceBackKernel, 5, sizeof(cl_mem), (void *)&outSeq1D);
		err |= clSetKernelArg(hTraceBackKernel, 6, sizeof(cl_mem), (void *)&outSeq2D);	
		err |= clSetKernelArg(hTraceBackKernel, 7, sizeof(cl_mem), (void *)&maxInfoD);
		err |= clSetKernelArg(hTraceBackKernel, 8, sizeof(int), (void *)&mfThreadNum);
		
		size_t tbGlobalSize[1] = {1};
		size_t tbLocalSize[1]  = {1};
		err = clEnqueueNDRangeKernel(commands, hTraceBackKernel, 1, NULL, tbGlobalSize,
									 tbLocalSize, 0, NULL, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "SWAT Kernels", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		CHKERR(err, "Launch kernel trace back error");
		clFinish(commands);
		//record time
		timerEnd();
		strTime.traceBackTime += elapsedTime();

		//record time
		timerStart();
		//copy matrix score structure back
		err = clEnqueueReadBuffer(commands, maxInfoD, CL_FALSE, 0, sizeof(MAX_INFO),
								  maxInfo, 0, 0, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "SWAT Max Info Copy", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		CHKERR(err, "Read maxInfo buffer error!");

		int maxOutputLen = rowNum + columnNum - 2;
		err  = clEnqueueReadBuffer(commands, outSeq1D, CL_FALSE, 0, maxOutputLen * sizeof(cl_char),
								   outSeq1, 0, 0, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "SWAT Sequence Copy", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		err != clEnqueueReadBuffer(commands, outSeq2D, CL_FALSE, 0, maxOutputLen * sizeof(cl_char),
								   outSeq2, 0, 0, &ocdTempEvent);
                clFinish(commands);
                START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "SWAT Sequence Copy", ocdTempTimer)
                END_TIMER(ocdTempTimer)
		CHKERR(err, "Read output sequence error!");
		//record time
		clFinish(commands);
		gettimeofday(&t2, NULL);
		timerEnd();
		strTime.copyTimeDeviceToHost += elapsedTime();

		//call the print function to print the match result
		printf("============================================================\n");
		printf("Sequence pair %d:\n", subSequenceNo);
		int nlength = maxInfo->noutputlen;
		PrintAlignment(outSeq1, outSeq2, nlength, CHAR_PER_LINE, openPenalty, extensionPenalty);
		printf("Max alignment score (on device) is %.1f\n", maxInfo->fmaxscore);
		//obtain max alignment score on host
		//err = clEnqueueReadBuffer(commands, nGapDistD, CL_TRUE, 0, sizeof(cl_float) * DPMatrixSize,
		//						  nGapDist, 0, 0, 0);
		//printf("Max alignment score (on host) is %.1f\n", maxScore(nGapDist, DPMatrixSize));

		printf("openPenalty = %.1f, extensionPenalty = %.1f\n", openPenalty, extensionPenalty);
		printf("Input sequence size, querySize: %d, subSequenceSize: %d\n", 
				querySize, subSequenceSize);

		printf("Max position, seq1 = %d, seq2 = %d\n", maxInfo->nposi, maxInfo->nposj);
	}
	tmpTime = 1000.0 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000.0;
	pfile = fopen("../kernelTime.txt", "at");
	fprintf(pfile, "verOpencl4:\t%.3f\n", tmpTime);
	fclose(pfile);

	//print time
	printTime_toStandardOutput();
	printTime_toFile();

	fclose(pDBLenFile);
	fclose(pDBDataFile);

	clReleaseKernel(hMatchStringKernel);
	clReleaseKernel(hTraceBackKernel);
	clReleaseKernel(hSetZeroKernel);

	delete allSequences;
	clReleaseMemObject(seq1D);
	clReleaseMemObject(seq2D);

	delete outSeq1;
	delete outSeq2;
	clReleaseMemObject(outSeq1D);
	clReleaseMemObject(outSeq2D);

	delete threadNum;
	clReleaseMemObject(threadNumD);
	delete diffPos;
	clReleaseMemObject(diffPosD);

	delete pathFlag;
	delete extFlag;
	delete nGapDist;
	delete hGapDist;
	delete vGapDist;
	clReleaseMemObject(pathFlagD);
	clReleaseMemObject(extFlagD);
	clReleaseMemObject(nGapDistD);
	clReleaseMemObject(hGapDistD);
	clReleaseMemObject(vGapDistD);

	delete maxInfo;
	clReleaseMemObject(maxInfoD);

	free(cSourceCL);

	clReleaseMemObject(blosum62D);
	clReleaseMemObject(mutexMem);

	clReleaseProgram(hProgram);
	clReleaseCommandQueue(commands);
	clReleaseContext(context);
	ocd_finalize();
	return 0;
}


