#include "global.h"
#include "functions.h"
#include "timeRec.h"
#include "dwarf_cuda.h"
#include "dwarf_timer.h"
#define AOCL_ALIGNMENT 64
#include "swat_kernel.cuh"

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


int main(int argc, char ** argv)
{
	dwarf_cuda_init(&argc, &argv);

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

	//check to make sure the device supports this block count (CL_DEVICE_MAX_COMPUTE_UNITS), then scale threads appropriately
	cudaDeviceProp swatProp; int swatDev;
	CUDA_CHECK(cudaGetDevice(&swatDev));
	CUDA_CHECK(cudaGetDeviceProperties(&swatProp, swatDev));
	unsigned int devBlockNum = swatProp.multiProcessorCount;
	if (devBlockNum == MIN(blockNum, devBlockNum)) {
		printf("Scaling blocks from %d to %d to fit on device\n",\
				blockNum, devBlockNum);
		blockNum = devBlockNum;
	}
	mfThreadNum = blockNum * blockSize;
	int i;

	sprintf(queryFilePathName, "%s", argv[1]);
	sprintf(dbDataFilePathName, "%s.data", argv[2]);
	sprintf(dbLenFilePathName, "%s.loc", argv[2]);

	char *allSequences, *querySequence, *subSequence;
	char *seq1, *seq2;
	char *seq1D, *seq2D;

	//allSequences = new char[2 * (MAX_LEN)];
    posix_memalign((void **) &allSequences, AOCL_ALIGNMENT, sizeof(char *)*2*MAX_LEN);
    
	if (allSequences == NULL)
	{
		printf("Allocate sequence buffer error!\n");
		return 1;
	}
	querySequence = allSequences;

	CUDA_CHECK(cudaMalloc((void**)&seq1D, sizeof(char) * MAX_LEN));
	CUDA_CHECK(cudaMalloc((void**)&seq2D, sizeof(char) * MAX_LEN));

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

	char *outSeq1D, *outSeq2D;
	CUDA_CHECK(cudaMalloc((void**)&outSeq1D, sizeof(char) * MAX_LEN * 2));
	CUDA_CHECK(cudaMalloc((void**)&outSeq2D, sizeof(char) * MAX_LEN * 2));

	//allocate thread number per launch and 
	//location difference information
	int *threadNum, *diffPos;
	//threadNum = new int[2 * MAX_LEN];
	//diffPos = new int[2 * MAX_LEN];
    posix_memalign((void **) &diffPos, AOCL_ALIGNMENT, sizeof(int *)*2*MAX_LEN);
    posix_memalign((void **) &threadNum, AOCL_ALIGNMENT, sizeof(int *)*2*MAX_LEN);
	if (threadNum == NULL ||
			diffPos == NULL)
	{
		printf("Allocate location buffer on host error!\n");
		return 1;
	}

	int *threadNumD, *diffPosD;
	CUDA_CHECK(cudaMalloc((void**)&threadNumD, sizeof(int) * (2 * MAX_LEN)));
	CUDA_CHECK(cudaMalloc((void**)&diffPosD, sizeof(int) * (2 * MAX_LEN)));

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

	char *pathFlagD, *extFlagD;
	float *nGapDistD, *hGapDistD, *vGapDistD;
	CUDA_CHECK(cudaMalloc((void**)&pathFlagD, sizeof(char) * maxElemNum));
	CUDA_CHECK(cudaMalloc((void**)&extFlagD, sizeof(char) * maxElemNum));
	CUDA_CHECK(cudaMalloc((void**)&nGapDistD, sizeof(float) * maxElemNum));
	CUDA_CHECK(cudaMalloc((void**)&hGapDistD, sizeof(float) * maxElemNum));
	CUDA_CHECK(cudaMalloc((void**)&vGapDistD, sizeof(float) * maxElemNum));
    
	//Allocate the MAX INFO structure
	MAX_INFO *maxInfo;
	maxInfo = new MAX_INFO[1];
	if (maxInfo == NULL)
	{
		printf("Alloate maxInfo on host error!\n");
		return 1;
	}

	MAX_INFO *maxInfoD;
	CUDA_CHECK(cudaMalloc((void**)&maxInfoD, sizeof(MAX_INFO) * mfThreadNum));

	//allocate the distance table
	float *blosum62D;
	int nblosumHeight = 23;
	CUDA_CHECK(cudaMalloc((void**)&blosum62D, sizeof(float) * nblosumWidth * nblosumHeight));
    float *blosum62_1d = NULL; 
    posix_memalign((void **) &blosum62_1d, AOCL_ALIGNMENT, sizeof(float *)*23*23);
    for(int i =0; i< 23; i++) 
        for(int j=0;j<23;j++)
            blosum62_1d[23*i + j] = blosum62[i][j];

	CUDA_CHECK(cudaMemcpy(blosum62D, blosum62[0], nblosumWidth * nblosumHeight * sizeof(float), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());
	int *mutexMem;
	CUDA_CHECK(cudaMalloc((void**)&mutexMem, sizeof(int)));

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

	/* common timer: repetitions over the whole database scan; only the device work of each sequence (zeroing, uploads,
	 * matrix fill, trace-back, downloads) is accumulated into the timed region; file reads, preprocessing and printing are not. */
	dwarf_init("swat", "dynamic-programming", "cuda");
	dwarf_set_size(querySize, 0);
	int rep;
	for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	int quiet = (rep != dwarf_first_rep());   /* print the alignments once */
	fseek(pDBLenFile, 0, SEEK_SET); fseek(pDBDataFile, 0, SEEK_SET);
	fread(&subSequenceNum, sizeof(int), 1, pDBLenFile);
	dwarf_acc_reset();
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

		dwarf_acc_begin();
		//record time
		timerStart();

		//use a kernel to initialize the matrix (free-size: DWARF_LS threads/block, the kernel guards index<arraySize)
		arraySize = DPMatrixSize * sizeof(char);
		setZero<<<DWARF_NBLOCKS(arraySize), DWARF_LS>>>((char*)pathFlagD, arraySize);
		CUDA_SYNC_LAUNCH();
		setZero<<<DWARF_NBLOCKS(arraySize), DWARF_LS>>>((char*)extFlagD, arraySize);
		CUDA_SYNC_LAUNCH();

		arraySize = matrixIniNum * sizeof(float);
		setZero<<<DWARF_NBLOCKS(arraySize), DWARF_LS>>>((char*)nGapDistD, arraySize);
		CUDA_SYNC_LAUNCH();
		setZero<<<DWARF_NBLOCKS(arraySize), DWARF_LS>>>((char*)hGapDistD, arraySize);
		CUDA_SYNC_LAUNCH();
		setZero<<<DWARF_NBLOCKS(arraySize), DWARF_LS>>>((char*)vGapDistD, arraySize);
		CUDA_SYNC_LAUNCH();

		arraySize = sizeof(MAX_INFO) * mfThreadNum;
		setZero<<<DWARF_NBLOCKS(arraySize), DWARF_LS>>>((char*)maxInfoD, arraySize);
		CUDA_SYNC_LAUNCH();

		arraySize = sizeof(int);
		setZero<<<DWARF_NBLOCKS(arraySize), DWARF_LS>>>((char*)mutexMem, arraySize);
		CUDA_SYNC_LAUNCH();

		//copy input sequences to device
#define SWAT_COPY(dst, src, bytes, kind) do { CUDA_CHECK(cudaMemcpy(dst, src, bytes, kind)); CUDA_CHECK(cudaDeviceSynchronize()); } while (0)
		SWAT_COPY(seq1D, seq1, (rowNum - 1) * sizeof(char), cudaMemcpyHostToDevice);
		SWAT_COPY(seq2D, seq2, (columnNum - 1) * sizeof(char), cudaMemcpyHostToDevice);
		SWAT_COPY(diffPosD, diffPos, launchNum * sizeof(int), cudaMemcpyHostToDevice);
		SWAT_COPY(threadNumD, threadNum, launchNum * sizeof(int), cudaMemcpyHostToDevice);

		//record time
		timerEnd();
		strTime.copyTimeHostToDevice += elapsedTime();

		//record time
		timerStart();

		/* ALGORITHMIC shape: blockNum x blockSize = 14 x 64 resident workers (per-worker max tracking, lock-based global barrier) */
		MatchStringGPUSync<<<(unsigned int)(mfThreadNum / blockSize), (unsigned int)blockSize>>>(pathFlagD, extFlagD, nGapDistD, hGapDistD, vGapDistD,
				diffPosD, threadNumD, rowNum, columnNum, seq1D, seq2D, nblosumWidth, openPenalty, extensionPenalty, maxInfoD, blosum62D, mutexMem);
		CUDA_SYNC_LAUNCH();

		//record time
		timerEnd();
		strTime.matrixFillingTime += elapsedTime();

		//record time
		timerStart();
		trace_back2<<<1, 1>>>(pathFlagD, extFlagD, diffPosD, seq1D, seq2D, outSeq1D, outSeq2D, maxInfoD, (int)mfThreadNum);
		CUDA_SYNC_LAUNCH();
		//record time
		timerEnd();
		strTime.traceBackTime += elapsedTime();

		//record time
		timerStart();
		//copy matrix score structure back
		SWAT_COPY(maxInfo, maxInfoD, sizeof(MAX_INFO), cudaMemcpyDeviceToHost);

		int maxOutputLen = rowNum + columnNum - 2;
		SWAT_COPY(outSeq1, outSeq1D, maxOutputLen * sizeof(char), cudaMemcpyDeviceToHost);
		SWAT_COPY(outSeq2, outSeq2D, maxOutputLen * sizeof(char), cudaMemcpyDeviceToHost);
		//record time
		CUDA_CHECK(cudaDeviceSynchronize());
		dwarf_acc_end();
		gettimeofday(&t2, NULL);
		timerEnd();
		strTime.copyTimeDeviceToHost += elapsedTime();

		//call the print function to print the match result
		if(!quiet) {
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
	}
	dwarf_acc_commit(rep);
	}
	tmpTime = 1000.0 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000.0;
	/* [FIX-OCL-19] removed: appended to ../kernelTime.txt through an unchecked fopen (NULL dereference when the parent directory
	 * is not writable) and polluted the parent directory of whoever ran the program. */

	//print time
	printTime_toStandardOutput();
	float swat_best = maxInfo->fmaxscore;   /* read before maxInfo is freed below */
	printTime_toFile();

	fclose(pDBLenFile);
	fclose(pDBDataFile);


	free(allSequences);
	cudaFree(seq1D);
	cudaFree(seq2D);

	delete outSeq1;
	delete outSeq2;
	cudaFree(outSeq1D);
	cudaFree(outSeq2D);

	free(threadNum);
	cudaFree(threadNumD);
	free(diffPos);
	cudaFree(diffPosD);

	delete pathFlag;
	delete extFlag;
	delete nGapDist;
	delete hGapDist;
	delete vGapDist;
	cudaFree(pathFlagD);
	cudaFree(extFlagD);
	cudaFree(nGapDistD);
	cudaFree(hGapDistD);
	cudaFree(vGapDistD);

	delete maxInfo;
	cudaFree(maxInfoD);

	cudaFree(blosum62D);
	cudaFree(mutexMem);

	/* [FIX-OCL-20] validity check (none before): the best score of the last pair must be positive; alignment text is the
	 * reference for comparisons across implementations. */
	{
		int ok = swat_best > 0.0f;
		printf("%s: best score of the last pair %.1f\n", ok ? "SUCCESS" : "FAILURE", swat_best);
		return dwarf_report(ok);
	}
}

