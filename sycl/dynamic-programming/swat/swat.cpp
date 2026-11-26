#include "global.h"
#include "dwarf_sycl.hpp"
#include "functions.h"
#include "timeRec.h"
#define AOCL_ALIGNMENT 64

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


/* ---- kernels, translated from kernels.cl ---- */
#define PATH_END_K 0

/* __barrier_opencl_lock_based: group barrier, one atomic arrival per work-group, spin until all groups of this
 * wavefront arrived, group barrier. OpenCL 1.2 uses atom_add + mem_fence; here the same two atomics carry
 * acquire/release ordering at device scope. As in the OpenCL kernel, all work-groups must be resident at once. */
static inline void barrier_lock_based(sycl::nd_item<1> it, int localID, int goalValue, int *g_mutex)
{
	sycl::group_barrier(it.get_group());
	sycl::atomic_fence(sycl::memory_order::release, sycl::memory_scope::device);
	if (localID == 0)
	{
		sycl::atomic_ref<int, sycl::memory_order::acq_rel, sycl::memory_scope::device, sycl::access::address_space::global_space> m(*g_mutex);
		m.fetch_add(1);
		while (m.fetch_add(0) < goalValue)
		{ }
	}
	sycl::atomic_fence(sycl::memory_order::acquire, sycl::memory_scope::device);
	sycl::group_barrier(it.get_group());
}

static void MatchStringGPUSync(sycl::queue &q, size_t global, size_t local, char *pathFlag, char *extFlag, float *nGapDist, float *hGapDist,
		float *vGapDist, int *diffPos, int *threadNum, int rowNum, int columnNum, char *seq1, char *seq2, int blosumWidth,
		float openPenalty, float extensionPenalty, MAX_INFO *maxInfo, float *blosum62D, int *mutexMem)
{
	q.parallel_for<class swat_match>(sycl::nd_range<1>(global, local), [=](sycl::nd_item<1> it) {
		int npos, ntablepos, tid;
		int npreposngap, npreposhgap, npreposvgap;
		int nLocalID = it.get_local_id(0);
		int blockNum = it.get_group_range(0);
		int blockSize = it.get_local_range(0);
		int totalThreadNum = blockSize * blockNum;
		int threadid = it.get_global_id(0);

		int launchNo;
		int launchNum = rowNum + columnNum - 1;
		int indexi1 = -1;
		int indexj1 = 0;
		int indexi, indexj;
		int startPos = 2 * COALESCED_OFFSET;
		int noffset = 0;

		float fdist;
		float fdistngap, fdisthgap, fdistvgap;
		float ext_dist;
		float fmaxdist;

		//__local float temp[5000];
		//temp[0] = 0.0f;

		for (launchNo = 2; launchNo < launchNum; launchNo++)
		{
			if (launchNo < rowNum)
			{
				indexi1++;
			}
			else if (launchNo == rowNum)
			{
				indexi1++;
				noffset = 1;
			}
			else
			{
				indexj1++;
			}

			for (tid = threadid; tid < threadNum[launchNo]; tid += totalThreadNum)
			{
				indexi = indexi1 - tid;
				indexj = indexj1 + tid;

				npos = startPos + tid;

				npreposhgap = npos - diffPos[launchNo];
				npreposvgap = npreposhgap - 1;
				npreposngap = npreposvgap - diffPos[launchNo - 1];

				ntablepos = seq1[indexi] * blosumWidth + seq2[indexj];
				//ntablepos = seq1C[indexi] * blosumWidth + seq2C[indexj];
				fdist = blosum62D[ntablepos];

				fmaxdist = nGapDist[npreposngap];
				fdistngap = fmaxdist + fdist;

				ext_dist  = hGapDist[npreposhgap] - extensionPenalty;
				fdisthgap = nGapDist[npreposhgap] - openPenalty;

				if (fdisthgap <= ext_dist) 
				{
					fdisthgap = ext_dist;
					extFlag[npreposhgap] = 1;
				}

				ext_dist  = vGapDist[npreposvgap] - extensionPenalty;
				fdistvgap = nGapDist[npreposvgap] - openPenalty;

				if (fdistvgap <= ext_dist)
				{
					fdistvgap = ext_dist;
					pathFlag[npreposvgap] += 8;
				}

				fdistngap = (fdistngap < 0.0f) ? 0.0f : fdistngap;
				fdisthgap = (fdisthgap < 0.0f) ? 0.0f : fdisthgap;
				fdistvgap = (fdistvgap < 0.0f) ? 0.0f : fdistvgap;

				hGapDist[npos] = fdisthgap;
				vGapDist[npos] = fdistvgap;

				//priority 00, 01, 10
				if (fdistngap >= fdisthgap && fdistngap >= fdistvgap)
				{
					fmaxdist = fdistngap;
					pathFlag[npos] = 2;
				}
				else if (fdisthgap >= fdistngap && fdisthgap >= fdistvgap)
				{
					fmaxdist = fdisthgap;
					pathFlag[npos] = 1;
				}
				else //fdistvgap >= fdistngap && fdistvgap >= fdisthgap
				{
					fmaxdist = fdistvgap;
					pathFlag[npos] = 3;
				}

				nGapDist[npos] = fmaxdist;

				//Here, the maximum match distance is 0, which means
				//previous alignment is useless
				if (fmaxdist <= 0.00000001f)
				{
					pathFlag[npos] = PATH_END_K;
				}

				if (maxInfo[threadid].fmaxscore < fmaxdist)
				{
					maxInfo[threadid].nposi = indexi + 1;
					maxInfo[threadid].nposj = indexj + 1;
					maxInfo[threadid].nmaxpos = npos;
					maxInfo[threadid].fmaxscore = fmaxdist;
				}
			}

			//GPU synchronization
			barrier_lock_based(it, nLocalID, (launchNo - 1) * blockNum, mutexMem);
			startPos += diffPos[launchNo + 1] + noffset;
		}
}).wait_and_throw();
}

static void trace_back2(sycl::queue &q, char *str_npathflagp, char *str_nExtFlagp, int *ndiffpos, char *instr1D, char *instr2D,
		char *outstr1, char *outstr2, MAX_INFO *strMaxInfop, int mfThreadNum)
{
	q.parallel_for<class swat_traceback>(sycl::nd_range<1>(1, 1), [=](sycl::nd_item<1>) {
		int i, j;
		int npos, maxPos, nlen;
		int npathflag;
		int nlaunchno;
		float maxScore;

		maxPos = 0;
		maxScore = strMaxInfop[0].fmaxscore;
		for (i = 1; i < mfThreadNum; i++)
		{
			if (maxScore < strMaxInfop[i].fmaxscore)
			{
				maxPos = i;
				maxScore = strMaxInfop[i].fmaxscore;
			}
		}

		npos = strMaxInfop[maxPos].nmaxpos;
		npathflag = str_npathflagp[npos] & 0x3;
		nlen = 0;

		i = strMaxInfop[maxPos].nposi;
		j = strMaxInfop[maxPos].nposj;
		nlaunchno = i + j;

		while (1)
		{
			if (npathflag == 3)
			{
				outstr1[nlen] = 23;
				outstr2[nlen] = instr2D[j - 1];
				nlen++;
				j--;

				//position in the transformed matrix
				npos = npos - ndiffpos[nlaunchno] - 1;
				nlaunchno--;
			}
			else if (npathflag == 1)
			{
				outstr1[nlen] = instr1D[i - 1];
				outstr2[nlen] = 23;
				nlen++;
				i--;

				//position in the transformed matrix
				npos = npos - ndiffpos[nlaunchno];
				nlaunchno--;
			}
			else if (npathflag == 2)
			{
				outstr1[nlen] = instr1D[i - 1];
				outstr2[nlen] = instr2D[j - 1];
				nlen++;
				i--;
				j--;

				//position in the transformed matrix
				npos = npos - ndiffpos[nlaunchno] - ndiffpos[nlaunchno - 1] - 1;
				nlaunchno = nlaunchno - 2;
			}
			else
			{
				//printf("npathflag = %d, npos = %d\n", npathflag, npos);
				//printf("find path error!\n");
				return;
			}

			//only if it is not an extension gap, will the path direction change
			//otherwise, back on the same direction
			int nExtFlag = str_npathflagp[npos] / 4;
			if (npathflag == 3 && (nExtFlag == 2 || nExtFlag == 3))
			{
				npathflag = 3;
			}
			//else if (npathflag == 1 && (nExtFlag == 1 || nExtFlag == 3))
			else if (npathflag == 1 && str_nExtFlagp[npos] == 1)
			{
				npathflag = 1;
			}
			else
			{
				npathflag = str_npathflagp[npos] & 0x3;
			}

			if (i == 0 || j == 0)
			{
				break;
			}

			if (npathflag == PATH_END_K)
			{
				break;
			}
		}

		i--;
		j--;

		while(i >= 0)
		{
			outstr1[nlen] = instr1D[i];
			outstr2[nlen] = 23;
			nlen++;
			i--;
		}

		while(j >= 0)
		{
			outstr1[nlen] = 23;
			outstr2[nlen] = instr2D[j];
			nlen++; 
			j--;
		}

		strMaxInfop[0] = strMaxInfop[maxPos];
		strMaxInfop[0].noutputlen = nlen;

		return;
}).wait_and_throw();
}

static void setZero(sycl::queue &q, void *a_, int arraySize)
{
	char *a = (char *)a_;
	q.parallel_for<class swat_setzero>(sycl::nd_range<1>(DWARF_ROUNDUP(arraySize), DWARF_LS), [=](sycl::nd_item<1> it) {
		unsigned int index = it.get_global_id(0);
		if (index < arraySize)
		{
			a[index] = 0;
		}
}).wait_and_throw();
}

int main(int argc, char ** argv)
{
	dsycl::parse_args(&argc, &argv);

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

	sycl::queue q = dsycl::make_queue();

	//check to make sure the device supports this block count
	//then scale threads appropriately
	unsigned int devBlockNum = q.get_device().get_info<sycl::info::device::max_compute_units>();
	if (devBlockNum == MIN(blockNum, devBlockNum)) {
		printf("Scaling blocks from %d to %d to fit on device\n",\
				blockNum, devBlockNum);
		blockNum = devBlockNum;
	}
	mfThreadNum = blockNum * blockSize;

	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */
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

	seq1D = dsycl::alloc<char>(q, MAX_LEN);
	seq2D = dsycl::alloc<char>(q, MAX_LEN);

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

	char *outSeq1D = dsycl::alloc<char>(q, MAX_LEN * 2), *outSeq2D = dsycl::alloc<char>(q, MAX_LEN * 2);

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

	int *threadNumD = dsycl::alloc<int>(q, 2 * MAX_LEN), *diffPosD = dsycl::alloc<int>(q, 2 * MAX_LEN);

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

	char *pathFlagD = dsycl::alloc<char>(q, maxElemNum), *extFlagD = dsycl::alloc<char>(q, maxElemNum);
	float *nGapDistD = dsycl::alloc<float>(q, maxElemNum), *hGapDistD = dsycl::alloc<float>(q, maxElemNum), *vGapDistD = dsycl::alloc<float>(q, maxElemNum);

	//Allocate the MAX INFO structure
	MAX_INFO *maxInfo;
	maxInfo = new MAX_INFO[1];
	if (maxInfo == NULL)
	{
		printf("Alloate maxInfo on host error!\n");
		return 1;
	}

	MAX_INFO *maxInfoD = dsycl::alloc<MAX_INFO>(q, mfThreadNum);

	//allocate the distance table
	int nblosumHeight = 23;
	float *blosum62D = dsycl::alloc<float>(q, nblosumWidth * nblosumHeight);
	dsycl::write(q, blosum62D, blosum62[0], nblosumWidth * nblosumHeight * sizeof(float));
	int *mutexMem = dsycl::alloc<int>(q, 1);

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
	dwarf_init("swat", "dynamic-programming", "sycl");
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

		//use a kernel to initialize the matrix ([FIX-OCL-34] free-size setZero: global=roundup(n), local=DWARF_LS)
		arraySize = DPMatrixSize * sizeof(char);
		setZero(q, pathFlagD, arraySize);
		setZero(q, extFlagD, arraySize);
		arraySize = matrixIniNum * sizeof(float);
		setZero(q, nGapDistD, arraySize);
		setZero(q, hGapDistD, arraySize);
		setZero(q, vGapDistD, arraySize);
		arraySize = sizeof(MAX_INFO) * mfThreadNum;
		setZero(q, maxInfoD, arraySize);
		arraySize = sizeof(int);
		setZero(q, mutexMem, arraySize);

		//copy input sequences to device
		dsycl::write(q, seq1D, seq1, (rowNum - 1) * sizeof(char));
		dsycl::write(q, seq2D, seq2, (columnNum - 1) * sizeof(char));
		dsycl::write(q, diffPosD, diffPos, launchNum * sizeof(int));
		dsycl::write(q, threadNumD, threadNum, launchNum * sizeof(int));

		//record time
		timerEnd();
		strTime.copyTimeHostToDevice += elapsedTime();

		//record time
		timerStart();

		/* [FIX-OCL-34] ALGORITHMIC: MatchStringGPUSync keeps blockNum x blockSize = 14 x 64 resident workers (per-worker max tracking, lock-based global barrier) */
		MatchStringGPUSync(q, mfThreadNum, blockSize, pathFlagD, extFlagD, nGapDistD, hGapDistD, vGapDistD, diffPosD, threadNumD, rowNum, columnNum,
				seq1D, seq2D, nblosumWidth, openPenalty, extensionPenalty, maxInfoD, blosum62D, mutexMem);

		//record time
		timerEnd();
		strTime.matrixFillingTime += elapsedTime();

		//record time
		timerStart();
		trace_back2(q, pathFlagD, extFlagD, diffPosD, seq1D, seq2D, outSeq1D, outSeq2D, maxInfoD, mfThreadNum);
		//record time
		timerEnd();
		strTime.traceBackTime += elapsedTime();

		//record time
		timerStart();
		//copy matrix score structure back
		dsycl::read(q, maxInfo, maxInfoD, sizeof(MAX_INFO));
		int maxOutputLen = rowNum + columnNum - 2;
		dsycl::read(q, outSeq1, outSeq1D, maxOutputLen * sizeof(char));
		dsycl::read(q, outSeq2, outSeq2D, maxOutputLen * sizeof(char));
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


	delete allSequences;
	dsycl::release(q, seq1D);
	dsycl::release(q, seq2D);

	delete outSeq1;
	delete outSeq2;
	dsycl::release(q, outSeq1D);
	dsycl::release(q, outSeq2D);

	delete threadNum;
	dsycl::release(q, threadNumD);
	delete diffPos;
	dsycl::release(q, diffPosD);

	delete pathFlag;
	delete extFlag;
	delete nGapDist;
	delete hGapDist;
	delete vGapDist;
	dsycl::release(q, pathFlagD);
	dsycl::release(q, extFlagD);
	dsycl::release(q, nGapDistD);
	dsycl::release(q, hGapDistD);
	dsycl::release(q, vGapDistD);

	delete maxInfo;
	dsycl::release(q, maxInfoD);

	dsycl::release(q, blosum62D);
	dsycl::release(q, mutexMem);

	/* [FIX-OCL-20] validity check (none before): the best score of the last pair must be positive; alignment text is the
	 * reference for comparisons across implementations. */
	{
		int ok = swat_best > 0.0f;
		printf("%s: best score of the last pair %.1f\n", ok ? "SUCCESS" : "FAILURE", swat_best);
		return dwarf_report(ok);
	}
}

