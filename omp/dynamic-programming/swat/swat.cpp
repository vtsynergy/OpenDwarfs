/* Smith-Waterman (dynamic-programming dwarf) -- OpenMP target-offload translation of
 * opencl/dynamic-programming/swat/swat.cpp + kernels.cl (affine-gap local alignment of a query against
 * a sequence database; "verOpenCL4": GPU lock-based synchronisation).
 *
 * Same host flow per database sequence: preprocessing (anti-diagonal bookkeeping), zero the DP arrays,
 * upload the sequences / diagonal tables, matrix fill, single-work-item trace-back, download max info +
 * aligned strings, print the alignment. Same kernels (MatchStringGPUSync, trace_back2, setZero), same
 * data layout and the same per-thread maximum tracking (mfThreadNum = blockNum*blockSize = 14*64 workers,
 * strided over each anti-diagonal), so ties are resolved to the same cell as in the OpenCL code.
 *
 * The one structural change: the OpenCL kernel keeps all work-groups resident and synchronises the
 * anti-diagonals with a hand-written lock-based global barrier (atomic counter + spinning). OpenMP has no
 * cross-team barrier inside a target region, so the anti-diagonal loop moves to the host and each
 * anti-diagonal is one kernel launch -- the kernel boundary is the barrier. The per-diagonal state
 * (indexi1, indexj1, noffset, startPos) that the kernel kept in registers is carried on the host.
 *
 * CLI as the OpenCL app:  swat <query> <db> [openPenalty extensionPenalty blocks]
 * Timed region: the device work of every database sequence (zeroing, uploads, matrix fill, trace-back,
 * downloads); file reading, preprocessing and printing are excluded.
 * Not carried over: the append to ../kernelTime.txt (a debugging side effect of the OpenCL app).
 */
#include "global.h"
#include "functions.h"
#include "timeRec.h"
#include <omp.h>
#include "dwarf_timer.h"
#include "dwarf_launch.h" /* [FIX-OMP-13] pinned launch config */
/* [FIX-OMP-13] k_match is ALGORITHMIC: 14 teams x 64 threads (= OpenCL blockNum x blockSize); only on the GPU build */
#if defined(DWARF_USE_TARGET) && DWARF_USE_TARGET
#define SWAT_MATCH_LAUNCH(total) num_teams((total) / 64) thread_limit(64) num_threads(64)
#else
#define SWAT_MATCH_LAUNCH(total)
#endif

#define MIN(a, b) (a < b ? a : b)

/* ---- kernel setZero (byte-wise) ---- */
static void k_setzero(char *a, int arraySize, size_t cap) {
/* [FIX-OMP-13] free-size setZero: pinned 256 x ceil(n/256) */
#pragma omp target teams distribute parallel for DWARF_LAUNCH(arraySize) map(present, alloc : a[0 : cap])
  for (int index = 0; index < arraySize; index++) a[index] = 0;
}

/* ---- kernel MatchStringGPUSync: ONE anti-diagonal (launchNo) for all workers ---- */
static void k_match(char *pathFlag, char *extFlag, float *nGapDist, float *hGapDist, float *vGapDist, int *diffPos,
                    int *threadNum, char *seq1, char *seq2, int blosumWidth, float openPenalty, float extensionPenalty,
                    MAX_INFO *maxInfo, float *blosum62D, int launchNo, int indexi1, int indexj1, int startPos,
                    int totalThreadNum, size_t mElems, size_t tabLen) {
#pragma omp target teams distribute parallel for SWAT_MATCH_LAUNCH(totalThreadNum) map(present, alloc : pathFlag[0 : mElems], extFlag[0 : mElems], nGapDist[0 : mElems], hGapDist[0 : mElems], vGapDist[0 : mElems], diffPos[0 : tabLen], threadNum[0 : tabLen], seq1[0 : MAX_LEN], seq2[0 : MAX_LEN], maxInfo[0 : totalThreadNum], blosum62D[0 : 23 * 23])
  for (int threadid = 0; threadid < totalThreadNum; threadid++) {
    for (int tid = threadid; tid < threadNum[launchNo]; tid += totalThreadNum) {
      int indexi = indexi1 - tid;
      int indexj = indexj1 + tid;
      int npos = startPos + tid;
      int npreposhgap = npos - diffPos[launchNo];
      int npreposvgap = npreposhgap - 1;
      int npreposngap = npreposvgap - diffPos[launchNo - 1];
      int ntablepos = seq1[indexi] * blosumWidth + seq2[indexj];
      float fdist = blosum62D[ntablepos];
      float fmaxdist = nGapDist[npreposngap];
      float fdistngap = fmaxdist + fdist;
      float ext_dist = hGapDist[npreposhgap] - extensionPenalty;
      float fdisthgap = nGapDist[npreposhgap] - openPenalty;
      if (fdisthgap <= ext_dist) {
        fdisthgap = ext_dist;
        extFlag[npreposhgap] = 1;
      }
      ext_dist = vGapDist[npreposvgap] - extensionPenalty;
      float fdistvgap = nGapDist[npreposvgap] - openPenalty;
      if (fdistvgap <= ext_dist) {
        fdistvgap = ext_dist;
        pathFlag[npreposvgap] += 8;
      }
      fdistngap = (fdistngap < 0.0f) ? 0.0f : fdistngap;
      fdisthgap = (fdisthgap < 0.0f) ? 0.0f : fdisthgap;
      fdistvgap = (fdistvgap < 0.0f) ? 0.0f : fdistvgap;
      hGapDist[npos] = fdisthgap;
      vGapDist[npos] = fdistvgap;
      // priority 00, 01, 10
      if (fdistngap >= fdisthgap && fdistngap >= fdistvgap) {
        fmaxdist = fdistngap;
        pathFlag[npos] = 2;
      } else if (fdisthgap >= fdistngap && fdisthgap >= fdistvgap) {
        fmaxdist = fdisthgap;
        pathFlag[npos] = 1;
      } else {  // fdistvgap >= fdistngap && fdistvgap >= fdisthgap
        fmaxdist = fdistvgap;
        pathFlag[npos] = 3;
      }
      nGapDist[npos] = fmaxdist;
      // the maximum match distance is 0: previous alignment is useless
      if (fmaxdist <= 0.00000001f) pathFlag[npos] = PATH_END;
      if (maxInfo[threadid].fmaxscore < fmaxdist) {
        maxInfo[threadid].nposi = indexi + 1;
        maxInfo[threadid].nposj = indexj + 1;
        maxInfo[threadid].nmaxpos = npos;
        maxInfo[threadid].fmaxscore = fmaxdist;
      }
    }
  }
}

/* ---- kernel trace_back2: a single work-item ---- */
static void k_traceback(char *str_npathflagp, char *str_nExtFlagp, int *ndiffpos, char *instr1D, char *instr2D,
                        char *outstr1, char *outstr2, MAX_INFO *strMaxInfop, int mfThreadNum, size_t mElems, size_t tabLen) {
#pragma omp target map(present, alloc : str_npathflagp[0 : mElems], str_nExtFlagp[0 : mElems], ndiffpos[0 : tabLen], instr1D[0 : MAX_LEN], instr2D[0 : MAX_LEN], outstr1[0 : 2 * MAX_LEN], outstr2[0 : 2 * MAX_LEN], strMaxInfop[0 : mfThreadNum])
  {
    int i, j;
    int npos, maxPos, nlen;
    int npathflag;
    int nlaunchno;
    float maxScore;
    maxPos = 0;
    maxScore = strMaxInfop[0].fmaxscore;
    for (i = 1; i < mfThreadNum; i++) {
      if (maxScore < strMaxInfop[i].fmaxscore) {
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
    int done = 0;
    while (1) {
      if (npathflag == 3) {
        outstr1[nlen] = 23;
        outstr2[nlen] = instr2D[j - 1];
        nlen++;
        j--;
        npos = npos - ndiffpos[nlaunchno] - 1;
        nlaunchno--;
      } else if (npathflag == 1) {
        outstr1[nlen] = instr1D[i - 1];
        outstr2[nlen] = 23;
        nlen++;
        i--;
        npos = npos - ndiffpos[nlaunchno];
        nlaunchno--;
      } else if (npathflag == 2) {
        outstr1[nlen] = instr1D[i - 1];
        outstr2[nlen] = instr2D[j - 1];
        nlen++;
        i--;
        j--;
        npos = npos - ndiffpos[nlaunchno] - ndiffpos[nlaunchno - 1] - 1;
        nlaunchno = nlaunchno - 2;
      } else {
        done = 1;  // the OpenCL kernel returns here (without writing maxInfo[0].noutputlen)
        break;
      }
      // only if it is not an extension gap, will the path direction change; otherwise, back on the same direction
      int nExtFlag = str_npathflagp[npos] / 4;
      if (npathflag == 3 && (nExtFlag == 2 || nExtFlag == 3)) {
        npathflag = 3;
      } else if (npathflag == 1 && str_nExtFlagp[npos] == 1) {
        npathflag = 1;
      } else {
        npathflag = str_npathflagp[npos] & 0x3;
      }
      if (i == 0 || j == 0) break;
      if (npathflag == PATH_END) break;
    }
    if (!done) {
      i--;
      j--;
      while (i >= 0) { outstr1[nlen] = instr1D[i]; outstr2[nlen] = 23; nlen++; i--; }
      while (j >= 0) { outstr1[nlen] = 23; outstr2[nlen] = instr2D[j]; nlen++; j--; }
      strMaxInfop[0] = strMaxInfop[maxPos];
      strMaxInfop[0].noutputlen = nlen;
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 3) {
    printf("Calculate similarities between two strings.\n");
    printf("Maximum length of each string is: %d\n", MAX_LEN);
    printf("Usage: %s query database\n", argv[0]);
    printf("or: %s query database [openPenalty extensionPenalty block#]\n", argv[0]);
    printf("openPenalty (5.0), extensionPenalty (0.5)\n");
    return 1;
  }
  char queryFilePathName[255], dbDataFilePathName[255], dbLenFilePathName[255];
  int querySize, subSequenceNum, subSequenceSize;
  float openPenalty, extensionPenalty;
  int nblosumWidth = 23;
  size_t blockSize = 64;
  size_t mfThreadNum;
  int blockNum = 14;
  memset(&strTime, 0, sizeof(STRUCT_TIME));
  timerStart();

  openPenalty = 5.0f;
  extensionPenalty = 0.5;
  if (argc == 6) {
    openPenalty = atof(argv[3]);
    extensionPenalty = atof(argv[4]);
    blockNum = atoi(argv[5]);
  }
  mfThreadNum = blockNum * blockSize;

  sprintf(queryFilePathName, "%s", argv[1]);
  sprintf(dbDataFilePathName, "%s.data", argv[2]);
  sprintf(dbLenFilePathName, "%s.loc", argv[2]);

  char *allSequences = (char *)malloc(sizeof(char *) * 2 * MAX_LEN), *querySequence, *subSequence;
  char *seq1, *seq2;
  querySequence = allSequences;
  querySize = readQuerySequence(queryFilePathName, querySequence);
  if (querySize <= 0 || querySize > MAX_LEN) {
    printf("Query size %d is out of range (0, %d)\n", MAX_LEN, querySize);
    return 1;
  }
  encoding(querySequence, querySize);
  subSequence = allSequences + querySize;

  char *outSeq1 = new char[2 * MAX_LEN], *outSeq2 = new char[2 * MAX_LEN];
  int *threadNum = (int *)malloc(sizeof(int *) * 2 * MAX_LEN), *diffPos = (int *)malloc(sizeof(int *) * 2 * MAX_LEN);
  const size_t tabLen = 2 * MAX_LEN;
  int maxElemNum = (MAX_LEN + 1) * (MAX_LEN + 1);
  size_t mElems = (size_t)maxElemNum;
  char *pathFlag = (char *)calloc(mElems, 1), *extFlag = (char *)calloc(mElems, 1);
  float *nGapDist = (float *)calloc(mElems, sizeof(float)), *hGapDist = (float *)calloc(mElems, sizeof(float)),
        *vGapDist = (float *)calloc(mElems, sizeof(float));
  MAX_INFO *maxInfo = new MAX_INFO[mfThreadNum > 1 ? mfThreadNum : 1];
  char *seq1D = (char *)calloc(MAX_LEN, 1), *seq2D = (char *)calloc(MAX_LEN, 1);
  float *blosum62D = (float *)malloc(sizeof(float) * 23 * 23);
  for (int i = 0; i < 23; i++)
    for (int j = 0; j < 23; j++) blosum62D[23 * i + j] = blosum62[i][j];

  /* device buffers (the OpenCL clCreateBuffer calls): host arrays carry the mapping */
  int mf = (int)mfThreadNum;
  char *outSeq1D = outSeq1, *outSeq2D = outSeq2;
#pragma omp target enter data map(alloc : pathFlag[0 : mElems], extFlag[0 : mElems], nGapDist[0 : mElems], hGapDist[0 : mElems], vGapDist[0 : mElems], diffPos[0 : tabLen], threadNum[0 : tabLen], seq1D[0 : MAX_LEN], seq2D[0 : MAX_LEN], outSeq1D[0 : 2 * MAX_LEN], outSeq2D[0 : 2 * MAX_LEN], maxInfo[0 : mf]) map(to : blosum62D[0 : 23 * 23])

  pDBDataFile = fopen(dbDataFilePathName, "rb");
  if (pDBDataFile == NULL) { printf("DB data file %s open error!\n", dbDataFilePathName); return 1; }
  pDBLenFile = fopen(dbLenFilePathName, "rb");
  if (pDBLenFile == NULL) { printf("DB length file %s open error!\n", dbLenFilePathName); return 1; }
  timerEnd();
  strTime.iniTime = elapsedTime();

  dwarf_init("swat", "dynamic-programming", "omp");
  dwarf_set_size(querySize, 0);
  int rowNum, columnNum, matrixIniNum, DPMatrixSize, launchNum;
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    int quiet = (rep != dwarf_first_rep());  // print alignments once
    fseek(pDBLenFile, 0, SEEK_SET);
    fseek(pDBDataFile, 0, SEEK_SET);
    if (fread(&subSequenceNum, sizeof(int), 1, pDBLenFile) != 1) return 1;
    dwarf_acc_reset();
    for (int subSequenceNo = 0; subSequenceNo < subSequenceNum; subSequenceNo++) {
      timerStart();
      if (fread(&subSequenceSize, sizeof(int), 1, pDBLenFile) != 1) break;
      if (subSequenceSize <= 0 || subSequenceSize > MAX_LEN) {
        printf("Size %d of bubject sequence %d is out of range!\n", subSequenceSize, subSequenceNo);
        break;
      }
      if (fread(subSequence, sizeof(char), subSequenceSize, pDBDataFile) != (size_t)subSequenceSize) break;
      if (subSequenceSize > querySize) {
        seq1 = subSequence; seq2 = querySequence; rowNum = subSequenceSize + 1; columnNum = querySize + 1;
      } else {
        seq1 = querySequence; seq2 = subSequence; rowNum = querySize + 1; columnNum = subSequenceSize + 1;
      }
      launchNum = rowNum + columnNum - 1;
      DPMatrixSize = preProcessing(rowNum, columnNum, threadNum, diffPos, matrixIniNum);
      timerEnd();
      strTime.preprocessingTime += elapsedTime();

      /* ---- device work (timed): zero, upload, matrix fill, trace-back, download ---- */
      dwarf_acc_begin();
      timerStart();
      k_setzero(pathFlag, DPMatrixSize * (int)sizeof(char), mElems);
      k_setzero(extFlag, DPMatrixSize * (int)sizeof(char), mElems);
      k_setzero((char *)nGapDist, matrixIniNum * (int)sizeof(float), mElems * sizeof(float));
      k_setzero((char *)hGapDist, matrixIniNum * (int)sizeof(float), mElems * sizeof(float));
      k_setzero((char *)vGapDist, matrixIniNum * (int)sizeof(float), mElems * sizeof(float));
      k_setzero((char *)maxInfo, (int)(sizeof(MAX_INFO) * mfThreadNum), sizeof(MAX_INFO) * mfThreadNum);
      memcpy(seq1D, seq1, rowNum - 1);
      memcpy(seq2D, seq2, columnNum - 1);
#pragma omp target update to(seq1D[0 : rowNum - 1], seq2D[0 : columnNum - 1], diffPos[0 : launchNum], threadNum[0 : launchNum])
      timerEnd();
      strTime.copyTimeHostToDevice += elapsedTime();

      timerStart();
      {
        int indexi1 = -1, indexj1 = 0, noffset = 0, startPos = 2 * COALESCED_OFFSET;
        for (int launchNo = 2; launchNo < launchNum; launchNo++) {
          if (launchNo < rowNum) indexi1++;
          else if (launchNo == rowNum) { indexi1++; noffset = 1; }
          else indexj1++;
          k_match(pathFlag, extFlag, nGapDist, hGapDist, vGapDist, diffPos, threadNum, seq1D, seq2D, nblosumWidth,
                  openPenalty, extensionPenalty, maxInfo, blosum62D, launchNo, indexi1, indexj1, startPos, (int)mfThreadNum,
                  mElems, tabLen);
          startPos += diffPos[launchNo + 1] + noffset;
        }
      }
      timerEnd();
      strTime.matrixFillingTime += elapsedTime();

      timerStart();
      k_traceback(pathFlag, extFlag, diffPos, seq1D, seq2D, outSeq1D, outSeq2D, maxInfo, (int)mfThreadNum, mElems, tabLen);
      timerEnd();
      strTime.traceBackTime += elapsedTime();

      timerStart();
      int maxOutputLen = rowNum + columnNum - 2;
#pragma omp target update from(maxInfo[0 : 1], outSeq1D[0 : maxOutputLen], outSeq2D[0 : maxOutputLen])
      timerEnd();
      strTime.copyTimeDeviceToHost += elapsedTime();
      dwarf_acc_end();

      if (!quiet) {
        printf("============================================================\n");
        printf("Sequence pair %d:\n", subSequenceNo);
        int nlength = maxInfo[0].noutputlen;
        PrintAlignment(outSeq1, outSeq2, nlength, CHAR_PER_LINE, openPenalty, extensionPenalty);
        printf("Max alignment score (on device) is %.1f\n", maxInfo[0].fmaxscore);
        printf("openPenalty = %.1f, extensionPenalty = %.1f\n", openPenalty, extensionPenalty);
        printf("Input sequence size, querySize: %d, subSequenceSize: %d\n", querySize, subSequenceSize);
        printf("Max position, seq1 = %d, seq2 = %d\n", maxInfo[0].nposi, maxInfo[0].nposj);
      }
    }
    dwarf_acc_commit(rep);
  }
  printTime_toStandardOutput();
  fclose(pDBLenFile);
  fclose(pDBDataFile);
#pragma omp target exit data map(delete : pathFlag[0 : mElems], extFlag[0 : mElems], nGapDist[0 : mElems], hGapDist[0 : mElems], vGapDist[0 : mElems], diffPos[0 : tabLen], threadNum[0 : tabLen], seq1D[0 : MAX_LEN], seq2D[0 : MAX_LEN], outSeq1D[0 : 2 * MAX_LEN], outSeq2D[0 : 2 * MAX_LEN], maxInfo[0 : mf], blosum62D[0 : 23 * 23])

  /* validation: the best score must be positive for a non-empty homologous search; fidelity to the OpenCL
   * app is checked by comparing the printed alignments with the OpenCL run */
  int ok = maxInfo[0].fmaxscore > 0.0f;
  printf("%s: best score of the last pair %.1f\n", ok ? "SUCCESS" : "FAILURE", maxInfo[0].fmaxscore);
  return dwarf_report(ok);
}
