#ifndef __FINE_GRAINED_FUNCTIONS_H__
#define __FINE_GRAINED_FUNCTIONS_H__
#include "global.h"

void trace_back(int*, int*, char*, char*, char*, char*, MAX_INFO*);

void print_maxinfo(MAX_INFO *maxinfo, int pairNum);
void print_seq(char *str, int seqSize);
void PrintAlignment(char *, char *, int, int, float, float);

int preProcessing(int rowNum,
				   int columnNum,
				   int *threadNum,
				   int *diffPos,
				   int& matrixIniElem);
float maxScore(float *scoreArray, int arraySize);

int readQuerySequence(char*, char*);
void encoding(char *, int &);
short char2index(char);

void copyScoringMatrixToConstant();

#endif
