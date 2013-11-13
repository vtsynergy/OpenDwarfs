#ifndef __FINE_GRAINDED_GLOBAL_H__
#define __FINE_GRAINDED_GLOBAL_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define MAX_LEN 5200
#define CHAR_PER_LINE 50
#define PATH_END 0 //when PATH_END module 4, the result should be 3
#define COALESCED_OFFSET 32

typedef struct {
	int npathflag;
    float fngapdist;
	float fhgapdist;
	float fvgapdist;
} match_info;

typedef struct {
	int nposi, nposj;
	int nmaxpos;
	float fmaxscore;
	int noutputlen;
}   MAX_INFO;

extern FILE *pDBLenFile;
extern FILE *pDBDataFile;
extern float blosum62[23][23];
extern char amino_acids[24];

#endif

