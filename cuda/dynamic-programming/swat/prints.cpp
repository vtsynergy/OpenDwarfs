#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "functions.h"

void print_maxinfo(MAX_INFO *maxinfo, int pairNum)
{
	int i;
	for (i = 0; i < pairNum; i++)
	{
		printf("PairNo %d:\n", i);
		printf("\tnposi = %d\n", maxinfo[i].nposi);
		printf("\tnposj = %d\n", maxinfo[i].nposj);
		printf("\tnmaxpos = %d\n", maxinfo[i].nmaxpos);
		printf("\tfmaxscore = %.1f\n", maxinfo[i].fmaxscore);
		printf("\tnoutoutlen = %d\n", maxinfo[i].noutputlen);
	}
}

void print_seq(char *str, int seqSize)
{
	int i;
	printf("Sequence size = %d\n", seqSize);
	for (i = 0; i < seqSize; i++)
	{
		printf("%c", amino_acids[str[i]]);
	}
	printf("\n");
	
	return;
}

void PrintAlignment(char *outstr1,
                 	char *outstr2,
					int noutLen,
					int nCharsPerLine,
					float fopen_pa,
					float fext_pa)
{
	int i, j, npos;
	int nstart1, nstart2;
	int nstartmiss;
	char *interstring;
	int ninterpos;
	int nflag;

    interstring = new char[noutLen];
	if (interstring == NULL)
	{
		printf("Allocate interstring error!\n");
		return;
	}

    ninterpos = 0;
	nflag = 0;
	nstartmiss = 0;
	nstart1 = 0;
	nstart2 = 0;

	for (i = 0; i < noutLen; i++)
	{
		if (nflag == 0)
		{
			if (outstr1[noutLen - 1 - i] != outstr2[noutLen - 1 - i])
			{
				if (outstr1[noutLen - 1 - i] != 23) //23 is the code of '-'
					nstart1++;

				if (outstr2[noutLen - 1 - i] != 23) //23 is the code of '-'
					nstart2++;
				nstartmiss++;
				continue;
			}
			else
			{
				nflag = 1;
			}
		}

		if (outstr1[noutLen - 1 - i] == 23 || outstr2[noutLen - 1 - i] == 23)
		{
			interstring[ninterpos++] = ' ';
		}
		else if (outstr1[noutLen - 1 - i] != outstr2[noutLen - 1 - i])
		{
			interstring[ninterpos++] = '.';
		}
		else
		{
			interstring[ninterpos++] = '|';
		}
	}

	printf("Start position: Sequence1: %d\n", nstart1); 
	printf("                Sequence2: %d\n\n", nstart2);

	int nOutLine = (noutLen - nstartmiss)/nCharsPerLine;
	for (i = 0; i < nOutLine; i++)
	{
		for (j = 0; j < nCharsPerLine; j++)
		{
			npos = i * nCharsPerLine + nstartmiss + j;
			printf("%c", amino_acids[outstr1[noutLen - 1 - npos]]);
			//printf("%2d", outstr1[noutLen - 1 - npos]);
		}
		printf("\n");
		for (j = 0; j < nCharsPerLine; j++)
		{
			npos = i * nCharsPerLine + j;
			printf("%c", interstring[npos]);
		}
		printf("\n");
		for (j = 0; j < nCharsPerLine; j++)
		{
			npos = i * nCharsPerLine + nstartmiss + j;
			printf("%c", amino_acids[outstr2[noutLen - 1 - npos]]);
			//printf("%2d", outstr2[noutLen - 1 - npos]);
		}
		printf("\n\n");
	}

	int nRemainingCh = noutLen - nstartmiss - nOutLine * nCharsPerLine;
	if (nRemainingCh > 0)
	{
		for (j = 0; j < nRemainingCh; j++)
		{
			printf("%c", amino_acids[outstr1[nRemainingCh - 1 - j]]);
			//printf("%2d", outstr1[nRemainingCh - 1 - j]);
		}
		printf("\n");
		for (j = 0; j < nRemainingCh; j++)
		{
			npos = nOutLine * nCharsPerLine + j;
			printf("%c", interstring[npos]);
		}
		printf("\n");
		for (j = 0; j < nRemainingCh; j++)
		{
			printf("%c", amino_acids[outstr2[nRemainingCh - 1 - j]]);
			//printf("%2d", outstr2[nRemainingCh - 1 - j]);
		}
		printf("\n\n");
	}

	//calculate the match score
	if (1)
	{
		float fscore = 0.0f;
		int npreIsgap = 0;
		int indextmp1, indextmp2;
		for (i = 0; i < noutLen - nstartmiss; i++)
		{
			if (interstring[i] == '|' || interstring[i] == '.')
			{
				npos = i + nstartmiss;
				indextmp1 = outstr1[noutLen - 1 - npos];
				indextmp2 = outstr2[noutLen - 1 - npos];
				fscore = fscore + blosum62[indextmp1][indextmp2];
				npreIsgap = 0;
			}
			else
			{
				if (npreIsgap == 0)
				{
					fscore = fscore - fopen_pa;
				}
				else
				{
					fscore = fscore - fext_pa;
				}
				npreIsgap = 1;
			}
		}
		printf("In print result function, match score = %.2f\n", fscore);
		printf("output string length = %d\n", noutLen - nstartmiss);
	}

	return;
}


