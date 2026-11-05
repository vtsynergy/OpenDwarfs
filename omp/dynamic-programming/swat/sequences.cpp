#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "functions.h"

FILE *pDBLenFile = NULL;
FILE *pDBDataFile = NULL;

int readQuerySequence(char *queryFileName, char *querySequence)
{
	int i;
	char *strTemp;
	FILE *pQueryFile;
	int nSeqLen;

	strTemp = new char[10000];
	if (strTemp == NULL)
	{
		printf("In function readQuerySequence: allocate strTemp error!\n");
		return -1;
	}

	pQueryFile = fopen(queryFileName, "rt");
	if (pQueryFile == NULL)
	{
		printf("In function readQuerySequence: query file %s open error!\n", queryFileName);
		return -1;
	}

	nSeqLen = 0;
	while (!feof(pQueryFile))
	{
		*strTemp = '\0';
		fgets(strTemp, 10000, pQueryFile);
		if (*strTemp == '>')
		{
			continue;
		}

		memcpy(querySequence + nSeqLen, strTemp, strlen(strTemp));
		nSeqLen += (strlen(strTemp));
	}

	fclose(pQueryFile);
	delete strTemp;

	return nSeqLen;
}

void encoding(char *seq, int& nsize)
{
	int i;
	int nOutLoc = 0;
	char code;
	for (i = 0; i < nsize; i++)
	{
		code = char2index(seq[i]);
		if (code >= 0)
		{
			seq[nOutLoc] = code;
			nOutLoc++;
		}
	}

	nsize = nOutLoc;
	return;
}

short char2index(char inch)
{
	int result;
	if(inch >= 65 && inch <= 73) //'A' --> 'I'
	{
		result = inch - 65;
	}
	else if (inch >= 75 && inch <= 78) //'K' --> 'N'
	{
		result = inch - 66;
	}
	else if (inch >= 80 && inch <= 84) //'P' --> 'T'
	{
		result = inch - 67;
	}
	else if (inch >= 86 && inch <= 90) //'V' --> 'Z'
	{
		result = inch - 68;
	}
	else if (inch >= 97 && inch <= 105)
	{
		result = inch - 97;
	}
	else if (inch >= 107 && inch <= 110)
	{
		result = inch - 98;
	}
	else if (inch >= 112 && inch <= 116)
	{
		result = inch - 99;
	}
	else if (inch >= 118 && inch <= 122)
	{
		result = inch - 100;
	}
	else
	{
		return -1;
	}

	return result;
}


