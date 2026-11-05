#include "global.h"
#include "functions.h"

int preProcessing(int rowNum,
		int columnNum,
		int *threadNum,
		int *diffPos,
		int& matrixIniElem)
{
	int launchNum = rowNum + columnNum - 1;
	int launchNo;
	int preColumnNum, tempDiffPos,lastEndFlagPos;
	int threadPerLaunch;
	int coalescedOffset = COALESCED_OFFSET;
	int startPos;
	int nOffset;

	nOffset = 0;
	startPos = coalescedOffset;
	threadPerLaunch = 1;
	preColumnNum = coalescedOffset - 1;

	for (launchNo = 1; launchNo < launchNum; launchNo++)
	{
		if (launchNo < rowNum)
		{
			if (launchNo < columnNum)
			{
				threadPerLaunch++;
				threadNum[launchNo] = threadPerLaunch - 2;
			}
			else
			{
				threadNum[launchNo] = threadPerLaunch - 1;
			}
			lastEndFlagPos = startPos;
		}
		else
		{
			threadPerLaunch--;
			threadNum[launchNo] = threadPerLaunch;
		}

		if (launchNo > rowNum)
		{
			nOffset = 1;
		}

		if ((rowNum == columnNum) && (launchNo == rowNum - 1))
		{
			tempDiffPos = ((threadPerLaunch + coalescedOffset - 2) / coalescedOffset) * coalescedOffset;
		}
		else
		{
			tempDiffPos = ((threadPerLaunch + coalescedOffset - 1) / coalescedOffset) * coalescedOffset;
		}
		startPos += tempDiffPos;

		if (launchNo <= rowNum)
		{
			diffPos[launchNo] = preColumnNum;
		}
		else
		{
			diffPos[launchNo] = preColumnNum - nOffset;
		}
		preColumnNum = tempDiffPos;
	}

	matrixIniElem = lastEndFlagPos + rowNum;

	return startPos;
}

float maxScore(float *scoreArray, int arraySize)
{
	float fMaxScore = 0.0f;
	int i;
	for (i = 0; i < arraySize; i++)
	{
		if (fMaxScore < scoreArray[i])
		{
			fMaxScore = scoreArray[i];
		}
	}

	return fMaxScore;
}
