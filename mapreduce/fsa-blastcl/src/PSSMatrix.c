// PSSMatrix.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for creating a Position Specific Substitution Matrix that contains,
// for each residue in the query (each position) a score for aligning that
// position against each possible character

#include "blast.h"

unsigned char* PSSMatrix_packedRightMatches = NULL;
unsigned char* PSSMatrix_packedLeftMatches = NULL;
int2* PSSMatrix_packedRightMatchScores = NULL;
int2* PSSMatrix_packedLeftMatchScores = NULL;
char* PSSMatrix_packedScore = NULL;

// Create a PSSM for the given query sequence and score matrix.
// The PSSMatrix will have length(query) columns and 25 rows
struct PSSMatrix PSSMatrix_create(struct scoreMatrix scoreMatrix, char* query)
{
	struct PSSMatrix PSSMatrix;
	int4 columnCount;
	unsigned char code, code2, bestCode;
    int4 packedByte, matches;

    if (encoding_alphabetType == encoding_nucleotide && parameters_strands == 3)
	{
    	PSSMatrix.length = 2 * strlen(query);
        PSSMatrix.strandLength = strlen(query);
    }
    else if (encoding_alphabetType == encoding_nucleotide && parameters_strands == 2)
	{
    	PSSMatrix.length = strlen(query);
        PSSMatrix.strandLength = 0;
    }
    else
    {
		PSSMatrix.strandLength = PSSMatrix.length = strlen(query);
	}

	// For each character in the query sequence
	columnCount = 0;
	while (columnCount < strlen(query))
	{
    		// Convert the character to upper case
		query[columnCount] = toupper(query[columnCount]);	
		columnCount++;
	}
	
	// Declare memory for use by PSSMatrix
	PSSMatrix.matrix = (int2**)global_malloc(sizeof(int2*) * (PSSMatrix.length + 2));

    // Declare memory to store encoded query
    PSSMatrix.queryCodes = (unsigned char*)global_malloc(sizeof(unsigned char) * (PSSMatrix.length));

    // Declare memory to store the highest scoring regular letter for each position
    PSSMatrix.bestMatchCodes = (unsigned char*)global_malloc(sizeof(unsigned char) * (PSSMatrix.length));

    // The first column will be a sentinal column
	PSSMatrix.matrix++;

	// Set initial maximum, minimum values
	PSSMatrix.highestValue = scoreMatrix.highestValue;
	PSSMatrix.lowestValue = scoreMatrix.lowestValue;

	// For each character in the query sequence
	columnCount = 0;
	while (columnCount < PSSMatrix.length)
	{
		// Get its code
        if (columnCount < PSSMatrix.strandLength)
        {
        	// Positive strand
			code = encoding_getCode(query[columnCount]);
		}
        else
        {
        	// Negative strand
			code = encoding_getComplement(encoding_getCode(query[PSSMatrix.length - columnCount - 1]));
        }

		// Use column from scoreMatrix for PSSMatrix
		PSSMatrix.matrix[columnCount] = scoreMatrix.matrix[code];

        // Add code to encoded query
        PSSMatrix.queryCodes[columnCount] = code;

        // For each possible regular letter
		bestCode = 0;
        code2 = 0;
		while (code2 < encoding_numRegularLetters)
        {
        	// Determine which scores highest at this position
            // (typically the same as the query letter for this position)
            if (PSSMatrix.matrix[columnCount][code2] > PSSMatrix.matrix[columnCount][bestCode])
            {
            	bestCode = code2;
            }
            code2++;
        }

        // Set best code at this position
		if (encoding_alphabetType == encoding_nucleotide && code >= encoding_numRegularLetters)
        {
        	PSSMatrix.bestMatchCodes[columnCount] = encoding_randomEncodedLetter(code);
        }
        else
        {
            PSSMatrix.bestMatchCodes[columnCount] = bestCode;
        }

		columnCount++;
	}

    // If this is a nucleotide sequence store byte-packed code information
    if (encoding_alphabetType == encoding_nucleotide)
    {
        // Declare memory to store xor codes; each 2-bit code copied 4 times in the byte
        PSSMatrix.xorCodes = (unsigned char*)global_malloc(sizeof(unsigned char) * (PSSMatrix.length));

        // For each nucleotide character in the query sequence
        columnCount = 0;
        while (columnCount < PSSMatrix.length)
        {
        	code = PSSMatrix.queryCodes[columnCount];

            // Set xorCode which is the 2-bit code copied 4 times in the byte
            PSSMatrix.xorCodes[columnCount] = (code << 6) | code << 4 | code << 2 | code;

            columnCount++;
        }

        // Declare memory to store byte-packed query
        PSSMatrix.bytePackedCodes = (unsigned char*)global_malloc(sizeof(unsigned char) * (PSSMatrix.length + 5)) + 4;

        // For each position in the sequence
        columnCount = 0;
        while (columnCount < PSSMatrix.length)
        {
        	// Store packed version of next 4 characters in query
        	if (columnCount + 3 < PSSMatrix.length)
        		PSSMatrix.bytePackedCodes[columnCount]
                	= encoding_bytePack(PSSMatrix.bestMatchCodes + columnCount);
            // Of the remaining < 4 characters
            else
				PSSMatrix.bytePackedCodes[columnCount] = encoding_bytePackRemaining(
                	PSSMatrix.bestMatchCodes + columnCount, PSSMatrix.length - columnCount);

            columnCount++;
		}

        // Store DUMMY packed value before and after sequence
		PSSMatrix.bytePackedCodes[-4] = 0;
		PSSMatrix.bytePackedCodes[PSSMatrix.length] = 0;

		// Store three before-sequence codes
        PSSMatrix.bytePackedCodes[-3] = encoding_bytePackBeginning(PSSMatrix.bestMatchCodes, 1);
        PSSMatrix.bytePackedCodes[-2] = encoding_bytePackBeginning(PSSMatrix.bestMatchCodes, 2);
        PSSMatrix.bytePackedCodes[-1] = encoding_bytePackBeginning(PSSMatrix.bestMatchCodes, 3);

		PSSMatrix_packedRightMatches = (unsigned char*)global_malloc(sizeof(unsigned char) * 256);
		PSSMatrix_packedLeftMatches = (unsigned char*)global_malloc(sizeof(unsigned char) * 256);
		PSSMatrix_packedRightMatchScores = (int2*)global_malloc(sizeof(int2) * 256);
		PSSMatrix_packedLeftMatchScores = (int2*)global_malloc(sizeof(int2) * 256);
		PSSMatrix_packedScore = (unsigned char*)global_malloc(sizeof(unsigned char) * 256);

        // For each packed byte
        packedByte = 0;
        while (packedByte < 256)
        {
            // For each possible XORed value determine number of matches going
            // from left to RIGHT
            matches = 0;
            if (!((packedByte >> 6) & 0x3))
            {
            	matches++;
                if (!((packedByte >> 4) & 0x3))
                {
                	matches++;
                    if (!((packedByte >> 2) & 0x3))
                    {
                        matches++;
                        if (!(packedByte & 0x3))
                        {
                            matches++;
                        }
                    }
				}
            }
            PSSMatrix_packedRightMatches[packedByte] = matches;
            PSSMatrix_packedRightMatchScores[packedByte] = matches * parameters_matchScore;

            // The same going from right to LEFT
            matches = 0;
            if (!(packedByte & 0x3))
            {
            	matches++;
                if (!((packedByte >> 2) & 0x3))
                {
                	matches++;
                    if (!((packedByte >> 4) & 0x3))
                    {
                        matches++;
                        if (!((packedByte >> 6) & 0x3))
                        {
                            matches++;
                        }
                    }
				}
            }
            PSSMatrix_packedLeftMatches[packedByte] = matches;
            PSSMatrix_packedLeftMatchScores[packedByte] = matches * parameters_matchScore;

            // Determine the score for this XORed pair
            matches = 0;
            if (!((packedByte >> 6) & 0x3))
	           	matches += parameters_matchScore;
			else
            	matches += parameters_mismatchScore;

            if (!((packedByte >> 4) & 0x3))
	           	matches += parameters_matchScore;
			else
            	matches += parameters_mismatchScore;

            if (!((packedByte >> 2) & 0x3))
	           	matches += parameters_matchScore;
			else
            	matches += parameters_mismatchScore;

            if (!(packedByte & 0x3))
	           	matches += parameters_matchScore;
			else
            	matches += parameters_mismatchScore;

            PSSMatrix_packedScore[packedByte] = matches;

            packedByte++;
        }
	}

    // Make columns flanking the query sentinal columns
	PSSMatrix.matrix[-1] = scoreMatrix.matrix[encoding_sentinalCode];
	PSSMatrix.matrix[PSSMatrix.length] = scoreMatrix.matrix[encoding_sentinalCode];

	return PSSMatrix;
}

//Transform the PSSMatrix to a layout compatible with GPU
struct PSSMatrixFP PSSMatrixFP_transform(struct PSSMatrix *PSSMatrixptr)
{
	struct PSSMatrixFP PSSMatrixFP;
	int columnCount;
	int i, j;
	int offset;
	int2 *FP;
	int2 **Orig;

	PSSMatrixFP.length = PSSMatrixptr->length;
	PSSMatrixFP.strandLength = PSSMatrixptr->strandLength;
	PSSMatrixFP.highestValue = PSSMatrixptr->highestValue;
	PSSMatrixFP.lowestValue = PSSMatrixptr->lowestValue;

	// Declare memory for use by PSSMatrix
	PSSMatrixFP.matrix = (cl_int2 *)global_malloc(sizeof(cl_int2) * (PSSMatrixFP.length + 2) * encoding_numCodes);
	// Copy matrix from PSSMatrix to PSSMatrixFP
	columnCount = 0;
	memcpy((char *)PSSMatrixFP.matrix, (char *)PSSMatrixptr->matrix[-1], sizeof(int2) * encoding_numCodes); 
	PSSMatrixFP.matrix += encoding_numCodes;
	while (columnCount < PSSMatrixFP.length)
	{
		memcpy((char *)&PSSMatrixFP.matrix[columnCount * encoding_numCodes], (char *)PSSMatrixptr->matrix[columnCount], sizeof(int2) * encoding_numCodes);
		columnCount++;
	}

	memcpy((char *)&PSSMatrixFP.matrix[PSSMatrixFP.length * encoding_numCodes], (char *)PSSMatrixptr->matrix[PSSMatrixFP.length], sizeof(int2) * encoding_numCodes);

	// Declare memory to store encoded query
	PSSMatrixFP.queryCodes = (unsigned char *)global_malloc(sizeof(unsigned char) * (PSSMatrixFP.length));
	// Copy queryCodes from PSSMatrix to PSSMatrixFP
	memcpy(PSSMatrixFP.queryCodes, PSSMatrixptr->queryCodes, sizeof(unsigned char) * (PSSMatrixFP.length));

	// Declare memory to store the highest scoring regular letter for each position
	PSSMatrixFP.bestMatchCodes = (unsigned char*)global_malloc(sizeof(unsigned char) * (PSSMatrixFP.length));
	//Copy bestMatchCode from PSSMatrix to PSSMatrixFP
	memcpy(PSSMatrixFP.bestMatchCodes, PSSMatrixptr->bestMatchCodes, sizeof(unsigned char) * (PSSMatrixFP.length));


	return PSSMatrixFP;

	offset = 0;
	Orig = PSSMatrixptr->matrix + offset;
	FP = PSSMatrixFP.matrix + offset * encoding_numCodes;
	
	for (j = 0; j < PSSMatrixFP.length; j++)
	{
		for (i = 0; i < encoding_numCodes + 10; i++)
		{
			if ((*Orig)[i] != FP[i])
			{
				printf("i = %d, %d %d\n", i, (*Orig)[i], FP[i]);
			}
		}
	
		Orig++;
		FP += encoding_numCodes;
		printf("---------------------------------\n");
	}


	return PSSMatrixFP;
}

// Calculate the start of strand for the given query offset
uint4 PSSMatrix_strandStart(struct PSSMatrix PSSMatrix, uint4 queryOffset)
{
	if (queryOffset > PSSMatrix.strandLength)
    	return PSSMatrix.strandLength;
    else
    	return 0;
}

// Calculate the end of strand for the given query offset
uint4 PSSMatrix_strandEnd(struct PSSMatrix PSSMatrix, uint4 queryOffset)
{
	if (queryOffset > PSSMatrix.strandLength)
    	return PSSMatrix.length;
    else
    	return PSSMatrix.strandLength;
}

// Returns a PSSMatrix with the first "amount" entries removed and the length shortened
// accordingly
struct PSSMatrix PSSMatrix_chop(struct PSSMatrix PSSMatrix, int4 amount)
{
	struct PSSMatrix chopped;

	chopped.matrix = PSSMatrix.matrix + amount;
	chopped.queryCodes = PSSMatrix.queryCodes + amount;
	chopped.bytePackedCodes = PSSMatrix.bytePackedCodes + amount;
    chopped.xorCodes = PSSMatrix.xorCodes + amount;
	chopped.length = PSSMatrix.length - amount;
	chopped.highestValue = PSSMatrix.highestValue;
	chopped.lowestValue = PSSMatrix.lowestValue;
    chopped.strandLength = PSSMatrix.strandLength - amount;

    if (chopped.strandLength < 0)
    	chopped.strandLength = 0;

	return chopped;
}

// Returns a PSSMatrix which is a reversed copy
struct PSSMatrix PSSMatrix_reverse(struct PSSMatrix PSSMatrix)
{
	struct PSSMatrix reversed;
	uint4 queryPosition = 0;

	reversed.matrix = (int2**)global_malloc(sizeof(int2*) * (PSSMatrix.length + 2));
    reversed.queryCodes = (unsigned char*)global_malloc(sizeof(unsigned char) * (PSSMatrix.length));
	reversed.length = PSSMatrix.length;
	reversed.strandLength = PSSMatrix.strandLength;
	reversed.highestValue = PSSMatrix.highestValue;
	reversed.lowestValue = PSSMatrix.lowestValue;
    reversed.bestMatchCodes = NULL;
    reversed.bytePackedCodes = NULL;
    reversed.xorCodes = NULL;

    // Make and reverse matrix and query codes
    while (queryPosition < reversed.length)
    {
    	reversed.matrix[queryPosition]
        	= PSSMatrix.matrix[reversed.length - queryPosition - 1];
    	reversed.queryCodes[queryPosition]
        	= PSSMatrix.queryCodes[reversed.length - queryPosition - 1];
        queryPosition++;
    }

    return reversed;
}

// Print the contents of the PSSMatrix
void PSSMatrix_print(struct PSSMatrix PSSMatrix)
{
	int4 x, y;

	// Iterate through each row
	y = 0;
	while (y < encoding_numCodes)
	{
		// For each cell in the row
		x = 0;
		while (x < PSSMatrix.length)
		{
			// Print value
			if (PSSMatrix.matrix[x][y] == constants_sentinalScore)
				printf(" X ");
			else if (PSSMatrix.matrix[x][y] >= 0 && PSSMatrix.matrix[x][y] <= 9)
				printf(" %d ", PSSMatrix.matrix[x][y]);
			else
				printf("%d ", PSSMatrix.matrix[x][y]);
			x++;
		}
		printf("\n");
		y++;
	}
}

// Free memory used by the matrix
void PSSMatrix_free(struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP)
{
	// Do not iterate through each column and free it, since PSSMatrix shares memory
	// with scoreMatrix, which will free the columns
	PSSMatrix.matrix--;
	free(PSSMatrix.matrix);
    free(PSSMatrix.queryCodes);
    free(PSSMatrix.bestMatchCodes);

	//Shucai
	PSSMatrixFP.matrix -= encoding_numCodes;
	free(PSSMatrixFP.matrix);
	free(PSSMatrixFP.queryCodes);
	free(PSSMatrixFP.bestMatchCodes);

    if (encoding_alphabetType == encoding_nucleotide)
    {
        free(PSSMatrix.bytePackedCodes - 4);
        free(PSSMatrix.xorCodes);
        free(PSSMatrix_packedRightMatches);
        free(PSSMatrix_packedLeftMatches);
        free(PSSMatrix_packedScore);
		free(PSSMatrix_packedRightMatchScores);
		free(PSSMatrix_packedLeftMatchScores);
	}
}

// Free a copy of the PSSMatrix
void PSSMatrix_freeCopy(struct PSSMatrix PSSMatrix)
{
	free(PSSMatrix.matrix);
    free(PSSMatrix.queryCodes);
}
