// scoreMatrix.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for loading a scoring matrix (ie. BLOSUM62) from a text description
// on disk and constructing an in-memory matrix

#include <errno.h>
#include <string.h>
#include "blast.h"

// Load the score matrix (eg. BLOSUM62) from disk and return contents in an array
// 25 by 25 for the 25 possible amino acids (actually 20 plus 3 wilds, 1 unknown,
// and a sentinal code which scores poorly, and flanks sequences)
struct scoreMatrix scoreMatrix_load(char* filename)
{
	FILE* matrixFile;
	int4 MAXLINELENGTH = 8096;
	int4 lineNumber = 0;
	int4 tokenCount;
	char line[MAXLINELENGTH];
	char *token, *tempAddress;
	unsigned char columnHeadings[24];
	unsigned char rowHeadings[24];
	int2 value;
    struct scoreMatrix scoreMatrix;
	int4 x, y;

	scoreMatrix.highestValue = 0;
	scoreMatrix.lowestValue = 0;

	// Declare memory used by scoreMatrix
	scoreMatrix.matrix = (int2**)global_malloc(sizeof(int2*) * encoding_numCodes
                       + sizeof(int2) * encoding_numCodes * encoding_numCodes);
    tempAddress = (char*)scoreMatrix.matrix;
	x = 0;
	while (x < encoding_numCodes)
	{
		scoreMatrix.matrix[x] = (int2*)(tempAddress + sizeof(int2*) * encoding_numCodes +
                                             sizeof(int2) * encoding_numCodes * x);
        // Initialize the score matrix, by setting all values to sentinal score
		y = 0;
		while (y < encoding_numCodes)
		{
			scoreMatrix.matrix[x][y] = constants_sentinalScore;
			y++;
		}
        x++;
	}

    // Open file for reading
	if ((matrixFile = fopen(filename, "r")) == NULL)
	{
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error opening matrix file %s for reading\n", filename);
		exit(-1);
	}

	// Read each line in turn
	while (fgets(line, MAXLINELENGTH, matrixFile) != NULL)
	{
		// Check we didn't max out the buffer
		if (strlen(line) >= MAXLINELENGTH - 1)
		{
	        fprintf(stderr, "%s\n", strerror(errno));
			fprintf(stderr, "Error reading file %s: maximum line length %d exceeded\n",
			        filename, MAXLINELENGTH);
			exit(-1);
		}
		// Check not a comment or blank line
		if (line[0] != '\0' && line[0] != '#' && lineNumber < 25)
		{
			// Read each of the space seperated tokens from the line
			tokenCount = 0;
			token = strtok(line, " \n");
			while (token != NULL && tokenCount < 25)
			{
				// First line - tokens are column headings
				if (lineNumber == 0)
				{
					columnHeadings[tokenCount] = token[0];
				}
				// Subsequent lines, first token is row heading
				else if (tokenCount == 0) 
				{
					rowHeadings[lineNumber - 1] = token[0];
				}
				// Subsequent lines, subsequent tokens are array values
				else
				{
                	// Get integer value of token
                	value = atoi(token);

                    // Add to scoring matrix
                    scoreMatrix.matrix[encoding_getCode(rowHeadings[lineNumber - 1])]
                                      [encoding_getCode(columnHeadings[tokenCount - 1])] = value;

                    // Determine the highest and lowest values in the matrix
                    if (value > scoreMatrix.highestValue)
                    {
                        scoreMatrix.highestValue = value;
                    }
                    if (value < scoreMatrix.lowestValue)
                    {
                        scoreMatrix.lowestValue = value;
                    }
				}

				token = strtok(NULL, " \n");
				tokenCount++;
			}

			lineNumber++;
		}
	}

    fclose(matrixFile);

    // For cells in the score matrix that did not recieve a score, use 1 and
    // lowestValue instead
	x = 0;
	while (x < encoding_numCodes)
	{
		y = 0;
		while (y < encoding_numCodes)
		{
        	if (scoreMatrix.matrix[x][y] == constants_sentinalScore)
            {
            	if (x == y)
                {
                	scoreMatrix.matrix[x][y] = 1;
                }
                else
                {
                	scoreMatrix.matrix[x][y] = scoreMatrix.lowestValue;
                }
            }
			y++;
		}
        x++;
	}

    // Every letter scores well against the wildcard
	x = 0;
	while (x < encoding_numCodes)
	{
		scoreMatrix.matrix[x][encoding_aaStartWildcards] = 1;
		scoreMatrix.matrix[encoding_aaStartWildcards][x] = 1;
		x++;
	}

    // Every letter scores poorly against the sentinal code
	x = 0;
	while (x < encoding_numCodes)
	{
		scoreMatrix.matrix[x][encoding_sentinalCode] = constants_sentinalScore;
		scoreMatrix.matrix[encoding_sentinalCode][x] = constants_sentinalScore;
		x++;
	}

    // Process wildcard scores
    y = 0;
    while (y < wildcards_numClusterWildcards)
    {
        x = 0;
        while (x < encoding_numLetters)
        {
			scoreMatrix.matrix[y + encoding_aaStartWildcards][x]
            	= wildcards_clusterWildcards[y].scoreMatrixRow[x];
			scoreMatrix.matrix[x][y + encoding_aaStartWildcards]
            	= wildcards_clusterWildcards[y].scoreMatrixRow[x];
            x++;
        }
        y++;
	}

    // Calculate average match score for two residues
    scoreMatrix.averageMatchScore = 0;
    y = 0;
    while (y < encoding_numLetters)
    {
        x = 0;
        while (x < encoding_numLetters)
        {
        	scoreMatrix.averageMatchScore += scoreMatrix.matrix[y][x] * Robinson_prob[x] * Robinson_prob[y];
            x++;
        }
        y++;
	}
    scoreMatrix.averageMatchScore /= 1000000;

    return scoreMatrix;
}

// Create a nucleotide scoring matrix use match and mismatch penalties
struct scoreMatrix scoreMatrix_create(int2 match, int2 mismatch)
{
    struct scoreMatrix scoreMatrix;
	char *tempAddress;
	int4 x, y, numXcodes, numYcodes, count, xCode, yCode, total;

	scoreMatrix.highestValue = match;
	scoreMatrix.lowestValue = mismatch;

	// Declare memory used by scoreMatrix
	scoreMatrix.matrix = (int2**)global_malloc(sizeof(int2*) * encoding_numCodes
                       + sizeof(int2) * encoding_numCodes * encoding_numCodes);
    tempAddress = (char*)scoreMatrix.matrix;

    // For each row in the matrix
    x = 0;
	while (x < encoding_numCodes)
	{
    	// Initialize memory
		scoreMatrix.matrix[x] = (int2*)(tempAddress + sizeof(int2*) * encoding_numCodes +
                                             sizeof(int2) * encoding_numCodes * x);

        // For each column, determine value
		y = 0;
		while (y < encoding_numCodes)
		{
        	// If either is the sentinal code, use sentinal score
        	if (x == encoding_sentinalCode || y == encoding_sentinalCode)
            {
				scoreMatrix.matrix[x][y] = constants_sentinalScore;
            }
            // If both characters are wilds, calculate score for match
            else if (x >= encoding_numRegularLetters && y >= encoding_numRegularLetters)
            {
            	// For each possible letter for x
                total = 0; count = 0;
                xCode = 0;
                numXcodes = encoding_wildcards[x].numCodes;
                while (xCode < numXcodes)
                {
                	// For each possible letter for y
                    yCode = 0;
                    numYcodes = encoding_wildcards[y].numCodes;
                    while (yCode < numYcodes)
                    {
						// If the letters match
                        if (encoding_wildcards[x].replacementCodes[xCode] ==
                            encoding_wildcards[y].replacementCodes[yCode])
                        {
                        	count++;
                        }
                    	total++;
                        yCode++;
                    }
                    xCode++;
                }

                // Calculate frequency of the letters matching and probably score
                scoreMatrix.matrix[x][y] = (int4)ceilf(((float)match * (float)count / (float)total) +
                                           ((float)mismatch * (float)(total - count) / (float)total));
            }
            // If the characters match
            else if (x == y)
            {
				scoreMatrix.matrix[x][y] = match;
            }
            // Mismatch
            else
            {
				scoreMatrix.matrix[x][y] = mismatch;

                // If y is in x's list of ambigious codes
                count = 0;
                numXcodes = encoding_wildcards[x].numCodes;
                while (count < numXcodes)
                {
                	if (encoding_wildcards[x].replacementCodes[count] == y)
                    {
                    	// Give score based on probability of a match
						scoreMatrix.matrix[x][y] = (int4)ceilf(((float)match / (float)numXcodes) +
                        	((float)mismatch * (float)(numXcodes - 1) / (float)numXcodes));
                    }
                    count++;
                }

                // Similarly if x is in y's list of ambigious codes
                count = 0;
                numYcodes = encoding_wildcards[y].numCodes;
                while (count < numYcodes)
                {
                	if (encoding_wildcards[y].replacementCodes[count] == x)
                    {
                    	// Give score based on probability of a match
						scoreMatrix.matrix[x][y] = (int4)ceilf(((float)match / (float)numYcodes) +
                        	((float)mismatch * (float)(numYcodes - 1) / (float)numYcodes));
                    }
                    count++;
                }
            }

			y++;
		}
        x++;
	}

    return scoreMatrix;
}

// Print the contents of the score matrix
void scoreMatrix_print(struct scoreMatrix scoreMatrix)
{
	int4 x, y;

    // Print column headers
    printf("   ");
	y = 0;
	while (y < encoding_numCodes)
	{
    	printf(" %c ", encoding_getLetter(y));
        y++;
	}
    printf("\n");

    // Iterate through each row
	x = 0;
	while (x < encoding_numCodes)
	{
    	// Print row header
    	printf(" %c ", encoding_getLetter(x));

		// For each cell in the row
		y = 0;
		while (y < encoding_numCodes)
		{
			// Print value
			if (scoreMatrix.matrix[x][y] == constants_sentinalScore)
				printf(" X ");
			else if (scoreMatrix.matrix[x][y] >= 0 && scoreMatrix.matrix[x][y] <= 9)
				printf(" %d ", scoreMatrix.matrix[x][y]);
			else
				printf("%d ", scoreMatrix.matrix[x][y]);
			y++;
		}
		printf("\n");
		x++;
	}
}

// Free memory used by the score matrix
void scoreMatrix_free(struct scoreMatrix scoreMatrix)
{
	free(scoreMatrix.matrix);
}

// Given a score matrix in the format used by cafe, converts this to the
// struct used by BLAST
#ifdef CAFEMODE
struct scoreMatrix scoreMatrix_convertCafe(int4 **cafeMatrix, int4 cafeMatrixSize)
{
	FILE* matrixFile;
	int4 MAXLINELENGTH = 8096;
	int4 lineNumber = 0;
	int4 tokenCount;
	char line[MAXLINELENGTH];
	char *token, *tempAddress;
	unsigned char columnHeadings[24];
	unsigned char rowHeadings[24];
	int2 value;
    struct scoreMatrix scoreMatrix;
	int4 x, y;

	scoreMatrix.highestValue = 0;
	scoreMatrix.lowestValue = 0;

	// Declare memory used by scoreMatrix
	scoreMatrix.matrix = (int2**)global_malloc(sizeof(int2*) * encoding_numCodes
                       + sizeof(int2) * encoding_numCodes * encoding_numCodes);
    tempAddress = (char*)scoreMatrix.matrix;

    x = 0;
	while (x < encoding_numCodes)
	{
		scoreMatrix.matrix[x] = (int2*)(tempAddress + sizeof(int2*) * encoding_numCodes +
                                             sizeof(int2) * encoding_numCodes * x);
        // Initialize the score matrix, by setting all values to sentinal score
		y = 0;
		while (y < encoding_numCodes)
		{
			scoreMatrix.matrix[x][y] = constants_sentinalScore;
			y++;
		}
        x++;
	}

    // Go through the CAFE scoring matrix
    x = 0;
    while (x < cafeMatrixSize)
    {
		y = 0;
        while (y < cafeMatrixSize)
        {
        	// If both current row and column represents a valid amino acid
        	if (encoding_getCode('A' + x) != encoding_unknownCode &&
                encoding_getCode('A' + y) != encoding_unknownCode)
            {
            	// Copy value to BLAST scoring matrix
				scoreMatrix.matrix[encoding_getCode('A' + x)][encoding_getCode('A' + y)]
                	= cafeMatrix[x][y];

                // Update highest and lowest values
                if (cafeMatrix[x][y] > scoreMatrix.highestValue)
                	scoreMatrix.highestValue = cafeMatrix[x][y];

                if (cafeMatrix[x][y] < scoreMatrix.lowestValue)
                	scoreMatrix.lowestValue = cafeMatrix[x][y];
            }
            y++;
        }
        x++;
    }

    // For cells in the score matrix that did not recieve a score, use 1 and
    // lowestValue instead
	x = 0;
	while (x < encoding_numCodes)
	{
		y = 0;
		while (y < encoding_numCodes)
		{
        	if (scoreMatrix.matrix[x][y] == constants_sentinalScore)
            {
            	if (x == y)
                {
                	scoreMatrix.matrix[x][y] = 1;
                }
                else
                {
                	scoreMatrix.matrix[x][y] = scoreMatrix.lowestValue;
                }
            }
			y++;
		}
        x++;
	}

    // Every letter scores poorly against the sentinal code
	x = 0;
	while (x < encoding_numCodes)
	{
		scoreMatrix.matrix[x][encoding_sentinalCode] = constants_sentinalScore;
		scoreMatrix.matrix[encoding_sentinalCode][x] = constants_sentinalScore;
		x++;
	}

    return scoreMatrix;
}
#endif

