// wildcards.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for encoding representives using wildcards

#include "blast.h"

struct scoreMatrix wildcards_scoreMatrix;
int wildcards_printWilds = 0;
float wildcards_scoringConstant = 0.2;
int2 wildcards_scoreMatrixRow[encoding_aaStartWildcards];
uint4 *wildcards_wildPositions, *wildcards_wildCount, wildcards_numWilds, wildcards_numWildCodes;

struct clusterWildcard wildcards_clusterWildcards[wildcards_numClusterWildcards] =
{
    {"LVIFM",82961,{3,-1,-3,-2,3,-2,-1,-2,-3,-2,3,-2,-3,-1,2,0,4,-2,-1,-1,-3,-2,-1,-4},-0.443727},
    {"GEKRQH",141476,{-2,-1,1,0,-2,3,-1,4,0,-1,-3,3,0,3,-3,-1,-1,3,-3,-2,0,2,-1,-4},-0.069308},
    {"AVTIX",4195410,{1,2,-1,0,3,-1,3,-1,-2,-1,3,-2,-2,-1,-1,-1,0,-2,-1,-3,-2,-1,-1,-4},-0.062484},
    {"SETKDN",4584,{-2,0,0,2,-2,2,2,2,4,-1,-2,0,4,1,-3,-2,-1,0,-2,-3,3,1,-1,-4},0.170094},
    {"LVTPRFYMHCW",1034833,{1,-1,-2,-1,1,-1,1,-1,-2,2,0,1,-1,0,2,3,2,2,2,3,-2,-1,-1,-4},0.233131},
    {"AGSDPH",131854,{-2,3,3,3,-1,0,0,0,2,2,-2,-1,1,0,-2,-2,-1,1,-1,-3,1,0,0,-4},0.311481},
    {"LAGSVETKDPIRNQFYMHCWBZXU",16777215,{1,1,1,2,1,2,1,1,1,1,1,1,2,2,1,1,1,2,1,1,1,1,-1,-3},1.267173},
};

// Initialize wildCodes in clusterWildcards structure
void wildcards_initialize()
{
	uint4 count, wildCode;
	unsigned char *letters, code;
    struct clusterWildcard *clusterWildcard;

    // For each of the cluster wildcards
    count = 0;
    while (count < wildcards_numClusterWildcards)
    {
    	wildCode = 0;
		clusterWildcard = wildcards_clusterWildcards + count;

        // For each code
        letters = clusterWildcard->letters;
        while (*letters != '\0')
        {
        	// Set bit in wildcode
			code = encoding_getCode(*letters);
            setbit(wildCode, code);
        	letters++;
        }

        clusterWildcard->wildCode = wildCode;

        count++;
    }
}

// Count the number of bits in code that are set
uint4 wildcards_countBits(uint4 code)
{
	uint4 count = 0, numBits = 0;

    while (count < encoding_numLetters)
    {
    	if (getbit(code, count))
        	numBits++;
    	count++;
    }

    return numBits;
}

// Compare the number of bits set in two wilds
int wildcards_compareWilds(const void* wild1, const void* wild2)
{
	const struct wild *w1, *w2;
	uint4 num1, num2;

	w1 = (struct wild*)wild1;
	w2 = (struct wild*)wild2;

    num1 = wildcards_countBits(w1->code);
    num2 = wildcards_countBits(w2->code);

	if (num1 > num2)
	{
		return 1;
	}
	else if (num1 < num2)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

// Build a scoring row for wildcard candidate and calculate average score
float wildcards_scoreCandidates(struct wild* originalWildCandidates, uint4 numWildCandidates,
                                  struct wild* wilds, uint4 numWilds, float defaultWildscore)
{
	struct wild* wildSubset, wildCandidate, *pastSubsets = NULL, *wildCandidates;
    uint4 sizeWildSubset, sizePastSubsets = 0;
	float score, totalSaving = 0;
    uint4 numOccurences, candidateNum = 0;

    wildCandidates = (struct wild*)global_malloc(sizeof(struct wild) * numWildCandidates);
    memcpy(wildCandidates, originalWildCandidates, sizeof(struct wild) * numWildCandidates);

    // Sort wild candidates
//	qsort(wildCandidates, numWildCandidates, sizeof(struct wild), wildcards_compareWilds);

    // For each wild candidate from shortest to longest
    while (candidateNum < numWildCandidates)
    {
    	wildCandidate = wildCandidates[candidateNum];
//        wildcards_printWildcard(wildCandidate.code);

        // Get subset of candidate wildcard
        wildSubset = wildcards_getSubset(wildCandidate, wilds, numWilds, &sizeWildSubset, &numOccurences);

        // Remove already process wildcards from subset
        wildcards_removeSubset(wildSubset, &sizeWildSubset, pastSubsets, sizePastSubsets, &numOccurences);

//        printf("wildSubset=%d\n", sizeWildSubset); fflush(stdout);
        score = wildcards_averageResidueWildMatch(wildCandidate, wildSubset, sizeWildSubset);
//        printf("Average score=%f occurences=%d defScore=%f\n", score, numOccurences, defaultWildscore);
//        printf("Total saving: %f\n", (defaultWildscore - score) * numOccurences);
        totalSaving += (defaultWildscore - score) * numOccurences;

        // Add to past subsets
        pastSubsets = wildcards_joinSubset(pastSubsets, &sizePastSubsets,
                                             wildSubset, sizeWildSubset);

		// Update cluster wildcards
		memcpy(wildcards_clusterWildcards[candidateNum].scoreMatrixRow, wildcards_scoreMatrixRow,
               sizeof(int2) * encoding_numLetters);

		wildcards_clusterWildcards[candidateNum].averageScore = score;
		wildcards_clusterWildcards[candidateNum].wildCode = wildCandidate.code;

        free(wildSubset);

        candidateNum++;
	}

/*    printf("Processed wildcards:\n");
    count = 0;
    while (count < sizePastSubsets)
    {
		wildcards_printWildcard(pastSubsets[count].code);
    	count++;
    }*/

    free(pastSubsets);
    free(wildCandidates);

    return totalSaving;
}

// Compare the average score of a pair of cluster wildcards
int4 wildcards_compareClusterWildcards(const void* wild1,
                                       const void* wild2)
{
	const struct clusterWildcard *w1, *w2;

	w1 = (struct clusterWildcard*)wild1;
	w2 = (struct clusterWildcard*)wild2;

	if (w1->averageScore > w2->averageScore)
		return 1;
	else if (w1->averageScore < w2->averageScore)
		return -1;
	else
        return 0;
}

// Remove from set1 any wildcards in set2
void wildcards_removeSubset(struct wild* set1, uint4 *size1,
                            struct wild* set2, uint4 size2, uint4* numOccurences)
{
	uint4 newSize1, count1, count2, match;

    // For each wild in set1
    *numOccurences = 0;
    count1 = 0; newSize1 = 0;
    while (count1 < *size1)
    {
    	// Determine if also occurs in set2
    	match = 0; count2 = 0;
        while (count2 < size2)
        {
			if (set1[count1].code == set2[count2].code)
            	match = 1;
        	count2++;
        }

        // Only keep if does not appear in set2
        if (!match)
        {
        	set1[newSize1] = set1[count1];
            newSize1++;

            *numOccurences += set1[count1].count;
        }

    	count1++;
    }

    *size1 = newSize1;
}

// Join two sets together
struct wild* wildcards_joinSubset(struct wild* set1, uint4 *size1, struct wild* set2, uint4 size2)
{
	set1 = (struct wild*)global_realloc(set1, sizeof(struct wild) * (*size1 + size2));
    memcpy(set1 + *size1, set2, size2 * sizeof(struct wild));
	(*size1) += size2;
	return set1;
}

// Given a list of wildcard
struct wild* wildcards_getSubset(struct wild wildCandidate, struct wild* wilds,
                                 uint4 numWilds, uint4* sizeWildSubset, uint4* numOccurences)
{
	struct wild *wildSubset;
	uint4 count = 0;

    // Calculate size of subset
    *sizeWildSubset = 0;
    while (count < numWilds)
    {
        if ((wilds[count].code & wildCandidate.code) == wilds[count].code)
			(*sizeWildSubset)++;
		count++;
    }

    // Build subset
	wildSubset = (struct wild*)global_malloc(sizeof(struct wild) * (*sizeWildSubset));
	count = 0; *sizeWildSubset = 0;
    while (count < numWilds)
    {
        if ((wilds[count].code & wildCandidate.code) == wilds[count].code)
        {
//			printf("%4d: ", wilds[count].count);
//		    wildcards_printWildcard(wilds[count].code);

            wildSubset[*sizeWildSubset] = wilds[count];
			(*sizeWildSubset)++;
		}
        count++;
    }

    // Calculate total number of occurences of wildcard candidate
    count = 0;
    *numOccurences = 0;
    while (count < *sizeWildSubset)
    {
        *numOccurences += wildSubset[count].count;
        count++;
    }

    return wildSubset;
}

// Calculate the average score for aligning a residue to the wildcard candidate
float wildcards_averageResidueWildMatch(struct wild wildCandidate, struct wild* wilds,
                                        uint4 numWilds)
{
	unsigned char code;
	float score, averageScore = 0;

    code = 0;
    while (code < encoding_numLetters)
    {
    	// Get score for aligning residue 'code' with wildcard
        score = wildcards_scoreResidueWildMatch(wildCandidate, wilds, numWilds, code);
        averageScore += score * (Robinson_prob[code] / 1000.0);

        wildcards_scoreMatrixRow[code] = ceil(score - 0.5);

        code++;
    }

    return averageScore;
}

// Calculate the average score for aligning residue 'code' with given wildcard candidate
float wildcards_scoreResidueWildMatch(struct wild wildCandidate, struct wild* wilds,
                                        uint4 numWilds, uint4 code)
{
	int4 wildCount = 0;
    int4 wildResidue, bestScore;
    float averageScore, totalScore = 0, totalCount = 0;

    // For each wild
    while (wildCount < numWilds)
    {
		// Determine the highest scoring matching residue in this wild
        bestScore = constants_sentinalScore;
		wildResidue = 0;
        while (wildResidue < encoding_numLetters)
        {
			if (getbit(wilds[wildCount].code, wildResidue) &&
                wildcards_scoreMatrix.matrix[code][wildResidue] > bestScore)
            {
            	bestScore = wildcards_scoreMatrix.matrix[code][wildResidue];
            }

        	wildResidue++;
        }
//		wildcards_printWildcard(wilds[wildCount].code);
//        printf("[%d,%d]\n", wilds[wildCount].count, bestScore);

		totalScore += ((int4)wilds[wildCount].count * bestScore);
    	totalCount += wilds[wildCount].count;

        wildCount++;
    }

    averageScore = (float)totalScore / (float)totalCount;

    // If code is included in wildcandidate, adjust score
    if (getbit(wildCandidate.code, code))
    {
    	averageScore = wildcards_scoringConstant * wildcards_scoreMatrix.matrix[code][code]
                     + (1 - wildcards_scoringConstant) * averageScore;
    }

    return averageScore;
}

//Initialize for counting occurences of wildcards
void wildcards_initializeCountOccurences(uint4 longestSequenceLength)
{
	uint4 position, wildCode;

    // Initialize array for recording wildcard position edits
    wildcards_wildPositions = (uint4*)global_malloc(sizeof(uint4) * longestSequenceLength);
	position = 0;
    while (position < longestSequenceLength)
    {
    	wildcards_wildPositions[position] = 0;
        position++;
    }

    // Initialize array for count wilding occurences
    wildcards_numWildCodes = ceil(pow(2, encoding_numLetters));
    wildcards_wildCount = (uint4*)global_malloc(sizeof(uint4) * wildcards_numWildCodes);
    wildCode = 0;
    while (wildCode < wildcards_numWildCodes)
    {
		wildcards_wildCount[wildCode] = 0;
    	wildCode++;
    }

    wildcards_numWilds = 0;
}

// Process set of children and add to occurence counts
void wildcards_countOccurences(struct child* children, uint4 numChildren, uint4 sequenceLength)
{
	uint4 childNum, count, position;
    struct edit* edits;
    struct child* child;
    uint4 numEdits;

    // For each child
    childNum = 0;
    while (childNum < numChildren)
    {
        child = children + childNum;

        edits = child->edits;
        numEdits = child->numEdits;

        // For each edit
        count = 0;
        while (count < child->numEdits)
        {
            // Set bit for this code appearing at this position
            position = edits->position + child->regionStart;
            setbit(wildcards_wildPositions[position], edits->code);

            edits++;
            count++;
        }

        childNum++;
    }

    // Increment count for each wildcard code and clear position counts
    position = 0;
    while (position < sequenceLength)
    {
        if (wildcards_wildPositions[position])
        {
            if (!wildcards_wildCount[wildcards_wildPositions[position]])
                wildcards_numWilds++;
            wildcards_wildCount[wildcards_wildPositions[position]]++;

            wildcards_wildPositions[position] = 0;
        }
        position++;
    }
}

// Get final list of number of occurences of each wild
struct wild* wildcards_getOccurences(uint4 *numWilds)
{
	uint4 count, wildCode;
	struct wild* wilds;

    // Construct list of wilds
    wilds = (struct wild*)global_malloc(sizeof(struct wild) * wildcards_numWilds);
    count = 0; wildCode = 0;
    while (wildCode < wildcards_numWildCodes)
    {
		if (wildcards_wildCount[wildCode])
        {
			wilds[count].code = wildCode;
			wilds[count].count = wildcards_wildCount[wildCode];
//			printf("(%d) %4d: ", count, wilds[count].count);
//            wildcards_printWildcard(wildCode);
            count++;
        }
    	wildCode++;
    }

    *numWilds = wildcards_numWilds;

    free(wildcards_wildPositions);
    free(wildcards_wildCount);

    return wilds;
}

// Read in a set of wildcards
void wildcards_readWildcards(char* filename)
{
    FILE* file;
	uint4 count = 0, letterCount;
    unsigned char code = 0;
    char* letters;

	// Open file for reading
	if ((file = fopen(filename, "r")) == NULL)
	{
		// If this fails, simply return
        return;
	}

    letters = (char*)global_malloc(encoding_aaStartWildcards + 1);

	// For each wildcard in order of average score
    while (count < wildcards_numClusterWildcards)
    {
        // Read letters
        letterCount = 0;
        while (1)
        {
        	fscanf(file, "%c", letters + letterCount);
			if (letters[letterCount] == ',')
            	break;
			letterCount++;
        }
        letters[letterCount] = '\0';
		wildcards_clusterWildcards[count].letters = letters;

        // Read details
		fscanf(file, "%d,[", &(wildcards_clusterWildcards[count].wildCode));

        // Print scores
        code = 0;
        while (code < encoding_aaStartWildcards)
        {
	        fscanf(file, "%d", &(wildcards_clusterWildcards[count].scoreMatrixRow[code]));
        	code++;

            if (code < encoding_aaStartWildcards)
            	fscanf(file, ",");
        }

		fscanf(file, "],%f\n", &(wildcards_clusterWildcards[count].averageScore));

    	count++;
    }

    free(letters);
    fclose(file);
}

// Print out the final wildcards
void wildcards_outputWildcards(char* filename)
{
	uint4 count = 0;
    unsigned char code = 0;
    FILE* file;

	// Open file for writing
	if ((file = fopen(filename, "w")) == NULL)
	{
		fprintf(stderr, "Error opening file %s for writing\n", filename);
		exit(-1);
	}

    // Sort wildcards in increasing order of average/expected score
//	qsort(wildcards_clusterWildcards, wildcards_numClusterWildcards,
//          sizeof(struct clusterWildcard), wildcards_compareClusterWildcards);

	// For each wildcard in order of average score
    while (count < wildcards_numClusterWildcards)
    {
        // Print letters
        code = 0;
        while (code < encoding_numLetters)
        {
            if (getbit(wildcards_clusterWildcards[count].wildCode, code))
            {
                fprintf(file, "%c", encoding_getLetter(code));
            }
            code++;
        }

        // Print details
		fprintf(file, ",%d,[", wildcards_clusterWildcards[count].wildCode);

        // Print scores
        code = 0;
        while (code < encoding_aaStartWildcards)
        {
	        fprintf(file, "%d", wildcards_clusterWildcards[count].scoreMatrixRow[code]);
        	code++;

            if (code < encoding_aaStartWildcards)
            	fprintf(file, ",");
        }

		fprintf(file, "],%f\n", wildcards_clusterWildcards[count].averageScore);

    	count++;
    }

    fclose(file);
}

// Print a wildcard character set
void wildcards_printWildcard(uint4 wildcard)
{
	uint4 code = 0;

    while (code < encoding_numLetters)
    {
    	if (getbit(wildcard, code))
        {
        	printf("%c", encoding_getLetter(code));
        }
        code++;
    }
    printf("\n");
}
