// readdbApp.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Tool for choosing wildcards

#include "blast.h"

// Print matrix showing amino acid co-occurence until wildcards
void chooseWilds_printOccurenceMatrix(struct wild* wilds, uint4 numWilds);

int4 main(int4 argc, char* argv[])
{
	unsigned char *filename, *readdb_address, *sequence, code, *wildcardsFilename;
	uint4 descriptionStart = 0, descriptionLength = 0, sequenceLength;
	uint4 encodedLength, numChildren, count;
	char *description;
    struct child* children, *child;
    uint4 candidateNum, change, childNum;
    uint4 numWilds = 0;
	struct wild *wilds, defaultWild, *candidates, bestNewCandidate;
	struct wild *wildSubset, *newCandidates, *bestNewCandidates;
    uint4 sizeWildSubset, numOccurences, numCandidates;
	float defaultWildscore, candidatesScore, bestScore;

	// User must provide FASTA format file at command line
	if (argc < 4)
	{
		fprintf(stderr, "Useage: chooseWilds <database> <Wildcard score constant> <Wildcards output file>\n");
		exit(-1);
	}
	filename = argv[1];
    wildcards_scoringConstant = atof(argv[2]);
    wildcardsFilename = argv[3];

	readdb_open(filename);

    printf("Number of clusters = %u\n", readdb_numberOfClusters);
    printf("Number of sequences = %u\n", readdb_numberOfSequences);
    printf("Number of volumes = %u\n", readdb_numberOfVolumes);
	printf("Total number of letters = %llu\n", readdb_numberOfLetters);
	printf("Length of longest sequence = %u\n", readdb_longestSequenceLength);
	printf("Alphabet type = %s\n", encoding_alphabetTypes[readdb_dbAlphabetType]);

	// Initialize codes array
	encoding_initialize(readdb_dbAlphabetType);

    // Load score matrix
    parameters_findScoringMatrix();
    wildcards_scoreMatrix = scoreMatrix_load(parameters_scoringMatrixPath);

    // Count occurences of each wildcard set
    wildcards_initializeCountOccurences(readdb_longestSequenceLength);
    do
    {
        // Read each sequence in the collection
        while (readdb_readSequence(&sequence, &sequenceLength, &descriptionStart,
                                   &descriptionLength, &encodedLength))
        {
        	// If a protein sequence cluster
            if (encoding_alphabetType == encoding_protein && sequenceLength + 2 != encodedLength)
            {
                // Get the children
                children = readdb_getChildren(sequence, sequenceLength, encodedLength,
                                              descriptionStart, &numChildren);

				// Add to list of occurences
                wildcards_countOccurences(children, numChildren, sequenceLength);

                childNum = 0;
                while (childNum < numChildren)
                {
                    free(children[childNum].edits);
                    free(children[childNum].sequence - 1);
                    childNum++;
                }

                free(children);
            }
        }
	}
    while (readdb_nextVolume());

    // Get final list of number of occurences of each wild
    wilds = wildcards_getOccurences(&numWilds);

    chooseWilds_printOccurenceMatrix(wilds, numWilds);

    // Build default wildcard
	defaultWild.code = 0;
    defaultWild.count = 0;
    code = 0;
    while (code < encoding_numLetters)
    {
    	setbit(defaultWild.code, code);
        code++;
    }

    // Get average score for default wildcard
    wildSubset = wildcards_getSubset(defaultWild, wilds, numWilds, &sizeWildSubset, &numOccurences);
    defaultWildscore = wildcards_averageResidueWildMatch(defaultWild, wildSubset, sizeWildSubset);
	printf("defaultWildScore=%f occurences=%d\n", defaultWildscore, numOccurences);

    // Build up list of wildcard candidates
	candidates = (struct wild*)global_malloc(sizeof(struct wild) * wildcards_numClusterWildcards);
    numCandidates = 0;
    while (numCandidates < wildcards_numClusterWildcards - 1)
    {
    	// Explore each possible option to add to list of candidates
    	count = 0;
        bestScore = 0;
        while (count < numWilds)
        {
//        	printf("set pos %d to ", numCandidates);
//			wildcards_printWildcard(wilds[count].code);
			candidates[numCandidates] = wilds[count];

            // Score a set of candidates
            candidatesScore = wildcards_scoreCandidates(candidates, numCandidates + 1,
                                                        wilds, numWilds, defaultWildscore);
//            printf("Candidates saving=%f\n", candidatesScore);
			if (candidatesScore > bestScore)
            {
            	bestScore = candidatesScore;
                bestNewCandidate = wilds[count];
            }

            count++;
        }

        printf("Score=%f Best new candidate (%d): ", bestScore, numCandidates);
		wildcards_printWildcard(bestNewCandidate.code);
		candidates[numCandidates] = bestNewCandidate;

		numCandidates++;
    }

    newCandidates = (struct wild*)global_malloc(sizeof(struct wild) * wildcards_numClusterWildcards);
    bestNewCandidates = (struct wild*)global_malloc(sizeof(struct wild) * wildcards_numClusterWildcards);

    // Perform hill climbing; consider changing each position
    change = 1;
    while (change)
    {
    	change = 0;
		candidateNum = 0;
        bestScore = 0;
        while (candidateNum < numCandidates)
        {
        	// Start with current candidates
			memcpy(newCandidates, candidates, sizeof(struct wild) * wildcards_numClusterWildcards - 1);

            // Change current position to every possible candidate
            count = 0;
            while (count < numWilds)
            {
                newCandidates[candidateNum] = wilds[count];

                // Score a possible new set of candidates
                candidatesScore = wildcards_scoreCandidates(newCandidates, numCandidates,
                                                            wilds, numWilds, defaultWildscore);

				// Check if best new candidates
                if (candidatesScore > bestScore)
                {
                    bestScore = candidatesScore;
                    memcpy(bestNewCandidates, newCandidates, sizeof(struct wild) *
                           wildcards_numClusterWildcards - 1);
                }

                count++;
            }

        	candidateNum++;
        }

        // Update candidates
        if (bestScore > wildcards_scoreCandidates(candidates, numCandidates,
                                                  wilds, numWilds, defaultWildscore))
		{
        	printf("New bestScore=%f\n", bestScore);
        	memcpy(candidates, bestNewCandidates, sizeof(struct wild) *
                   wildcards_numClusterWildcards - 1);
			change = 1;
        }

		candidateNum = 0;
        while (candidateNum < numCandidates)
        {
			wildcards_printWildcard(candidates[candidateNum].code);
        	candidateNum++;
		}
    }

    // Print out final set of clusters with default wild added
	candidates[numCandidates] = defaultWild; numCandidates++;
	wildcards_scoreCandidates(candidates, numCandidates, wilds, numWilds, defaultWildscore);
	wildcards_outputWildcards(wildcardsFilename);

    printf("%d sequences read.\n", readdb_numberOfSequences);
	fflush(stdout);

    free(candidates);
    free(newCandidates);
    free(bestNewCandidates);

	return 0;
}

// Print matrix showing amino acid co-occurence until wildcards
void chooseWilds_printOccurenceMatrix(struct wild* wilds, uint4 numWilds)
{
	uint4 count, aa1, aa2, code;
	uint4 counts[encoding_numLetters][encoding_numLetters];

    aa1 = 0;
    while (aa1 < encoding_numLetters)
    {
        aa2 = 0;
        while (aa2 < encoding_numLetters)
        {
            counts[aa1][aa2] = 0;
            aa2++;
        }
        aa1++;
    }

    count = 0;
    while (count < numWilds)
    {
		code = wilds[count].code;

        aa1 = 0;
        while (aa1 < encoding_numLetters)
        {
            aa2 = 0;
            while (aa2 < encoding_numLetters)
            {
            	if (getbit(code, aa1) && getbit(code, aa2))
                {
					counts[aa1][aa2] += wilds[count].count;
                }
                aa2++;
            }
        	aa1++;
        }

//        printf("%d: ", wilds[count].count);
//    	wildcards_printWildcard(wilds[count].code);
    	count++;
    }

    aa1 = 0;
    while (aa1 < encoding_numLetters)
    {
    	printf("%c ", encoding_getLetter(aa1));

        aa2 = 0;
        while (aa2 < aa1)
        {
            printf("%6d ", counts[aa1][aa2]);
            aa2++;
        }
        printf("\n");
        aa1++;
    }

    printf("  ");
    aa2 = 0;
    while (aa2 < encoding_numLetters)
    {
        printf("     %c ", encoding_getLetter(aa2));
        aa2++;
    }
    printf("\n");
}

