// cluster.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Find near identical sequences in a collection

#include "blast.h"

#include "vbyte2.h"
#include "readindex.h"
#include "vec.h"
#include "vbyte2.c"
#include "readindex.c"
#include "vec.c"

#define decoWordLength 30
#define numIterations 3
#define decoQantum 9
#define maxListLength 100
#define wildcardSumWindow 100
#define percentMatchesRequired 15
float cluster_averageMatchScore = 0;
float *cluster_averageWildcodeScores = NULL;

struct diagonal
{
	int4 relativeOffset;
    uint4 matchSequence;
    uint4 numMatches;
};

struct sequenceMatch
{
	uint4 sequence1, sequence2;
    int4 relativeOffset;
    uint4 score;
};

struct sequence
{
	uint4 number;
    uint4 length;
	unsigned char* sequence;
    struct parent* parent;
    uint4 regionStart;
	uint4 descriptionLength;
    uint4 descriptionLocation;
    char* description;
};

struct parent
{
    uint4 number;
	unsigned char* sequence;
    uint4 length;
    uint4 numChildren;
    struct sequence** children;
    struct wildCode* wildCodes;
    uint4 numWildcards;
};

struct wildCode
{
	uint4 wildCode;
    uint4 position;
};

float wildcardsThreshold = 0.25;
uint4 cluster_numMergeDiagonals;
struct memBlocks* cluster_parents;
uint4 cluster_numParents = 0;
uint4 **cluster_PSSM;
int4 cluster_PSSMstart, cluster_PSSMend, cluster_PSSMlength;

void cluster_printSequence(uint4 whitespace, unsigned char* sequence, uint4 length);
// Process sequence matches
void cluster_processSequenceMatches(struct memSingleBlock* sequenceMatches, struct sequence* sequences);
// Simple approach to clustering sequences
void cluster_simpleClusterSequences(struct sequence* sequences, uint4 numberOfSequences,
                                    uint4 numberOfLetters);
void cluster_clusterSequences(struct sequence* sequences, uint4 numSequences,
                              uint4 numberOfLetters, char* filename);
int4 cluster_compareScore(const void* sequenceMatch1, const void* sequenceMatch2);
// Process a list of word locations and insert/update diagonal entries
void cluster_processList(struct memSingleBlock *diagonals, struct wordLocation* wordLocations,
                         uint4 numWordLocations);
// Print a cluster
void cluster_printCluster(struct parent* parent);
// Calculate a score for a sequence match
uint4 cluster_score(struct sequence* seq1, struct sequence* seq2,
                    int4 relativeOffset, float wildscoreAllowed);
// Calculate the number of characters overlap between two sequences
uint4 cluster_overlap(struct sequence* sequence1, struct sequence* sequence2, int4 relativeOffset);
// Create a new cluster where the parent has just one child
void cluster_newParent(struct parent* newParent, struct sequence* child);
// Add a child sequence to an existing parent/cluster
void cluster_addChild(struct parent* parent, struct sequence* child, int4 relativeOffset);
// Merge two clusters together into one with a common parent
void cluster_mergeParents(struct parent* parent1, struct parent* parent2, int4 relativeOffset);
// Report the number of matching entries in two word location lists
uint4 cluster_checkMerge(struct wordLocation* wordLocations1, uint4 numWordLocations1,
                         struct wordLocation* wordLocations2, uint4 numWordLocations2,
                         int4* bestDiagonal);
// Merge two lists in the hashtable and return the intersection
struct wordLocation* cluster_mergeLists(struct wordLocation* wordLocations1, uint4 *numWordLocations1,
                     struct wordLocation* wordLocations2, uint4 *numWordLocations2,
                     int4 mergeRelativeOffset, uint* lengthNewList);
// Print a word list
void cluster_printList(struct wordLocation* list, uint4 length);
// Given a list of sequences containing a word, perform a multiple alignment and find
// the optimal cluster
int cluster_buildCluster(struct wordLocation* wordLocations, uint4 numWordLocations,
                         struct sequence* sequences);
// Remove any locations from the list that refer to child sequences
uint4 cluster_removeChildren(struct wordLocation* wordLocations, uint4 numWordLocations,
                             struct sequence* sequences);
// Calculate total number of bytes saved by clustering
int4 cluster_calculateSavings();
// Write the clusters to disk
void cluster_writeClusters(char* filename, struct sequence* sequences);
// Get edits for this child relative to the parent sequence
struct edit* cluster_getEdits(struct parent* parent, struct sequence* child, uint4 *numEdits);
// Build a PSSM from a list of sequences containing a word
int cluster_buildPSSM(struct wordLocation* wordLocations, uint4 numWordLocations,
                      struct sequence* sequences);
// Calculate the average best score for aligning to a set of residues defined by given wildcode
float cluster_averageWildcodeScore(uint4 wildcode);
// Free the PSSM
void cluster_freePSSM();
// Perform SPEX algorithm and cluster sequences
void cluster_spexClusterSequences(struct sequence* sequences, uint4 numberOfSequences,
                                  uint4 numberOfLetters);

int4 main(int4 argc, char* argv[])
{
	unsigned char *filename, *sequence;
	uint4 descriptionStart = 0, descriptionLength = 0, sequenceLength;
	uint4 sequenceNumber = 0, encodedLength;
	struct sequence* sequences;

	// User must provide FASTA format file at command line
	if (argc < 2)
	{
		fprintf(stderr, "Useage: cluster <Database>\n");
		exit(-1);
	}
	filename = argv[1];

    // Open formatted file for reading
    readdb_open(filename);

	printf("Number of sequences = %u\n", readdb_numberOfSequences);
	printf("Total number of letters = %llu\n", readdb_numberOfLetters);
	printf("Length of longest sequence = %u\n", readdb_longestSequenceLength);
	printf("Alphabet type = %s\n", encoding_alphabetTypes[readdb_dbAlphabetType]);

    if (readdb_dbAlphabetType != encoding_protein)
    {
    	fprintf(stderr, "Error: Only protein collection clustering currently supported\n");
		exit(-1);
    }

    if (readdb_numberOfSequences != readdb_numberOfClusters)
    {
    	fprintf(stderr, "Error clustering %s: collection has already been clustered\n", filename);
        exit(-1);
    }

    sequences = (struct sequence*)global_malloc(sizeof(struct sequence) * readdb_numberOfSequences);

	// Initialize codes array
	encoding_initialize(readdb_dbAlphabetType);

    // Initialize cluster wildcards
	wildcards_initialize();

    if (argc == 5)
    {
        // Read in a set of wildcards
        wildcards_readWildcards(argv[4]);
    }

    // Load score matrix
    parameters_findScoringMatrix();
    wildcards_scoreMatrix = scoreMatrix_load(parameters_scoringMatrixPath);

    // Calculate average score for matching two residues
    cluster_averageMatchScore = wildcards_scoreMatrix.averageMatchScore;

    // Move through index file reading sequence codes and converting back to ASCII sequences
    sequenceNumber = 0;
	while (readdb_readSequence(&sequence, &sequenceLength, &descriptionStart,
                               &descriptionLength, &encodedLength))
	{
        // Record sequence, sequence length, number, no parent
        sequences[sequenceNumber].sequence = sequence;
        sequences[sequenceNumber].length = sequenceLength;
        sequences[sequenceNumber].number = sequenceNumber;
        sequences[sequenceNumber].parent = NULL;
        sequences[sequenceNumber].descriptionLength = descriptionLength;
        sequences[sequenceNumber].descriptionLocation = descriptionStart;

		sequenceNumber++;
	}

    if (sequenceNumber != readdb_numberOfSequences)
    {
    	fprintf(stderr, "Error: only %d sequences read\n", sequenceNumber);
        exit(-1);
    }

    printf("%d sequences read.\n", sequenceNumber);
	fflush(stdout);

    // Initialize array for storing parents
	cluster_parents = memBlocks_initialize(sizeof(struct parent), 10000);

    // Simple approach to clustering sequences
//    cluster_simpleClusterSequences(sequences, sequenceNumber, readdb_numberOfLetters);

	// Perform SPEX algorithm and cluster sequences
    cluster_spexClusterSequences(sequences, sequenceNumber, readdb_numberOfLetters);

    // Cluster sequences in the collection
//    cluster_clusterSequences(sequences, sequenceNumber, readdb_numberOfLetters, filename);

    printf("Total bytes saved=%d\n", cluster_calculateSavings()); fflush(stdout);

//    printf("END Total malloced=%s\n", global_int4toString(global_totalMalloc)); fflush(stdout);

    // Load sequence descriptions from disk
//    cluster_loadDescriptions(sequences);

	// Write clusters to disk
    cluster_writeClusters(filename, sequences);

    // Close reading formatted database
    readdb_close();

    free(cluster_averageWildcodeScores);
    free(sequences);
    memBlocks_free(cluster_parents);
    encoding_free();
    global_free();
    scoreMatrix_free(wildcards_scoreMatrix);
    parameters_free();
    return 0;
}

// Write the clusters to disk
void cluster_writeClusters(char* filename, struct sequence* sequences)
{
    struct parent* parent;
	char *newFilename, *wildcardsOutFilename, *renamedFilename;
    uint4 count, sequenceNumber, sequenceLength, descriptionLength, descriptionLocation;
	unsigned char *sequence, *description;
    uint4 childNum, numChildren, descriptionTotalLength, numWilds;
	struct child* children;
	struct wild* wilds, *wildcards;

    printf("Writing clusters to disk..."); fflush(stdout);

	// Construct sequence and description filenames
	newFilename = (char*)global_malloc(strlen(filename) + 12);
	sprintf(newFilename, "%s_clustered", filename);

    // Initialize writing to formatted database
    writedb_initialize(newFilename, readdb_dbAlphabetType);

    // Initialize wildcard occurence counting
    wildcards_initializeCountOccurences(readdb_longestSequenceLength);

    // For each parent
    memBlocks_resetCurrent(cluster_parents);
	while ((parent = memBlocks_getCurrent(cluster_parents)) != NULL)
    {
    	if (parent->children != NULL && parent->numChildren > 1)
        {
            sequenceLength = parent->length;
            sequence = parent->sequence;
            numChildren = parent->numChildren;

            children = (struct child*)global_malloc(sizeof(struct child) * numChildren);

            // For each child
            childNum = 0;
            descriptionTotalLength = 0;
            while (childNum < numChildren)
            {
                // Load description from old collection
                descriptionLocation = parent->children[childNum]->descriptionLocation;
                descriptionLength = parent->children[childNum]->descriptionLength;
                description = descriptions_getDescription(descriptionLocation, descriptionLength);
//                description = parent->children[childNum]->description;

				children[childNum].description = description;
				children[childNum].descriptionLength = descriptionLength;
				children[childNum].regionStart = parent->children[childNum]->regionStart;
                children[childNum].length = parent->children[childNum]->length;
                children[childNum].edits = cluster_getEdits(parent, parent->children[childNum],
                                                            &(children[childNum].numEdits));
                childNum++;
            }

            // Add sequence to the formatted collection
            writedb_addSequence(sequence, sequenceLength, NULL, 0, NULL,
                                0, children, numChildren);

			// Count final wildcard occurences
            wildcards_countOccurences(children, numChildren, sequenceLength);

            // Free memory used for children and edits
            childNum = 0;
            while (childNum < numChildren)
            {
            	free(children[childNum].edits);
                free(children[childNum].description);
                childNum++;
            }
            free(children);

            // Print
//            if (numChildren > 1)
//    			cluster_printCluster(parent);
		}

		// Free the parent
        free(parent->sequence);
        free(parent->children);
        free(parent->wildCodes);
    }

    printf("done.\nWriting remaining sequences to disk..."); fflush(stdout);

    // For each sequence without a parent or only child
    sequenceNumber = 0;
	while (sequenceNumber < readdb_numberOfSequences)
	{
    	if (sequences[sequenceNumber].parent == NULL ||
            sequences[sequenceNumber].parent->numChildren == 1)
        {
            sequenceLength = sequences[sequenceNumber].length;
            sequence = sequences[sequenceNumber].sequence;
            descriptionLength = sequences[sequenceNumber].descriptionLength;
            descriptionLocation = sequences[sequenceNumber].descriptionLocation;

            // Load description from old collection
//            description = sequences[sequenceNumber].description;
            description = descriptions_getDescription(descriptionLocation, descriptionLength);

            // Add sequence to the formatted collection
            writedb_addSequence(sequence, sequenceLength, description,
                                descriptionLength, NULL, 0, NULL, 0);

			free(description);
		}

        // Free the sequence
//        free(sequences[sequenceNumber].sequence);

        sequenceNumber++;
	}

    printf("done.\n"); fflush(stdout);

    // Get list of wild occurences
	wilds = wildcards_getOccurences(&numWilds);

    // Score the wildcards using occurences
    wildcards = (struct wild*)global_malloc(sizeof(struct wild) * wildcards_numClusterWildcards); 
    count = 0;
    while (count < wildcards_numClusterWildcards)
    {
		wildcards[count].code = wildcards_clusterWildcards[count].wildCode;
        count++;
	}
	wildcards_scoreCandidates(wildcards, wildcards_numClusterWildcards, wilds, numWilds, 0);
    free(wilds);

    // Output wildcards scoring info
    wildcardsOutFilename = (char*)global_malloc(strlen(newFilename) + 12);
	sprintf(wildcardsOutFilename, "%s.wildcards", newFilename);
    wildcards_outputWildcards(wildcardsOutFilename);

    // Finalize writing to the new formatted collection
    writedb_close();

    // Rename clustered collection
    renamedFilename = (char*)global_malloc(strlen(filename) + 20);
	sprintf(renamedFilename, "%s.wildcards", filename);
    if (rename(wildcardsOutFilename, renamedFilename))
    	fprintf(stderr, "Error renaming file %s to %s\n", wildcardsOutFilename, renamedFilename);

	sprintf(renamedFilename, "%s.sequences", filename);
    if (rename(writedb_sequenceFilename, renamedFilename))
    	fprintf(stderr, "Error renaming file %s to %s\n", writedb_sequenceFilename, renamedFilename);

	sprintf(renamedFilename, "%s.descriptions", filename);
    if (rename(writedb_descriptionsFilename, renamedFilename))
    	fprintf(stderr, "Error renaming file %s to %s\n", writedb_descriptionsFilename, renamedFilename);

	sprintf(renamedFilename, "%s.data", filename);
    if (rename(writedb_dataFilename, renamedFilename))
    	fprintf(stderr, "Error renaming file %s to %s\n", writedb_dataFilename, renamedFilename);
}

// Get edits for this child relative to the parent sequence
struct edit* cluster_getEdits(struct parent* parent, struct sequence* child, uint4 *numEdits)
{
	uint4 count = 0;
	struct edit* edits;

    *numEdits = 0;
    edits = (struct edit*)global_malloc(sizeof(struct edit) * child->length);

    while (count < child->length)
    {
    	if (parent->sequence[count + child->regionStart] != child->sequence[count])
        {
        	edits[*numEdits].position = count;
            edits[*numEdits].code = child->sequence[count];
            (*numEdits)++;
        }
		count++;
    }

    return edits;
}

// Compare the two sequence lengths
int4 cluster_compareSequenceLengths(const void* sequence1, const void* sequence2)
{
	const struct sequence *a1, *a2;

	a1 = (struct sequence*)sequence1;
	a2 = (struct sequence*)sequence2;

	if (a1->length > a2->length)
	{
		return 1;
	}
	else if (a1->length < a2->length)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

void cluster_spexClusterSequences(struct sequence* sequences, uint4 numberOfSequences,
                                  uint4 numberOfLetters)
{
   	Hashcounter *hashCounter1 = NULL, *hashCounter2 = NULL;
    uint4 words, duplicates, iteration = 0, seed, wordLength;
    uint4 sequenceNumber, sequenceSize, offset, suboffset;
	unsigned char* sequence, *word;
    uint4 numHashtablePasses = 50, hashtablePass = 0, passhash;
	uint4 hashCounterSize = 1, sinceChange, submatches;
    int4 key, numWordLocations;
    struct wordLocation* wordLocations;
    int4 relativeOffset;
    uint4 sequence1, sequence2, score;
	struct memSingleBlock* diagonals, *sequenceDiagonals;
    struct diagonal* diagonalMatch;
    struct sequenceMatch *sequenceMatch;
    struct memSingleBlock* sequenceMatches;
	struct sequence* seq1, *seq2;
    int4 time;
    struct finalList* finalLists;
    uint4 hashValue1, hashValue2, position, lookup;

    // Calculate size of hashcounter
    while (pow(2,hashCounterSize) < numberOfLetters / decoQantum)
    {
    	hashCounterSize++;
    }
    hashCounter2 = hashcounter_new(hashCounterSize, 0);

    wordLength = decoWordLength - ((numIterations - 1) * decoQantum);

    while (iteration < numIterations)
    {
    	time = -clock();
    	words = 0; duplicates = 0;
        sinceChange = 0;

        // For each sequence
        sequenceNumber = 0;
        while (sequenceNumber < numberOfSequences)
        {
            sequence = sequences[sequenceNumber].sequence;
            sequenceSize = sequences[sequenceNumber].length;

            // For each sequence longer enough that hasn't been excluded
            if (sequenceSize >= wordLength)
            {
                // For each word
                offset = 0;
                while (offset < sequenceSize - wordLength)
                {
                    // If word has occured before
                    word = sequence + offset;
                    hashcounter_hash(hashCounter2, word, wordLength, hashValue2, position);
                    lookup = hashcounter_single(hashCounter2, hashValue2);

                    if (lookup)
                        sinceChange = 0;

                    if (lookup == 1)
                    {
                    	// If it has only occured once, now it has occured twice
                    	hashcounter_insert(hashCounter2, hashValue2);
                    }
                    else
                    {
						submatches = 0;
                        // First iteration
                        if (hashCounter1 != NULL)
                        {
                        	suboffset = offset;
                        	while (suboffset <= offset + decoQantum)
                            {
                                word = sequence + suboffset;
                                hashcounter_hash(hashCounter1, word, wordLength - decoQantum, hashValue1, position);
                            	if (hashcounter_multiple(hashCounter1, hashValue1))
                                	submatches++;
								suboffset++;
                            }
                        }
                        if (sinceChange >= decoQantum && (hashCounter1 == NULL || submatches >= 2))
                        {
							hashcounter_insert(hashCounter2, hashValue2);
                            sinceChange = 0;
                        }
					}

                    sinceChange++;
                    offset++;
                }
			}

            sequenceNumber++;
        }

        if (hashCounter1)
        	hashcounter_free(hashCounter1);
		hashCounter1 = hashCounter2;
    	seed = rand();

        if (iteration != numIterations - 1)
        	hashCounter2 = hashcounter_new(hashCounterSize, 0);

        printf("Iteration %d WordLength=%d\n", iteration, wordLength); fflush(stdout);
        hashcounter_count(hashCounter1, stdout); fflush(stdout);
//        printf("Total malloced=%s\n", global_int4toString(global_totalMalloc)); fflush(stdout);

        time += clock();
        printf("Iteration time=%.2f secs\n", (float)time / CLOCKS_PER_SEC); fflush(stdout);

        iteration++;
        wordLength += decoQantum;
	}

    // Temporary buffer for storing word locations
    wordLocations = (struct wordLocation*)global_malloc(sizeof(struct wordLocation) * numberOfSequences);

    // Initialize diagonals
    diagonals = (struct memSingleBlock*)global_malloc(sizeof(struct memSingleBlock)
                                                      * numberOfSequences);
	sequenceNumber = 0;
    while (sequenceNumber < numberOfSequences)
    {
        memSingleBlock_initializeExisting(diagonals + sequenceNumber, sizeof(struct diagonal), 1);
    	sequenceNumber++;
    }

    printf("Initialized diagonals\n"); fflush(stdout);

    wordLength -= decoQantum;
    seed = rand();

    printf("Constructing and processing postings lists"); fflush(stdout);
    // Process collection a section at a time
    while (hashtablePass < numHashtablePasses)
    {
    	time = -clock();
        postings_initialize();

        // For each sequence
        sequenceNumber = 0;
        while (sequenceNumber < numberOfSequences)
        {
            sequence = sequences[sequenceNumber].sequence;
            sequenceSize = sequences[sequenceNumber].length;

            // For each sequence long enough that hasn't been excluded
            if (sequenceSize >= wordLength)
            {
            	// Create initial passhash
            	passhash = 0;
				offset = 0;
                while (offset < wordLength)
                {
					passhash += sequence[offset];
                	offset++;
                }

                // For each word
                offset = 0;
                while (offset < sequenceSize - wordLength)
                {
                	// If to be processed this pass
                    if (passhash % numHashtablePasses == hashtablePass)
                    {
                        // Insert into hash table
                        word = sequence + offset;
                        hashcounter_hash(hashCounter1, word, wordLength, hashValue1, position);
                        if (hashcounter_multiple(hashCounter1, hashValue1))
                        {
                        	postings_addEntry(sequence + offset, wordLength, sequenceNumber, offset);
                        }
					}

                    // Update passhash
                    passhash -= sequence[offset];
                    passhash += sequence[offset + wordLength];
                    offset++;
                }
            }

            sequenceNumber++;
        }

        // Report the load on the hashtable
//        printf("Hashtable pass %d\n", hashtablePass);
//        postings_print();
//        printf("Total malloced=%s\n", global_int4toString(global_totalMalloc));
		time += clock();
//		printf("Hashtable construction time=%f\n", (float)time / CLOCKS_PER_SEC); fflush(stdout);
		time = -clock();

        // Sort hash table entries in order of size
		finalLists = postings_getSortedLists();

        // Start with the longest entry
        key = postings_numLists - 1;
        while (key >= 0)
        {
            // Get each entry in descending order of size
			numWordLocations = postings_decodeList(finalLists[key], wordLocations);

            if (numWordLocations > maxListLength)
            {
//                printf("Was %d\n", numWordLocations); fflush(stdout);

                // Remove any locations refering to already clustered sequences
                numWordLocations = cluster_removeChildren(wordLocations, numWordLocations, sequences);

                if (numWordLocations >= maxListLength &&
                    !cluster_buildPSSM(wordLocations, numWordLocations, sequences))
                {
                    while (numWordLocations >= maxListLength)
                    {
                        // Build cluster from longest list of word locations
                        cluster_buildCluster(wordLocations, numWordLocations, sequences);

                        //	cluster_printList(wordLocations, numWordLocations);

                        numWordLocations = cluster_removeChildren(wordLocations, numWordLocations, sequences);
//                        printf("Now %d\n", numWordLocations); fflush(stdout);
                    }

                    // Free the PSSM
                    cluster_freePSSM();
                }
            }

            // Process a list of word locations and insert/update diagonal entries
            if (numWordLocations > 1)
                cluster_processList(diagonals, wordLocations, numWordLocations);

            key--;
        }

        free(finalLists);

		time += clock();
//		printf("Hashtable processing time=%f\n", (float)time / CLOCKS_PER_SEC); fflush(stdout);
		printf("."); fflush(stdout);

        hashtablePass++;
    }

	printf("done.\n"); fflush(stdout);
    free(wordLocations);

    // Free hash counter
    hashcounter_free(hashCounter1);

    printf("Identifying high-scoring sequence pairs..."); fflush(stdout);

    sequenceMatches = memSingleBlock_initialize(sizeof(struct sequenceMatch), numberOfSequences);

    // For each sequence
    sequenceNumber = 0;
    while (sequenceNumber < numberOfSequences)
    {
        // Go through list of matching diagonals
    	sequenceDiagonals = diagonals + sequenceNumber;
        memSingleBlock_resetCurrent(sequenceDiagonals);
        while ((diagonalMatch = memSingleBlock_getCurrent(sequenceDiagonals)) != NULL)
        {
        	sequence1 = sequenceNumber;
            sequence2 = diagonalMatch->matchSequence;
            relativeOffset = diagonalMatch->relativeOffset;
            seq1 = sequences + sequence1;
            seq2 = sequences + sequence2;

            // Calculate rough score based on number of matches, length of overlap region
            score = 100 * diagonalMatch->numMatches /
                    cluster_overlap(seq1, seq2, relativeOffset) + decoWordLength;

            if (score > percentMatchesRequired)// && relativeOffset == 0 && seq1->length == seq2->length)
            {
            	// Compare the two sequences in detail
                score = cluster_score(seq1, seq2, relativeOffset, wildcardsThreshold);
                if (score > 0)
                {
//                    printf("Sequences %d,%d lengths %d,%d offset %d matches=%d score=%d\n",
//                           sequence1, sequence2, sequences[sequence1].length, sequences[sequence2].length,
//                           relativeOffset, diagonalMatch->numMatches, score);

                    // Record high-scoring match between pair of sequences
                    sequenceMatch = memSingleBlock_newEntry(sequenceMatches);
                    sequenceMatch->sequence1 = sequence1;
                    sequenceMatch->sequence2 = sequence2;
                    sequenceMatch->relativeOffset = relativeOffset;
                    sequenceMatch->score = score;
                }
            }

        }
    	sequenceNumber++;
    }

    printf("done. (%d pairs)\n", sequenceMatches->numEntries); fflush(stdout);

    // Free memory used to store diagonal matches
	sequenceNumber = 0;
    while (sequenceNumber < numberOfSequences)
    {
        free(diagonals[sequenceNumber].block);
    	sequenceNumber++;
    }
    free(diagonals);

    // Process sequence matches
    cluster_processSequenceMatches(sequenceMatches, sequences);

    memSingleBlock_free(sequenceMatches);
}

// Simple approach to clustering sequences
void cluster_simpleClusterSequences(struct sequence* sequences, uint4 numberOfSequences,
                                    uint4 numberOfLetters)
{
    uint4 sequenceLength = 0, sequenceCount, startGroup, endGroup, numWordLocations;
    struct sequence* sequence, *seq1, *seq2;
	struct wordLocation* wordLocations;
    uint4 count1, count2, score;
    struct memSingleBlock* sequenceMatches;
    struct sequenceMatch *sequenceMatch;

    // Temporary buffer for storing word locations
    wordLocations = (struct wordLocation*)global_malloc(sizeof(struct wordLocation) * numberOfSequences);

    // Sort sequences by length
    qsort(sequences, numberOfSequences, sizeof(struct sequence),
          cluster_compareSequenceLengths);

    sequenceMatches = memSingleBlock_initialize(sizeof(struct sequenceMatch), numberOfSequences);

    // For each sequence in order of length
    sequenceCount = 0;
    while (sequenceCount < numberOfSequences)
    {
        sequence = sequences + sequenceCount;
        sequenceLength = sequence->length;
        startGroup = sequenceCount;

        // Find end of group of sequences with same length
        while (sequenceCount < numberOfSequences)
        {
            sequence = sequences + sequenceCount;

            if (sequence->length != sequenceLength)
            {
                break;
            }
            sequenceCount++;
		}
        endGroup = sequenceCount;

        numWordLocations = 0;
        while (numWordLocations < endGroup - startGroup)
        {
			wordLocations[numWordLocations].sequenceNumber = startGroup + numWordLocations;
			wordLocations[numWordLocations].offset = 0;
        	numWordLocations++;
        }

        if (numWordLocations >= maxListLength &&
            !cluster_buildPSSM(wordLocations, numWordLocations, sequences))
		{
            while (numWordLocations >= maxListLength)
            {
                // Build cluster from longest list of word locations
                cluster_buildCluster(wordLocations, numWordLocations, sequences);

                // cluster_printList(wordLocations, numWordLocations);

                numWordLocations = cluster_removeChildren(wordLocations, numWordLocations, sequences);
//                printf("Now %d\n", numWordLocations); fflush(stdout);
            }

            // Free the PSSM
            cluster_freePSSM();
		}

        // For each pair of remaining sequences
		count1 = 0;
        while (count1 < numWordLocations)
        {
            count2 = count1 + 1;
            while (count2 < numWordLocations)
            {
				seq1 = sequences + wordLocations[count1].sequenceNumber;
				seq2 = sequences + wordLocations[count2].sequenceNumber;

                // Compare the two sequences in detail
                score = cluster_score(seq1, seq2, 0, wildcardsThreshold);

                if (score > 0)
                {
                    // Record high-scoring match between pair of sequences
                    sequenceMatch = memSingleBlock_newEntry(sequenceMatches);
                    sequenceMatch->sequence1 = wordLocations[count1].sequenceNumber;
                    sequenceMatch->sequence2 = wordLocations[count2].sequenceNumber;
                    sequenceMatch->relativeOffset = 0;
                    sequenceMatch->score = score;
                }

                count2++;
			}
            count1++;
		}
        printf("length=%d num=%d\n", sequence->length, endGroup - startGroup);
    }

    // Process sequence matches
    cluster_processSequenceMatches(sequenceMatches, sequences);

    memSingleBlock_free(sequenceMatches);
}

// Cluster sequences in the collection
void cluster_clusterSequences(struct sequence* sequences, uint4 numberOfSequences,
                              uint4 numberOfLetters, char* filename)
{
	struct wordLocation* wordLocations;
	uint4 numWordLocations;
    int4 relativeOffset, sequenceNumber;
    uint4 sequence1, sequence2, score;
	struct memSingleBlock* diagonals, *sequenceDiagonals;
    struct diagonal* diagonalMatch;
    struct sequenceMatch *sequenceMatch;
    struct memSingleBlock* sequenceMatches;
	struct sequence* seq1, *seq2;
	struct index_scanner *index_scanner;
	struct listinfo* listinfo = NULL;
    uint4 offsetCount, sequenceCount;

    // Initialize diagonals
    diagonals = (struct memSingleBlock*)global_malloc(sizeof(struct memSingleBlock)
                                                      * numberOfSequences);
	sequenceNumber = 0;
    while (sequenceNumber < numberOfSequences)
    {
        memSingleBlock_initializeExisting(diagonals + sequenceNumber, sizeof(struct diagonal), 1);
    	sequenceNumber++;
    }

    printf("Initialized diagonals\n"); fflush(stdout);

    // Temporary buffer for storing word locations
    wordLocations = (struct wordLocation*)global_malloc(sizeof(struct wordLocation) * numberOfSequences);

    index_scanner = open_index(filename);

    // First process the longest lists
    while ((listinfo = get_next_biglist(index_scanner, listinfo)) != NULL &&
           listinfo->doc_count >= maxListLength)
	{
    	printf("doc_count=%ld\n", listinfo->doc_count); fflush(stdout);

    	// For each sequence
        numWordLocations = 0;
    	sequenceCount = 0;
        while (sequenceCount < listinfo->doc_count)
        {
			// Only consider first occurence in that sequence
			offsetCount = 0;

            wordLocations[numWordLocations].sequenceNumber
                = listinfo->doc_numbers[sequenceCount] - 1;
            wordLocations[numWordLocations].offset
                = listinfo->phrase_offsets[sequenceCount][offsetCount] - 1;

            numWordLocations++;
            offsetCount++;

        	sequenceCount++;
        }

    	printf("Was %d\n", numWordLocations); fflush(stdout);

        // Remove any locations refering to already clustered sequences
        numWordLocations = cluster_removeChildren(wordLocations, numWordLocations, sequences);

        if (numWordLocations >= maxListLength &&
            !cluster_buildPSSM(wordLocations, numWordLocations, sequences))
		{
            while (numWordLocations >= maxListLength)
            {
                // Build cluster from longest list of word locations
                cluster_buildCluster(wordLocations, numWordLocations, sequences);

                //	cluster_printList(wordLocations, numWordLocations);

                numWordLocations = cluster_removeChildren(wordLocations, numWordLocations, sequences);
                printf("Now %d\n", numWordLocations); fflush(stdout);
            }

            // Free the PSSM
            cluster_freePSSM();
		}

        // Process a list of word locations and insert/update diagonal entries
        cluster_processList(diagonals, wordLocations, numWordLocations);
    }

    // Process remaining shorter lists
	index_scanner = reset_index(index_scanner);
    while ((listinfo = get_next_list(index_scanner, listinfo)) != NULL)
	{
    	if (listinfo->doc_count < maxListLength)
        {
            // For each sequence
            numWordLocations = 0;
            sequenceCount = 0;
            while (sequenceCount < listinfo->doc_count)
            {
                // Only consider first occurence in that sequence
                offsetCount = 0;

                wordLocations[numWordLocations].sequenceNumber
                    = listinfo->doc_numbers[sequenceCount] - 1;
                wordLocations[numWordLocations].offset
                    = listinfo->phrase_offsets[sequenceCount][offsetCount] - 1;

                offsetCount++;
                numWordLocations++;

                sequenceCount++;
            }

//            printf("cluster_processList with %d entries\n", numWordLocations);
            // Process a list of word locations and insert/update diagonal entries
            cluster_processList(diagonals, wordLocations, numWordLocations);
        }
    }
    close_index(index_scanner);

    free(wordLocations);

    printf("END STAGE 2\n"); fflush(stdout);

    sequenceMatches = memSingleBlock_initialize(sizeof(struct sequenceMatch), numberOfSequences);

    printf("\n");
    // For each sequence
    sequenceNumber = 0;
    while (sequenceNumber < numberOfSequences)
    {
        // Go through list of matching diagonals
    	sequenceDiagonals = diagonals + sequenceNumber;
        memSingleBlock_resetCurrent(sequenceDiagonals);
        while ((diagonalMatch = memSingleBlock_getCurrent(sequenceDiagonals)) != NULL)
        {
        	sequence1 = sequenceNumber;
            sequence2 = diagonalMatch->matchSequence;
            relativeOffset = diagonalMatch->relativeOffset;
            seq1 = sequences + sequence1;
            seq2 = sequences + sequence2;

            // Calculate rough score based on number of matches, length of overlap region
            score = 100 * diagonalMatch->numMatches /
                    cluster_overlap(seq1, seq2, relativeOffset) + decoWordLength;

            if (score > percentMatchesRequired)// && relativeOffset == 0 && seq1->length == seq2->length)
            {
            	// Compare the two sequences in detail
                score = cluster_score(seq1, seq2, relativeOffset, wildcardsThreshold);
                if (score > 0)
                {
//                    printf("Sequences %d,%d lengths %d,%d offset %d matches=%d score=%d\n",
//                           sequence1, sequence2, sequences[sequence1].length, sequences[sequence2].length,
//                           relativeOffset, diagonalMatch->numMatches, score);

                    // Record high-scoring match between pair of sequences
                    sequenceMatch = memSingleBlock_newEntry(sequenceMatches);
                    sequenceMatch->sequence1 = sequence1;
                    sequenceMatch->sequence2 = sequence2;
                    sequenceMatch->relativeOffset = relativeOffset;
                    sequenceMatch->score = score;
                }
            }

        }
    	sequenceNumber++;
    }

    printf("END STAGE 3 numSequenceMatches=%d\n", sequenceMatches->numEntries); fflush(stdout);

    // Free memory used to store diagonal matches
	sequenceNumber = 0;
    while (sequenceNumber < numberOfSequences)
    {
        free(diagonals[sequenceNumber].block);
    	sequenceNumber++;
    }
    free(diagonals);

    // Process sequence matches
    cluster_processSequenceMatches(sequenceMatches, sequences);

    memSingleBlock_free(sequenceMatches);
}

// Process sequence matches
void cluster_processSequenceMatches(struct memSingleBlock* sequenceMatches, struct sequence* sequences)
{
    struct sequenceMatch *sequenceMatch;
	struct sequence *seq1, *seq2;
    struct parent* parent;
    int4 score, temp;

    printf("Performing hierarchical clustering..."); fflush(stdout);

    // Sort the matches by score
    qsort(sequenceMatches->block, sequenceMatches->numEntries, sizeof(struct sequenceMatch),
          cluster_compareScore);

    // For each match
	memSingleBlock_resetCurrent(sequenceMatches);
    while ((sequenceMatch = memSingleBlock_getCurrent(sequenceMatches)) != NULL)
    {
    	// If second sequence is a parent, reverse sequences
    	if (sequences[sequenceMatch->sequence2].parent != NULL)
        {
			temp = sequenceMatch->sequence1;
            sequenceMatch->sequence1 = sequenceMatch->sequence2;
            sequenceMatch->sequence2 = temp;

            sequenceMatch->relativeOffset *= -1;
        }

//    	printf("sequence match %d,%d\n", sequenceMatch->sequence1, sequenceMatch->sequence2);

        seq1 = sequences + sequenceMatch->sequence1;
        seq2 = sequences + sequenceMatch->sequence2;

        //printf("[%d,%d] ", seq1->number, seq2->number);
        // Case 1 - neither sequence has a parent
		if (seq1->parent == NULL && seq2->parent == NULL)
        {
        //	printf("New cluster\n"); fflush(stdout);

//        	printf("New parent\n");
            parent = (struct parent*)memBlocks_newEntry(cluster_parents);
            cluster_newParent(parent, seq1);
            cluster_addChild(parent, seq2, sequenceMatch->relativeOffset);
        }
        // Case 2 - sequence1 has a parent, sequence2 does not
		else if (seq1->parent != NULL && seq2->parent == NULL)
        {
        //	printf("Add to %d\n", seq1->parent->numChildren); fflush(stdout);

            // Check score for aligning sequence to parent
        	sequenceMatch->relativeOffset += seq1->regionStart;
            score = cluster_score(seq1, seq2, sequenceMatch->relativeOffset, wildcardsThreshold);
//			printf("[%d][%d,%d,%d]\n", score, seq1->parent->length, seq2->length, sequenceMatch->relativeOffset);

            if (score > 0)
            {
//                printf("Add to existing\n");
                cluster_addChild(seq1->parent, seq2, sequenceMatch->relativeOffset);
          //  	printf("!"); fflush(stdout);
            }
        }
        // Case 3 - both sequences have parents
		else if (seq1->parent != seq2->parent)
        {
        //	printf("Merge %d and %d\n", seq1->parent->numChildren, seq2->parent->numChildren); fflush(stdout);

        	sequenceMatch->relativeOffset += seq1->regionStart;
        	sequenceMatch->relativeOffset -= seq2->regionStart;
            score = cluster_score(seq1, seq2, sequenceMatch->relativeOffset, wildcardsThreshold);
            if (score > 0)
            {
                cluster_mergeParents(seq1->parent, seq2->parent, sequenceMatch->relativeOffset);
			}
        }
    }

    printf("done.\n"); fflush(stdout);
}

// Process a list of word locations and insert/update diagonal entries
void cluster_processList(struct memSingleBlock *diagonals, struct wordLocation* wordLocations,
                         uint4 numWordLocations)
{
	uint4 count1, count2, numNewDiagonalMatches = 0;
    struct wordLocation location1, location2;
	int4 relativeOffset;
    char matchFound;
    struct memSingleBlock *sequenceDiagonals;
    struct diagonal* diagonalMatch;

    // For each pair of locations
    numNewDiagonalMatches = 0;
    count1 = 0;
    while (count1 < numWordLocations)
    {
        count2 = count1 + 1;
        while (count2 < numWordLocations)
        {
            // location1 < location2
            if (wordLocations[count1].sequenceNumber < wordLocations[count2].sequenceNumber)
            {
                location1 = wordLocations[count1];
                location2 = wordLocations[count2];
            }
            else
            {
                location2 = wordLocations[count1];
                location1 = wordLocations[count2];
            }

            // If the locations don't both refer to the same sequence
            if (location1.sequenceNumber != location2.sequenceNumber)
            {
                // Search list of matching diagonals for this sequence for one that is the
                // same diagonal and matching sequence as the new match
                relativeOffset = location1.offset - location2.offset;
                matchFound = 0;
                sequenceDiagonals = diagonals + location1.sequenceNumber;
                memSingleBlock_resetCurrent(sequenceDiagonals);
                while ((diagonalMatch = memSingleBlock_getCurrent(sequenceDiagonals)) != NULL)
                {
                    if (diagonalMatch->relativeOffset == relativeOffset &&
                        diagonalMatch->matchSequence == location2.sequenceNumber)
                    {
                        // If we found an identinal entry, increment number of matches
                        matchFound = 1;
                        diagonalMatch->numMatches++;
                        break;
                    }
                }

                // Otherwise create a new entry, ONLY if neither sequence1 or 2 belong to a cluster
                if (!matchFound)// && sequences[wordLocations[count1].sequenceNumber].parent == NULL
                                // && sequences[wordLocations[count2].sequenceNumber].parent == NULL)
                {
                    diagonalMatch = memSingleBlock_newEntry(sequenceDiagonals);
                    diagonalMatch->relativeOffset = relativeOffset;
                    diagonalMatch->matchSequence = location2.sequenceNumber;
                    diagonalMatch->numMatches = 1;
                    numNewDiagonalMatches++;
                }
            }

            count2++;
        }
        count1++;
    }

//    printf("listLength=%d numNewDiagonalMatches=%d\n", numWordLocations,
//           numNewDiagonalMatches);
}

// Build a PSSM from a list of sequences containing a word
int cluster_buildPSSM(struct wordLocation* wordLocations, uint4 numWordLocations,
                      struct sequence* sequences)
{
	int4 start, end, column;
	uint4 locationNum;
    unsigned char code;
    struct sequence sequence;

    cluster_PSSMstart = 0;
	cluster_PSSMend = 0;

    // Calculate the start and end of consensus region
    locationNum = 0;
    while (locationNum < numWordLocations)
    {
    	// Check sequences do not already have parents
    	if (sequences[wordLocations[locationNum].sequenceNumber].parent == NULL)
        {
            start = -wordLocations[locationNum].offset;
            if (start < cluster_PSSMstart)
                cluster_PSSMstart = -wordLocations[locationNum].offset;

            end = sequences[wordLocations[locationNum].sequenceNumber].length
                - wordLocations[locationNum].offset;
            if (end > cluster_PSSMend)
                cluster_PSSMend = end;
		}

        locationNum++;
    }

    if (cluster_PSSMstart == 0 && cluster_PSSMend == 0)
    	return 1;

//    printf("consensus start=%d end=%d\n", cluster_PSSMstart, cluster_PSSMend);

    // Initialize the PSSM
    cluster_PSSMlength = cluster_PSSMend - cluster_PSSMstart;
    cluster_PSSM = (uint4**)global_malloc(sizeof(uint4*) * cluster_PSSMlength);
    cluster_PSSM -= cluster_PSSMstart;
	column = cluster_PSSMstart;
    while (column < cluster_PSSMend)
    {
    	// Initialize each column to zero values
    	cluster_PSSM[column] = (uint4*)global_malloc(sizeof(uint4) * encoding_numCodes);
        code = 0;
        while (code < encoding_numCodes)
        {
			cluster_PSSM[column][code] = 0;
        	code++;
        }
    	column++;
    }

    // Calculate residue frequencies at each column
    locationNum = 0;
    while (locationNum < numWordLocations)
    {
    	if (sequences[wordLocations[locationNum].sequenceNumber].parent == NULL)
        {
            // For each sequence
            sequence = sequences[wordLocations[locationNum].sequenceNumber];

            // Find start and end
            start = -wordLocations[locationNum].offset;
            end = sequence.length - wordLocations[locationNum].offset;

            // Process each residue in the sequence and update count in PSSM
            column = 0;
            while (column < sequence.length)
            {
//            	printf("[%d]", sequence.sequence[column]); fflush(stdout);
                cluster_PSSM[column + start][sequence.sequence[column]]++;
                column++;
            }
		}

        locationNum++;
    }

    return 0;
}

// Given a list of sequences containing a word, perform a multiple alignment and find
// the optimal cluster
int cluster_buildCluster(struct wordLocation* wordLocations, uint4 numWordLocations,
                         struct sequence* sequences)
{
	uint4 *toCluster, numToCluster, *toCluster2, numToCluster2, *toCluster3;
	int4 start, end, column;
    int4 score, bestScore, bestLocation;
    struct sequence sequence;
	uint4 locationNum;
    struct parent* parent;
    struct wordLocation firstMember;

//    printf("cluster_buildCluster["); fflush(stdout);

    // Compare sequences to the PSSM
    bestScore = 0; bestLocation = 0;
    locationNum = 0;
    while (locationNum < numWordLocations)
    {
    	if (sequences[wordLocations[locationNum].sequenceNumber].parent == NULL)
        {
            // For each sequence
            sequence = sequences[wordLocations[locationNum].sequenceNumber];

            // Find start and end
            start = -wordLocations[locationNum].offset;
            end = sequence.length - wordLocations[locationNum].offset;

            // Score sequence against PSSM
            score = 0;
            column = start;
            while (column < end)
            {
                score += cluster_PSSM[column][sequence.sequence[column - start]];
                column++;
            }

            // Check if this sequence is most similiar to PSSM
            score = 100 * score / sequence.length;
            if (score > bestScore)
            {
                bestLocation = locationNum;
                bestScore = score;
            }

//            printf("[%3d] ", score);
//            cluster_printSequence(start - cluster_PSSMstart, sequence.sequence, sequence.length);
		}

        locationNum++;
    }

    // Closest to the consensus becomes the first member of the new cluster
    parent = (struct parent*)memBlocks_newEntry(cluster_parents);
    cluster_newParent(parent, sequences + wordLocations[bestLocation].sequenceNumber);
    firstMember = wordLocations[bestLocation];

    // Remove first member from PSSM
    sequence = sequences[wordLocations[bestLocation].sequenceNumber];
    start = -wordLocations[bestLocation].offset;
    end = sequence.length - wordLocations[bestLocation].offset;
    column = start;
    while (column < end)
    {
        cluster_PSSM[column][sequence.sequence[column - start]]--;
        column++;
    }

	toCluster = (uint4*)global_malloc(sizeof(uint4) * numWordLocations);
    toCluster2 = (uint4*)global_malloc(sizeof(uint4) * numWordLocations);
    numToCluster = 0; numToCluster2 = 0;

    // Initially consider every sequence
    locationNum = 0;
    while (locationNum < numWordLocations)
    {
    	// If not already a member of a cluster
    	if (sequences[wordLocations[locationNum].sequenceNumber].parent == NULL)
        {
        	// Add to list to align to new cluster's parent
        	toCluster[numToCluster] = locationNum;
            numToCluster++;
		}
        locationNum++;
	}

//    printf("!"); fflush(stdout);

    // Repeat until no more sequences are similiar enough to the parent
    while (numToCluster > 0)
	{
        // For each sequence still worth considering
        bestScore = 0; bestLocation = 0;
//        printf("%d to cluster\n", numToCluster);
        while (numToCluster > 0)
        {
            numToCluster--;
            locationNum = toCluster[numToCluster];

    		if (sequences[wordLocations[locationNum].sequenceNumber].parent == NULL)
			{
                // Get sequence, start and end
                sequence = sequences[wordLocations[locationNum].sequenceNumber];
                start = -wordLocations[locationNum].offset;
                end = sequence.length - wordLocations[locationNum].offset;

                // Calculate score for aligning sequence to parent of cluster
                score = cluster_score(sequences + firstMember.sequenceNumber, &sequence,
                        firstMember.offset - wordLocations[locationNum].offset,
                        wildcardsThreshold);

                // If above cutoff
                if (score > 0)
                {
                    // Add to list
                    toCluster2[numToCluster2] = locationNum;
                    numToCluster2++;

                    // Calculate highest scoring match
                    if (score > bestScore)
                    {
                        bestScore = score;
                        bestLocation = locationNum;
                    }

//                    cluster_printSequence(start - cluster_PSSMstart, sequence.sequence, sequence.length);
//                    printf("Score=%d\n", score);
                }
        	}
        }

        // If we have matches to our close-to-consensus sequence
        if (bestScore > 0)
        {
            // Add the closest match to the new cluster
            cluster_addChild(parent, sequences + wordLocations[bestLocation].sequenceNumber,
                             firstMember.offset - wordLocations[bestLocation].offset +
                             sequences[firstMember.sequenceNumber].regionStart);

            // Remove new member from PSSM
            sequence = sequences[wordLocations[bestLocation].sequenceNumber];
            start = -wordLocations[bestLocation].offset;
            end = sequence.length - wordLocations[bestLocation].offset;
            column = start;
            while (column < end)
            {
                cluster_PSSM[column][sequence.sequence[column - start]]--;
                column++;
            }

            // Switch clusters so processing new list of close matches to parent
            toCluster3 = toCluster;
            toCluster = toCluster2;
            toCluster2 = toCluster3;
            numToCluster = numToCluster2;
            numToCluster2 = 0;
        }

//    	cluster_printCluster(parent);
	}

//    cluster_printCluster(parent);

	free(toCluster);
    free(toCluster2);

//    printf("]\n"); fflush(stdout);

    return 0;
}

// Free the PSSM
void cluster_freePSSM()
{
	int4 column;

    column = cluster_PSSMstart;
    while (column < cluster_PSSMend)
    {
    	free(cluster_PSSM[column]);
    	column++;
    }
    cluster_PSSM += cluster_PSSMstart;
    free(cluster_PSSM);
}

// Remove any duplicate locations in the list that refer to children from the same parent;
// only keep the first child of that parent
uint4 cluster_removeChildren(struct wordLocation* wordLocations, uint4 numWordLocations,
                             struct sequence* sequences)
{
	struct wordLocation *wordLocations2;
    uint4 numWordLocations2, locationNum;
    struct sequence sequence;

    wordLocations2 = (struct wordLocation*)global_malloc(sizeof(struct wordLocation) * numWordLocations);
    numWordLocations2 = 0;
    locationNum = 0;
    while (locationNum < numWordLocations)
    {
    	// If not a member of a cluster
        sequence = sequences[wordLocations[locationNum].sequenceNumber];
    	if (sequence.parent == NULL)
        {
        	// Add to new list of word locations
			wordLocations2[numWordLocations2] = wordLocations[locationNum];
            numWordLocations2++;
        }
        locationNum++;
	}

    // Return new pruned list of word locations
	memcpy(wordLocations, wordLocations2, sizeof(struct wordLocation) * numWordLocations2);
    free(wordLocations2);
    return numWordLocations2;
}

void cluster_printSequence(uint4 whitespace, unsigned char* sequence, uint4 length)
{
    // Print whitespace before sequence begins
    while (whitespace > 0)
    {
        printf(" ");
        whitespace--;
    }

    print_singleSequence(sequence, length); printf("\n");
}

// Print a word list
void cluster_printList(struct wordLocation* list, uint4 length)
{
	uint4 count = 0;

    while (count < length)
    {
    	printf("[%d,%d] ", list[count].sequenceNumber, list[count].offset);
        count++;
    }
    printf(".\n");
}

// Calculate total number of bytes saved by clustering
int4 cluster_calculateSavings()
{
	int4 count, childNum, numCharsSaved, bytesSaved = 0, numWilds;
    struct sequence* child;
	struct parent* parent;

    // For each parent
    memBlocks_resetCurrent(cluster_parents);
	while ((parent = memBlocks_getCurrent(cluster_parents)) != NULL)
    {
    	count = 0; numWilds = 0;
        while (count < parent->length)
        {
			if (parent->sequence[count] >= encoding_aaStartWildcards)
				numWilds++;
            count++;
        }

    	// If it has children
    	if (parent->children != NULL)
        {
            // For each child calculate number of characters saved
            numCharsSaved = 0;
            childNum = 0;
            while (childNum < parent->numChildren)
            {
                child = parent->children[childNum];

                // CHECK: no child has two parents
                if (child->parent != parent)
                {
                	printf("Error sequence %d has multiple parents\n", child->number);
                    exit(-1);
                }

                numCharsSaved += child->length;
                childNum++;
            }

            // Minus cost of parent
            bytesSaved += numCharsSaved - parent->length;

            if (numCharsSaved < parent->length)
            {
				printf("Error parent length %d with %d children saves negative bytes",
                        parent->length, parent->numChildren);
                exit(-1);
            }
        }
    }

	return bytesSaved;
}

// Print a cluster
void cluster_printCluster(struct parent* parent)
{
	int4 count, childNum, numCharsSaved;
    struct sequence* child;

    printf("### Parent %d with %d children ###\nParent    ", parent->number, parent->numChildren);
	print_singleSequence(parent->sequence, parent->length); printf("\n");

    // For each child
    numCharsSaved = 0;
    childNum = 0;
    while (childNum < parent->numChildren)
    {
    	child = parent->children[childNum];

        printf("[%7d] ", child->number);
    	// Print spaces to align child to parent
        count = 0;
        while (count < child->regionStart)
        {
            printf(" ");
            count++;
        }

        // Print child
        print_singleSequence(child->sequence, child->length); printf("\n");
//       	printf("Child Length=%d start=%d\n", child->length, child->regionStart);

        numCharsSaved += child->length;
        childNum++;
    }
    printf("%d bytes saved\n", numCharsSaved - parent->length);
    printf("\n");
}

// Calculate the number of characters overlap between two sequences
uint4 cluster_overlap(struct sequence* sequence1, struct sequence* sequence2, int4 relativeOffset)
{
	uint4 start1, end1;

    // Find start of match region
    if (relativeOffset > 0)
        start1 = relativeOffset;
    else
        start1 = 0;

    // Find end of match region
    if (sequence1->length < sequence2->length + relativeOffset)
        end1 = sequence1->length;
    else
        end1 = sequence2->length + relativeOffset;

    return end1 - start1;
}

// Create a new cluster where the parent has just one child
void cluster_newParent(struct parent* parent, struct sequence* child)
{
	// Make copy of the child
    parent->number = cluster_numParents; cluster_numParents++;
	parent->length = child->length;
	parent->sequence = (char*)global_malloc(parent->length);
    memcpy(parent->sequence, child->sequence, child->length);

    // Initialize list of wildcards
    parent->wildCodes = NULL;
    parent->numWildcards = 0;

    // Link parent to child, and vice versa
    parent->numChildren = 1;
    parent->children = (struct sequence**)global_malloc(sizeof(struct sequence*));
    parent->children[0] = child;
    child->regionStart = 0;
    child->parent = parent;
}

// Update wildcodes and change position by "change"
void cluster_updateWildcodes(struct parent* parent, int4 change)
{
	uint4 count = 0;
    while (count < parent->numWildcards)
    {
		parent->wildCodes[count].position += change;
    	count++;
    }
}

// Get the parent's wildcode for position
uint4 cluster_getWildcode(struct parent* parent, uint4 position)
{
	uint4 count = 0;

    while (count < parent->numWildcards)
    {
		if (parent->wildCodes[count].position == position)
        	return parent->wildCodes[count].wildCode;
    	count++;
    }

    fprintf(stderr, "Internal error cluster_getWildcode(%d, %d)\n", parent->number, position);
    exit(-1);
}

// Update the parent's wildcode at position
void cluster_updateWildcode(struct parent* parent, uint4 position, uint4 wildCode)
{
	uint4 count;

	count = 0;
    while (count < parent->numWildcards)
    {
		if (parent->wildCodes[count].position == position)
        {
        	parent->wildCodes[count].wildCode = wildCode;
        	return;
		}
        count++;
    }

    // If position not found, add wildcode
	parent->numWildcards++;
	parent->wildCodes = (struct wildCode*)global_realloc(parent->wildCodes, parent->numWildcards
                      * sizeof(struct wildCode));

    parent->wildCodes[parent->numWildcards - 1].position = position;
    parent->wildCodes[parent->numWildcards - 1].wildCode = wildCode;
}

// Merge two clusters together into one with a common parent
void cluster_mergeParents(struct parent* parent1, struct parent* parent2, int4 relativeOffset)
{
	int4 childNum;
    struct sequence* child;
	struct parent* temp;
	uint4 wildcodeCount = 0;

    // Arrange so that parent1 comes before parent2
    if (relativeOffset < 0)
    {
    	temp = parent1;
        parent1 = parent2;
        parent2 = temp;
		relativeOffset = -relativeOffset;
	}

    // If parent2 is longer than parent1
    if (parent2->length + relativeOffset > parent1->length)
    {
    	// Add end of parent2 onto parent1
		parent1->sequence = (char*)global_realloc(parent1->sequence, parent2->length + relativeOffset);
		memcpy(parent1->sequence + parent1->length, parent2->sequence + parent1->length - relativeOffset,
               parent2->length - parent1->length + relativeOffset);
        parent1->length = parent2->length + relativeOffset;

        while (wildcodeCount < parent2->numWildcards)
        {
			cluster_updateWildcode(parent1, parent2->wildCodes[wildcodeCount].position + relativeOffset,
                                   parent2->wildCodes[wildcodeCount].wildCode);
			wildcodeCount++;
        }
    }

    // For each child in parent2
    childNum = 0;
    while (childNum < parent2->numChildren)
    {
    	// Add child to parent1
    	child = parent2->children[childNum];
		cluster_addChild(parent1, child, relativeOffset + child->regionStart);

        childNum++;
    }

    // parent2 no longer has children
    free(parent2->children);
    parent2->children = NULL;
    free(parent2->sequence);
	parent2->sequence = NULL;
    parent2->length = 0;
	parent2->numChildren = 0;
}

// Add a child sequence to an existing parent/cluster
void cluster_addChild(struct parent* parent, struct sequence* child, int4 relativeOffset)
{
	int4 childNum, count, wildCode, wildcardCount;
	unsigned char* newParentSequence;

    if (relativeOffset < 0)
    {
    	// Child starts before parent, add new sequence at beginning of parent
		newParentSequence = (char*)global_malloc(parent->length - relativeOffset);
        memcpy(newParentSequence - relativeOffset, parent->sequence, parent->length);
		memcpy(newParentSequence, child->sequence, -relativeOffset);
        free(parent->sequence);
        parent->sequence = newParentSequence;
        parent->length -= relativeOffset;

        // Update children coordinates
        childNum = 0;
        while (childNum < parent->numChildren)
        {
        	parent->children[childNum]->regionStart -= relativeOffset;
        	childNum++;
        }

        // Update wildcodes
        cluster_updateWildcodes(parent, -relativeOffset);

        relativeOffset = 0;
    }

    // If child is longer than parent
    if (child->length + relativeOffset > parent->length)
    {
    	// Add end of new sequence onto parent
		parent->sequence = (char*)global_realloc(parent->sequence, child->length + relativeOffset);
		memcpy(parent->sequence + parent->length, child->sequence + parent->length - relativeOffset,
               child->length - parent->length + relativeOffset);
        parent->length = child->length + relativeOffset;
    }

//    printf("[%d:%d:%d]\n", child->length, parent->length, relativeOffset);
    // For region that child and parent have in common
    count = 0;
    while (count < child->length)
    {
//    	printf("[%d/%d]", count + relativeOffset, parent->length); fflush(stdout);
    	// If characters don't match
		if (child->sequence[count] != parent->sequence[count + relativeOffset])
        {
//        	printf("%d: %d,%d\n", count + relativeOffset, child->sequence[count],
//                                  parent->sequence[count + relativeOffset]); fflush(stdout);

        	if (parent->sequence[count + relativeOffset] >= encoding_aaStartWildcards)
            {
                // Get parent wildcode
                wildCode = cluster_getWildcode(parent, count + relativeOffset);
			}
            else
            {
            	// Add parent character to new wildcode
	            wildCode = 0;
            	setbit(wildCode, parent->sequence[count + relativeOffset]);
            }

            // If child residue will change wildcode
            if (!getbit(wildCode, child->sequence[count]))
            {
                // Add child character to wildcode
                setbit(wildCode, child->sequence[count]);
				cluster_updateWildcode(parent, count + relativeOffset, wildCode);

//                printf("set %d pos %d = %d\n", parent->number, count + relativeOffset, wildCode);

                // Check for a match in set of wildcards
                wildcardCount = 0;
                while (wildcardCount < wildcards_numClusterWildcards)
                {
                    if ((wildcards_clusterWildcards[wildcardCount].wildCode & wildCode) == wildCode)
                    {
//                    	printf("Match: ");
//    		            wildcards_printWildcard(wildcards_clusterWildcards[wildcardCount].wildCode);
                        break;
                    }
                    wildcardCount++;
                }

                if (wildcardCount >= wildcards_numClusterWildcards)
                {
                    printf("Error %d >= %d\n", wildcardCount, wildcards_numClusterWildcards);
                    wildcards_printWildcard(wildCode);
                    exit(-1);
                }

                // Insert wildcard character into parent
                parent->sequence[count + relativeOffset] = encoding_aaStartWildcards + wildcardCount;
            }

//            printf("Inserted wildNum %d into parent\n", wildcardCount);
        }
    	count++;
    }

    // Add new child
    child->regionStart = relativeOffset;
    child->parent = parent;
    parent->numChildren++;
	parent->children = (struct sequence**)global_realloc(parent->children,
                       sizeof(struct sequence*) * parent->numChildren);
	parent->children[parent->numChildren - 1] = child;

//	printf("Added %d: \n", child->number);
//    cluster_printCluster(parent);
}

// Compare the two sequence match's scores
int4 cluster_compareScore(const void* sequenceMatch1, const void* sequenceMatch2)
{
	const struct sequenceMatch *a1, *a2;

	a1 = (struct sequenceMatch*)sequenceMatch1;
	a2 = (struct sequenceMatch*)sequenceMatch2;

	if (a1->score > a2->score)
	{
		return -1;
	}
	else if (a1->score < a2->score)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Calculate a score for a sequence match
uint4 cluster_score(struct sequence* seq1, struct sequence* seq2,
                    int4 relativeOffset, float wildscoreAllowed)
{
	int4 pos1, pos2, start1, end1, start2, end2;
    float sumWildscores = 0, windowWildscoreAllowed, sumWindowWildscores = 0;
    uint4 wildCode, wildcardCount;
    unsigned char* sequence1, *sequence2;
    int4 length1, length2;
    struct parent* parent1 = NULL, *parent2 = NULL;

    // Check if sequence1 has a parent
    if (seq1->parent == NULL)
    {
    	sequence1 = seq1->sequence;
        length1 = seq1->length;
    }
    else
    {
    	parent1 = seq1->parent;
    	sequence1 = seq1->parent->sequence;
        length1 = seq1->parent->length;
        relativeOffset += seq1->regionStart;
    }

    // Check if sequence2 has a parent
    if (seq2->parent == NULL)
    {
    	sequence2 = seq2->sequence;
        length2 = seq2->length;
    }
    else
    {
    	parent2 = seq2->parent;
    	sequence2 = seq2->parent->sequence;
        length2 = seq2->parent->length;
        relativeOffset -= seq2->regionStart;
    }

    // Find start of match region
    if (relativeOffset > 0)
    {
        start1 = relativeOffset;
        start2 = 0;
    }
    else
    {
        start1 = 0;
        start2 = -relativeOffset;
    }

    // Find end of match region
    if (length1 < length2 + relativeOffset)
    {
        end1 = length1;
        end2 = length1 - relativeOffset;
    }
    else
    {
        end1 = length2 + relativeOffset;
        end2 = length2;
    }

    // Parents do not overlap, return zero
    if (relativeOffset > length1 || -relativeOffset > length2)
    {
//    	printf("ro=%d\n", relativeOffset);
//    	printf("[%d,%d,%d]\n", start1, end1, length1);
//    	printf("[%d,%d,%d]\n", start2, end2, length2);
        return 0;
    }

	pos1 = start1;
    pos2 = start2;

    windowWildscoreAllowed = wildscoreAllowed * wildcardSumWindow;
    wildscoreAllowed *= (end1 - start1);

    // For each position where the two sequences overlap
    while (pos1 < end1)
    {
    	// If the sequences differ
		if (sequence1[pos1] != sequence2[pos2])
        {
            // Build a wildcode for this position
            wildCode = 0;

//            printf("%c != %c p1=%p p2=%p\n", encoding_getLetter(sequence1[pos1]),
//                   encoding_getLetter(sequence2[pos2]), parent1, parent2);

			// Sequence1 char is included in parent2's wildcard
			if (parent1 == NULL && parent2 != NULL &&
                sequence2[pos2] >= encoding_aaStartWildcards &&
                getbit(wildcards_clusterWildcards[sequence2[pos2] -
                encoding_aaStartWildcards].wildCode, sequence1[pos1]))
			{
            	wildCode = wildcards_clusterWildcards[sequence2[pos2] -
                           encoding_aaStartWildcards].wildCode;
            }
			// Sequence2 char is included in parent1's wildcard
			else if (parent2 == NULL && parent1 != NULL &&
                sequence1[pos1] >= encoding_aaStartWildcards &&
                getbit(wildcards_clusterWildcards[sequence1[pos1] -
                encoding_aaStartWildcards].wildCode, sequence2[pos2]))
			{
            	wildCode = wildcards_clusterWildcards[sequence1[pos1] -
                           encoding_aaStartWildcards].wildCode;
            }
			else
            {
                // If we are processing a single sequence (1)
                if (parent1 == NULL || sequence1[pos1] < encoding_aaStartWildcards)
                {
                    setbit(wildCode, sequence1[pos1]);
                }
                // Else we are processing a parent (1)
                else
                {
                    wildCode |= cluster_getWildcode(parent1, pos1);
                }

                // If we a processing a single sequence (2)
                if (parent2 == NULL || sequence2[pos2] < encoding_aaStartWildcards)
                {
                    setbit(wildCode, sequence2[pos2]);
                }
                // Else we are processing a parent (2)
                else
                {
                    wildCode |= cluster_getWildcode(parent2, pos2);
                }
			}

            //wildCodeAverageScore = cluster_averageWildcodeScore(wildCode);
            //wildCodeAverageScore = 0;
            //wildcards_printWildcard(wildCode);
			//printf("Average score=%f\n", wildCodeAverageScore);

            // Check for a match in set of wildcards
            wildcardCount = 0;
            while (wildcardCount < wildcards_numClusterWildcards)
            {
                if ((wildcards_clusterWildcards[wildcardCount].wildCode & wildCode) == wildCode)
                {
//                	printf("Match: ");
//		            wildcards_printWildcard(wildcards_clusterWildcards[wildcardCount].wildCode);
                    break;
                }
                wildcardCount++;
            }

            sumWildscores += wildcards_clusterWildcards[wildcardCount].averageScore
                           - cluster_averageMatchScore;
            sumWindowWildscores += wildcards_clusterWildcards[wildcardCount].averageScore
                                 - cluster_averageMatchScore;

/*            // Add wild's average score to tally
            if (wildcards_clusterWildcards[wildcardCount].averageScore > 0)
            {
                sumWildscores += wildcards_clusterWildcards[wildcardCount].averageScore;
                sumWindowWildscores += wildcards_clusterWildcards[wildcardCount].averageScore;
			}*/

            // Stop when more than allowed total average score of wilds
            if (sumWildscores > wildscoreAllowed)
            	return 0;

            if (sumWindowWildscores > windowWildscoreAllowed)
            	return 0;
        }

        // Reset the window sum at the end of each window
        if (pos1 % wildcardSumWindow == 0)
			sumWindowWildscores = 0;

        pos1++;
        pos2++;
    }

    // Return score
//    return (end1 - start1) * (wildscoreAllowed - sumWildscores) / wildscoreAllowed;
    return 2000 - (1000 * sumWildscores / (end1 - start1));
}

// Calculate the average best score for aligning to a set of residues defined by given wildcode
float cluster_averageWildcodeScore(uint4 wildCode)
{
	unsigned char code, wildResidue;
	float averageScore = 0;
    int4 bestScore, tempCode;

    // Initialize array to hold precomputed answers
    if (cluster_averageWildcodeScores == NULL)
    {
    	cluster_averageWildcodeScores = (float*)global_malloc(sizeof(float)
                                      * ceil(pow(2, encoding_numLetters)));
		tempCode = 0;
		while (tempCode < ceil(pow(2, encoding_numLetters)))
        {
			cluster_averageWildcodeScores[tempCode] = 0;
        	tempCode++;
        }
    }

    if (cluster_averageWildcodeScores[wildCode] == 0)
	{
        // For each residue
        code = 0;
        while (code < encoding_numLetters)
        {
            // For each residue in wildcode
            bestScore = constants_sentinalScore;
            wildResidue = 0;
            while (wildResidue < encoding_numLetters)
            {
                if (getbit(wildCode, wildResidue))
                {
                    if (wildcards_scoreMatrix.matrix[code][wildResidue] > bestScore)
                        bestScore = wildcards_scoreMatrix.matrix[code][wildResidue];
                }
                wildResidue++;
            }

    //        printf("averageScore=%f\n", (float)bestScore);
            averageScore += (float)bestScore * (Robinson_prob[code] / 1000.0);

            code++;
        }

        cluster_averageWildcodeScores[wildCode] = averageScore;
	}

    return cluster_averageWildcodeScores[wildCode];
}
