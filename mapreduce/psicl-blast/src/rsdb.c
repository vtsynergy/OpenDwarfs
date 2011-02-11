// rsdb.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Remove near identical sequences in a collection

#include "blast.h"

#include "vbyte2.h"
#include "readindex.h"
#include "vec.h"
#include "vbyte2.c"
#include "readindex.c"
#include "vec.c"
#include "identityAlign.c"

uint4 decoWordLength = 30;
uint4 numIterations = 3;
uint4 decoQantum = 9;
uint4 maxListLength = 20;

#define numHashtablePasses 50
#define wildcardSumWindow 100
#define percentMatchesRequired 1
#define rsdb_threshold 90
float rsdb_averageMatchScore = 0;
float *rsdb_averageWildcodeScores = NULL;

struct longList
{
	struct wordLocation* list;
    uint4 length;
};

struct sequenceLists
{
	struct longList* longList;
	struct sequenceLists* next;
};

struct diagonal
{
	int4 relativeOffset;
    uint4 matchSequence;
//    uint4 numMatches;
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
};

float wildcardsThreshold = 0.25;
uint4 rsdb_numMergeDiagonals;
struct memBlocks* rsdb_parents;
uint4 rsdb_numParents = 0;
uint4 **rsdb_PSSM;
int4 rsdb_PSSMstart, rsdb_PSSMend, rsdb_PSSMlength;
uint4 rsdb_accumulatorsSize = 0;

// Compare two sequences' description locations
int4 rsdb_compareDescriptionLocations(const void* sequence1, const void* sequence2);
// Compare sequence lengths
int4 rsdb_compareSequenceLengths(const void* sequence1, const void* sequence2);
void rsdb_printSequence(uint4 whitespace, unsigned char* sequence, uint4 length);
// Process sequence matches
void rsdb_processSequenceMatches(struct memSingleBlock* sequenceMatches, struct sequence* sequences);
// Process a list of word locations and insert/update diagonal entries
void rsdb_processList(struct memSingleBlock *diagonals, struct wordLocation* wordLocations,
                         uint4 numWordLocations);
// Print a cluster
void rsdb_printCluster(struct parent* parent);
// Calculate a score for a sequence match
uint4 rsdb_score(struct sequence* seq1, struct sequence* seq2,
                    int4 relativeOffset, float wildscoreAllowed);
// Create a new cluster where the parent has just one child
void rsdb_newParent(struct parent* newParent, struct sequence* child);
// Add a child sequence to an existing parent/cluster
void rsdb_addChild(struct parent* parent, struct sequence* child, int4 relativeOffset);
// Merge two lists in the hashtable and return the intersection
struct wordLocation* rsdb_mergeLists(struct wordLocation* wordLocations1, uint4 *numWordLocations1,
                     struct wordLocation* wordLocations2, uint4 *numWordLocations2,
                     int4 mergeRelativeOffset, uint* lengthNewList);
// Print a word list
void rsdb_printList(struct wordLocation* list, uint4 length);
// Calculate total number of bytes saved by clustering
int4 rsdb_calculateSavings();
// Write the clusters to disk
void rsdb_writeClusters(char* filename, struct sequence* sequences);
// Free the PSSM
void rsdb_freePSSM();
// Perform SPEX algorithm and cluster sequences
void rsdb_spexClusterSequences(struct sequence* sequences, uint4 numberOfSequences,
                                  uint4 numberOfLetters);
// Update accumulators by adding match between two locations
void rsdb_updateDiagonalMatches(struct memSingleBlock *diagonals, struct wordLocation location1,
                                   struct wordLocation location2);

int4 main(int4 argc, char* argv[])
{
	unsigned char *filename, *sequence;
	uint4 descriptionStart = 0, descriptionLength = 0, sequenceLength;
	uint4 sequenceNumber = 0, encodedLength;
	struct sequence* sequences;

	// User must provide FASTA format file at command line
	if (argc != 2 && argc != 6)
	{
		fprintf(stderr, "Useage: rsdb <database> [decoWordLength] [numIterations] [decoQantum] [maxListLength]\n");
		exit(-1);
	}
	filename = argv[1];

    if (argc == 6)
    {
        decoWordLength = atoi(argv[2]);
        numIterations = atoi(argv[3]);
        decoQantum = atoi(argv[4]);
        maxListLength = atoi(argv[5]);
	}

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
    rsdb_averageMatchScore = wildcards_scoreMatrix.averageMatchScore;

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

    // Sort sequences by length
    qsort(sequences, readdb_numberOfSequences, sizeof(struct sequence),
          rsdb_compareSequenceLengths);

    // Update sequence numbers
    sequenceNumber = 0;
	while (sequenceNumber < readdb_numberOfSequences)
	{
    	sequences[sequenceNumber].number = sequenceNumber; 
    	sequenceNumber++;
	}

    // Initialize array for storing parents
	rsdb_parents = memBlocks_initialize(sizeof(struct parent), 10000);

    // Simple approach to clustering sequences
//    rsdb_simpleClusterSequences(sequences, sequenceNumber, readdb_numberOfLetters);

	// Perform SPEX algorithm and cluster sequences
    rsdb_spexClusterSequences(sequences, sequenceNumber, readdb_numberOfLetters);

    // Cluster sequences in the collection
//    rsdb_clusterSequences(sequences, sequenceNumber, readdb_numberOfLetters, filename);

    printf("Total bytes saved=%d\n", rsdb_calculateSavings()); fflush(stdout);

//    printf("END Total malloced=%s\n", global_int4toString(global_totalMalloc)); fflush(stdout);

    // Load sequence descriptions from disk
//    rsdb_loadDescriptions(sequences);

	// Write clusters to disk
    rsdb_writeClusters(filename, sequences);

    // Close reading formatted database
    readdb_close();

    free(rsdb_averageWildcodeScores);
    free(sequences);
    memBlocks_free(rsdb_parents);
    encoding_free();
    global_free();
    scoreMatrix_free(wildcards_scoreMatrix);
    parameters_free();
    identityAlign_free();
    return 0;
}

// Write the clusters to disk
void rsdb_writeClusters(char* filename, struct sequence* sequences)
{
    struct parent* parent;
	char *newFilename, *wildcardsOutFilename, *renamedFilename;
    uint4 count, sequenceNumber, sequenceLength, descriptionLength, descriptionLocation;
	unsigned char *sequence, *description;
    uint4 childNum, numChildren, descriptionTotalLength, numWilds;
	struct child* children;
	struct wild* wilds, *wildcards;

	// Construct sequence and description filenames
	newFilename = (char*)global_malloc(strlen(filename) + 12);
	sprintf(newFilename, "%s_rsdb", filename);

    // Initialize writing to formatted database
    writedb_initialize(newFilename, readdb_dbAlphabetType);

    printf("Writing representative sequences to disk..."); fflush(stdout);

    // Sort sequences back to original order
    qsort(sequences, readdb_numberOfSequences, sizeof(struct sequence),
          rsdb_compareDescriptionLocations);

    // For each sequence without a parent or only child
    sequenceNumber = 0;
	while (sequenceNumber < readdb_numberOfSequences)
	{
    	if (sequences[sequenceNumber].parent == NULL ||
            sequences[sequenceNumber].parent->numChildren == 1 ||
            sequences[sequenceNumber].parent->sequence == sequences[sequenceNumber].sequence)
        {
            sequenceLength = sequences[sequenceNumber].length;
            if (sequenceLength > 10)
            {
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
        }

        // Free the sequence
//        free(sequences[sequenceNumber].sequence);

        sequenceNumber++;
	}

    // For each parent
    memBlocks_resetCurrent(rsdb_parents);
	while ((parent = memBlocks_getCurrent(rsdb_parents)) != NULL)
    {
		// Free the parent
        free(parent->children);
    }

    printf("done.\n"); fflush(stdout);

    printf("%llu letters written\n", writedb_numberOfLetters);

    // Finalize writing to the new formatted collection
    writedb_close();

/*    // Rename clustered collection
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
    	fprintf(stderr, "Error renaming file %s to %s\n", writedb_dataFilename, renamedFilename);*/
}

// Compare two sequences' description locations
int4 rsdb_compareDescriptionLocations(const void* sequence1, const void* sequence2)
{
	const struct sequence *a1, *a2;

	a1 = (struct sequence*)sequence1;
	a2 = (struct sequence*)sequence2;

	if (a1->descriptionLocation > a2->descriptionLocation)
	{
		return -1;
	}
	else if (a1->descriptionLocation < a2->descriptionLocation)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Compare the two sequence lengths
int4 rsdb_compareSequenceLengths(const void* sequence1, const void* sequence2)
{
	const struct sequence *a1, *a2;

	a1 = (struct sequence*)sequence1;
	a2 = (struct sequence*)sequence2;

	if (a1->length > a2->length)
	{
		return -1;
	}
	else if (a1->length < a2->length)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Compare the length of two postings lists
int4 rsdb_compareListLengths(const void* list1, const void* list2)
{
	const struct finalList *a1, *a2;

	a1 = (struct finalList*)list1;
	a2 = (struct finalList*)list2;

	if (a1->numEntries > a2->numEntries)
	{
		return -1;
	}
	else if (a1->numEntries < a2->numEntries)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Perform slotted SPEX algorithm
Hashcounter *rsdb_slottedSpex(struct sequence* sequences, uint4 numberOfSequences,
                              uint4 numberOfLetters)
{
   	Hashcounter *hashCounter1 = NULL, *hashCounter2 = NULL;
	uint4 hashCounterSize = 1;
    uint4 words, duplicates, iteration = 0, seed, wordLength;
	uint4 sinceChange, submatches, lookup;
    uint4 sequenceNumber, sequenceSize, offset, suboffset;
	unsigned char* sequence, *word;
	uint4 time;
    uint4 hashValue1, hashValue2, position;

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

    return hashCounter1;
}

// Perform winnowing SPEX algorithm
Hashcounter *rsdb_winnowingSpex(struct sequence* sequences, uint4 numberOfSequences,
                                   uint4 numberOfLetters)
{
   	Hashcounter *hashCounter1 = NULL, *hashCounter2 = NULL;
	uint4 hashCounterSize = 1;
    uint4 words, duplicates, iteration = 0, seed, wordLength;
	uint4 sinceChange, submatches, lookup;
    uint4 sequenceNumber, sequenceSize, offset, suboffset;
	unsigned char* sequence, *word;
	uint4 time;
    uint4 hashValue1, hashValue2, position, smallestPosition;
	uint4 windowHashes[decoQantum], inserted[decoQantum];

//    int debug = 0;

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

//        if (debug)
//        printf("Iteration %d wordLength %d\n", iteration, wordLength);

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

                // For first few words
                while (offset < decoQantum)
                {
                	// Just add hash to the window
                    word = sequence + offset;
                    hashcounter_hash(hashCounter2, word, wordLength, hashValue2, position);
					windowHashes[offset] = hashValue2;
                    inserted[offset] = 0;

//                    if (debug)
//                    printf("New hash [%d] %u\n", offset, hashValue2);

                    offset++;
                }

                // For remaining words
                while (offset < sequenceSize - wordLength)
                {
                	// If not the first iteration
                    submatches = 0;
                	if (hashCounter1 != NULL)
     				{
                    	// Count number of sub-chunks that appear more than once
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

                    // If this word may appear more than once in the collection
                    if (hashCounter1 == NULL || submatches >= 2)
                    {
						// Hash the word and add to window
                        word = sequence + offset;
                        hashcounter_hash(hashCounter2, word, wordLength, hashValue2, position);
                        windowHashes[offset % decoQantum] = hashValue2;
                        inserted[offset % decoQantum] = 0;

//	                    if (debug)
//                        printf("New hash [%d] %u (submatches=%d) winPos=%d\n", offset, hashValue2, submatches, offset % decoQantum);

                        // Search window for smallest hash value
                        position = 0; smallestPosition = 0;
                        while (position < decoQantum)
                        {
//		                    if (debug)
//							printf("Old hash %u [%d/%d]\n", windowHashes[position], position, decoQantum);

                            if (windowHashes[position] < windowHashes[smallestPosition])
                            	smallestPosition = position;

                        	position++;
                        }

                        // If no smaller value, insert word into hash counter
                        if (!inserted[smallestPosition])
                        {
                        	inserted[smallestPosition] = 1;
//		                    if (debug)
//                            {
//                                printf("Inserted %d! ", windowHashes[smallestPosition]);
//                                print_singleSequence(sequence + offset, wordLength);
//                                printf("\n");
//							}
                            hashcounter_insert(hashCounter2, hashValue2);
						}
                    }
                    else
                    {
//	                    if (debug)
//                    	printf("Null hash\n");
                    	// Otherwise insert placeholder into window
                    	windowHashes[offset % decoQantum] = hashCounterSize;
                    }

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

    return hashCounter1;
}

// Perform SPEX algorithm and cluster sequences
void rsdb_spexClusterSequences(struct sequence* sequences, uint4 numberOfSequences,
                                  uint4 numberOfLetters)
{
   	Hashcounter *hashCounter;
    uint4 seed, wordLength, position, hashValue;
    uint4 sequenceNumber, sequenceSize, offset;
	unsigned char* sequence, *word;
    uint4 hashtablePass = 0, passhash;
    int4 key, numWordLocations;
    struct wordLocation* wordLocations, *longListWordLocations;
    int4 relativeOffset;
    uint4 sequence1, sequence2, score, location1, location2;
	struct memSingleBlock* diagonals, *sequenceDiagonals;
    struct diagonal* diagonalMatch;
    struct sequenceMatch *sequenceMatch;
    struct memSingleBlock* sequenceMatches;
	struct sequence* seq1, *seq2;
    uint4 time, numPairs;
    struct parent* parent;
    struct finalList* finalLists;
    struct longList *longList;
    struct sequenceLists **sequencesLists, *sequenceLists;
    struct memBlocks *blockLongLists, *blockSequenceLists;
    uint4 newListLength;

    // Perform SPEX algorithm
//    hashCounter = rsdb_winnowingSpex(sequences, numberOfSequences, numberOfLetters);
    hashCounter = rsdb_slottedSpex(sequences, numberOfSequences, numberOfLetters);

    // Memory to store long lists
    blockLongLists = memBlocks_initialize(sizeof(struct longList), 10000);
    blockSequenceLists = memBlocks_initialize(sizeof(struct sequenceLists), 10000);

    // Temporary buffer for storing word locations
    wordLocations = (struct wordLocation*)global_malloc(sizeof(struct wordLocation) * numberOfSequences);

    // Initialize struct to hold long lists
    sequencesLists = (struct sequenceLists**)global_malloc(sizeof(struct sequenceLists*) * numberOfSequences);

    // Initialize diagonals
    diagonals = (struct memSingleBlock*)global_malloc(sizeof(struct memSingleBlock)
                                                      * numberOfSequences);
	sequenceNumber = 0;
    while (sequenceNumber < numberOfSequences)
    {
        memSingleBlock_initializeExisting(diagonals + sequenceNumber, sizeof(struct diagonal), 1);
        sequencesLists[sequenceNumber] = NULL;
    	sequenceNumber++;
    }

    printf("Initialized diagonals\n"); fflush(stdout);

    wordLength = decoWordLength;
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
                        hashcounter_hash(hashCounter, word, wordLength, hashValue, position);
                        if (hashcounter_multiple(hashCounter, hashValue))
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
        printf("Hashtable pass %d\n", hashtablePass);
        postings_print();
        printf("Total malloced=%s\n", global_int4toString(global_totalMalloc));
		time += clock();
		printf("Hashtable construction time=%f\n", (float)time / CLOCKS_PER_SEC); fflush(stdout);
		time = -clock();

        // Sort hash table entries in order of size
		finalLists = postings_getSortedLists();

        // Start with the longest list
        key = postings_numLists - 1;
        while (key >= 0)
        {
            // If list is too long to process
            if (finalLists[key].numEntries > maxListLength)
            {
            	// Decode the long list
            	longListWordLocations = (struct wordLocation*)global_malloc(
                                        sizeof(struct wordLocation) * finalLists[key].numEntries);
                numWordLocations = postings_decodeList(finalLists[key], longListWordLocations);

                longList = (struct longList*)memBlocks_newEntry(blockLongLists);
                longList->list = longListWordLocations;
                longList->length = numWordLocations;

                rsdb_accumulatorsSize += sizeof(struct longList);
				rsdb_accumulatorsSize += sizeof(struct wordLocation) * finalLists[key].numEntries;

                // For each sequence in the list
                while (numWordLocations > 0)
                {
                	// Add a reference to this long list
                	numWordLocations--;
                    sequenceLists = (struct sequenceLists*)memBlocks_newEntry(blockSequenceLists);
                    sequenceLists->next = sequencesLists[longListWordLocations[numWordLocations].sequenceNumber];
                    sequenceLists->longList = longList;
					sequencesLists[longListWordLocations[numWordLocations].sequenceNumber] = sequenceLists;

                    rsdb_accumulatorsSize += sizeof(struct sequenceLists);
                }
            }
            else
            {
                // Get each list
                numWordLocations = postings_decodeList(finalLists[key], wordLocations);

                // Process a list of word locations and insert/update diagonal entries
                rsdb_processList(diagonals, wordLocations, numWordLocations);
			}

            key--;
        }

		time += clock();
		printf("Hashtable processing time=%f\n", (float)time / CLOCKS_PER_SEC); fflush(stdout);
		printf("."); fflush(stdout);

        printf("Accumulators size=%.1f Mb\n", (float)rsdb_accumulatorsSize / 1024.0 / 1024.0);

        free(finalLists);

        hashtablePass++;
    }

	printf("done.\n"); fflush(stdout);

    // Free hash counter
    hashcounter_free(hashCounter);

    free(wordLocations);

    printf("Identifying high-scoring sequence pairs / hierarchical clustering..."); fflush(stdout);
    time = -clock();

    sequenceMatches = memSingleBlock_initialize(sizeof(struct sequenceMatch), numberOfSequences);

    // For each sequence
    sequenceNumber = 0;
    numPairs = 0;
    while (sequenceNumber < numberOfSequences)
    {
        sequence1 = sequenceNumber;
        seq1 = sequences + sequence1;

        // Don't process sequences that already have parents
        if (seq1->parent == NULL)
        {
            // For each long list this sequence appears in
            sequenceLists = sequencesLists[sequenceNumber];
            while (sequenceLists != NULL)
            {
            	longList = sequenceLists->longList;

            	// Get the long list
                wordLocations = longList->list;
                numWordLocations = longList->length;

//				printf("Seq %d listLength %d\n", sequenceNumber, numWordLocations);

                // Find entry for sequence 1
                location1 = 0;
                while (location1 < numWordLocations)
                {
                	if (wordLocations[location1].sequenceNumber == sequence1)
                    	break;

                    location1++;
                }

                // For every entry after that
                location2 = location1 + 1;
                newListLength = location2;
                while (location2 < numWordLocations)
                {
                    if (sequences[wordLocations[location2].sequenceNumber].parent == NULL)
                    {
                        // Update accumulators
                        rsdb_updateDiagonalMatches(diagonals, wordLocations[location1],
                                                      wordLocations[location2]);

                        // Remove clustered sequences from the list
                        wordLocations[newListLength] = wordLocations[location2];

						newListLength++;
                    }

                    location2++;
                }

//                printf("Length Now %d\n", newListLength);
                longList->length = newListLength;

                sequenceLists = sequenceLists->next;
            }

            // Go through list of matching diagonals
            sequenceDiagonals = diagonals + sequenceNumber;
            memSingleBlock_resetCurrent(sequenceDiagonals);
            while ((diagonalMatch = memSingleBlock_getCurrent(sequenceDiagonals)) != NULL)
            {
                sequence2 = diagonalMatch->matchSequence;
                relativeOffset = diagonalMatch->relativeOffset;
                seq2 = sequences + sequence2;

/*
                if (minimum(seq1->length, seq2->length) > 10)
                {
                    // Calculate rough score based on number of matches, length of overlap region
                    score = 100 * diagonalMatch->numMatches * decoWordLength /
                            minimum(seq1->length, seq2->length);

//                    score = rsdb_score(seq1, seq2, relativeOffset, wildcardsThreshold);
//                    if (score > 0)
                        printf("%d,%d matches=%d w=%d q=%d i=%d seqLength=%d\n",
                               seq1->number, seq2->number, diagonalMatch->numMatches, decoWordLength,
                               decoQantum, numIterations, minimum(seq1->length, seq2->length), score);
				}
*/

//                if (score >= percentMatchesRequired)
                {
                    numPairs++;

                    // Case 1 - neither sequence has a parent
                    if (seq1->parent == NULL && seq2->parent == NULL)
                    {
                        // Compare the two sequences in detail
                        score = rsdb_score(seq1, seq2, relativeOffset, wildcardsThreshold);

                        if (score > 0)
                        {
                //        	printf("New cluster\n"); fflush(stdout);
                            parent = (struct parent*)memBlocks_newEntry(rsdb_parents);
                            rsdb_newParent(parent, seq1);
                            rsdb_addChild(parent, seq2, relativeOffset);
                        }
                    }
                    // Case 2 - sequence1 has a parent, sequence2 does not
                    else if (seq1->parent != NULL && seq2->parent == NULL)
                    {
            //        	printf("Add to %d\n", seq1->parent->numChildren); fflush(stdout);

                        // Check score for aligning sequence to parent
                        parent = seq1->parent;
                        relativeOffset += seq1->regionStart;

                        score = rsdb_score(parent->children[0], seq2, relativeOffset, wildcardsThreshold);

                        if (score > 0)
                        {
                            rsdb_addChild(seq1->parent, seq2, relativeOffset);
                        }
                    }
                }
            }
		}

        free(diagonals[sequenceNumber].block);

     	sequenceNumber++;
    }

    printf("done. (%d pairs)\n", numPairs); fflush(stdout);
    time += clock();
    printf("Total time=%f\n", (float)time / CLOCKS_PER_SEC); fflush(stdout);

    // Free long lists
    memBlocks_resetCurrent(blockLongLists);
    while ((longList = memBlocks_getCurrent(blockLongLists)) != NULL)
    {
    	free(longList->list);
    }
	memBlocks_free(blockSequenceLists);
	memBlocks_free(blockLongLists);

	free(sequencesLists);
    free(diagonals);

    memSingleBlock_free(sequenceMatches);
}

// Process a list of word locations and insert/update diagonal entries
void rsdb_processList(struct memSingleBlock *diagonals, struct wordLocation* wordLocations,
                         uint4 numWordLocations)
{
	uint4 count1, count2;
    struct wordLocation location1, location2;

    // For each pair of locations
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

            // Update accumulators by adding match between two locations
            rsdb_updateDiagonalMatches(diagonals, location1, location2);

            count2++;
        }
        count1++;
    }

//    printf("listLength=%d numNewDiagonalMatches=%d\n", numWordLocations,
//           numNewDiagonalMatches);
}

// Update accumulators by adding match between two locations
void rsdb_updateDiagonalMatches(struct memSingleBlock *diagonals, struct wordLocation location1,
                                   struct wordLocation location2)
{
    struct memSingleBlock *sequenceDiagonals;
    struct diagonal* diagonalMatch;
	int4 relativeOffset;
    char matchFound;

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
//                diagonalMatch->numMatches++;
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
//            diagonalMatch->numMatches = 1;

            rsdb_accumulatorsSize += sizeof(struct diagonal);
        }
    }
}

void rsdb_printSequence(uint4 whitespace, unsigned char* sequence, uint4 length)
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
void rsdb_printList(struct wordLocation* list, uint4 length)
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
int4 rsdb_calculateSavings()
{
	int4 count, childNum, numCharsSaved, bytesSaved = 0, numWilds;
    struct sequence* child;
	struct parent* parent;

    // For each parent
    memBlocks_resetCurrent(rsdb_parents);
	while ((parent = memBlocks_getCurrent(rsdb_parents)) != NULL)
    {
    	numWilds = 0;

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
void rsdb_printCluster(struct parent* parent)
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

// Create a new cluster where the parent has just one child
void rsdb_newParent(struct parent* parent, struct sequence* child)
{
	// Make copy of the child
    parent->number = rsdb_numParents; rsdb_numParents++;
	parent->length = child->length;
	parent->sequence = child->sequence;

    // Link parent to child, and vice versa
    parent->numChildren = 1;
    parent->children = (struct sequence**)global_malloc(sizeof(struct sequence*));
    parent->children[0] = child;
    child->regionStart = 0;
    child->parent = parent;
}

// Add a child sequence to an existing parent/cluster
void rsdb_addChild(struct parent* parent, struct sequence* child, int4 relativeOffset)
{
    // If child is longer than parent
    if (child->length > parent->length)
    	return;

    // Add new child
    child->regionStart = 0;
//    child->length = parent->length;
    child->parent = parent;
//    child->sequence = parent->sequence;
    parent->numChildren++;
	parent->children = (struct sequence**)global_realloc(parent->children,
                       sizeof(struct sequence*) * parent->numChildren);
	parent->children[parent->numChildren - 1] = child;

//	printf("Added %d: \n", child->number);
//    rsdb_printCluster(parent);
}

// Calculate a score for a sequence match
uint4 rsdb_score(struct sequence* seq1, struct sequence* seq2,
                    int4 relativeOffset, float wildscoreAllowed)
{
    uint4 score, targetScore, identity;
    unsigned char* sequence1, *sequence2;
    int4 length1, length2, shorterLength;

    sequence1 = seq1->sequence;
    length1 = seq1->length;
    sequence2 = seq2->sequence;
    length2 = seq2->length;

    if (length1 < length2)
    	shorterLength = length1;
	else
    	shorterLength = length2;

    targetScore = (int)((float)rsdb_threshold * (float)shorterLength / 100.0);
    score = identityAlign_score(sequence2, length2, sequence1, length1, -relativeOffset,
                                targetScore);

	identity = 100.0 * (float)score / (float)shorterLength;

    if (identity < rsdb_threshold)
    	return 0;

//	printf("score/target=%d,%d\n", score, targetScore);

//    printf("lengths=%d,%d\n", length1, length2);
//    printf("%d:%d\n", (end1 - start1), score);

    // Return score
    return identity;
}

