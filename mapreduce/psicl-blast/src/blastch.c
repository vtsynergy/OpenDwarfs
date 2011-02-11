// blast.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Main code for blast

#include "blast.h"

//void blast_search(char* searchDbFile, struct PSSMatrix PSSMatrix, char* query);
//Shucai
void blast_search(char* searchDbFile, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, char* query);

int4 main(int4 argc, char* argv[])
{
	char *query, *queryDescription;
    unsigned char queryAlphabetType, previousQueryAlphabetType = 10;
	struct scoreMatrix scoreMatrix;
	struct PSSMatrix PSSMatrix;
	// Shucai
	struct PSSMatrixFP PSSMatrixFP;

    #ifdef SSEARCH
    parameters_ssearch = 1;
    #endif

	// Process command line arguments
	parameters_processArguments(argc, argv);

    if (parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
	{
        if (parameters_ssearch)
            printf("FSA-SSEARCH 1.05\n\n");
        else
            printf("FSA-BLAST 1.05\n\n");
	}

    // Read the first sequence from FASTA file (the query)
	readFasta_open(parameters_queryFile);
    if (!(readFasta_readSequence()))
	{
		fprintf(stderr, "Error reading query from FASTA file %s\n", argv[1]);
		exit(-1);
	}

    do
    {
    	// Initialize global variables
    	global_initialize();
        
        // Make copy of the sequence
        query = (char*)global_malloc(sizeof(char) * readFasta_sequenceLength + 1);
        strcpy(query, readFasta_sequenceBuffer);

        // Make copy of the description
        queryDescription = (char*)global_malloc(sizeof(char) * (readFasta_descriptionLength + 1));
        blast_queryDescription = (char*)global_malloc(sizeof(char) * (readFasta_descriptionLength + 1));
        strcpy(queryDescription, readFasta_descriptionBuffer);
        strcpy(blast_queryDescription, readFasta_descriptionBuffer);

        // Determine the alphabet type of the query
        queryAlphabetType = encoding_determineAlphabetType(query, strlen(query));

        // If not the same alphabet type as previous query, abort
        if (previousQueryAlphabetType < 10 && previousQueryAlphabetType != queryAlphabetType)
        {
            fprintf(stderr, "Error: Processing sequence %s\n", query);
            fprintf(stderr, "Error: Unable to process a mix of both protein and nucleotide queries\n");
            fflush(stderr);
            exit(-1);
        }
        previousQueryAlphabetType = queryAlphabetType;

        // Initialize encoding
        encoding_initialize(queryAlphabetType);

        // Filter the query using DUST or SEG
        if (parameters_filterEnabled)
        {
            if (queryAlphabetType == encoding_protein)
                seg_segSequence(query);
            else
                dust_dustSequence(query);
        }

        // Load parameter defaults based on query alphabet type
        parameters_loadDefaults(queryAlphabetType);

        queryDescription = print_formatDescription(queryDescription, 7, 0, 70);
	    if (parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
		{
            printf("Query= %s\n", queryDescription);
            printf("         (%u letters)\n\n", strlen(query));
		}

        // Open sequence data file and read information
        readdb_open(parameters_subjectDatabaseFile);

        // If a nucleotide alphabet
        if (queryAlphabetType == encoding_nucleotide)
        {
            // Create a nucleotide scoring matrix use match and mismatch penalties
            scoreMatrix = scoreMatrix_create(parameters_matchScore, parameters_mismatchScore);
    //		scoreMatrix_print(scoreMatrix);

            // Create the PSSMatrix
            PSSMatrix = PSSMatrix_create(scoreMatrix, query);
    //		PSSMatrix_print(PSSMatrix);

            nucleotideLookup_build(PSSMatrix, parameters_wordTableBytes);
    //		nucleotideLookup_print();
        }
        // If a protein alphabet
        else
        {
            // Load the scoring matrix (eg. BLOSUM)
            scoreMatrix = scoreMatrix_load(parameters_scoringMatrixPath);
    //		scoreMatrix_print(scoreMatrix);

            // Create the PSSMatrix
            PSSMatrix = PSSMatrix_create(scoreMatrix, query);
    //		PSSMatrix_print(PSSMatrix);
	
			//Shucai
			//Transform the PSSMatrix layout to make it compatible with GPU
			PSSMatrixFP = PSSMatrixFP_transform(&PSSMatrix);

            // Use query sequence to build the word lookup FSA structure
            if (readdb_numberOfSequences != readdb_numberOfClusters)
                wordLookupDFA_build(PSSMatrix, encoding_sentinalCode, parameters_wordSize);
    		else
            	wordLookupDFA_build(PSSMatrix, encoding_numRegularLetters, parameters_wordSize);
    //Shucai
//			wordLookupDFA_print();
        }

		//Shucai
        //blast_search(parameters_subjectDatabaseFile, PSSMatrix, query);
		blast_search(parameters_subjectDatabaseFile, PSSMatrix, PSSMatrixFP, query);

/*	    if (parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
		{
            printf("Prep=%f\nSearch=%f\nCopySubjects=%f\nFastScore=%f\n",
            (float)blast_prepTime / CLOCKS_PER_SEC, (float)blast_searchTime / CLOCKS_PER_SEC,
            (float)blast_copyTime / CLOCKS_PER_SEC, (float)blast_semiGappedScoreTime / CLOCKS_PER_SEC);

            printf("GappedScore=%f\nUnpack=%f\nGappedExtend=%f\nFinalize=%f\n",
            (float)blast_gappedScoreTime / CLOCKS_PER_SEC, (float)blast_unpackTime / CLOCKS_PER_SEC,
            (float)blast_gappedExtendTime / CLOCKS_PER_SEC, (float)blast_finalizeTime / CLOCKS_PER_SEC);
		}*/

        // Free score matrix, and PSSMatrix columes at the same time
        scoreMatrix_free(scoreMatrix);
        wordLookupDFA_free();
        nucleotideLookup_free();
	    encoding_free();
        free(query); free(queryDescription);
	}
    while (readFasta_readSequence());

	// close FASTA reader
	readFasta_close();

    // Free all global data
    global_free();
    semiGappedScoring_free();
    oldSemiGappedScoring_free();
    oldGappedScoring_free();
    gappedScoring_free();
    nuGappedScoring_free();
    bytepackGappedScoring_free();
    fasterBytepackGappedScoring_free();
	gappedExtension_free();
	fasterGappedExtension_free();
    parameters_free();

	return 0;
}

//Shucai
//void blast_search(char* searchDbFile, struct PSSMatrix PSSMatrix, char* query)
void blast_search(char* searchDbFile, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, char* query)
{
	char *indexFilename;
	int4 tickFrequency;

    // Construct sequence filename
	indexFilename = (char*)global_malloc(strlen(searchDbFile) + 9);
	sprintf(indexFilename, "%s.index", searchDbFile);

	// Check if index file exists. If not, disable use of index
/*	if ((indexFile = fopen(indexFilename, "r")) != NULL)
    	fclose(indexFile);
    else*/
    parameters_useIndex = 0;

    // Check that alphabet type of query and database match
    if (encoding_alphabetType == encoding_protein && readdb_dbAlphabetType == encoding_nucleotide)
    {
		fprintf(stderr, "Error: database %s contains nucleotide sequences\n", searchDbFile);
		fprintf(stderr,
        "Error: searching a nucleotide database with a protein query is not supported\n\n");
		exit(-1);
    }
    if (encoding_alphabetType == encoding_nucleotide && readdb_dbAlphabetType == encoding_protein)
    {
		fprintf(stderr, "Error: database %s contains protein sequences\n", searchDbFile);
		fprintf(stderr,
        "Error: searching a protein database with a nucleotide query is not supported\n\n");
		exit(-1);
    }

	// Determine tick frequence
	tickFrequency = ceil((float)readdb_numberOfSequences / 50.0);

	// Initialize BLAST statistics (calculate log(2), log(K), nominal drop-offs, etc.)
	statistics_initialize(PSSMatrix, readdb_numberOfLetters, readdb_numberOfSequences);

	// Determine the minimum gapped nominal score required for reporting the alignment
	blast_gappedNominalCutoff = statistics_gappedEvalue2nominal(parameters_cutoff);

    // Determine the minimum/maximum semi-gapped scores to achieve cutoff
	blast_nominalR1cutoff = ceil((float)blast_gappedNominalCutoff * parameters_semiGappedR1);
	blast_nominalR2cutoff = ceil((float)blast_gappedNominalCutoff * parameters_semiGappedR2);

	// Determine the minimum ungapped nominal score required to trigger gapping
	if (encoding_alphabetType == encoding_protein)
    {
        blast_ungappedNominalTrigger
            = statistics_ungappedNormalized2nominal(parameters_ungappedNormalizedTrigger);
	}
    else
    {
        blast_ungappedNominalTrigger
        	= statistics_ungappedNucleotideTrigger(PSSMatrix);
	}

	// Gapping trigger cannot be greater than final cutoff
	if (blast_ungappedNominalTrigger > blast_gappedNominalCutoff)
		blast_ungappedNominalTrigger = blast_gappedNominalCutoff;

	// Initialize collections of alignments
	alignments_initialize();

    // Initialize collections of ungapped extensions
    ungappedExtension_initialize();

    if (parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
    {
        printf("Database: %s\n", searchDbFile);
        printf("           %s sequences;", global_int4toString(readdb_numberOfSequences));
        printf(" %s total letters\n\n", global_int8toString(readdb_numberOfLetters));

        #ifndef VERBOSE
        printf("Searching...");
        fflush(stdout);
        #endif
	}

    // Initialize the hitMatrix
	hitMatrix_initialize(PSSMatrix.length, readdb_longestSequenceLength, readdb_sequences);

    blast_prepTime = clock();
	blast_searchTime = -clock();

    while (1)
    {
        // If ssearch mode
        if (parameters_ssearch)
        {
            // Nucleotide search
            if (encoding_alphabetType == encoding_nucleotide)
            {
                search_nucleotideSsearch(PSSMatrix, readdb_sequenceData,
                                        readdb_numVolumeSequences, tickFrequency);
            }
            // Protein Search
            else
            {
                search_proteinSsearch(PSSMatrix, readdb_sequenceData,
                                     readdb_numVolumeSequences, tickFrequency);
            }
        }
        else
        {
            // Nucleotide search
            if (encoding_alphabetType == encoding_nucleotide)
            {
                if (parameters_useIndex)
                {
//                    search_nucleotideIndex(PSSMatrix, readdb_sequenceData,
//                                          readdb_numVolumeSequences, tickFrequency);
                }
                else if (parameters_wordExtraBytes > 0)
                {
                    if (nucleotideLookup_largeTable)
                    {
                        // TODO: add support for long word and large table
                    }
                    else
                    {
                        search_nucleotide_longWord(PSSMatrix, readdb_sequenceData,
                                                  readdb_numVolumeSequences, tickFrequency);
					}
                }
                else
                {
                    if (nucleotideLookup_largeTable)
                    {
                        search_nucleotide_largeTable(PSSMatrix, readdb_sequenceData,
                                                    readdb_numVolumeSequences, tickFrequency);
                    }
                    else
                    {
                        search_nucleotide(PSSMatrix, readdb_sequenceData,
                                         readdb_numVolumeSequences, tickFrequency);
                    }
                }
            }
            // Protein Search
            else
            {
                // Only one hit required to trigger ungapped extension
                if (parameters_oneHitTrigger)
                {
					//Shucai
					//search_protein1hit(PSSMatrix, readdb_sequenceData,
					//				  readdb_numVolumeSequences, tickFrequency);
                    search_protein1hit(PSSMatrix, PSSMatrixFP, readdb_sequenceData,
                                      readdb_numVolumeSequences, tickFrequency);

                }
                // Two hits to trigger an ungapped extensions
                else
                {	
					//Shucai
                    //search_protein2hit(PSSMatrix, readdb_sequenceData,
                    //                  readdb_numVolumeSequences, tickFrequency);
                	search_protein2hit(PSSMatrix, PSSMatrixFP, readdb_sequenceData,
									  readdb_numVolumeSequences, tickFrequency);
				}
            }
        }

        if (readdb_volume + 1 < readdb_numberOfVolumes)
        {
            #ifndef NO_STAGE3
            // Before loading next volume, perform initial semi-gapped or bytepacked alignment
            // on high-scoring ungapped extensions in this volume
            blast_searchTime += clock();
            blast_semiGappedScoreTime -= clock();
            alignments_findGoodAlignments(PSSMatrix, PSSMatrixFP);
            blast_semiGappedScoreTime += clock();
            blast_searchTime -= clock();
            #endif

            // Copy subject sequences from good alignments into memory
            blast_searchTime += clock();
            blast_copyTime -= clock();
            alignments_loadSubjectsIntoMemory(PSSMatrix);
            blast_copyTime += clock();
            blast_searchTime -= clock();

            // Load the next volume
            readdb_nextVolume();

            // Re-initialize the hitMatrix
            hitMatrix_reinitialize(PSSMatrix.length, readdb_longestSequenceLength, readdb_sequences);
        }
        else
		{
        	break;
		}
	}

    if (parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
    {
        #ifndef VERBOSE
        printf("done.\n\n\n\n");
        fflush(stdout);
        #endif
	}

    blast_searchTime += clock();

//	blast_compareScorings(PSSMatrix);
//	exit(0);
	
    #ifndef NO_STAGE3
    if (!parameters_ssearch)
    {	
		
    	// Perform semi-gapped / bytepacked alignment to find good alignments
		blast_semiGappedScoreTime -= clock();
		//Shucai
    	alignments_findGoodAlignments(PSSMatrix, PSSMatrixFP);
    	blast_semiGappedScoreTime += clock();

        // Perform gapped alignment to find final alignments
        blast_gappedScoreTime -= clock();
        alignments_findFinalAlignments(PSSMatrix);
        blast_gappedScoreTime += clock();
	}
    #endif

    #ifndef NO_STAGE4
    blast_gappedExtendTime -= clock();
    // Read the final alignment subject descriptions
    alignments_getFinalAlignmentDescriptions();

    if (!parameters_ssearch)
    {
        // Find traceback information
		alignments_getTracebacks(PSSMatrix, PSSMatrixFP);
	}
    blast_gappedExtendTime += clock();
    #endif

    blast_finalizeTime -= clock();

	// Print alignments
    if (alignments_finalAlignments->numEntries == 0 && parameters_outputType != parameters_xml
        && parameters_outputType != parameters_tabular)
    {
    	printf("\n ***** No hits found ******\n");
    }
    else
    {
        #ifndef NO_STAGE4
        if (parameters_outputType == parameters_xml)
        {
        	print_XMLheader(query, PSSMatrix);
            print_gappedAlignmentsFull(query, PSSMatrix);
        	print_XMLfooter();
        }
        else if (parameters_outputType == parameters_tabular)
        {
            print_gappedAlignmentsFull(query, PSSMatrix);
        }
        else
        {
            print_gappedAlignmentsBrief();
            print_gappedAlignmentsFull(query, PSSMatrix);
		}
        #endif
    }

    if (parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
    {
        if (readdb_numberOfVolumes > 0)
            printf("  Database: %s  (%d volumes)\n", searchDbFile, readdb_numberOfVolumes);
        else
            printf("  Database: %s\n", searchDbFile);
    //    printf("    Posted date:  Apr 5, 2004  5:12 PM\n");
        printf("  Number of letters in database: %s\n", global_int8toString(statistics_databaseSize));
        printf("  Number of sequences in database:  %u\n", readdb_numberOfSequences);

        printf("\nLambda     K      H     (ungapped)");
        printf("\n %.3f     %.3f  %.3f", statistics_ungappedLambda, statistics_ungappedK,
                                         statistics_ungappedH);
        printf("\n\nLambda     K      H     (gapped)");
        printf("\n %.3f     %.3f  %.3f", statistics_gappedParams.lambda, statistics_gappedParams.K,
                                         statistics_gappedParams.H);
        printf("\n\n\nMatrix: %s", parameters_scoringMatrix);
        printf("\nGap Penalties: Existence: %d, Extension: %d", parameters_startGap,
                                                                parameters_extendGap);
        if ((parameters_semiGappedScoring || parameters_bytepackedScoring) && !parameters_ssearch)
            printf("\nSemi-Gapped Gap Penalties: Existence: %d, Extension: %d",
                   parameters_semiGappedStartGap, parameters_semiGappedExtendGap);
        if (!parameters_ssearch)
            printf("\nNumber of Hits to DB: %s", global_int4toString(blast_numHits));
        printf("\nNumber of Sequences: %u\n", readdb_numberOfSequences);
        if (!parameters_ssearch)
        {
            printf("Number of extensions: %u\n", blast_numUngappedExtensions);
            printf("Number of successful extensions: %u\n", blast_numTriggerExtensions);
            printf("Number of sequences with successful extensions: %u\n", blast_numTriggerSequences);
        }
        if ((parameters_semiGappedScoring || parameters_bytepackedScoring || parameters_tableScoring)
            && !parameters_ssearch)
            printf("Number of sequences with semi-gapped score above cutoff: %u\n",
                   blast_numGoodAlignments);
        printf("Number of sequences better than %g: %u\n",
               parameters_cutoff, alignments_finalAlignments->numEntries);
        if (!parameters_ssearch)
        {
            if (parameters_semiGappedScoring || parameters_bytepackedScoring || parameters_tableScoring)
                printf("Number of HSP's that attempted semi-gapping: %u\n", blast_numSemiGapped);
            printf("Number of HSP's that attempted gapping: %u\n", blast_numGapped);
            printf("Number of HSP's contained and not gapped: %u\n", blast_numExtensionsPruned);
            printf("Number of HSP's succeeded/attempted join: %u/%u\n",
                blast_numSuccessfullyJoined, blast_numAttemptedJoin);
        }
        if (blast_numExpandedSequences)
            printf("Number of cluster members recreated = %d\n", blast_numExpandedSequences);
        printf("Total subject bytes copied/unpacked = %d/%d\n", blast_totalCopied, blast_totalUnpacked);
        printf("length of query: %u\n", statistics_querySize);
        printf("length of database: %s\n", global_int8toString(statistics_databaseSize));
        printf("effective HSP length: %u\n", statistics_lengthAdjust);
        printf("effective length of query: %u\n", statistics_effectiveQuerySize);
        printf("effective length of database: %s\n",
            global_int8toString(statistics_effectiveDatabaseSize));
        printf("effective search space: %llu\n", statistics_searchSpaceSize);
        printf("effective search space used: %llu\n", statistics_searchSpaceSize);

        if (encoding_alphabetType == encoding_protein)
        {
            printf("T: %d\n", parameters_T);
            printf("A: %d\n", parameters_A);
        }
        printf("X1: %d\n", statistics_ungappedNominalDropoff);
        printf("X2: %d\n", statistics_gappedNominalDropoff);
        printf("X3: %d\n", statistics_gappedFinalNominalDropoff);
        printf("S1: %d\n", blast_ungappedNominalTrigger);
        printf("S2: %d\n", blast_gappedNominalCutoff);
        if (blast_dynamicGappedNominalCutoff > 0)
            printf("S3: %d\n", blast_dynamicGappedNominalCutoff);
        printf("F2: %d\n", blast_nominalR1cutoff);
        if (blast_dynamicNominalR1cutoff > 0)
            printf("F3: %d\n", blast_dynamicNominalR1cutoff);

//    	printf("Total malloced=%s\n", global_int4toString(global_totalMalloc));
	}

	// Free memory used by hitMatrix, PSSMatrix, alignments and sequence filename
	hitMatrix_free();
	alignments_free();
	PSSMatrix_free(PSSMatrix, PSSMatrixFP);
    readdb_close();

	blast_finalizeTime += clock();
}


