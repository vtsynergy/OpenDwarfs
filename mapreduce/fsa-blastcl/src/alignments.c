// alignment.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Functions for processing and rescoring a collection of high-scoring alignments
// alignments_alignments contains all ungapped alignments scoring above S1
// alignments_goodAlignments contains semi-gapped alignments scoring above F2
// alignments_finalAlignments contains gapped alignments scoring above S2

#include "blast.h"

struct memBlocks* alignments_alignments;
struct memSingleBlock* alignments_goodAlignments;
struct memSingleBlock* alignments_finalAlignments;
struct alignment* alignments_currentAlignment;
//uint4 alignments_numClusters;


struct memBlocks* alignments_volumeAlignments[constants_maxNumVolumes];
uint4 alignments_numVolumes;
int alignments_finalAlignmentsSorted;
uint4 alignments_numClusters = 0;

// Expand a cluster and perform stages 1 and 2 of search on the children
int alignments_expandCluster(struct alignment* alignment, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP);

// Initialize array storing pointers to alignments
void alignments_initialize()
{
	// Initialize alignments, good alignments and final alignments blocks
    alignments_alignments = memBlocks_initialize(sizeof(struct alignment),
                            constants_initialAllocAlignments);
	alignments_goodAlignments = memSingleBlock_initialize(sizeof(struct finalAlignment),
                                constants_initialAllocGoodAlignments);
	alignments_finalAlignments = memSingleBlock_initialize(sizeof(struct finalAlignment),
                                 constants_initialAllocFinalAlignments);

    // Initialize code for processing regions of a subject
    //Shucai
	//unpack_initialize(blast_numGoodExtensions);
    unpack_initialize();

    alignments_currentAlignment = NULL;
	alignments_numVolumes = 0;
    alignments_finalAlignmentsSorted = 0;
}

// Use the next available slot in the alignment object array
void alignments_createNew(uint4 descriptionLocation, uint4 descriptionLength,
                          unsigned char* subject, int4 subjectLength,
                          int4 encodedLength)
{
	struct alignment* alignment;

	// Get slot for new alignment
	alignment = (struct alignment*)memBlocks_newEntry(alignments_alignments);

	// Create/initialize contents
	alignment->ungappedExtensions = NULL;
	alignment->gappedExtensions = NULL;
	alignment->descriptionLocation = descriptionLocation;
    alignment->descriptionLength = descriptionLength;
	alignment->subject = subject;
	alignment->subjectLength = subjectLength;
	alignment->joinChecked = 0;
    alignment->encodedLength = encodedLength;
    alignment->inMemorySubject = 0;
    alignment->unpackRegions = NULL;
	alignment->numUnpackRegions = 0;
    alignment->edits = NULL;
    alignment->cluster = 0;

    // Record pointer to wilcard edits if there are any for this subject
    if (encoding_alphabetType == encoding_nucleotide)
    {
        alignment->edits = subject + ((alignment->subjectLength + 3) / 4);
        if (alignment->edits == subject + encodedLength)
        {
            alignment->edits = NULL;
		}
    }

	alignments_currentAlignment = alignment;
    blast_numTriggerSequences++;
}

// Move contents of list of good alignments, to list of final alignments
void alignments_moveGoodToFinal()
{
	struct memSingleBlock* tempAlignments;

    // Switch good and final
	tempAlignments = alignments_goodAlignments;
	alignments_goodAlignments = alignments_finalAlignments;
	alignments_finalAlignments = tempAlignments;

    // Reset number of good
    alignments_goodAlignments->numEntries = 0;
}

// Add an ungapped extension at the end of the list
void alignments_addUngappedExtensionAtEnd(struct alignment* alignment,
                                          struct ungappedExtension* newExtension)
{
	struct ungappedExtension *currentExtension;

    newExtension->next = NULL;

    // List is empty
    if (alignment->ungappedExtensions == NULL)
	{
        alignment->ungappedExtensions = newExtension;
        return;
	}

	// Start at beginning of ungapped extensions and traverse to end
	currentExtension = alignment->ungappedExtensions;

    while (currentExtension->next != NULL)
    {
    	currentExtension = currentExtension->next;
    }

    // Add at end
    currentExtension->next = newExtension;
}

// Add an ungapped extension to an alignment's list of ungapped extensions, which
// are ordered highest nominal score first, lowest score last
void alignments_addUngappedExtension(struct ungappedExtension* newExtension)
{
	struct ungappedExtension *currentExtension, *previousExtension;

	// Start at beginning of ungapped extensions (one with highest score)
	currentExtension = alignments_currentAlignment->ungappedExtensions;

	// If there are none, add first/sole extension
	if (currentExtension == NULL)
	{
		alignments_currentAlignment->ungappedExtensions = newExtension;
	}
	else
	{
		previousExtension = NULL;

		// Else move through list of existing extensions until we either reach the
		// end or reach one with a score less than newExtension
		while (currentExtension != NULL &&
		      (currentExtension->nominalScore > newExtension->nominalScore))
		{
			previousExtension = currentExtension;
			currentExtension = currentExtension->next;
		}

		if (previousExtension == NULL)
		{
			// This is the highest scoring alignment, insert at front of
			// the queue
			alignments_currentAlignment->ungappedExtensions = newExtension;
			newExtension->next = currentExtension;
		}
		else
		{
			// Insert between higher and lower scoring extensions
			previousExtension->next = newExtension;
			newExtension->next = currentExtension;
		}
	}
}

// Add a high-scoring gapped extension to this alignment
void alignments_addGappedExtension(struct alignment* alignment, struct gappedExtension* newExtension)
{
	struct gappedExtension *currentExtension, *previousExtension;

	// If this is the first high-scoring gapped extension for this alignment
	if (alignment->gappedExtensions == NULL)
	{
		// Make this the first gapped extension in the alignment
		alignment->gappedExtensions = newExtension;
	}
	else
	{
		// Start at beginning of list of gapped extensions (one with highest score)
		currentExtension = alignment->gappedExtensions;
		previousExtension = NULL;

		// Move through list of existing extensions until we either reach the
		// end or reach one with a score less than newExtension
		while (currentExtension != NULL &&
		      (currentExtension->nominalScore > newExtension->nominalScore))
		{
			previousExtension = currentExtension;
			currentExtension = currentExtension->next;
		}

		if (previousExtension == NULL)
		{
			// This is the highest scoring extension, insert at front of the queue
			alignment->gappedExtensions = newExtension;
			newExtension->next = currentExtension;
		}
		else
		{
			// Insert between higher and lower scoring extensions
			previousExtension->next = newExtension;
			newExtension->next = currentExtension;
		}
	}
}

// Remove given extension from the list of gapped extensions associated
// with the given alignment
void alignments_removeGappedExtension(struct alignment* alignment,
                                      struct gappedExtension* toRemove)
{
	struct gappedExtension *currentExtension, *previousExtension;

	currentExtension = alignment->gappedExtensions;
	previousExtension = NULL;

	// For each extension in the list
	while (currentExtension != NULL)
	{
		// If matches extension to remove
		if (currentExtension == toRemove)
		{
			// Remove it
			if (previousExtension == NULL)
			{
				// It was the at the start of the list
				alignment->gappedExtensions = currentExtension->next;
				free(currentExtension->trace.traceCodes);
				free(currentExtension);
				return;
			}
			else
			{
				// It was mid or end-list
				previousExtension->next = currentExtension->next;
				free(currentExtension->trace.traceCodes);
				free(currentExtension);
				return;
			}
		}
		previousExtension = currentExtension;
		currentExtension = currentExtension->next;
	}
}

// Given an alignment with a list of gapped extensions, remove the lowest scoring of
// any pair of extensions that share either the same start or end point
void alignments_pruneOverlappingExtensions(struct alignment* alignment)
{
	struct gappedExtension *extension, *extension2;

	// For each gapped extension in the alignment
	extension = alignment->gappedExtensions;
	while (extension != NULL)
	{
		// For each extension after that one
		extension2 = extension->next;
		while (extension2 != NULL)
		{
			// If the two extensions start or end in the same place
			if ((extension->trace.queryStart == extension2->trace.queryStart &&
			     extension->trace.subjectStart == extension2->trace.subjectStart) ||
			    (extension->queryEnd == extension2->queryEnd &&
			     extension->subjectEnd == extension2->subjectEnd))
			{
				// Remove the lower scoring one
				if (extension->nominalScore > extension2->nominalScore)
				{
					// Remove extension two
					alignments_removeGappedExtension(alignment, extension2);
					// Go back to start of extensions following extension one
					extension2 = extension;
				}
				else
				{
					// Remove extension one
					alignments_removeGappedExtension(alignment, extension);
					// Go back to start
					extension = alignment->gappedExtensions;
					extension2 = extension;
				}
			}

			if (extension2 != NULL)
				extension2 = extension2->next;
		}
		if (extension != NULL)
			extension = extension->next;
	}
}

// Returns true if extension2 is contained in the rectangular region defined by extension2
int alignments_contains(struct ungappedExtension* extension1, struct ungappedExtension* extension2)
{
    if (extension2->start.queryOffset >= extension1->start.queryOffset &&
        extension2->end.queryOffset <= extension1->end.queryOffset &&
        extension2->start.subjectOffset >= extension1->start.subjectOffset &&
        extension2->end.subjectOffset <= extension1->end.subjectOffset &&
        extension2->nominalScore <= extension1->nominalScore &&
        extension1 != extension2)
	{
    	return 1;
    }
    else
    {
    	return 0;
    }
}

// Check for ungapped extensions that were mistakenly pruned. If one is found, add a non-deleted
// copy to the end of list of ungapped extensions
void alignments_unpruneRegion(struct alignment* alignment, struct ungappedExtension* oldUngappedExtension,
	struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP)
{
	struct ungappedExtension *currentExtension, *previousExtension, *restoreExtension;

    #ifdef VERBOSE
    if (alignment->descriptionLocation == parameters_verboseDloc)
    {
        printf("unpruneRegion dloc=%d oldRegion:", alignment->descriptionLocation);
		printf("Pruning old: "); ungappedExtension_print(oldUngappedExtension);
        printf("                      newRegion:");
		printf("Pruning new: "); ungappedExtension_print(ungappedExtension);
	}
    #endif

	currentExtension = alignment->ungappedExtensions;
	previousExtension = NULL;
	while (currentExtension != NULL)
	{
		// If old table-driven alignment contains extension but new gapped alignment
        // does not, undelete it
		if (currentExtension->status == ungappedExtension_DELETED &&
            alignments_contains(oldUngappedExtension, currentExtension) &&
            !alignments_contains(ungappedExtension, currentExtension))
		{
        	restoreExtension = currentExtension;

            #ifdef VERBOSE
            if (alignment->descriptionLocation == parameters_verboseDloc)
            {
                printf("Restoring: ");
                ungappedExtension_print(restoreExtension);
			}
            #endif

            // Removing from front of list
			if (previousExtension == NULL)
			{
				alignment->ungappedExtensions = currentExtension->next;
				currentExtension = alignment->ungappedExtensions;
			}
			// Removing from the middle of the list
			else
			{
				previousExtension->next = currentExtension->next;
				currentExtension = previousExtension->next;
			}

            // Add the restored extension to end of the list and change status
            alignments_addUngappedExtensionAtEnd(alignment, restoreExtension);
            restoreExtension->status = ungappedExtension_UNGAPPED;

            // Find the seed
            ungappedExtension_findSeed(restoreExtension, PSSMatrix, PSSMatrixFP, alignment->subject);
        }
        else
        {
			// Advance to next ungapped extension
			previousExtension = currentExtension;
			currentExtension = currentExtension->next;
        }
    }
}

// Remove any ungapped extensions in this alignment that are in the area covered
// by the given gapped scored extension
void alignments_pruneRegion(struct alignment* alignment, struct ungappedExtension* ungappedExtension)
{
	struct ungappedExtension *currentExtension;

    #ifdef VERBOSE
    if (alignment->descriptionLocation == parameters_verboseDloc)
        printf ("pruneRegion ungappedExtension query=%d to %d subject=%d to %d score=%d dloc=%d status=%d\n",
                ungappedExtension->start.queryOffset,
                ungappedExtension->end.queryOffset,
                ungappedExtension->start.subjectOffset,
                ungappedExtension->end.subjectOffset,
                ungappedExtension->nominalScore, alignment->descriptionLocation,
                ungappedExtension->status);
	#endif

	currentExtension = alignment->ungappedExtensions;
	while (currentExtension != NULL)
	{
		// If currentExtension start or end is contained within area bound by ungappedExtension
		// (after gappedScoring ungappedExtension), then remove it
		if (currentExtension->status != ungappedExtension_DELETED &&
            alignments_contains(ungappedExtension, currentExtension))
		{
            #ifdef VERBOSE
                if (alignment->descriptionLocation == parameters_verboseDloc)
                printf("REMOVING %d query %d to %d subject %d to %d\n",
                currentExtension->nominalScore, currentExtension->start.queryOffset,
                currentExtension->end.queryOffset,
                currentExtension->start.subjectOffset,
                currentExtension->end.subjectOffset);
            #endif

            blast_numExtensionsPruned++;

			currentExtension->status = ungappedExtension_DELETED;
		}

        // Advance to next ungapped extension
        currentExtension = currentExtension->next;
	}
}

// Add the current alignment (which contains at least one gapped extension
// with semi-gapped score above semi-gapped cutoff) to to-be-sorted list of good alignments
void alignments_addGoodAlignment(int4 highestNominalScore, struct alignment* alignment)
{
	struct finalAlignment* goodAlignment;
    // Get a new good alignment entry
    goodAlignment = (struct finalAlignment *)memSingleBlock_newEntry(alignments_goodAlignments);

	// Insert alignment information into the new good alignment
	goodAlignment->highestNominalScore = highestNominalScore;
	goodAlignment->alignment = alignment;
}

// Add the current alignment (which contains at least one gapped extension
// scoring above cutoff) to to-be-sorted list of final alignments
void alignments_addFinalAlignment(int4 highestNominalScore, struct alignment* alignment)
{
	struct finalAlignment* finalAlignment;

    // Get a new final alignment entry
    finalAlignment = (struct finalAlignment *)memSingleBlock_newEntry(alignments_finalAlignments);

	// Insert alignment information into the new final alignment
	finalAlignment->highestNominalScore = highestNominalScore;
	finalAlignment->alignment = alignment;
	finalAlignment->description = NULL;
}

// Compare the two alignments' highest scores. Return -1 if alignment1 < alignment2
// 1 if alignment1 > alignment2 and 0 if they are equal
int4 alignments_compareFinalAlignments(const void* alignment1,
                                      const void* alignment2)
{
	const struct finalAlignment *a1, *a2;

	a1 = (struct finalAlignment*)alignment1;
	a2 = (struct finalAlignment*)alignment2;

	if (a1->highestNominalScore > a2->highestNominalScore)
	{
		return -1;
	}
	else if (a1->highestNominalScore < a2->highestNominalScore)
	{
		return 1;
	}
	else
	{
    	// Resolve conflicts using subject length
        if (a1->alignment->subjectLength > a2->alignment->subjectLength)
        	return 1;
        else if (a1->alignment->subjectLength < a2->alignment->subjectLength)
        	return -1;
        else
        	return 0;
	}
}

// Sort the array of good alignments in order of score
void alignments_sortGoodAlignments()
{
	qsort(alignments_goodAlignments->block, alignments_goodAlignments->numEntries,
          sizeof(struct finalAlignment), alignments_compareFinalAlignments);
}

// Sort the array of final alignments in order of score
void alignments_sortFinalAlignments()
{
	qsort(alignments_finalAlignments->block, alignments_finalAlignments->numEntries,
          sizeof(struct finalAlignment), alignments_compareFinalAlignments);
}

// Select the region for this alignment and return the subject
unsigned char* alignments_selectRegion(struct alignment* alignment,
               struct ungappedExtension* ungappedExtension)
{
	struct unpackRegion* unpackRegion;

	if (alignment->unpackRegions == NULL)
    {
		return alignment->subject;
	}
    else
    {
    	unpackRegion = unpack_selectRegion(alignment->unpackRegions, alignment->numUnpackRegions,
                                           ungappedExtension->seed.subjectOffset);

//        printf("Select region start=%d end=%d/%d seed=%d\n", unpackRegion->startOffset,
//        unpackRegion->endOffset, alignment->subjectLength, ungappedExtension->seed.subjectOffset); fflush(stdout);

    	return unpackRegion->subject;
	}
}

// Check other ungapped extensions that have been scored to see if we can join
// it to this one to get a higher scoring HSP
void alignments_checkForJoin(struct alignment* alignment, struct ungappedExtension* extension1,
                             struct PSSMatrix PSSMatrix)
{
	struct ungappedExtension* extension2;
	int4 queryDistance, subjectDistance, newScore;
	unsigned char* subject;

	subject = alignments_selectRegion(alignment, extension1);

	extension2 = alignment->ungappedExtensions;
	blast_dloc = alignment->descriptionLocation;

    // For each extension that has already been gapped scored
	while (extension2 != NULL)
	{
		// Check extension2 has been scored, and that the combined scores
		// of extension 1&2 could possibly exceed the cutoff
		if (extension2 != extension1 && extension2->status != ungappedExtension_DELETED)
		{
			// If extension2 comes after extension1, determine distance between them
			if (extension1->start.queryOffset < extension2->start.queryOffset)
			{
				queryDistance = extension2->start.queryOffset - extension1->end.queryOffset;
				subjectDistance = extension2->start.subjectOffset - extension1->end.subjectOffset;
			}
			// Else extension1 comes after extension2
			else
			{
				queryDistance = extension1->start.queryOffset - extension2->end.queryOffset;
				subjectDistance = extension1->start.subjectOffset - extension2->end.subjectOffset;
			}

            #ifdef VERBOSE
            if (parameters_verboseDloc == alignment->descriptionLocation)
            {
            	printf("Check for join:\n");
				ungappedExtension_print(extension1);
				ungappedExtension_print(extension2);
			}
            #endif

			// Quick check the distance between HSPs is not too great
			if ((queryDistance >= 0 || subjectDistance >= 0) &&
			    minimum(abs(queryDistance), abs(subjectDistance)) + abs(queryDistance - subjectDistance)
                * parameters_extendGap <= 2 * statistics_gappedFinalNominalDropoff)
			{
                // Perform gappedScoring with higher dropoff to try and bridge the gap
                newScore
                = gappedScoring_score(extension1, PSSMatrix, alignment->subjectLength,
                                      subject, statistics_gappedFinalNominalDropoff);

                #ifdef VERBOSE
                if (parameters_verboseDloc == alignment->descriptionLocation)
                {
                	printf("Performed join. New score=%d\n", newScore);
                }
                #endif

                blast_numAttemptedJoin++;
				// Check if we successfully joined two HSPs
				if (newScore > extension1->nominalScore)
				{
					extension1->nominalScore = newScore;
					extension1->status = ungappedExtension_JOINED;
					extension2->status = ungappedExtension_JOINED;
					blast_numSuccessfullyJoined++;
				}

				return;
			}
		}
		extension2 = extension2->next;
	}
}

// Compare the two alignments' description locations scores.
int4 alignments_compareAlignmentDescriptionLocations(const void* alignment1, const void* alignment2)
{
	const struct finalAlignment *a1, *a2;

	a1 = (struct finalAlignment*)alignment1;
	a2 = (struct finalAlignment*)alignment2;

	if (a1->alignment->descriptionLocation > a2->alignment->descriptionLocation)
	{
		return 1;
	}
	else if (a1->alignment->descriptionLocation < a2->alignment->descriptionLocation)
	{
		return -1;
	}
	else
	{
		return 0;
    }
}

// Get the subject sequence descriptions for all of the final alignments
void alignments_getFinalAlignmentDescriptions()
{
	struct finalAlignment* finalAlignment;

    // Sort descriptions in order of description location (to speed up disk access)
	qsort(alignments_finalAlignments->block, alignments_finalAlignments->numEntries,
          sizeof(struct finalAlignment), alignments_compareAlignmentDescriptionLocations);

    // For each alignment, read the description
    memSingleBlock_resetCurrent(alignments_finalAlignments);
    while ((finalAlignment = (struct finalAlignment *)memSingleBlock_getCurrent(alignments_finalAlignments)) != NULL)
    {
		finalAlignment->description
        	= descriptions_getDescription(finalAlignment->alignment->descriptionLocation,
                                          finalAlignment->alignment->descriptionLength);
    }

    // Re-sort the final alignments by score
	qsort(alignments_finalAlignments->block, alignments_finalAlignments->numEntries,
          sizeof(struct finalAlignment), alignments_compareFinalAlignments);
}

// Print list of final alignments
void alignments_printFinalAlignments()
{
	int4 count = 0;
	struct alignment* alignment;
	struct finalAlignment* finalAlignment;

	while (count < alignments_finalAlignments->numEntries)
	{
		finalAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_finalAlignments, count);
        alignment = finalAlignment->alignment;

		printf("%3d) Score=%d Dloc=%d\n", count, finalAlignment->highestNominalScore, alignment->descriptionLocation);
		count++;
	}
}

// Print list of good alignments
void alignments_printGoodAlignments()
{
	int4 count = 0;
	struct alignment* alignment;
	struct finalAlignment* goodAlignment;

	while (count < alignments_goodAlignments->numEntries)
	{
		goodAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_goodAlignments, count);
        alignment = goodAlignment->alignment;

		printf("%3d) Score=%d Dloc=%d\n", count, goodAlignment->highestNominalScore, alignment->descriptionLocation);
		count++;
	}
}

// Free memory used by an alignment to store tracebacks
void alignment_freeAlignment(struct alignment* alignment)
{
	struct gappedExtension *gappedExtension, *nextGappedExtension;

    // Free each gapped extension
    gappedExtension = alignment->gappedExtensions;
    while (gappedExtension != NULL)
    {
        nextGappedExtension = gappedExtension->next;
        // Free the trace codes first
        free(gappedExtension->trace.traceCodes);
        // Free the gapped extension
        free(gappedExtension);
        gappedExtension = nextGappedExtension;
    }

    alignment->gappedExtensions = NULL;

    // Free any in-memory subject sequences
    if (alignment->inMemorySubject)
    {
        free(alignment->subject - 1); alignment->subject = NULL;
        free(alignment->edits); alignment->edits = NULL;
    }
}

// Free memory used to store alignments
void alignments_free()
{
	struct alignment* alignment;
	struct finalAlignment* finalAlignment;
	uint4 volumeCount;

    memBlocks_resetCurrent(alignments_alignments);

    // Free unused block of alignments
	memBlocks_free(alignments_alignments);

    // For each volume
	volumeCount = 0;
	while (volumeCount < alignments_numVolumes)
	{
		alignments_alignments = alignments_volumeAlignments[volumeCount];

		// For each alignment
        memBlocks_resetCurrent(alignments_alignments);
		while ((alignment = (struct alignment*)memBlocks_getCurrent(alignments_alignments)) != NULL)
		{
			alignment_freeAlignment(alignment);
		}

		// Free blocks of alignments
		memBlocks_free(alignments_alignments);

		volumeCount++;
	}

    // Free memory used by ungapped extensions
	memBlocks_free(ungappedExtension_extensions);

    // Free good alignments
	memSingleBlock_free(alignments_goodAlignments);

    // For each final alignment, free the description
    memSingleBlock_resetCurrent(alignments_finalAlignments);
    while ((finalAlignment = (struct finalAlignment *)memSingleBlock_getCurrent(alignments_finalAlignments)) != NULL)
    {
		free(finalAlignment->description);
    }

    // Free final alignments
	memSingleBlock_free(alignments_finalAlignments);

	// Free all the unpacked regions
	unpack_free();
}

// Perform initial scoring of all ungapped extensions to find "good" alignments that may
// score above the cutoff
void alignments_findGoodAlignments(struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP)
{
	struct alignment* alignment;
	struct ungappedExtension* ungappedExtension;
	int4 bestScore, numExtensions, hasChildren;

	//debug=========================
	int count = 0;
	//==============================
	
    // For each alignment
    memBlocks_resetCurrent(alignments_alignments);
	while ((alignment = (struct alignment*)memBlocks_getCurrent(alignments_alignments)) != NULL)
	{
    	bestScore = 0;
		blast_dloc = alignment->descriptionLocation;

		//debug=======================
		//printf("count = %d\n", count++);
		//===============================

        // Record if subject has children
        if (encoding_alphabetType == encoding_protein &&
            alignment->encodedLength > alignment->subjectLength + 2)
        	hasChildren = 1;
		else
        	hasChildren = 0;

    	// For each ungapped extension (in descending order of score)
        numExtensions = 0;
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL)
        {
        	if (ungappedExtension->status != ungappedExtension_DELETED)
            {
                // Find the seed
                ungappedExtension_findSeed(ungappedExtension, PSSMatrix, PSSMatrixFP, alignment->subject);

                // Byte-packed scoring
                if (parameters_bytepackedScoring)
                {
                    ungappedExtension->nominalScore
                        = bytepackGappedScoring_score(ungappedExtension, PSSMatrix, alignment->subjectLength,
                           alignment->subject, statistics_gappedNominalDropoff);

                    // Mark as semigapped
                    ungappedExtension->status = ungappedExtension_SEMIGAPPED;
                }
                // Table driven scoring
                else if (parameters_tableScoring)
                {
                    ungappedExtension->nominalScore
                        = tableGappedScoring_score(ungappedExtension, PSSMatrix, alignment->subjectLength,
                           alignment->subject, statistics_gappedNominalDropoff);

                    // Mark as semigapped
                    ungappedExtension->status = ungappedExtension_SEMIGAPPED;
                }
                // Semi-gapped scoring
                else if (parameters_semiGappedScoring)
                {
                    ungappedExtension->nominalScore
                        = semiGappedScoring_score(ungappedExtension, PSSMatrix, alignment->subjectLength,
                                                  alignment->subject, statistics_gappedNominalDropoff);

                    // Mark as semigapped
                    ungappedExtension->status = ungappedExtension_SEMIGAPPED;
                }
                // Regular gapped scoring
                else
                {
                    ungappedExtension->nominalScore
                        = gappedScoring_score(ungappedExtension, PSSMatrix, alignment->subjectLength,
                                              alignment->subject, statistics_gappedNominalDropoff);

                    // Mark as gapped
                    ungappedExtension->status = ungappedExtension_GAPPED;
                }
                // If alignment scores above R1 cutoff
                if (ungappedExtension->nominalScore >= blast_nominalR1cutoff)
                {
                    if (hasChildren)
                    {
                    	// Subject has children so perform stage1 and 2 on children
                    	alignments_expandCluster(alignment, PSSMatrix, PSSMatrixFP);
						bestScore = 0;
                        break;
					}
                    else if (ungappedExtension->nominalScore > bestScore)
                    {
                        // Update best score for the alignment
                        bestScore = ungappedExtension->nominalScore;
                    }
                }
                else
                {
                    // Else mark it as deleted
                    ungappedExtension->status = ungappedExtension_DELETED;
                }

                // Remove any ungapped extensions in this alignment that are in the area covered
                // by the gapped scoring just performed
                alignments_pruneRegion(alignment, ungappedExtension);

                numExtensions++;
		}
            ungappedExtension = ungappedExtension->next;
	}
        // If this alignment contains gapped extensions that could score above cutoff
        if (bestScore >= blast_nominalR1cutoff)
        {
            // If a single sequence add to list of "good" alignments
            alignments_addGoodAlignment(bestScore, alignment);

            blast_numGoodExtensions += numExtensions;
            blast_numGoodAlignments++;
        }
    }

	// Record point to list of alignments for this volume
	alignments_volumeAlignments[alignments_numVolumes] = alignments_alignments;
	alignments_numVolumes++;

	// Construct new list for next volume (if there is one)
    alignments_alignments = memBlocks_initialize(sizeof(struct alignment),
                            constants_initialAllocAlignments);
}

// Copy subject sequences associated with good alignments into memory before
// the volume is closed
void alignments_loadSubjectsIntoMemory(struct PSSMatrix PSSMatrix)
{
	struct finalAlignment *goodAlignment, *finalAlignment;
	struct alignment* alignment;
    int4 minimumDynamicCutoff, position, totalLengths = 0, totalCopied = 0;

    // First dismiss good alignments that can't make it into the top N
	if (parameters_numDisplayAlignments > 0 &&
	    alignments_goodAlignments->numEntries > parameters_numDisplayAlignments)
    {
        // Sort good alignments by score
        alignments_sortGoodAlignments();

        // Get the Nth alignment
        goodAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_goodAlignments, parameters_numDisplayAlignments);

        // Calculate the minimum possible dynamic cutoff
        minimumDynamicCutoff = goodAlignment->highestNominalScore / parameters_semiGappedR2
                             * parameters_semiGappedR1;

		// Move through good alignments
        position = parameters_numDisplayAlignments;
        while (position < alignments_goodAlignments->numEntries)
        {
        	// Stop when score is below minimum dynamic cutoff
            goodAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_goodAlignments, position);
            if (goodAlignment->highestNominalScore < minimumDynamicCutoff)
			{
//            	printf("Changed numGoodAlignments from %d to %d\n",
//					alignments_goodAlignments->numEntries, position); fflush(stdout);
				alignments_goodAlignments->numEntries = position;
            	break;
			}
            position++;
		}
	}

    // For each good alignment (for BLAST)
    position = 0;
    while (position < alignments_goodAlignments->numEntries)
    {
        goodAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_goodAlignments, position);
        alignment = goodAlignment->alignment;

        // Load subject sequence into memory if not already loaded
		if (alignment->inMemorySubject == 0)
        {
        	totalCopied += unpack_loadSubject(PSSMatrix, alignment);
	        totalLengths += alignment->encodedLength;
        }

		position++;
	}

    // For each final alignment (for protein SSEARCH)
    position = 0;
    while (position < alignments_finalAlignments->numEntries)
    {
        finalAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_finalAlignments, position);
        alignment = finalAlignment->alignment;

        // Load subject sequence into memory if not already loaded
		if (alignment->inMemorySubject == 0)
            unpack_loadSubject(PSSMatrix, alignment);

		position++;
	}

//    printf("alignments_loadSubjectsIntoMemory=%d/%d\n", totalCopied, totalLengths);
}

// Perform regular gapped scoring on an extension
void alignments_regularGappedAlignment(struct PSSMatrix PSSMatrix,
	struct ungappedExtension* ungappedExtension, struct alignment* alignment)
{
	int4 newScore;
	unsigned char* subject;

	subject = alignments_selectRegion(alignment, ungappedExtension);

//    printf("[%d] subject=%p", subject[ungappedExtension->seed.subjectOffset / 4], subject); fflush(stdout);

	blast_dloc = alignment->descriptionLocation;

    newScore
        = gappedScoring_score(ungappedExtension, PSSMatrix,
          alignment->subjectLength, subject,
          statistics_gappedNominalDropoff);

    #ifdef VERBOSE
    if (blast_dloc == parameters_verboseDloc)
    	printf("Was %d Now %d\n", ungappedExtension->nominalScore, newScore);
    #endif

    // If the extension's score has dropped considerably
    if (newScore * parameters_semiGappedR2 < ungappedExtension->nominalScore)
    {
        // Rescore with larger dropoff
        newScore
        = gappedScoring_score(ungappedExtension, PSSMatrix, alignment->subjectLength,
                              subject, statistics_gappedFinalNominalDropoff);

        #ifdef VERBOSE
        if (blast_dloc == parameters_verboseDloc)
            printf("Rescore now %d\n", newScore);
        #endif
    }

    ungappedExtension->nominalScore = newScore;

    if (ungappedExtension->nominalScore < blast_gappedNominalCutoff)
    {
        // If extension scores below cutoff, mark it as deleted
        ungappedExtension->status = ungappedExtension_DELETED;
    }
    else
    {
        // Else mark it as scored with regular gapped alignment
        ungappedExtension->status = ungappedExtension_GAPPED;
    }

    // Remove any ungapped extensions in this alignment that are in the area covered
    // by the gapped scoring just performed
    alignments_pruneRegion(alignment, ungappedExtension);
}


// Find the top N final alignments
void alignments_findTopFinalAlignments(struct PSSMatrix PSSMatrix)
{
    struct finalAlignment* goodAlignment, *finalAlignment;
	struct ungappedExtension* ungappedExtension;
    struct alignment* alignment;
	int4 position, bestScore;
    int4 maximumDynamicCutoff;

//    alignments_printGoodAlignments();

    // Get the Nth alignment
    goodAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_goodAlignments, parameters_numDisplayAlignments);

    // Calculate the maximum possibly dynamic cutoff
	maximumDynamicCutoff = goodAlignment->highestNominalScore / parameters_semiGappedR1;

    // Alignments above that maximum cutoff are definately final alignments
    position = 0;
    while (position < parameters_numDisplayAlignments)
    {
    	// Get each alignment in descending order of score
        goodAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_goodAlignments, position);
        alignment = goodAlignment->alignment;

		// Stop when score before maximum cutoff
		if (goodAlignment->highestNominalScore < maximumDynamicCutoff)
        	break;

		// Add to final alignments
		alignments_addFinalAlignment(goodAlignment->highestNominalScore, alignment);

        position++;
    }

//    alignments_printFinalAlignments();

    // For the rest of the top N alignments
    while (position < parameters_numDisplayAlignments)
    {
    	// Continue getting each alignment in descending order of score
        goodAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_goodAlignments, position);
        alignment = goodAlignment->alignment;

//        printf("Alignment score=%d\n", goodAlignment->highestNominalScore);

        // Perform regular gapped scoring
        bestScore = 0;
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL)
        {
        	// Find the highest scoring semigapped alignment
        	if (ungappedExtension->status != ungappedExtension_DELETED &&
                ungappedExtension->nominalScore == goodAlignment->highestNominalScore)
            {
                // Perform gapped alignment
                alignments_regularGappedAlignment(PSSMatrix, ungappedExtension, alignment);

                // Check for join to another extension
                alignments_checkForJoin(alignment, ungappedExtension, PSSMatrix);

                // Update the alignment's best score
                bestScore = ungappedExtension->nominalScore;
            }

            ungappedExtension = ungappedExtension->next;
        }

        // Update the alignment's score
        goodAlignment->highestNominalScore = bestScore;

		// Add to final alignments
		alignments_addFinalAlignment(goodAlignment->highestNominalScore, alignment);

        position++;
    }

//    alignments_printFinalAlignments();

    // Initialize dynamic cutoffs
    blast_dynamicGappedNominalCutoff = 0;
    blast_dynamicNominalR1cutoff = 0;

    // For the remaining good alignments
    while (position < alignments_goodAlignments->numEntries)
    {
    	// Continue getting each alignment in descending order of score
        goodAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_goodAlignments, position);
        alignment = goodAlignment->alignment;

        // Stop when below dynamic cutoff
        if (goodAlignment->highestNominalScore < blast_dynamicNominalR1cutoff)
        	break;

        // Perform regular gapped scoring
        bestScore = 0;
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL)
        {
        	// Find the highest scoring semigapped alignment
        	if (ungappedExtension->status != ungappedExtension_DELETED &&
                ungappedExtension->nominalScore == goodAlignment->highestNominalScore)
            {
                // Perform gapped alignment
                alignments_regularGappedAlignment(PSSMatrix, ungappedExtension, alignment);

                // Check for join to another extension
                alignments_checkForJoin(alignment, ungappedExtension, PSSMatrix);

                // Update the alignment's best score
                bestScore = ungappedExtension->nominalScore;
            }

            ungappedExtension = ungappedExtension->next;
        }

//        printf("Position=%d OldScore=%d NewScore=%d DynamicCutoff=%d FastCutoff=%d\n", position,
//               goodAlignment->highestNominalScore, bestScore, blast_dynamicGappedNominalCutoff,
//               blast_dynamicNominalR1cutoff);

        // Update the alignment's score
        goodAlignment->highestNominalScore = bestScore;

        // If it scores high enough to join the top N final alignments
		if (bestScore > blast_dynamicGappedNominalCutoff)

		// If it score above cutoff
//		if (bestScore > blast_gappedNominalCutoff)
        {
            // Add to final alignments
            alignments_addFinalAlignment(goodAlignment->highestNominalScore, alignment);

            // Re-sort final alignments and final minimum score of top N
            alignments_sortFinalAlignments();
            finalAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_finalAlignments,
                                                     parameters_numDisplayAlignments - 1);

            // Calculate dynamic cutoffs
            blast_dynamicGappedNominalCutoff = finalAlignment->highestNominalScore;
            blast_dynamicNominalR1cutoff = blast_dynamicGappedNominalCutoff * parameters_semiGappedR1;
		}

        position++;
    }
}

//debug==================================
void print_goodAlignmentInfo()
{
	struct finalAlignment *goodAlignment;
	struct alignment *alignment;
	struct ungappedExtension *ungappedExtension;
	int alignmentNo = 0;

	memSingleBlock_resetCurrent(alignments_goodAlignments);
	while ((goodAlignment = (struct finalAlignment *)memSingleBlock_getCurrent(alignments_goodAlignments)) != NULL)
	{
		alignment = goodAlignment->alignment;

		// For each ungapped extension that hasn't been deleted
		ungappedExtension = alignment->ungappedExtensions;
		while (ungappedExtension != NULL)
		{
			printf("%6d ", alignmentNo);
			ungappedExtension_print(ungappedExtension);

			ungappedExtension = ungappedExtension->next;
		}

		alignmentNo++;
	}
	return;
}
//================================

// Given a collection of good alignments, find the final top N alignments above cutoff
void alignments_findFinalAlignments(struct PSSMatrix PSSMatrix)
{
    struct finalAlignment* goodAlignment;
    struct alignment* alignment;
	struct ungappedExtension* ungappedExtension;
    int4 bestScore;

    // Sort good alignments by score
    // TODO: How long does this take?
    alignments_sortGoodAlignments();

	//debug======================
	//print_goodAlignmentInfo();
	//alignments_printGoodAlignments();
	//===========================
	
    // If we have more good alignments than we can display
    if (parameters_numDisplayAlignments > 0 &&
	    alignments_goodAlignments->numEntries > parameters_numDisplayAlignments)
	{
        goodAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_goodAlignments,
                                                parameters_numDisplayAlignments - 1);

		// If it is clear we will have more final alignments than we can display
		if (goodAlignment->highestNominalScore > blast_nominalR2cutoff)
        {
        	// Find the top N alignments
			alignments_findTopFinalAlignments(PSSMatrix);
            return;
        }
    }

    // Further score all good alignments to see if they score above cutoff
    memSingleBlock_resetCurrent(alignments_goodAlignments);
    while ((goodAlignment = (struct finalAlignment *)memSingleBlock_getCurrent(alignments_goodAlignments)) != NULL)
    {
		alignment = goodAlignment->alignment;

        bestScore = 0;
        // For each ungapped extension that hasn't been deleted
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL)
        {
            #ifdef VERBOSE
            if (parameters_verboseDloc == alignment->descriptionLocation)
				ungappedExtension_print(ungappedExtension);
            #endif

            if (ungappedExtension->status != ungappedExtension_DELETED)
            {
                // If it isn't clear if the extension is above cutoff
                if (ungappedExtension->nominalScore >= blast_nominalR1cutoff &&
                    ungappedExtension->nominalScore <= blast_nominalR2cutoff)
                {
                    // Perform regular gapped scoring
                    alignments_regularGappedAlignment(PSSMatrix, ungappedExtension, alignment);
                }

                // Update the alignment's best score
                if (ungappedExtension->nominalScore > bestScore)
                    bestScore = ungappedExtension->nominalScore;
            }
            ungappedExtension = ungappedExtension->next;
        }

        goodAlignment->highestNominalScore = bestScore;

        // If the alignment scores above cutoff, include in list of final alignments
        if (goodAlignment->highestNominalScore >= blast_gappedNominalCutoff)
        {
            alignments_addFinalAlignment(goodAlignment->highestNominalScore, alignment);
        }
    }
}

// Expand a cluster and perform stages 1 and 2 of search on the children
int alignments_expandCluster(struct alignment* alignment, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP)
{
	struct child* children, *child;
    uint4 numChildren, childNum, numNewAlignments = 0;
	struct sequenceData sequenceData;

    // Get children
    children = readdb_getChildren(alignment->subject, alignment->subjectLength,
               alignment->encodedLength, alignment->descriptionLocation, &numChildren);

//	printf("Alignment dloc=%d\n", alignment->descriptionLocation);

    if (numChildren == 0)
		return 0;

    // For each child
    childNum = 0;
    while (childNum < numChildren)
    {
        child = children + childNum;

        sequenceData.descriptionStart = child->descriptionLocation;
        sequenceData.descriptionLength = child->descriptionLength;
		sequenceData.sequenceLength = child->length;
		sequenceData.encodedLength = child->length;
		sequenceData.sequence = child->sequence;

        // Re-initialize the hitMatrix
        hitMatrix_reinitialize(PSSMatrix.length, readdb_longestSequenceLength, child->sequence);

        // Perform stage 1 and 2 of search
        if (parameters_oneHitTrigger)
        {
            search_protein1hit(PSSMatrix, PSSMatrixFP, &sequenceData, 1, constants_maxSignedInt);
		}
        else
        {
            search_protein2hit(PSSMatrix, PSSMatrixFP, &sequenceData, 1, constants_maxSignedInt);
        }

        if (alignments_currentAlignment != NULL)
        {
        	if (numNewAlignments == 0)
            	alignments_numClusters++;

        	numNewAlignments++;
        	alignments_currentAlignment->inMemorySubject = 1;
            alignments_currentAlignment->cluster = alignments_numClusters;
//            printf("Child %d start=%d length=%d dloc=%d\n", childNum, child->regionStart,
//                   child->length, child->descriptionLocation);
		}

        blast_numExpandedSequences++;
        childNum++;
    }

/*    if (numNewAlignments == 0)
    {
    	print_singleSequence(alignment->subject, alignment->subjectLength); printf("\n");
	}*/

    return 1;
}

// Get the tracebacks for all of the final alignments
void alignments_getTracebacks(struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP)
{
    struct finalAlignment* finalAlignment;
	struct gappedExtension* gappedExtension;
	struct ungappedExtension* ungappedExtension, *highestScoringExtension, oldUngappedExtension;
    struct alignment* alignment;
    struct unpackRegion* unpackRegion;
    int4 numProcessed = 0, numAboveCutoff = 0, repeatComputeTracebacks = 1;
	int4 alignmentCount;
	unsigned char* subject;

    // Sort final alignments by score
    alignments_sortFinalAlignments();

//    alignments_printFinalAlignments();

    // Only keep alignments above cutoff
    while (alignments_finalAlignments->numEntries > 0)
    {
    	// Get last alignment in list
		finalAlignment = (struct finalAlignment *)memSingleBlock_getLastEntry(alignments_finalAlignments);

        // Stop if above cutoff
        if (finalAlignment->highestNominalScore != 0)
        	break;

        // Otherwise remove it from the list
    	alignments_finalAlignments->numEntries--;
    }

    // For each alignment that is in the top numDisplayAlignments but not the top numDisplayTracebacks
    alignmentCount = parameters_numDisplayTracebacks;
    while (alignmentCount < parameters_numDisplayAlignments &&
           alignmentCount < alignments_finalAlignments->numEntries)
	{
    	finalAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_finalAlignments, alignmentCount);
        alignment = finalAlignment->alignment;
        blast_dloc = alignment->descriptionLocation;

        // Get the highest scoring ungapped extension
        ungappedExtension = alignment->ungappedExtensions;
        highestScoringExtension = ungappedExtension;
        while (ungappedExtension != NULL)
        {
        	if (ungappedExtension->nominalScore > highestScoringExtension->nominalScore &&
                ungappedExtension->status != ungappedExtension_DELETED)
            {
				highestScoringExtension = ungappedExtension;
            }
        	ungappedExtension = ungappedExtension->next;
		}

        if (highestScoringExtension != NULL)
        {
    		subject = alignments_selectRegion(alignment, highestScoringExtension);

            // Perform gapped scoring with higher dropoff
            highestScoringExtension->nominalScore
                = gappedScoring_score(highestScoringExtension, PSSMatrix, alignment->subjectLength,
                                      subject, statistics_gappedFinalNominalDropoff);

			finalAlignment->highestNominalScore = highestScoringExtension->nominalScore;
        }

//        printf("Rescore with larger dropoff num %d: Score=%d\n", alignmentCount, highestScoringExtension->nominalScore);
        alignmentCount++;
    }

	while (repeatComputeTracebacks)
    {
    	numAboveCutoff = 0; numProcessed = 0;
        memSingleBlock_resetCurrent(alignments_finalAlignments);
        while (((finalAlignment = (struct finalAlignment *)memSingleBlock_getCurrent(alignments_finalAlignments)) != NULL)
        && (numAboveCutoff < parameters_numDisplayTracebacks || parameters_numDisplayTracebacks == 0))
        {
            alignment = finalAlignment->alignment;
            blast_dloc = alignment->descriptionLocation;

            // If traceback haven't been computed for this alignment
            if (alignment->gappedExtensions == NULL)
            {
                // Unpack part or all of the subject (for nucleotide)
                blast_gappedExtendTime += clock();
                blast_unpackTime -= clock();
                unpack_unpackSubject(PSSMatrix, alignment);
                blast_unpackTime += clock();
                blast_gappedExtendTime -= clock();

                // For each ungapped extension that hasn't been deleted
                ungappedExtension = alignment->ungappedExtensions;
                while (ungappedExtension != NULL)
                {
                    // If extension scores above cutoff
                    if (ungappedExtension->status != ungappedExtension_DELETED)
                    {
                        // Make copy of ungapped extension
                    	oldUngappedExtension = *ungappedExtension;

                        // If subject and query are short enough and sequence does not have multiple
                        // unpack regions, use faster but less memory efficient gapped alignment with traceback
                        if (((uint8)PSSMatrix.length * (uint8)alignment->subjectLength <
                            (uint8)constants_maximumTracebackSize) && unpack_entireSubjectUnpacked(alignment))
                        {
                            gappedExtension
                                = fasterGappedExtension_build(ungappedExtension, PSSMatrix,
                                  alignment->subjectLength, alignment->unpackRegions[0].unpackedSubject,
                                  statistics_gappedFinalNominalDropoff);
                        }
                        // Otherwise use slower but more memory-efficient gapped alignment
                        else
                        {
                            unpackRegion = unpack_selectRegion(alignment->unpackRegions,
                                           alignment->numUnpackRegions, ungappedExtension->seed.subjectOffset);

                            gappedExtension
                                = gappedExtension_build(ungappedExtension, PSSMatrix, alignment->subjectLength,
                                      alignment->subject, unpackRegion, statistics_gappedFinalNominalDropoff);
                        }

                        // Calculate normalized score and e-value
                        gappedExtension_score(gappedExtension);

                        // Add it to the current alignment
                        if (gappedExtension->nominalScore >= blast_gappedNominalCutoff)
                            alignments_addGappedExtension(alignment, gappedExtension);
                        else
                        {
                            free(gappedExtension->trace.traceCodes);
                            free(gappedExtension);
                        }

                        // Check for ungapped extensions that were mistakenly pruned
                        if (parameters_tableScoring)
                        	alignments_unpruneRegion(alignment, &oldUngappedExtension,
                                                     ungappedExtension, PSSMatrix, PSSMatrixFP);

                        // Remove any ungapped extensions in this alignment that are in the area covered
                        // by the gapped scoring just performed
                        alignments_pruneRegion(alignment, ungappedExtension);
                   }
                   ungappedExtension = ungappedExtension->next;
                }

                // Finally prune extensions that share a start or end point
                alignments_pruneOverlappingExtensions(alignment);

//                printf("Was %d Now %d\n", finalAlignment->highestNominalScore, alignment->gappedExtensions->nominalScore);

                // Update final alignment's high-score to that of the first (and highest-scoring)
                // gapped extension in the list
                if (alignment->gappedExtensions != NULL)
                    finalAlignment->highestNominalScore = alignment->gappedExtensions->nominalScore;
                else
                    finalAlignment->highestNominalScore = 0;

            	numProcessed++;
//	            printf("Computed alignment score=%d\n", finalAlignment->highestNominalScore);
			}

            // Tally number of final alignments above cutoff
            if (finalAlignment->highestNominalScore > 0)
                numAboveCutoff++;
        }

//    	printf("Traceback alignments performed=%d\n", numProcessed);

        // Sort final alignments by score
        alignments_sortFinalAlignments();

//        printf("repeatComputeTracebacks:");

        // If the first numDisplayTracebacks alignments have traceback computed, stop
        numProcessed = 0; repeatComputeTracebacks = 0;
        memSingleBlock_resetCurrent(alignments_finalAlignments);
        while ((finalAlignment = (struct finalAlignment *)memSingleBlock_getCurrent(alignments_finalAlignments)) != NULL
               && numProcessed < parameters_numDisplayTracebacks)
        {
            alignment = finalAlignment->alignment;
            if (alignment->gappedExtensions == NULL && finalAlignment->highestNominalScore != 0)
            {
//            	printf("1");
            	repeatComputeTracebacks = 1;
				break;
            }
//            printf("0");

            numProcessed++;
		}

	}

    // Only keep top N alignments
    if (parameters_numDisplayAlignments != 0 &&
        alignments_finalAlignments->numEntries > parameters_numDisplayAlignments)
	{
		alignments_finalAlignments->numEntries = parameters_numDisplayAlignments;
    }

    // Only keep alignments above cutoff
    while (alignments_finalAlignments->numEntries > 0)
    {
    	// Get last alignment in list
		finalAlignment = (struct finalAlignment *)memSingleBlock_getLastEntry(alignments_finalAlignments);

        // Stop if above cutoff
        if (finalAlignment->highestNominalScore != 0)
        	break;

        // Otherwise remove it from the list
    	alignments_finalAlignments->numEntries--;
    }
}

// Returns true if the given alignment score is high enough for it to be
// one of the N displayed final alignments
int alignments_isFinalAlignment(uint4 score)
{
	struct finalAlignment* finalAlignment;
    struct alignment* alignment;

//    printf("%d,%d\n", parameters_numDisplayAlignments, alignments_finalAlignments->numEntries);

	if (parameters_numDisplayAlignments == 0 ||
        alignments_finalAlignments->numEntries < parameters_numDisplayAlignments)
    	return 1;

	if (!alignments_finalAlignmentsSorted)
    {
        // Sort final alignments by score
        alignments_sortFinalAlignments();
        alignments_finalAlignmentsSorted = 1;
    }

    // Get the N-1th final alignment
    finalAlignment = (struct finalAlignment *)memSingleBlock_getEntry(alignments_finalAlignments,
                                             parameters_numDisplayAlignments - 1);

	// If new alignment scores high enough to be in the top N
	if (score > finalAlignment->highestNominalScore)
    {
		// Final alignments will need to be resorted
		alignments_finalAlignmentsSorted = 0;

        // Free the lowest scoring final alignment
		alignment = finalAlignment->alignment;
		alignment_freeAlignment(alignment);
        alignments_finalAlignments->numEntries--;

        return 1;
    }
    else
    {
    	return 0;
    }
}

