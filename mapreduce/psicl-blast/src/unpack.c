// unpack.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for unpacking nucleotide sequences

#include "blast.h"

struct memBlocks* unpack_unpackRegions;
struct memBlocks* unpack_subjectRegions;

// Initialize region copying/unpacking
void unpack_initialize()
{
	unpack_unpackRegions = memBlocks_initialize(sizeof(struct unpackRegion),
                           constants_initialAllocUnpackRegions);

	unpack_subjectRegions = memBlocks_initialize(sizeof(struct unpackRegion),
                            constants_initialAllocUnpackRegions);
}

// Compare the start offset of two unpack regions
int4 unpack_compareUnpackRegions(const void* unpackRegion1, const void* unpackRegion2)
{
	const struct unpackRegion *u1, *u2;

	u1 = (struct unpackRegion*)unpackRegion1;
	u2 = (struct unpackRegion*)unpackRegion2;

	if (u1->startOffset > u2->startOffset)
	{
		return 1;
	}
	else if (u1->startOffset < u2->startOffset)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

// Return the unpack region in the given list that contains the given subject offset
struct unpackRegion* unpack_selectRegion(struct unpackRegion* unpackRegions, uint4 numUnpackRegions,
                                         uint4 subjectOffset)
{
    // For each region
	while (numUnpackRegions > 0)
    {
    	if (unpackRegions->startOffset <= subjectOffset && unpackRegions->endOffset > subjectOffset)
        {
			// Return it if seed falls within its range
            return unpackRegions;
        }

        unpackRegions++;
		numUnpackRegions--;
    }

    fprintf(stderr, "Error unpacking subject and performing gapped alignment\n");
    exit(-1);
}

// Extend the start of a region if necessary
void unpack_extendRegionStart(int4 position, struct unpackRegion* unpackRegion)
{
	unsigned char* newUnpackedSubject;
    int4 newRegionStart, newRegionEnd;

	if (position < unpackRegion->startOffset)
    {
    	// Extend the region start
    	newRegionStart = unpackRegion->startOffset - constants_unpackRegionExtend;
        if (newRegionStart < 0) newRegionStart = 0;
        newRegionEnd = unpackRegion->endOffset;

        // Make start of region a multiple of 4
        newRegionStart = (newRegionStart / 4) * 4;

        // Declare memory for the new region
    	newUnpackedSubject = (unsigned char*)global_malloc(sizeof(char) * (newRegionEnd - newRegionStart));
		newUnpackedSubject -= newRegionStart;

        // Copy unpacked subject from old region to new
        memcpy(newUnpackedSubject + unpackRegion->startOffset,
               unpackRegion->unpackedSubject + unpackRegion->startOffset,
               sizeof(char) * (unpackRegion->endOffset - unpackRegion->startOffset));

		// Free old subject
		unpackRegion->unpackedSubject += unpackRegion->startOffset;
        free(unpackRegion->unpackedSubject);

		// Unpack the new part of the region
		encoding_byteUnpackRegion(newUnpackedSubject + newRegionStart, unpackRegion->subject + (newRegionStart / 4),
                                  unpackRegion->startOffset - newRegionStart);

		unpackRegion->unpackedSubject = newUnpackedSubject;

        unpackRegion->startOffset = newRegionStart;
	}
}

// Extend the end of a region if necessary
void unpack_extendRegionEnd(int4 position, struct unpackRegion* unpackRegion)
{
	unsigned char* newUnpackedSubject;
    int4 newRegionStart, newRegionEnd;

//    printf("pos=%d region=%d,%d subjectLength=%d\n", position, unpackRegion->startOffset,
//           unpackRegion->endOffset, unpackRegion->subjectLength); fflush(stdout);

	if (position > unpackRegion->endOffset)
    {
    	// Extend the region end
    	newRegionStart = unpackRegion->startOffset;
        newRegionEnd = unpackRegion->endOffset + constants_unpackRegionExtend;
        if (newRegionEnd > unpackRegion->subjectLength) newRegionEnd = unpackRegion->subjectLength;

        // Realloc memory for the new region
		unpackRegion->unpackedSubject += unpackRegion->startOffset;
		newUnpackedSubject = (unsigned char*)global_realloc(unpackRegion->unpackedSubject,
                              sizeof(char) * (newRegionEnd - newRegionStart));
		newUnpackedSubject -= newRegionStart;

        // Round old end
        unpackRegion->endOffset = (unpackRegion->endOffset / 4) * 4;

		// Unpack the new part of the region
		encoding_byteUnpackRegion(newUnpackedSubject + unpackRegion->endOffset,
        	unpackRegion->subject + (unpackRegion->endOffset / 4),
            newRegionEnd - unpackRegion->endOffset);

		unpackRegion->unpackedSubject = newUnpackedSubject;

        unpackRegion->endOffset = newRegionEnd;
	}
}

// Returns a list of subject region start/end values that cover the area likely to be
// processed during gapped extension, or the maximum area to be processed if maximum=1
uint4 unpack_getRegions(struct PSSMatrix PSSMatrix, struct alignment* alignment, char maximum,
                        struct memBlocks* regions)
{
	struct ungappedExtension* ungappedExtension;
	int4 regionStart, regionEnd, subjectOffset, queryOffset, queryBefore, queryAfter;
    int4 maxRegionStart, maxRegionEnd;
	uint4 maxNumRegions = 0;
	struct unpackRegion *firstRegion = NULL, *currentRegion, *lastRegion, *nextRegion;

    // Count the number of ungapped extensions
    ungappedExtension = alignment->ungappedExtensions;
    while (ungappedExtension != NULL)
    {
        if (ungappedExtension->status != ungappedExtension_DELETED)
        {
        	maxNumRegions++;
        }
        ungappedExtension = ungappedExtension->next;
    }

    // Allocate memory for some regions
    lastRegion = firstRegion = memBlocks_newEntries(regions, maxNumRegions);

    // For each ungapped extension that hasn't been deleted
    ungappedExtension = alignment->ungappedExtensions;
    while (ungappedExtension != NULL)
    {
        if (ungappedExtension->status != ungappedExtension_DELETED)
        {
            // Calculate the maximum possible region to be covered by this gapped alignment
            queryOffset = ungappedExtension->seed.queryOffset;
            subjectOffset = ungappedExtension->seed.subjectOffset;
            queryBefore = queryOffset - PSSMatrix_strandStart(PSSMatrix, queryOffset);
            queryAfter = PSSMatrix_strandEnd(PSSMatrix, queryOffset) - queryOffset;

            maxRegionStart = subjectOffset - queryBefore - statistics_gappedFinalNominalDropoff
                        - (queryBefore * PSSMatrix.highestValue / parameters_extendGap) - 2;
            if (maxRegionStart < 0) maxRegionStart = 0;

            maxRegionEnd = subjectOffset + queryAfter + statistics_gappedFinalNominalDropoff +
                        (queryAfter * PSSMatrix.highestValue / parameters_extendGap) + 2;
            if (maxRegionEnd > alignment->subjectLength) maxRegionEnd = alignment->subjectLength;

            if (maximum)
            {
            	// Just use the maximum
                regionStart = maxRegionStart;
                regionEnd = maxRegionEnd;
			}
            else
            {
                // Calculate start and end of region likely to be covered by this gapped alignment
                regionStart = ungappedExtension->start.subjectOffset - constants_unpackRegionExtend;
                if (regionStart < 0) regionStart = 0;

                regionEnd = ungappedExtension->end.subjectOffset + constants_unpackRegionExtend + 1;
                if (regionEnd > alignment->subjectLength) regionEnd = alignment->subjectLength;

                // The region should be no greater than the maximum
                if (regionStart < maxRegionStart)
                	regionStart = maxRegionStart;

                if (regionEnd > maxRegionEnd)
                	regionEnd = maxRegionEnd;
            }

            // Make start of region a multiple of 4
            regionStart = (regionStart / 4) * 4;

            // Create new region
            lastRegion->startOffset = regionStart;
            lastRegion->endOffset = regionEnd;
            lastRegion++;
       }
       ungappedExtension = ungappedExtension->next;
    }

    memBlocks_returnUnused(regions, (maxNumRegions - (lastRegion - firstRegion)));

    // If there are multiple regions for this subject
    if (lastRegion - firstRegion > 1)
    {
        // Sort the regions in order of start position
        qsort(firstRegion, lastRegion - firstRegion,
              sizeof(struct unpackRegion), unpack_compareUnpackRegions);

        // For each region in order of starting position
        currentRegion = firstRegion;
        while (currentRegion < lastRegion)
        {
        	// If it has been marked as deleted
			if (currentRegion->startOffset < alignment->subjectLength)
            {
                // For each region after current region that intersects with it
                nextRegion = currentRegion + 1;
                while ((nextRegion < lastRegion) &&
                       (nextRegion->startOffset <= currentRegion->endOffset))
                {
                    if (nextRegion->endOffset <= currentRegion->endOffset)
                    {
                        // This region is contained by current region, mark as deleted
                        nextRegion->startOffset = alignment->subjectLength + 1;
                        nextRegion->endOffset = alignment->subjectLength + 1;
                    }
                    else
                    {
                        // This region partially overlaps with current region; extend currentRegion
                        currentRegion->endOffset = nextRegion->endOffset;
                        // Remove this region
                        nextRegion->startOffset = alignment->subjectLength + 1;
                        nextRegion->endOffset = alignment->subjectLength + 1;
                    }

                    nextRegion++;
                }
            }

        	currentRegion++;
        }

        // Sort the regions in order of start position
        qsort(firstRegion, lastRegion - firstRegion,
              sizeof(struct unpackRegion), unpack_compareUnpackRegions);

		// Remove deleted regions
        currentRegion = lastRegion - 1;
        while (currentRegion->startOffset > alignment->subjectLength)
        {
        	lastRegion = currentRegion;
        	currentRegion--;
            memBlocks_returnUnused(regions, 1);
        }
    }

/*    printf("Regions:\n");

    currentRegion = firstRegion;
    while (currentRegion < lastRegion)
    {
    	printf("[%p]", currentRegion);
        printf("Region from %d to %d\n", currentRegion->startOffset, currentRegion->endOffset); fflush(stdout);
        currentRegion++;
    }*/

    return lastRegion - firstRegion;
}

// Load a single subject into memory
int4 unpack_loadSubject(struct PSSMatrix PSSMatrix, struct alignment* alignment)
{
	uint4 totalCopied = 0;
	unsigned char *subject, *edits, *endEdits;
    struct unpackRegion *firstRegion = NULL, *lastRegion, *currentRegion;
    int4 numRegions, regionStart, regionEnd;

    // If protein search
    if (encoding_alphabetType == encoding_protein)
    {
        // Make copy of sequence
        subject = (unsigned char*)global_malloc(sizeof(unsigned char) * alignment->encodedLength);
        subject++;
        memcpy(subject - 1, alignment->subject - 1, alignment->encodedLength);
        alignment->subject = subject;

        blast_totalCopied += alignment->encodedLength;
    }
    // If a nucleotide search
    else
    {
		// Get a list of regions to copy
        numRegions = unpack_getRegions(PSSMatrix, alignment, 1, unpack_subjectRegions);
        lastRegion = memBlocks_getLastEntry(unpack_subjectRegions);
        lastRegion++;
        firstRegion = lastRegion - numRegions;

        #ifdef VERBOSE
        if (parameters_verboseDloc == alignment->descriptionLocation)
        {
        	printf("%d regions for subject\n", lastRegion - firstRegion); fflush(stdout);
		}
        #endif

        // Copy each region into memory
        currentRegion = firstRegion;
        while (currentRegion < lastRegion)
        {
            #ifdef VERBOSE
            if (parameters_verboseDloc == alignment->descriptionLocation)
            {
                printf("Load region %d to %d into memory\n", currentRegion->startOffset,
                       currentRegion->endOffset); fflush(stdout);
                fflush(stdout);
            }
            #endif

            regionStart = currentRegion->startOffset / 4;
            regionEnd = (currentRegion->endOffset + 3) / 4;

            currentRegion->unpackedSubject = NULL;
			currentRegion->subject = (unsigned char*)global_malloc(sizeof(unsigned char)
                                                   * (regionEnd - regionStart));

            totalCopied += regionEnd - regionStart;
			memcpy(currentRegion->subject, alignment->subject + regionStart, regionEnd - regionStart);
			currentRegion->subject -= regionStart;
			currentRegion->subjectLength = alignment->subjectLength;

        	blast_totalCopied += (regionEnd - regionStart);

            currentRegion++;
		}

        // Store new alignment regions
        alignment->unpackRegions = firstRegion;
        alignment->numUnpackRegions = lastRegion - firstRegion;

        // If there are edits for this subject
        if (alignment->edits != NULL)
        {
            edits = alignment->edits;
            endEdits = alignment->subject + alignment->encodedLength;

            // Make an in-memory copy of them
            alignment->edits = (unsigned char*)malloc(sizeof(char) * (endEdits - edits));
            memcpy(alignment->edits, edits, endEdits - edits);
        }

        alignment->subject = NULL;
    }

    alignment->inMemorySubject = 1;

    return totalCopied;

}

// Unpack entire or sections of a subject sequence before gapped alignment
void unpack_unpackSubject(struct PSSMatrix PSSMatrix, struct alignment* alignment)
{
    unsigned char *subject, *unpackedSubject, wildcard, *edits, *endEdits;
    uint4 wildcardPosition;
    struct unpackRegion *firstRegion = NULL, *lastRegion, *currentRegion, *unpackRegion;
    int4 regionStart, regionEnd, numRegions;

    // No need to unpack a protein subject, or already unpacked nucleotide subject
    if (parameters_ssearch || encoding_alphabetType == encoding_protein)
    {
    	// Just create a single region covering the entire sequence
        firstRegion = memBlocks_newEntry(unpack_unpackRegions);
        firstRegion->startOffset = 0;
        firstRegion->endOffset = alignment->subjectLength;
        firstRegion->subject = alignment->subject;
        firstRegion->unpackedSubject = alignment->subject;
        firstRegion->subjectLength = alignment->subjectLength;
        alignment->unpackRegions = firstRegion;
        alignment->numUnpackRegions = 1;
        return;
    }

    // Get the subject regions for this alignment
    numRegions = unpack_getRegions(PSSMatrix, alignment, 0, unpack_unpackRegions);
    lastRegion = memBlocks_getLastEntry(unpack_unpackRegions);
    lastRegion++;
    firstRegion = lastRegion - numRegions;

    // Sort the regions in order of start position
	qsort(firstRegion, lastRegion - firstRegion,
          sizeof(struct unpackRegion), unpack_compareUnpackRegions);

    // Unpack each region
    currentRegion = firstRegion;
    while (currentRegion < lastRegion)
    {
    	regionEnd = currentRegion->endOffset;
        regionStart = currentRegion->startOffset;

        #ifdef VERBOSE
        if (parameters_verboseDloc == alignment->descriptionLocation)
		{
	        printf("Unpack subject region %d to %d (length=%d)\n", regionStart, regionEnd,
                   alignment->subjectLength);
            fflush(stdout);
        }
		#endif

		// Get the subject region to be unpacked
        if (alignment->unpackRegions == NULL)
        {
        	subject = alignment->subject;
        }
        else
        {
            unpackRegion = unpack_selectRegion(alignment->unpackRegions, alignment->numUnpackRegions,
                                               regionStart);
			subject = unpackRegion->subject;
        }

		// Declare memory for the region
        unpackedSubject = (unsigned char*)global_malloc(sizeof(char) * (regionEnd - regionStart));

		// Unpack the region of interest
        encoding_byteUnpackRegion(unpackedSubject, subject + (regionStart / 4),
                                  regionEnd - regionStart);
        unpackedSubject -= regionStart;
		currentRegion->unpackedSubject = unpackedSubject;

        currentRegion->subject = subject;
        currentRegion->subjectLength = alignment->subjectLength;

        blast_totalUnpacked += (regionEnd - regionStart);

        currentRegion++;
    }

	currentRegion = firstRegion;

	// Get wildcard edits for the sequence
	edits = alignment->edits;
	endEdits = alignment->edits + alignment->encodedLength - ((alignment->subjectLength + 3) / 4);

    // If there are edits
    if (edits < endEdits)
    {
        // Read first wildcard
        wildcard = *edits;
        edits++;

        // Read its position
        vbyte_getVbyte(edits, &wildcardPosition);

        // For each region in order of position in the subject
        while (currentRegion < lastRegion)
        {
            // Skip past edits that are before current region
            while (edits < endEdits && wildcardPosition < currentRegion->startOffset)
            {
                // Read wildcard
                wildcard = *edits;
                edits++;

                // Read its position
                vbyte_getVbyte(edits, &wildcardPosition);
            }

            // Process edits that are in the current region
            while (edits < endEdits && wildcardPosition < currentRegion->endOffset)
            {
                // Insert wildcard into sequence
                currentRegion->unpackedSubject[wildcardPosition] = wildcard;

                // Read next wildcard
                wildcard = *edits;
                edits++;

                // Read its position
                vbyte_getVbyte(edits, &wildcardPosition);
            }

            // Advance to the next region
            currentRegion++;
        }
	}

    alignment->unpackRegions = firstRegion;
    alignment->numUnpackRegions = lastRegion - firstRegion;
}

// Return 1 if the entire subject has been unpacked, otherwise 0
int unpack_entireSubjectUnpacked(struct alignment* alignment)
{
	if (alignment->unpackRegions[0].startOffset == 0 &&
        alignment->unpackRegions[0].endOffset == alignment->subjectLength)
		return 1;
    else
    	return 0;
}

// Free memory used to store unpacked regions
void unpack_free()
{
	struct unpackRegion* region;

    // For each unpack region
	if (!parameters_ssearch && encoding_alphabetType == encoding_nucleotide)
	{
        memBlocks_resetCurrent(unpack_unpackRegions);
        while ((region = memBlocks_getCurrent(unpack_unpackRegions)) != NULL)
        {
            // Free the unpacked sequence
            free(region->unpackedSubject + region->startOffset);
        }
	}

    // For each copied subject region
	memBlocks_resetCurrent(unpack_subjectRegions);
	while ((region = memBlocks_getCurrent(unpack_subjectRegions)) != NULL)
    {
		// Free the subject
        free(region->subject + region->startOffset / 4);
    }

    memBlocks_free(unpack_unpackRegions);
    memBlocks_free(unpack_subjectRegions);
}
