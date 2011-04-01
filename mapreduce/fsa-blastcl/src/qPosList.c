// qPosList.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Given a collection of lists of query positions, create a compact array that
// contains the lists. This allows sharing of query positions between words
// with identical or similar lists.

#include "blast.h"

struct initialQPosList
{
	int2* queryPositions;
    int2 numQueryPositions;
    int4 codeword;
};

struct memSingleBlock* qPosList_qPosLists;
int4 qPosList_numQPosLists;
int4 qPosList_maxQPosLists;

struct initialQPosList* qPosList_initialQPosLists;
int4 qPosList_numInitialQPosLists;

int4 qPosList_compareList(const void* list1, const void* list2);
void qPosList_processList(int2* queryPositions, int2 numQueryPositions, int4 codeword);

// Initialize for the construction of a new query position list
void qPosList_initialize(int4 maxNumLists)
{
    struct memSingleBlock* list;

    qPosList_qPosLists = (struct memSingleBlock*)global_malloc(sizeof(struct memSingleBlock) * maxNumLists);
    qPosList_numQPosLists = 0;
    qPosList_maxQPosLists = 0;

    // Initialize lists of query positions
    while (qPosList_maxQPosLists < maxNumLists)
    {
        list = qPosList_qPosLists + qPosList_maxQPosLists;
        memSingleBlock_initializeExisting(list, sizeof(struct queryPosition), 10);
        qPosList_maxQPosLists++;
    }

    qPosList_initialQPosLists = (struct initialQPosList*)global_malloc(sizeof(struct initialQPosList) * maxNumLists);
    qPosList_numInitialQPosLists = 0;
}

// Reset between query position list constructions
void qPosList_reset()
{
    qPosList_numQPosLists = 0;
    qPosList_numInitialQPosLists = 0;
}

// Free the structure
void qPosList_free()
{
    struct memSingleBlock* list;

    // Free initial lists
	free(qPosList_initialQPosLists);

    // Free each final qPos list
    while (qPosList_maxQPosLists > 0)
    {
    	qPosList_maxQPosLists--;
        list = qPosList_qPosLists + qPosList_maxQPosLists;
        free(list->block);
    }

    free(qPosList_qPosLists);
}

// Print the list
void qPosList_print()
{
	int4 listCount = 0, queryPositionCount;
    struct memSingleBlock* list;
    struct queryPosition* queryPosition;

    // For each list
    while (listCount < qPosList_numQPosLists)
    {
    	list = qPosList_qPosLists + listCount;

        printf("%d) ", listCount); fflush(stdout);

        // Traverse in reverse order and print elements
    	queryPositionCount = list->numEntries;
        while (queryPositionCount > 0)
        {
        	queryPositionCount--;
            queryPosition = (struct queryPosition*)list->block + queryPositionCount;

            if (queryPosition->codewords != NULL)
            {
            	printf("%d* ", queryPosition->queryPosition); fflush(stdout);
            }
            else
            {
            	printf("%d ", queryPosition->queryPosition); fflush(stdout);
            }
        }
        printf("\n");

		listCount++;
    }
}

// Add a list to the initial collection of query positions lists
void qPosList_addList(int2* queryPositions, int2 numQueryPositions, int4 codeword)
{
	qPosList_initialQPosLists[qPosList_numInitialQPosLists].queryPositions = queryPositions;
	qPosList_initialQPosLists[qPosList_numInitialQPosLists].numQueryPositions = numQueryPositions;
	qPosList_initialQPosLists[qPosList_numInitialQPosLists].codeword = codeword;

    qPosList_numInitialQPosLists++;
}

// Compare the two initial query position list lengths.
// Return 1 if list1 length > list2 length
// -1 if list1 length < list2 length and 0 if they are equal
int4 qPosList_compareInitialList(const void* list1, const void* list2)
{
	struct initialQPosList *l1, *l2;

	l1 = (struct initialQPosList*)list1;
	l2 = (struct initialQPosList*)list2;

	if (l1->numQueryPositions < l2->numQueryPositions)
	{
		return 1;
	}
	else if (l1->numQueryPositions > l2->numQueryPositions)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

// Process the initial collection of query position lists in order from
// shortest to longest
void qPosList_processLists()
{
	int2* queryPositions;
    int2 numQueryPositions;
    int4 codeword;

    // Sort initial query position lists from shortest to longest
    qsort(qPosList_initialQPosLists, qPosList_numInitialQPosLists,
          sizeof(struct initialQPosList), qPosList_compareInitialList);

    // Add in sorted order
	while (qPosList_numInitialQPosLists > 0)
    {
    	qPosList_numInitialQPosLists--;

		queryPositions = qPosList_initialQPosLists[qPosList_numInitialQPosLists].queryPositions;
		numQueryPositions = qPosList_initialQPosLists[qPosList_numInitialQPosLists].numQueryPositions;
		codeword = qPosList_initialQPosLists[qPosList_numInitialQPosLists].codeword;

        if (numQueryPositions != 0)
        {
			qPosList_processList(queryPositions, numQueryPositions, codeword);
        }
    }
}

// Process a given query position list
void qPosList_processList(int2* queryPositions, int2 numQueryPositions, int4 codeword)
{
	int4 listCount = 0, queryPositionCount, subset, present;
    struct memSingleBlock* list;
    struct queryPosition* queryPosition = NULL;
	struct codeword* newCodeword;

    // Iterative through existing query positions lists (ordered from longest to shortest)
    while (listCount < qPosList_numQPosLists)
    {
        // Check for one that contains a subset of to-be-added query positions
        list = qPosList_qPosLists + listCount;

        // Start by assuming it is
        subset = 1;

        // Iterate through each query position in the current existing list
        memSingleBlock_resetCurrent(list);
        while ((queryPosition = memSingleBlock_getCurrent(list)) != NULL && subset)
        {
            // Iterate through each query position in the new list (which is sorted)
            queryPositionCount = 0;
            while (queryPositionCount < numQueryPositions)
            {
                // Found a match, break out and proceed to next position in current list
                if (queryPosition->queryPosition == queryPositions[queryPositionCount])
                {
                    break;
                }
                // The query position is not present in the new list, then existing list
                // is not a subset of the new one
                else if (queryPosition->queryPosition < queryPositions[queryPositionCount])
                {
                    subset = 0;
                    break;
                }
                // Otherwise keep going
                queryPositionCount++;
            }

            // If we got to the end of the list, and didn't find a match, not a subset
            if (queryPositionCount == numQueryPositions)
                subset = 0;

            // If the query positions in the existing list processed so far match all of
            // the positions in the new list
            if (list->currentEntry == numQueryPositions && subset)
            {
                // We have a match, starting here
                newCodeword = global_malloc(sizeof(struct codeword));
				newCodeword->codeword = codeword;
                newCodeword->next = queryPosition->codewords;
                queryPosition->codewords = newCodeword;

                return;
            }
        }

        if (subset)
        {
            // If this existing list is a subset of the new list then add the new/additional
            // query positions to the end of it
            queryPosition = memSingleBlock_getLastEntry(list);

            // Iterate through each query position in the new list
            while (numQueryPositions > 0)
            {
                numQueryPositions--;
                present = 0;

                // Check if present in the existing list
                memSingleBlock_resetCurrent(list);
                while ((queryPosition = memSingleBlock_getCurrent(list)) != NULL && subset)
                {
                    // Found it
                    if (queryPosition->queryPosition == queryPositions[numQueryPositions])
                    {
                        present = 1;
                        break;
                    }
                }

                // Not present - add to the existing list with a null reference codeword
                if (!present)
                {
                    queryPosition = memSingleBlock_newEntry(list);
                    queryPosition->queryPosition = queryPositions[numQueryPositions];
                    // No refering codeword for any of the positions except the last
                    queryPosition->codewords = NULL;
                }

                queryPositionCount++;
            }

            // Get the last, new query position
            queryPosition = memSingleBlock_getLastEntry(list);

            // Add reference codeword to the last query position (will become first)
            newCodeword = global_malloc(sizeof(struct codeword));
            newCodeword->next = NULL;
            newCodeword->codeword = codeword;
            queryPosition->codewords = newCodeword;

            // Re-sort the lists of query positions from longest to shortest
            qsort(qPosList_qPosLists, qPosList_numQPosLists,
                  sizeof(struct memSingleBlock), qPosList_compareList);

            return;
        }

        listCount++;
    }

    // Instead use a new list of query positions
    list = qPosList_qPosLists + qPosList_numQPosLists;
    list->numEntries = 0;
    qPosList_numQPosLists++;

    // And copy values into it
    while (numQueryPositions > 0)
    {
        numQueryPositions--;
        queryPosition = memSingleBlock_newEntry(list);
        queryPosition->queryPosition = queryPositions[numQueryPositions];
        // No refering codeword for any of the positions except the last
        queryPosition->codewords = NULL;
    }

    // Reference at the last query position (will become the first) to the
    // new query position list's codeword
    newCodeword = global_malloc(sizeof(struct codeword));
    newCodeword->next = NULL;
    newCodeword->codeword = codeword;
    queryPosition->codewords = newCodeword;

    // Sort the lists from longest to shortest
    qsort(qPosList_qPosLists, qPosList_numQPosLists,
          sizeof(struct memSingleBlock), qPosList_compareList);
}

// Compare the two query position list lengths. Return 1 if list1 length > list2 length
// -1 if list 1 length < list 2 length and 0 if they are equal
int4 qPosList_compareList(const void* list1, const void* list2)
{
	const struct memSingleBlock *l1, *l2;

	l1 = (struct memSingleBlock*)list1;
	l2 = (struct memSingleBlock*)list2;

	if (l1->numEntries < l2->numEntries)
	{
		return 1;
	}
	else if (l1->numEntries > l2->numEntries)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

