#ifndef _wordLookupDFA_
#define _wordLookupDFA_

#ifdef __cplusplus
extern "C" {
#endif

struct group
{
	unsigned char* nextWords;
    struct group* nextGroups;
};

//Shucai
struct groupFP   //for parallelization
{
	int4 nextWords;
	int4 nextGroups;
};

extern uint2* wordLookupDFA_additionalQueryPositions;
extern int4 wordLookupDFA_numAdditionalQueryPositions;
//Shucai
extern uint4 wordLookupDFA_numExtPositions;
extern uint4 wordLookupDFA_numWords;
extern uint4 wordLookupDFA_numGroups;

extern struct group *wordLookupDFA_groups;
//Shucai
extern struct groupFP *wordLookupDFA_groupsFP;

extern int4 wordLookupDFA_numCodes, wordLookupDFA_wordLength, wordLookupDFA_numBlocks;

// Build the word-lookup structure
void wordLookupDFA_build(struct PSSMatrix PSSMatrix, int4 numCodes, int4 wordLength);

// Print4 the contents of the word lookup table
void wordLookupDFA_print();

// Free memory used by the word lookup table
void wordLookupDFA_free();

#ifdef __cplusplus
}
#endif

#endif

