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
struct __attribute__ ((aligned )) groupFP   //for parallelization
{
	cl_int nextWords;
	cl_int nextGroups;
};

extern ushort* wordLookupDFA_additionalQueryPositions;
extern int wordLookupDFA_numAdditionalQueryPositions;
//Shucai
extern uint wordLookupDFA_numExtPositions;
extern uint wordLookupDFA_numWords;
extern uint wordLookupDFA_numGroups;

extern struct group *wordLookupDFA_groups;
//Shucai
extern struct groupFP *wordLookupDFA_groupsFP;

extern int wordLookupDFA_numCodes, wordLookupDFA_wordLength, wordLookupDFA_numBlocks;

// Build the word-lookup structure
void wordLookupDFA_build(struct PSSMatrix PSSMatrix, int numCodes, int wordLength);

// Print the contents of the word lookup table
void wordLookupDFA_print();

// Free memory used by the word lookup table
void wordLookupDFA_free();

#ifdef __cplusplus
}
#endif

#endif

