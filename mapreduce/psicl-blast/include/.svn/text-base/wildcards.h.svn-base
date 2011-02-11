#ifndef _wildcards_
#define _wildcards_

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

struct clusterWildcard
{
	char* letters;
	uint4 wildCode;
	int2 scoreMatrixRow[encoding_aaStartWildcards];
    float averageScore;
};

struct wild
{
	uint4 code;
    uint4 count;
};

#define setbit(i, pos) (i |= (1 << pos))
#define getbit(i, pos) (i & (1 << pos))
#define wildcards_numClusterWildcards 7

extern struct scoreMatrix wildcards_scoreMatrix;
extern int wildcards_printWilds;
extern float wildcards_scoringConstant;
extern int2 wildcards_scoreMatrixRow[encoding_aaStartWildcards];

extern struct clusterWildcard wildcards_clusterWildcards[wildcards_numClusterWildcards];

// Initialize wildCodes in clusterWildcards structure
void wildcards_initialize();
// Print a wildcard character set
void wildcards_printWildcard(uint4 wildcard);
// Build a scoring row for wildcard candidate and calculate average score
float wildcards_scoreCandidates(struct wild* wildCandidates, uint4 numWildCandidates,
                                  struct wild* wilds, uint4 numWilds, float defaultWildscore);
// Given a list of wildcard
struct wild* wildcards_getSubset(struct wild wildCandidate, struct wild* wilds,
                                   uint4 numWilds, uint4* sizeWildSubset, uint4* numOccurences);
// Calculate the average score for aligning a residue to the wildcard candidate
float wildcards_averageResidueWildMatch(struct wild wildCandidate, struct wild* wilds,
                                          uint4 numWilds);
// Calculate the average score for aligning residue 'code' with given wildcard candidate
float wildcards_scoreResidueWildMatch(struct wild wildCandidate, struct wild* wilds,
                                        uint4 numWilds, uint4 code);
// Join two sets together
struct wild* wildcards_joinSubset(struct wild* set1, uint4 *size1, struct wild* set2, uint4 size2);

// Remove from set1 any wildcards in set2
void wildcards_removeSubset(struct wild* set1, uint4 *size1,
                              struct wild* set2, uint4 size2, uint4* numOccurences);
// Print out the final wildcards
void wildcards_outputWildcards(char* filename);
// Read in a set of wildcards
void wildcards_readWildcards(char* filename);

//Initialize for counting occurences of wildcards
void wildcards_initializeCountOccurences(uint4 longestSequenceLength);
// Process set of children and add to occurence counts
void wildcards_countOccurences(struct child* children, uint4 numChildren, uint4 sequenceLength);
// Get final list of number of occurences of each wild
struct wild* wildcards_getOccurences(uint4 *numWilds);

#ifdef __MY_EXTERN_C__
}
#endif

#endif
