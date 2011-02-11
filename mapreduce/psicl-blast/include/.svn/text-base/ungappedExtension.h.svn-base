
#ifndef _ungappedExtension_
#define _ungappedExtension_

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned char* ungappedExtension_subjectEndReached;
//Shucai
extern uint4 ungappedExtension_subjectEndReachedFP;

extern int4 ungappedExtension_bestScore;
extern struct memBlocks* ungappedExtension_extensions;

// Scoring status constants
#define ungappedExtension_DELETED 0
#define ungappedExtension_UNGAPPED 1
#define ungappedExtension_SEMIGAPPED 2
#define ungappedExtension_GAPPED 3
#define ungappedExtension_JOINED 4

// Initialize the creation of ungapped extensions
void ungappedExtension_initialize();

// Perform an ungapped extension between point4s queryStart,subjectStart and queryEnd,subjectEnd
// and extend in each direction until score drops below best score yet minus a dropoff parameter
//Shucai
//struct ungappedExtension* ungappedExtension_extend(int2** queryHit, unsigned char* subjectHit,
//struct ungappedExtension* ungappedExtension_extend(int2* queryHit, unsigned char* subjectHit,
struct ungappedExtension* ungappedExtension_extend(int4 queryoffset, unsigned char* subjectHit,
//	unsigned char* lastHit, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, unsigned char* subject);
    uint4 lastHitFP, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, unsigned char* subject, 
	unsigned char *startAddressFP);

// Perform an ungapped extension when the seed is only a single hit on the diagonal, rather
// than a pair of hits.
// Shucai
//struct ungappedExtension* ungappedExtension_oneHitExtend(int2** queryHit,
//struct ungappedExtension* ungappedExtension_oneHitExtend(int2* queryHit,
struct ungappedExtension* ungappedExtension_oneHitExtend(int4 queryoffset,
unsigned char* subjectHit, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, unsigned char* subject, unsigned char *startAddressFP);

// Perform one-hit seeded ungapped extension for nucleotide, 1 packed-byte at a time
struct ungappedExtension* ungappedExtension_nucleotideExtend(int4 queryHitOffset,
	int4 subjectHitOffset, struct PSSMatrix PSSMatrix, unsigned char* subject,
    uint4 subjectLength);

// Find seed point4 for an ungapped extension
extern inline void ungappedExtension_findSeed(struct ungappedExtension* ungappedExtension,
		                                      struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, unsigned char* subject);

void ungappedExtension_print(struct ungappedExtension* extension);

#ifdef __cplusplus
}
#endif

#endif
