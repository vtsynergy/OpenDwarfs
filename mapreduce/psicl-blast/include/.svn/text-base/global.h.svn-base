#ifndef _global_
#define _global_

#ifdef __cplusplus
extern "C" {
#endif

#define maximum(a,b) ((a > b) ? a : b)
#define minimum(a,b) ((a < b) ? a : b)

// BLAST datatypes

// This structure contains values for a scoring matrix such as BLOSUM62
struct scoreMatrix
{
	int2** matrix;
	int4 highestValue;
	int4 lowestValue;
    float averageMatchScore;
};

// This structure contains a PSSM (Position Specific Scoring Matrix) which has
// a column for each amino acid in original query sequence (of length "length")
// and 25 rows.
struct PSSMatrix
{
    int4 length;
    int4 strandLength;
	int4 highestValue;
	int4 lowestValue;
	int2** matrix;
	unsigned char* queryCodes;
    unsigned char* bestMatchCodes;
    unsigned char* bytePackedCodes;
	unsigned char* xorCodes;
};

// Shucai
// This structure is the PSSM structure used for parallelization on GPU
struct PSSMatrixFP
{
	int4 length;
	int4 strandLength;
	int4 highestValue;
	int4 lowestValue;
	int2* matrix;
    unsigned char* queryCodes;
	unsigned char* bestMatchCodes;
	unsigned char* bytePackedCodes;
	unsigned char* xorCodes;
};

//Shucai
//This stucture is for time recording
typedef struct strTimeRecord {
	uint4 iniTime;
	uint4 preProcessTime;
	uint4 dataCopyTimeH2D;
	uint4 searchTime;
	uint4 dataCopyTimeD2H;
	uint4 addUngappedExtensionTime;
	uint4 postProcessTime;
	uint4 hitUngappedExtTime;
	uint4 gappedAlignmentTime;
	uint4 finalAlignmentTime;
	uint4 totalTime;
} TIMERECORD;

// A query/subject coordinate pair
struct coordinate
{
    int4 queryOffset;
	int4 subjectOffset;
};

// Information about an ungapped extension
struct ungappedExtension
{
	struct coordinate start;
	struct coordinate end;
	struct coordinate seed;
	int4 nominalScore;
	//Shucai
	uint4 sequenceCount;
	uint4 tid;
	char status;
	struct ungappedExtension* next;
};

// A region either copied or unpacked by blast
struct unpackRegion
{
	int4 startOffset;
    int4 endOffset;
	unsigned char* unpackedSubject;
    unsigned char* subject;
    int4 subjectLength;
};

// Result of dynamic programming
struct dpResults
{
	struct coordinate best;
	uint4 bestScore;
	unsigned char** traceback;
};

// An alignment trace resulting from dynamic programming
struct trace
{
	uint4 length;
	uint4 queryStart;
	uint4 subjectStart;
	unsigned char* traceCodes;
};

// A gapped alignment
struct gappedExtension
{
	struct trace trace;
	int4 nominalScore;
	int4 queryEnd;
	int4 subjectEnd;
	float normalizedScore;
	double eValue;
	struct gappedExtension* next;
};

// Information about the alignment between the query and a subject
struct alignment
{
	int4 descriptionLocation;
    int4 descriptionLength;
	unsigned char* subject;
	int4 subjectLength;
	struct ungappedExtension* ungappedExtensions;
	struct gappedExtension* gappedExtensions;
    unsigned char* edits;
    int4 encodedLength;
	char joinChecked;
    char inMemorySubject;
    struct unpackRegion* unpackRegions;
    uint4 numUnpackRegions;
    uint4 cluster;
};

// A final alignment above cutoff
struct finalAlignment
{
	int4 highestNominalScore;
    char* description;
	struct alignment* alignment;
};

// Timing variables
extern int4 blast_prepTime, blast_searchTime;
extern int4 blast_gappedScoreTime, blast_gappedExtendTime, blast_finalizeTime;
extern int4 blast_semiGappedScoreTime, blast_copyTime, blast_unpackTime;

// BLAST statistics
extern uint4 blast_numHits;
extern uint4 blast_numUngappedExtensions, blast_numTriggerExtensions, blast_numTriggerSequences;
extern uint4 blast_numGapped;
extern uint4 blast_numSemiGapped;
extern uint4 blast_numExtensionsPruned;
extern uint4 blast_numAttemptedJoin, blast_numSuccessfullyJoined;
extern uint4 blast_numGoodAlignments;
extern uint4 blast_numGoodExtensions;
extern uint4 blast_totalUnpacked;
extern uint4 blast_totalCopied;
extern uint4 blast_numExpandedSequences;

// BLAST global variables
extern int4 blast_ungappedNominalTrigger;
extern int4 blast_gappedNominalCutoff;
extern int4 blast_gappedNominalCutoff, blast_nominalR1cutoff, blast_nominalR2cutoff;
extern int4 blast_dynamicGappedNominalCutoff, blast_dynamicNominalR1cutoff;
extern int4 blast_dloc;
extern char* blast_queryDescription;

//Shucai
extern TIMERECORD timeRecord;

// Initialize global variables
void global_initialize();
// Convert a 32-bit int4eger int4o a string with commas
char* global_int4toString(uint4 number);
// Convert a 64-bit int4eger int4o a string with commas
char* global_int8toString(uint8 number);

extern uint4 global_totalMalloc;

// Malloc new memory, check that malloc was successful
void* global_malloc(size_t size);

// Realloc memory, check that realloc was successful
void* global_realloc(void* ptr, size_t size);

// Free the global convert string
void global_free();

#ifdef __cplusplus
}
#endif

#endif

