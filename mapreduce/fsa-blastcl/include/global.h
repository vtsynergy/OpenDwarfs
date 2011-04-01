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
	short** matrix;
	int highestValue;
	int lowestValue;
    float averageMatchScore;
};

// This structure contains a PSSM (Position Specific Scoring Matrix) which has
// a column for each amino acid in original query sequence (of length "length")
// and 25 rows.
struct PSSMatrix
{
    int length;
    int strandLength;
	int highestValue;
	int lowestValue;
	short** matrix;
	unsigned char* queryCodes;
    unsigned char* bestMatchCodes;
    unsigned char* bytePackedCodes;
	unsigned char* xorCodes;
};

// Shucai
// This structure is the PSSM structure used for parallelization on GPU
struct __attribute__((aligned )) PSSMatrixFP 
{
	cl_int length;
	cl_int strandLength;
	cl_int highestValue;
	cl_int lowestValue;
	cl_short* matrix;
    unsigned char* queryCodes;
	unsigned char* bestMatchCodes;
	unsigned char* bytePackedCodes;
	unsigned char* xorCodes;
};

//Shucai
//This stucture is for time recording
typedef struct __attribute__((aligned )) strTimeRecord {
	cl_uint iniTime;
	cl_uint preProcessTime;
	cl_uint dataCopyTimeH2D;
	cl_uint searchTime;
	cl_uint dataCopyTimeD2H;
	cl_uint addUngappedExtensionTime;
	cl_uint postProcessTime;
	cl_uint hitUngappedExtTime;
	cl_uint gappedAlignmentTime;
	cl_uint finalAlignmentTime;
	cl_uint totalTime;
} TIMERECORD;

// A query/subject coordinate pair
struct __attribute__ ((aligned )) coordinate
{
    cl_int queryOffset;
	cl_int subjectOffset;
};

// Information about an ungapped extension
struct __attribute__ ((aligned )) ungappedExtension
{
	cl_int nominalScore;
	//Shucai
	cl_uint sequenceCount;
	cl_uint tid;
	struct ungappedExtension* next;
	struct coordinate start;
	struct coordinate end;
	struct coordinate seed;
	int status;
//	char pad[11];
};

// A region either copied or unpacked by blast
struct unpackRegion
{
	int startOffset;
    int endOffset;
	unsigned char* unpackedSubject;
    unsigned char* subject;
    int subjectLength;
};

// Result of dynamic programming
struct dpResults
{
	cl_uint bestScore;
	struct coordinate best;
	unsigned char** traceback;
};

// An alignment trace resulting from dynamic programming
struct trace
{
	uint length;
	uint queryStart;
	uint subjectStart;
	unsigned char* traceCodes;
};

// A gapped alignment
struct gappedExtension
{
	cl_int nominalScore;
	cl_int queryEnd;
	cl_int subjectEnd;
	float normalizedScore;
	double eValue;
	struct trace trace;
	struct gappedExtension* next;
};

// Information about the alignment between the query and a subject
struct alignment
{
	int descriptionLocation;
    int descriptionLength;
	unsigned char* subject;
	int subjectLength;
	struct ungappedExtension* ungappedExtensions;
	struct gappedExtension* gappedExtensions;
    unsigned char* edits;
    int encodedLength;
	char joinChecked;
    char inMemorySubject;
    struct unpackRegion* unpackRegions;
    uint numUnpackRegions;
    uint cluster;
};

// A final alignment above cutoff
struct finalAlignment
{
	int highestNominalScore;
    char* description;
	struct alignment* alignment;
};

// Timing variables
extern int blast_prepTime, blast_searchTime;
extern int blast_gappedScoreTime, blast_gappedExtendTime, blast_finalizeTime;
extern int blast_semiGappedScoreTime, blast_copyTime, blast_unpackTime;

// BLAST statistics
extern uint blast_numHits;
extern uint blast_numUngappedExtensions, blast_numTriggerExtensions, blast_numTriggerSequences;
extern uint blast_numGapped;
extern uint blast_numSemiGapped;
extern uint blast_numExtensionsPruned;
extern uint blast_numAttemptedJoin, blast_numSuccessfullyJoined;
extern uint blast_numGoodAlignments;
extern uint blast_numGoodExtensions;
extern uint blast_totalUnpacked;
extern uint blast_totalCopied;
extern uint blast_numExpandedSequences;

// BLAST global variables
extern int blast_ungappedNominalTrigger;
extern int blast_gappedNominalCutoff;
extern int blast_gappedNominalCutoff, blast_nominalR1cutoff, blast_nominalR2cutoff;
extern int blast_dynamicGappedNominalCutoff, blast_dynamicNominalR1cutoff;
extern int blast_dloc;
extern char* blast_queryDescription;

//Shucai
extern TIMERECORD timeRecord;

// Initialize global variables
void global_initialize();
// Convert a 32-bit integer into a string with commas
char* global_inttoString(uint number);
// Convert a 64-bit integer into a string with commas
char* global_int8toString(uint8 number);

extern uint global_totalMalloc;

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

