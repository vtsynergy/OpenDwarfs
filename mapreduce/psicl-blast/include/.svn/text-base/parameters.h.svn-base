#ifndef _parameters_
#define _parameters_

#ifdef __cplusplus
extern "C" {
#endif
//Shucai
extern int parameters_blockNum;
extern int parameters_threadNum;

extern char* parameters_queryFile;
extern char* parameters_subjectDatabaseFile;
extern char* parameters_hspindexFilename;
extern int4 parameters_numDisplayAlignments;
extern int4 parameters_numDisplayTracebacks;
extern int4 parameters_databaseSize;
extern int4 parameters_numberOfSequences;
extern char parameters_outputType;
#define parameters_pairwise 0
#define parameters_xml 7
#define parameters_tabular 8

extern char parameters_scoringSystem;
extern char parameters_semiGappedScoring;
extern char parameters_restrictedInsertionScoring;
extern char parameters_bytepackedScoring;
extern char parameters_tableScoring;
extern char parameters_oneHitTrigger;
extern char parameters_getDescriptions;
extern char parameters_ssearch;
extern char parameters_strands;
extern char parameters_useIndex;
extern char parameters_filterEnabled;
extern char parameters_allClusterMembers;

extern char* parameters_scoringMatrix;
extern char* parameters_scoringMatrixPath;
extern int4 parameters_mismatchScore;
extern int4 parameters_matchScore;

extern char parameters_wordSize;
extern int4 parameters_T;
extern int4 parameters_A;
extern char parameters_overlap;
extern int2 parameters_startGap;
extern int2 parameters_extendGap;
extern int2 parameters_openGap;

extern int4 parameters_wordTableBytes;
extern int4 parameters_wordTableLetters;
extern int4 parameters_wordExtraLetters;
extern int4 parameters_wordExtraBytes;

extern double parameters_cutoff;

extern float parameters_ungappedNormalizedTrigger;
extern float parameters_ungappedNormalizedDropoff;
extern float parameters_gappedNormalizedDropoff;
extern float parameters_gappedFinalNormalizedDropoff;

extern float parameters_semiGappedR1;
extern float parameters_semiGappedR2;
extern int4 parameters_semiGappedExtensionN;
extern int2 parameters_semiGappedStartGap;
extern int2 parameters_semiGappedExtendGap;
extern int2 parameters_semiGappedOpenGap;
extern int4 parameters_semiGappedDropoffIncrease;

extern int2 parameters_bytepackStartGap;
extern int2 parameters_bytepackOpenGap;
extern int2 parameters_bytepackExtendGap;
extern int2 parameters_bytepackOpen4Gap;
extern int2 parameters_bytepackExtend4Gap;
extern int4 parameters_bytepackDropoffDecrease;

extern int4 parameters_verboseDloc;

// Process command line arguments and set parameters value
void parameters_processArguments(int4 argc, char* argv[]);

// Given the alphabet type of the query, use the relevant parameter defaults
void parameters_loadDefaults(unsigned char alphabetType);

// Load a BLOSUM or PAM scoring matrix
void parameters_findScoringMatrix();

void parameters_free();

#ifdef __cplusplus
}
#endif

#endif
