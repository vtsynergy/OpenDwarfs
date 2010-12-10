#ifndef _PSSMatrix_
#define _PSSMatrix_

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned char* PSSMatrix_packedRightMatches;
extern unsigned char* PSSMatrix_packedLeftMatches;
extern int2* PSSMatrix_packedRightMatchScores;
extern int2* PSSMatrix_packedLeftMatchScores;
extern char* PSSMatrix_packedScore;

// Create a PSSM for the given query sequence and score matrix.
// The PSSMatrix will have length(query) columns and 25 rows
struct PSSMatrix PSSMatrix_create(struct scoreMatrix scoreMatrix, char* query);

//Shucai
// Transform the PSSMatrix to layout compatible with GPU
struct PSSMatrixFP PSSMatrixFP_transform(struct PSSMatrix *PSSMatrixptr);

// Returns a PSSMatrix with the first "amount" entries removed and the length shortened
// accordingly
struct PSSMatrix PSSMatrix_chop(struct PSSMatrix PSSMatrix, int4 amount);

// Returns a PSSMatrix which is a reversed copy
struct PSSMatrix PSSMatrix_reverse(struct PSSMatrix PSSMatrix);

// Calculate the start of strand for the given query offset
uint4 PSSMatrix_strandStart(struct PSSMatrix PSSMatrix, uint4 queryOffset);

// Calculate the end of strand for the given query offset
uint4 PSSMatrix_strandEnd(struct PSSMatrix PSSMatrix, uint4 queryOffset);

// Print4 the contents of the PSSMatrix
void PSSMatrix_print(struct PSSMatrix PSSMatrix);

// Free memory used by the matrix
void PSSMatrix_free(struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP);

// Free a copy of the PSSMatrix
void PSSMatrix_freeCopy(struct PSSMatrix PSSMatrix);

#ifdef __cplusplus
}
#endif

#endif

