#ifndef _scoreMatrix_
#define _scoreMatrix_

#ifdef __cplusplus
extern "C" {
#endif

// Load the score matrix (eg. BLOSUM62) from disk and return contents in an array
// 25 by 25 for the 25 possible amino acids (actually 20 plus 3 wilds, 1 unknown,
// and a sentinal code which scores poorly, and flanks sequences)
struct scoreMatrix scoreMatrix_load(char* filename);

// Create a nucleotide scoring matrix use match and mismatch penalties
struct scoreMatrix scoreMatrix_create(int2 match, int2 mismatch);

// Print4 the contents of the score matrix
void scoreMatrix_print(struct scoreMatrix scoreMatrix);

// Free memory used by the score matrix
void scoreMatrix_free(struct scoreMatrix scoreMatrix);

#ifdef __cplusplus
}
#endif

#endif

