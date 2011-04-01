#ifndef _hitMatrix_
#define _hitMatrix_

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned char** hitMatrix_furthest;

//Shucai
extern uint4 *hitMatrix_furthestFP;

extern int4 hitMatrix_queryLength;
extern uint4 hitMatrix_diagonalMask;

// Initialize the hitMatrix by declaring memory for the maximum number
// of diagonals required by a subject sequence
void hitMatrix_initialize(int4 queryLength, int4 maximumSubjectLength, unsigned char* startAddress);

// Reinitialize the hit matrix so all values point to start of new file
void hitMatrix_reinitialize(int4 queryLength, int4 maximumSubjectLength, unsigned char* startAddress);

// Free the matrix
void hitMatrix_free();
#ifdef __cplusplus
}
#endif

#endif

