#ifndef _nucleotideLookup_
#define _nucleotideLookup_

#ifdef __cplusplus
extern "C" {
#endif

extern int2* nucleotideLookup;
extern uint2* nucleotideLookup_additionalPositions;
extern uint4 nucleotideLookup_numAdditionalPositions;
extern uint4 *nucleotideLookup_bitLookup, nucleotideLookup_mask;
extern int4 nucleotideLookup_packedWordLength;

extern int4* nucleotideLookup_large;
extern uint4* nucleotideLookup_additionalPositions_large;
extern char nucleotideLookup_largeTable;

// Build the nucleotide lookup table
void nucleotideLookup_build(struct PSSMatrix PSSMatrix, int4 packedWordLength);

// Print4 the lookup table
void nucleotideLookup_print();

// Free the table
void nucleotideLookup_free();
#ifdef __cplusplus
}
#endif

#endif

