#ifndef _search_
#define _search_

#ifdef __cplusplus
extern "C" {
#endif

// Search a protein database using 1-hit extension mode
void search_protein1hitParallel(struct scoreMatrix *scoreMatrixp,
								struct PSSMatrixFP PSSMatrixFP, 
								struct sequenceData* sequenceData,
								uint4 numSequences, 
								uint4 tickFrequency);

void search_protein2hitParallel(struct scoreMatrix *scoreMatrixp,
								struct PSSMatrix PSSMatrix,
								struct PSSMatrixFP PSSMatrixFP, 
								struct sequenceData* sequenceData,
								uint4 numSequences, 
								uint4 tickFrequency);

void search_protein1hit(struct PSSMatrix PSSMatrix, 
						struct PSSMatrixFP PSSMatrixFP, 
						struct sequenceData* sequenceData,
						uint4 numSequences, 
						uint4 tickFrequency);

// Search a protein database using 2-hit extension mode
void search_protein2hit(struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, struct sequenceData* sequenceData,
                       uint4 numSequences, uint4 tickFrequency);
// Search a nucleotide database using 1-hit extension mode
void search_nucleotide(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                       uint4 numSequences, uint4 tickFrequency);
// Search a nucleotide database using 1-hit extension mode with large wordsize > 14
void search_nucleotide_longWord(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                                uint4 numSequences, uint4 tickFrequency);
// Search a nucleotide database using 1-hit extension mode, using a large word lookup table
// due to long query sequence
void search_nucleotide_largeTable(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                                  uint4 numSequences, uint4 tickFrequency);
void search_nucleotideIndex(unsigned char* address, unsigned char* endAddress, int4 tickFrequency,
                            struct PSSMatrix PSSMatrix);
// SSearch a protein database using Smith-Waterman algorithm
void search_proteinSsearch(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                           uint4 numSequences, uint4 tickFrequency);
// Ssearch a nucleotide database using Smith-waterman algorithm
void search_nucleotideSsearch(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                              uint4 numSequences, uint4 tickFrequency);
#ifdef __cplusplus
}
#endif

#endif
