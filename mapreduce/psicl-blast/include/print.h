#ifndef _print_
#define _print_

#ifdef __cplusplus
extern "C" {
#endif

// Format a sequence description so it wraps around lines nicely
char* print_formatDescription(char* description, int4 firstLineIndent, int4 remainingLinesIndent,
                              int4 maxLineLength);

// Print4 full gapped alignments
void print_gappedAlignmentsFull(char* query, struct PSSMatrix PSSMatrix);

// Print4 1 line description of gapped alignments
void print_gappedAlignmentsBrief();

// Print4 a single sequence
void print_singleSequence(unsigned char* sequence, int4 length);

// Convert an e-value to a string representation, using the same
// technique as NCBI-BLAST
char* print_eValue2String(double evalue);

// Header for XML output
void print_XMLheader(char* query, struct PSSMatrix PSSMatrix);

// Footer of XML output
void print_XMLfooter();

//Shucai
void print_seqData(struct sequenceData *, uint4);

void print_seqDataFP(struct sequenceDataFP *, uint4);

#ifdef __cplusplus
}
#endif

#endif
