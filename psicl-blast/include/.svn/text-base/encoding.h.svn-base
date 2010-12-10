#ifndef _encoding_
#define _encoding_

#ifdef __cplusplus
extern "C" {
#endif

struct wildcard
{
	unsigned char letter;
    unsigned char* replacements;
    unsigned char* replacementCodes;
    int4 numCodes;
    int4 replacementCounter;
};

struct wildcardEdit
{
	unsigned char code;
    int4 position;
};

extern struct packedByte* encoding_packedByteLookup;

extern unsigned char* encoding_codesArray;
extern unsigned char* encoding_lettersArray;
extern unsigned char encoding_unknownCode, encoding_sentinalCode;
extern unsigned char encoding_numLetters, encoding_numRegularLetters, encoding_numCodes;
extern unsigned char encoding_alphabetType;
extern struct wildcard* encoding_wildcards;

extern char* encoding_alphabetTypes[2];

#define encoding_nucleotide 1
#define encoding_protein 0
#define encoding_aaStartWildcards 24

// Initialize encoding functions
void encoding_initialize(unsigned char alphabetType);

// Given a nucleotide code, returns the code of complement
unsigned char encoding_getComplement(unsigned char code);

// Given an amino acid or nucleotide character, returns the code
unsigned char encoding_getCode(unsigned char letter);

// Given a code, returns the corresponding amino acid or nucleotide character, or '_' if not valid
unsigned char encoding_getLetter(unsigned char code);

// Given a sequence, determine the alphabet type
unsigned char encoding_determineAlphabetType(char* sequence, uint4 sequenceSize);

// Generate a random regular letter for a given wildcard
unsigned char encoding_randomEncodedLetter(unsigned char code);

// Insert wildcards back int4o the sequence
extern inline void encoding_insertWilds(unsigned char* subject, unsigned char* edits,
                                        unsigned char* endEdits);

// Extract a single character from a packed byte
#define encoding_extractBase(byte, bytePosition) ((byte >> (6 - (bytePosition * 2))) & 0x3)

// Unpack a sequence from byte-packed form
unsigned char* encoding_byteUnpack(unsigned char* bytePackedSequence, int4 sequenceLength);

// Byte pack fourth letters
extern inline unsigned char encoding_bytePack(unsigned char* sequence);

// Byte pack the last 1 to 4 characters in a sequence
extern inline unsigned char encoding_bytePackRemaining(unsigned char* sequence, int4 numLetters);

// Byte pack the first 1 to 4 characters in a sequence
extern inline unsigned char encoding_bytePackBeginning(unsigned char* sequence, int4 numLetters);

// Replace the wildcards in a protein or nucleotide sequence
int4 encoding_replaceWildcards(struct memSingleBlock* wildcardEdits, unsigned char* sequence,
                              int4 sequenceSize);

// Byte-pack a nucleotide sequence and return the packed size
int4 encoding_bytePackSequence(unsigned char* sequence, int4 sequenceSize);

// Given a byte-packed code, prints the 4 letters
void encoding_printLetters(unsigned char code, int4 numLetters);

// Given a sequence and its length, encode it
void encoding_encodeSequence(char* sequence, int4 sequenceSize, uint4 alphabetType);

// Free structures used to encoding/decode
void encoding_free();

// Unpack part of a sequence
unsigned char* encoding_byteUnpackRegion(unsigned char* subject, unsigned char* bytePackedSequence,
                                         int4 sequenceLength);
#ifdef __cplusplus
}
#endif

#endif
