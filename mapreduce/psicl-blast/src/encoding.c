// encoding.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Functions for converting ASCII amino acid and nucleotide characters and sequences
// to and from a byte encoded form

#include <ctype.h>
#include <stdlib.h>
#include "blast.h"

#ifdef __cplusplus
extern "C"
#endif

struct packedByte* encoding_packedByteLookup = NULL;

unsigned char* encoding_codesArray;
unsigned char* encoding_lettersArray;
unsigned char encoding_unknownCode, encoding_sentinalCode;
unsigned char encoding_numLetters, encoding_numRegularLetters, encoding_numCodes;
unsigned char encoding_alphabetType;
struct wildcard* encoding_wildcards;
unsigned char* encoding_complements;

char* encoding_alphabetTypes[2] = {"Protein", "Nucleotide"};

// Declare array which translates codes (bytes) back into amino acid characters
#define encoding_aaUnknownCode 23
#define encoding_aaSentinalCode 31
#define encoding_aaNumLetters 24
#define encoding_aaNumRegularLetters 20
#define encoding_aaNumCodes 32
unsigned char encoding_aaLettersArray[encoding_aaNumCodes] =
{'L','A','G','S','V','E','T','K','D','P','I','R','N',
'Q','F','Y','M','H','C','W','B','Z','X','U','1','2','3','4','5','6','#','!'};

#define encoding_nuUnknownCode 14
#define encoding_nuSentinalCode 15
#define encoding_nuNumLetters 15
#define encoding_nuNumRegularLetters 4
#define encoding_nuNumCodes 16
unsigned char encoding_nuLettersArray[encoding_nuNumCodes] =
{'A','C','G','T','N','B','D','H','K','M','R','S',
'V','W','Y','!'};
unsigned char encoding_nuComplementsArray[encoding_nuNumCodes] =
{'T','G','C','A','N','V','H','D','M','K','Y','S',
'B','W','R','!'};

struct packedByte
{
	unsigned char codes[4];
};

struct wildcard encoding_nuWildcards[encoding_nuNumCodes] =
{{'A', "A", NULL, 1},
 {'C', "C", NULL, 1},
 {'G', "G", NULL, 1},
 {'T', "T", NULL, 1},
 {'N', "ACGT", NULL, 4},
 {'B', "CGT", NULL, 3},
 {'D', "AGT", NULL, 3},
 {'H', "ACT", NULL, 3},
 {'K', "GT", NULL, 2},
 {'M', "AC", NULL, 2},
 {'R', "AG", NULL, 2},
 {'S', "CG", NULL, 2},
 {'V', "ACG", NULL, 3},
 {'W', "AT", NULL, 2},
 {'Y', "CT", NULL, 2},
 {'!', "", NULL, 0}};
struct wildcard encoding_aaWildcards[encoding_aaNumCodes] =
{{'L', "L", NULL, 1},
 {'A', "A", NULL, 1},
 {'G', "G", NULL, 1},
 {'S', "S", NULL, 1},
 {'V', "V", NULL, 1},
 {'E', "E", NULL, 1},
 {'T', "T", NULL, 1},
 {'K', "K", NULL, 1},
 {'D', "D", NULL, 1},
 {'P', "P", NULL, 1},
 {'I', "I", NULL, 1},
 {'R', "R", NULL, 1},
 {'N', "N", NULL, 1},
 {'Q', "Q", NULL, 1},
 {'F', "F", NULL, 1},
 {'Y', "Y", NULL, 1},
 {'M', "M", NULL, 1},
 {'H', "H", NULL, 1},
 {'C', "C", NULL, 1},
 {'W', "W", NULL, 1},
 {'B', "ND", NULL, 2},
 {'Z', "QE", NULL, 2},
 {'X', "LAGSVETKDPIRNQFYMHCW", NULL, 20},
 {'U', "LAGSVETKDPIRNQFYMHCW", NULL, 20},
 {'!', "", NULL, 0}};

// Initialize encoding functions
void encoding_initialize(unsigned char alphabetType)
{
	unsigned char letter = 0, code = 0, count;
    int4 packedByte;

    // Record alphabet type
    encoding_alphabetType = alphabetType;

    if (encoding_alphabetType == encoding_protein)
    {
        // Use amino acid alphabet
        encoding_lettersArray = encoding_aaLettersArray;
        encoding_unknownCode = encoding_aaUnknownCode;
        encoding_sentinalCode = encoding_aaSentinalCode;
        encoding_numLetters = encoding_aaNumLetters;
        encoding_numRegularLetters = encoding_aaNumRegularLetters;
        encoding_numCodes = encoding_aaNumCodes;
		encoding_wildcards = encoding_aaWildcards;
    }
    else
    {
        // Use nucleotide alphabet
        encoding_lettersArray = encoding_nuLettersArray;
        encoding_unknownCode = encoding_nuUnknownCode;
        encoding_sentinalCode = encoding_nuSentinalCode;
        encoding_numLetters = encoding_nuNumLetters;
        encoding_numRegularLetters = encoding_nuNumRegularLetters;
        encoding_numCodes = encoding_nuNumCodes;
		encoding_wildcards = encoding_nuWildcards;

        // Construct byte unpack table
		encoding_packedByteLookup = (struct packedByte*)global_malloc(sizeof(struct packedByte) * 256);

        // For each packed byte
        packedByte = 0;
        while (packedByte < 256)
        {
        	// Store unpacked codes
			encoding_packedByteLookup[packedByte].codes[0] = (packedByte >> 6) & 0x3;
			encoding_packedByteLookup[packedByte].codes[1] = (packedByte >> 4) & 0x3;
			encoding_packedByteLookup[packedByte].codes[2] = (packedByte >> 2) & 0x3;
			encoding_packedByteLookup[packedByte].codes[3] = packedByte & 0x3;

            packedByte++;
        }
    }

	// Create array for translating amino acid characters into the related byte code
	encoding_codesArray = (unsigned char*)global_malloc(256);

	// Iterate for each possible character value
	while (letter < 255)
	{
		// By default use unknown code (eg. U or N)
		encoding_codesArray[letter] = encoding_unknownCode;

		// If a valid code (found by looking up lettersArray)
		code = 0;
		while (code < encoding_numCodes)
		{
			if (toupper(letter) == encoding_lettersArray[code])
			{
				// Add to codes array
				encoding_codesArray[letter] = code;
			}
			code++;
		}
		letter++;
	}

    // Initialize array for converting wildcards to regular letters
    code = 0;
    while (code < encoding_numLetters)
    {
        // Declare memory for regular codes
        encoding_wildcards[code].replacementCodes = (char*)global_malloc(sizeof(char) *
                                                    encoding_wildcards[code].numCodes);

		// Initialize replacement counter
        encoding_wildcards[code].replacementCounter = 0;

        // For each regular letter
        count = 0;
        while (count < encoding_wildcards[code].numCodes)
        {
            // Convert to encoded form
            encoding_wildcards[code].replacementCodes[count]
                = encoding_codesArray[encoding_wildcards[code].replacements[count]];

            count++;
        }
        code++;
    }

    // Construct complements array
    if (encoding_alphabetType == encoding_nucleotide)
    {
    	encoding_complements = (unsigned char*)global_malloc(sizeof(unsigned char) * encoding_numLetters);

        // For each code
        code = 0;
        while (code < encoding_numLetters)
        {
        	// Store the code of the complement
			encoding_complements[code] = encoding_getCode(encoding_nuComplementsArray[code]);

        	code++;
        }
    }
}

// Given a nucleotide code, returns the code of complement
unsigned char encoding_getComplement(unsigned char code)
{
	return encoding_complements[code];
}

// Given an amino acid or nucleotide character, returns the code
unsigned char encoding_getCode(unsigned char letter)
{
	return encoding_codesArray[letter];
}

// Given a code, returns the corresponding amino acid or nucleotide character, or '_' if not valid
unsigned char encoding_getLetter(unsigned char code)
{
	if (code < encoding_numCodes)
		return encoding_lettersArray[code];
	else
		return '_';
}

// Given a sequence, determine the alphabet type
unsigned char encoding_determineAlphabetType(char* sequence, uint4 sequenceSize)
{
	int4 letterCount, alphabetCount, matches = 0, wildchars = 0;
    char letter;

    // Go through the sequence and check it for non-nucleic characters
    letterCount = 0;
    while (letterCount < sequenceSize)
    {
        letter = sequence[letterCount];

        // For each regular nucleic letter
        alphabetCount = 0;
        while (alphabetCount <= encoding_nuNumRegularLetters)
        {
            // Check for a match
            if (toupper(letter) == encoding_nuLettersArray[alphabetCount])
                matches++;

            alphabetCount++;
        }

        // For each wildcard character
        while (alphabetCount < encoding_nuNumLetters)
        {
            // Check for a match
            if (toupper(letter) == encoding_nuLettersArray[alphabetCount])
                wildchars++;

            alphabetCount++;
        }

        letterCount++;
    }

//    printf("[%d,%d,%d]\n", matches, wildchars, sequenceSize);

    // If not all nucleotide chars or less than 90% are G,A,T,C or N this is a protein query
    if (matches + wildchars < sequenceSize || matches < sequenceSize * 0.9)
    {
        return encoding_protein;
    }

    return encoding_nucleotide;
}

// Generate a random regular letter for a given wildcard
unsigned char encoding_randomEncodedLetter(unsigned char code)
{
    // For protein alphabet use robinson&robinson frequencies
	if (encoding_alphabetType == encoding_protein)
    {
        encoding_wildcards[code].replacementCounter++;
        if (encoding_wildcards[code].replacementCounter >= encoding_wildcards[code].numCodes)
			encoding_wildcards[code].replacementCounter = 0;

    	return encoding_wildcards[code].replacementCodes[encoding_wildcards[code].replacementCounter];
    }
    else
    {
		encoding_wildcards[code].replacementCounter++;
        if (encoding_wildcards[code].replacementCounter >= encoding_wildcards[code].numCodes)
			encoding_wildcards[code].replacementCounter = 0;

    	return encoding_wildcards[code].replacementCodes[encoding_wildcards[code].replacementCounter];
	}
}

// Insert wildcards back into the sequence
void encoding_insertWilds(unsigned char* subject, unsigned char* edits,
                          unsigned char* endEdits)
{
    uint4 wildcardPosition;
	unsigned char wildcard;

    // For each edit
    while (edits < endEdits)
    {
    	// Read wildcard
		wildcard = *edits;
        edits++;

        // Read its position
		vbyte_getVbyte(edits, &wildcardPosition);

        // Add wildcard
        subject[wildcardPosition] = wildcard;
    }
}

// Unpack part of a sequence
unsigned char* encoding_byteUnpackRegion(unsigned char* subject, unsigned char* bytePackedSequence,
                                         int4 sequenceLength)
{
	unsigned char *subjectPosition;
    int4 packedCount = 0;
    int4 numPackedBytes, numRemaining;

    subjectPosition = subject;

    // Calculate packed length and remainder
    numPackedBytes = sequenceLength / 4;
    numRemaining = sequenceLength % 4;

    // For each packed byte
	while (packedCount < numPackedBytes)
    {
    	// Copy the value to the subject sequence
		memcpy(subjectPosition, encoding_packedByteLookup[bytePackedSequence[packedCount]].codes,
               sizeof(char) * 4);

        subjectPosition+=4;
        packedCount++;
    }

    // Unpack last 1-3 letters
    if (numRemaining)
    {
		if (numRemaining == 1)
        {
        	*subjectPosition = (bytePackedSequence[packedCount] >> 6);
            packedCount++;
        }
        else if (numRemaining == 2)
        {
        	*subjectPosition = (bytePackedSequence[packedCount] >> 6);
            subjectPosition++;
        	*subjectPosition = (bytePackedSequence[packedCount] >> 4) & 0x3;
            packedCount++;
        }
        else if (numRemaining == 3)
        {
        	*subjectPosition = (bytePackedSequence[packedCount] >> 6);
            subjectPosition++;
        	*subjectPosition = (bytePackedSequence[packedCount] >> 4) & 0x3;
            subjectPosition++;
        	*subjectPosition = (bytePackedSequence[packedCount] >> 2) & 0x3;
            packedCount++;
        }
    }

    return subject;
}

// Unpack an entire sequence
unsigned char* encoding_byteUnpack(unsigned char* bytePackedSequence, int4 sequenceLength)
{
	unsigned char* subject;

    subject = global_malloc(sizeof(char) * sequenceLength);
	encoding_byteUnpackRegion(subject, bytePackedSequence, sequenceLength);

    return subject;
}

#define encoding_unpackBase1(packedByte) (packedByte >> 6)
#define encoding_unpackBase2(packedByte) ((packedByte >> 4) & 0x3)
#define encoding_unpackBase3(packedByte) ((packedByte >> 2) & 0x3)
#define encoding_unpackBase4(packedByte) (packedByte & 0x3)

// Byte pack fourth letters
unsigned char encoding_bytePack(unsigned char* sequence)
{
	return (*sequence << 6) | (*(sequence + 1) << 4) | (*(sequence + 2) << 2) | *(sequence + 3);
}

// Byte pack the last 1 to 4 characters in a sequence
unsigned char encoding_bytePackRemaining(unsigned char* sequence, int4 numLetters)
{
	if (numLetters == 1)
    {
		return (*sequence << 6);
	}
    else if (numLetters == 2)
    {
		return (*(sequence + 1) << 4) | (*sequence << 6);
	}
    else if (numLetters == 3)
    {
		return (*(sequence + 2) << 2) | (*(sequence + 1) << 4) | (*sequence << 6);
	}
    else // numLetters = 4
    {
		return *(sequence + 3) | (*(sequence + 2) << 2) | (*(sequence + 1) << 4) | (*sequence << 6);
    }
}

// Byte pack the first 1 to 4 characters in a sequence
unsigned char encoding_bytePackBeginning(unsigned char* sequence, int4 numLetters)
{
	if (numLetters == 1)
    {
		return *sequence;
	}
    else if (numLetters == 2)
    {
    	return *(sequence + 1) | (*sequence << 2);
	}
    else if (numLetters == 3)
    {
    	return *(sequence + 2) | (*(sequence + 1) << 2) | (*sequence << 4);
	}
    else // numLetters = 4
    {
		return *(sequence + 3) | (*(sequence + 2) << 2) | (*(sequence + 1) << 4) | (*sequence << 6);
    }
}

// Replace the wildcards in a protein or nucleotide sequence
int4 encoding_replaceWildcards(struct memSingleBlock* wildcardEdits, unsigned char* sequence,
                              int4 sequenceSize)
{
	struct wildcardEdit* wildcardEdit;
	int4 count = 0;

    // Initialize list of wildcard edits
    wildcardEdits->numEntries = 0;

    // First scan through sequence and replace wilds
	while (count < sequenceSize)
    {
    	// If a wild
        if (sequence[count] >= encoding_numRegularLetters)
        {
        	// Record location and code of wild
			wildcardEdit = memSingleBlock_newEntry(wildcardEdits);
            wildcardEdit->code = sequence[count];
            wildcardEdit->position = count;

            // Code replacement
            sequence[count] = encoding_randomEncodedLetter(sequence[count]);
        }
        count++;
	}

    return wildcardEdits->numEntries;
}

// Byte-pack a nucleotide sequence and return the packed size
int4 encoding_bytePackSequence(unsigned char* sequence, int4 sequenceSize)
{
	int4 count = 0;

    // For each group of four letters
    while (count < sequenceSize / 4)
    {
		sequence[count] = encoding_bytePack(sequence + count * 4);
        count++;
    }

    // Pack any remaining 1 to 3 letters
    if (sequenceSize % 4)
    {
		sequence[count] = encoding_bytePackRemaining(sequence + count * 4, sequenceSize % 4);
        count++;
    }

    return count;
}

// Given a byte-packed code, prints the 4 letters
void encoding_printLetters(unsigned char code, int4 numLetters)
{
	if (numLetters > 0)
		printf("%c", encoding_getLetter((code >> 6) & 0x3));
	if (numLetters > 1)
		printf("%c", encoding_getLetter((code >> 4) & 0x3));
	if (numLetters > 2)
		printf("%c", encoding_getLetter((code >> 2) & 0x3));
	if (numLetters > 3)
		printf("%c", encoding_getLetter(code & 0x3));
	fflush(stdout);
}

// Given a sequence and its length, encode it
void encoding_encodeSequence(char* sequence, int4 sequenceSize, uint4 alphabetType)
{
	int4 count;
    unsigned char code;

    // Convert letters to codes
    count = 0;
    while (count < sequenceSize)
    {
        code = encoding_codesArray[(unsigned char)(sequence[count])];
        sequence[count] = code;
        count++;
    }
}

// Free structures used to encoding/decode
void encoding_free()
{
	unsigned char code;

	free(encoding_codesArray);
	free(encoding_packedByteLookup);

    if (encoding_alphabetType == encoding_nucleotide)
    	free(encoding_complements);

    code = 0;
    while (code < encoding_numLetters)
    {
        // Declare memory for regular codes
        free(encoding_wildcards[code].replacementCodes);
		code++;
    }
}

#ifdef __cplusplus
}
#endif
