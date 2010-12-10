#ifndef _writedb_
#define _writedb_

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

extern char* writedb_filename;
extern FILE *writedb_sequenceFile, *writedb_descriptionsFile;
extern uint8 writedb_volumnSize;
extern char *writedb_sequenceFilename, *writedb_descriptionsFilename, *writedb_dataFilename;
extern uint4 writedb_maximumSequenceLength, writedb_alphabetType, writedb_minimumSequenceLength;
extern uint8 writedb_numberOfLetters;
extern uint4 writedb_volume, writedb_sequenceCount;

struct edit
{
	uint4 position;
    unsigned char code;
};

struct child
{
	unsigned char* sequence;
	char* description;
	uint4 descriptionLength;
    uint4 descriptionLocation;
    uint4 regionStart;
    uint4 length;
    struct edit* edits;
    uint4 numEdits;
};

struct sequenceData
{
	uint4 sequenceLength;
	uint4 descriptionStart;
	uint4 descriptionLength;
	uint4 encodedLength;
    unsigned char* sequence;
};

struct sequenceDataFP
{
	uint4 sequenceLength;
	uint4 descriptionStart;
	uint4 descriptionLength;
	uint4 encodedLength;
    uint4 offset;
};



// Initialize writing to formatted database
void writedb_initialize(char* filename, uint4 alphabetType);

// Add sequence to the formatted collection
void writedb_addSequence(unsigned char* sequence, uint4 sequenceLength, unsigned char* description,
                         uint4 descriptionLength, unsigned char* wildcards, uint4 wildcardsLength,
                         struct child* children, uint4 numChildren);

// Finalize writing to the formatted collection
void writedb_close();

#ifdef __MY_EXTERN_C__
}
#endif

#endif
