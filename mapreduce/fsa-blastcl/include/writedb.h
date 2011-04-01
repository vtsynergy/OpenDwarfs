#ifndef _writedb_
#define _writedb_

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

extern char* writedb_filename;
extern FILE *writedb_sequenceFile, *writedb_descriptionsFile;
extern uint8 writedb_volumnSize;
extern char *writedb_sequenceFilename, *writedb_descriptionsFilename, *writedb_dataFilename;
extern uint writedb_maximumSequenceLength, writedb_alphabetType, writedb_minimumSequenceLength;
extern uint8 writedb_numberOfLetters;
extern uint writedb_volume, writedb_sequenceCount;

struct edit
{
	uint position;
    unsigned char code;
};

struct child
{
	unsigned char* sequence;
	char* description;
	uint descriptionLength;
    uint descriptionLocation;
    uint regionStart;
    uint length;
    struct edit* edits;
    uint numEdits;
};

struct sequenceData
{
	uint sequenceLength;
	uint descriptionStart;
	uint descriptionLength;
	uint encodedLength;
    unsigned char* sequence;
};

struct __attribute__((aligned )) sequenceDataFP 
{
	cl_uint sequenceLength;
	cl_uint descriptionStart;
	cl_uint descriptionLength;
	cl_uint encodedLength;
    cl_uint offset;
};



// Initialize writing to formatted database
void writedb_initialize(char* filename, uint alphabetType);

// Add sequence to the formatted collection
void writedb_addSequence(unsigned char* sequence, uint sequenceLength, unsigned char* description,
                         uint descriptionLength, unsigned char* wildcards, uint wildcardsLength,
                         struct child* children, uint numChildren);

// Finalize writing to the formatted collection
void writedb_close();

#ifdef __MY_EXTERN_C__
}
#endif

#endif
