#ifndef _readdb_
#define _readdb_

#ifdef __cplusplus
extern "C" {
#endif

extern uint4 readdb_numberOfSequences, readdb_longestSequenceLength, readdb_dbAlphabetType;
extern uint8 readdb_numberOfLetters;
extern unsigned char *readdb_filename, *readdb_sequences;
extern uint4 readdb_fileSize, readdb_sequenceCount, readdb_descriptionStart;
extern uint4 readdb_volumeNumber, readdb_numberOfClusters, readdb_numberOfVolumes;
extern uint4 readdb_numVolumeSequences, readdb_volume;
extern struct sequenceData* readdb_sequenceData;

// Open formatted database for reading
void readdb_open(char* filename);

// Read a sequence and description information. Return 0 if end-of-collection.
int readdb_readSequence(unsigned char** sequence, uint4* sequenceLength, uint4* descriptionStart,
                        uint4* descriptionLength, uint4* encodedLength);

// Load the next volume
int readdb_nextVolume();

// Get the children
struct child* readdb_getChildren(unsigned char* sequence, uint4 sequenceLength, uint4 encodedLength,
                                 uint4 descriptionLocation, uint4* numChildren);

// Close the database for reading
void readdb_close();

#ifdef __cplusplus
}
#endif

#endif
