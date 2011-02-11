#ifndef _dust_
#define _dust_

#ifdef __cplusplus
extern "C" {
#endif

struct chunk
{
	int score;
    int start;
    int end;
};

struct maskRegion
{
	struct maskRegion* next;
	int from, to;
};

// Prototypes
void dust_dustSequence(char* originalSequence);
void dust_processChunk(int windowLength, unsigned char* sequence, int chunkStart, struct chunk* chunk);
void dust_processWindow(int windowLength, int windowStart, struct chunk* chunk, unsigned char* sequence);

#ifdef __cplusplus
}
#endif


#endif
