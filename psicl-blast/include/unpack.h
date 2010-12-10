#ifndef _unpack_
#define _unpack_

#ifdef __cplusplus
extern "C" {
#endif

// Initialize region copying/unpacking
void unpack_initialize();

// Return the unpack region in the given list that contains the given subject offset
struct unpackRegion* unpack_selectRegion(struct unpackRegion* unpackRegions, uint4 numUnpackRegions,
                                         uint4 subjectOffset);

// Extend the start of a region if necessary
void unpack_extendRegionStart(int4 position, struct unpackRegion* unpackRegion);

// Extend the end of a region if necessary
void unpack_extendRegionEnd(int4 position, struct unpackRegion* unpackRegion);

// Load a single subject into memory
int4 unpack_loadSubject(struct PSSMatrix PSSMatrix, struct alignment* alignment);

// Unpack entire of sections of a subject sequence before gapped alignment
void unpack_unpackSubject(struct PSSMatrix PSSMatrix, struct alignment* alignment);

// Return 1 if the entire subject has been unpacked, otherwise 0
int unpack_entireSubjectUnpacked(struct alignment* alignment);

// Free memory used to store unpacked regions
void unpack_free();

#ifdef __cplusplus
}
#endif

#endif
