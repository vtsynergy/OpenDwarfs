#ifndef _gappedExtension_
#define _gappedExtension_

#ifdef __cplusplus
extern "C" {
#endif

// Build a gapped extension with a trace and nominal score from the seed point4 of an ungapped
// extension using dynamic programming
struct gappedExtension* gappedExtension_build(struct ungappedExtension* ungappedExtension,
                        struct PSSMatrix PSSMatrix, int4 subjectSize, unsigned char* subject,
                        struct unpackRegion* unpackRegion, int4 dropoff);

// Debugging routine
void gappedExtension_printBeforeRow(int4* row, unsigned char* subject, unsigned char* rowDropoff,
                                    unsigned char* columnDropoff);
// Debugging routine
void gappedExtension_printAfterRow(int4* row, unsigned char* subject, unsigned char* rowDropoff,
                                   unsigned char* columnDropoff);

// Given a gapped extension with a nominal score, calculate the normalized score
// and E-Value
void gappedExtension_score(struct gappedExtension* gappedExtension);

// Given a gappedExtension and list of ungappedExtensions, prune the latter to
// remove those which overlap/int4ersect with the gappedExtension
void gappedExtension_prune(struct gappedExtension* gappedExtension,
                           struct ungappedExtension* ungappedExtension);

void gappedExtension_free();

#ifdef __cplusplus
}
#endif


#endif

