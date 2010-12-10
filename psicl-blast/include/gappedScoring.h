#ifndef _gappedScoring_
#define _gappedScoring_

#ifdef __cplusplus
extern "C" {
#endif

// Perform gapped alinment with restricted insertion
int4 gappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff);

void gappedScoring_free();

#ifdef __cplusplus
}
#endif

#endif

