#ifndef _oldGappedScoring_
#define _oldGappedScoring_

#ifdef __cplusplus
extern "C" {
#endif

int4 oldGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff);

void oldGappedScoring_free();

#ifdef __cplusplus
}
#endif

#endif

                        
