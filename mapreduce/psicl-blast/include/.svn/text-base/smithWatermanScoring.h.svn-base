#ifndef _smithWatermanScoring_
#define _smithWatermanScoring_

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

// Perform smith-waterman dynamic programming on the given query and subject sequences
struct dpResults smithWatermanScoring_score(struct PSSMatrix PSSMatrix, int4 subjectSize,
                                            unsigned char* subject);

// Reverse the query and subject sequences, and perform dynamic programming to find the
// start of the optimal alignment, rather than the end
struct dpResults smithWatermanScoring_scoreReverse(struct PSSMatrix PSSMatrix, int4 subjectSize,
                                                   unsigned char* subject, struct coordinate end);
#ifdef __MY_EXTERN_C__
}
#endif

#endif

