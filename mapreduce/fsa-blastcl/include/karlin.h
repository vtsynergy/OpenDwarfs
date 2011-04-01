#ifndef _karlin_
#define _karlin_

#ifdef __cplusplus
extern "C" {
#endif

extern float BlastKarlin_lambda;
extern float BlastKarlin_K;
extern float BlastKarlin_H;

void BlastKarlinBlkCalc(double* scoreProbabilities, int4 min, int4 max);

int4 BlastComputeLengthAdjustment(float K,
                             float logK,
                             float alpha_d_lambda,
                             float beta,
                             int4 query_length,
                             uint4 db_length,
                             int4 db_num_seqs,
                             int4 *length_adjustment);
#ifdef __cplusplus
}
#endif

#endif

