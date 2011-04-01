#ifndef _statistics_
#define _statistics_

#ifdef __cplusplus
extern "C" {
#endif

struct statisticalParameters
{
	char* matrix;
    int2 startGap;
    int2 extendGap;
    float lambda, K, H;
    float alpha, beta;
};

extern float statistics_log2;
extern float statistics_gappedLogK;

// Search space size metrics
extern int4 statistics_querySize;
extern int4 statistics_effectiveQuerySize;
extern uint8 statistics_databaseSize;
extern uint8 statistics_effectiveDatabaseSize;
extern uint8 statistics_searchSpaceSize;
extern int4 statistics_lengthAdjust;

// Ungapped statistical parameters
extern float statistics_ungappedLambda;
extern float statistics_ungappedH;
extern float statistics_ungappedK;
extern float statistics_ungappedLogK;

// Gapped statistical parameters
extern struct statisticalParameters statistics_gappedParams;

// Dropoff values
extern int4 statistics_ungappedNominalDropoff;
extern int4 statistics_gappedNominalDropoff;
extern int4 statistics_gappedFinalNominalDropoff;

// Initialize statistics by calculating some global parameters
void statistics_initialize(struct PSSMatrix PSSMatrix, uint8 databaseSize, int4 numberOfSequences);

// Convert a normalized score to a nominal score
int4 statistics_ungappedNormalized2nominal(float normalizedScore);

// Convert a gapped nominal score to a normalized score
float statistics_gappedNominal2normalized(int4 nominalScore);

// Calculate the evalue for a given gapped normalizedScore
double statistics_gappedCalculateEvalue(float normalizedScore);

// Given an evalue (such as a cutoff) calculate the minimum gapped nominal score needed to attain it
int4 statistics_gappedEvalue2nominal(double evalue);

// Calculate minimum nominal score required to trigger gapping for nucleotide searches
int4 statistics_ungappedNucleotideTrigger(struct PSSMatrix PSSMatrix);

#ifdef __cplusplus
}
#endif

#endif
