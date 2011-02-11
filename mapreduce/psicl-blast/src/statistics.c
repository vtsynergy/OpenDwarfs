// statistics.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Statistical functions related to BLAST (ie. converting nominal scores
// to normalized scores and to e-values)

#include <math.h>
#include "blast.h"

float statistics_log2;
float statistics_gappedLogK;

// Search space size metrics
int4 statistics_querySize;
int4 statistics_effectiveQuerySize;
uint8 statistics_databaseSize = 0;
uint8 statistics_effectiveDatabaseSize;
uint8 statistics_searchSpaceSize;
int4 statistics_lengthAdjust;
int4 statistics_numberOfSequences;

// Ungapped statistical parameters
float statistics_ungappedLambda;
float statistics_ungappedH;
float statistics_ungappedK;
float statistics_ungappedLogK;

// Gapped statistical parameters
struct statisticalParameters statistics_gappedParams;

// Dropoff values
int4 statistics_ungappedNominalDropoff;
int4 statistics_gappedNominalDropoff;
int4 statistics_gappedFinalNominalDropoff;

void statistics_calculateUngappedKarlinParameters(struct PSSMatrix PSSMatrix);
float statistics_calcLengthAdjust(int4 numberOfSequences);
float statistics_calcLengthAdjustNew(int4 numberOfSequences);
struct statisticalParameters statistics_lookupPrecomputedParams(
       char* matrix, int2 startGap, int2 extendGap);

// Initialize statistics by calculating some global parameters
void statistics_initialize(struct PSSMatrix PSSMatrix, uint8 databaseSize, int4 numberOfSequences)
{
	// Record the size (number of letters) of query and subject database
	if (PSSMatrix.strandLength == 0)
        statistics_querySize = PSSMatrix.length;
	else
        statistics_querySize = PSSMatrix.strandLength;

    // Use command-line parameter for database size if given
    if (parameters_databaseSize == 0)
		statistics_databaseSize = databaseSize;
    else
    	statistics_databaseSize = parameters_databaseSize;

    statistics_numberOfSequences = numberOfSequences;

    // Use command-line parameter for number of sequences if given
    if (parameters_numberOfSequences != 0)
    	numberOfSequences = parameters_numberOfSequences;

    // Calculate ungapped karlin-atschul parameters lambda, K and H for given
	// query matrix, and average residue compositions for the subject database
	statistics_calculateUngappedKarlinParameters(PSSMatrix);

    if (encoding_alphabetType == encoding_protein)
    {
        // Lookup precomputed gapped lambda, K, H, alpha and beta values based on
        // scoring matrix and gap penalties
        statistics_gappedParams = statistics_lookupPrecomputedParams(
            parameters_scoringMatrix, parameters_startGap, parameters_extendGap);
	}
    else
    {
    	// For nucleotide, use ungapped parameters for gapped alignment
    	statistics_gappedParams.lambda = statistics_ungappedLambda;
        statistics_gappedParams.K = statistics_ungappedK;
        statistics_gappedParams.H = statistics_ungappedH;
        statistics_gappedParams.alpha = 0;
        statistics_gappedParams.beta = 0;
    }

    // Log of K (for gapped and ungapped extension) and log of 2
	statistics_ungappedLogK = log(statistics_ungappedK);
	statistics_gappedLogK = log(statistics_gappedParams.K);
	statistics_log2 = log(2.0);

	// Calculate nominal dropoff for ungapped extension using normalized dropoff
	statistics_ungappedNominalDropoff = ceil(parameters_ungappedNormalizedDropoff *
	                                    statistics_log2 / statistics_ungappedLambda);

	// Same for initial and final gapped extension dropoffs
	statistics_gappedNominalDropoff = floor(parameters_gappedNormalizedDropoff *
	                                  statistics_log2 / statistics_gappedParams.lambda);
	statistics_gappedFinalNominalDropoff = floor(parameters_gappedFinalNormalizedDropoff *
	                                       statistics_log2 / statistics_gappedParams.lambda);

	// Calculate length adjust for query and subject sequences
    if (encoding_alphabetType == encoding_protein)
    {
        // Round lengthAdjust value down
        statistics_lengthAdjust = floor(statistics_calcLengthAdjustNew(numberOfSequences));
	}
    else
    {
	    statistics_lengthAdjust = floor(statistics_calcLengthAdjust(numberOfSequences));
    }

	// Using length adjustment, calculate effective query, database length and search space
    statistics_effectiveQuerySize
		= statistics_querySize - statistics_lengthAdjust;
	statistics_effectiveDatabaseSize
		= statistics_databaseSize - numberOfSequences * statistics_lengthAdjust;
	statistics_searchSpaceSize
		= statistics_effectiveQuerySize * statistics_effectiveDatabaseSize;
}

#define statistics_maximum(a,b) ((a > b) ? a : b)

// Calculate length adjustment (aka. "effective HSP length") by iteratively
// applying equation from "Local Alignment Statistics" Altschul,Gish
float statistics_calcLengthAdjustNew(int4 numberOfSequences)
{
	int4 lengthAdjust;

    BlastComputeLengthAdjustment(statistics_gappedParams.K,
                             statistics_gappedLogK,
                             statistics_gappedParams.alpha / statistics_gappedParams.lambda,
                             statistics_gappedParams.beta,
                             statistics_querySize,
                             statistics_databaseSize,
                             numberOfSequences,
                             &lengthAdjust);

	return lengthAdjust;
}

// Convert a normalized score to a nominal score
int4 statistics_ungappedNormalized2nominal(float normalizedScore)
{
	float nominalScore;

	nominalScore = (normalizedScore * statistics_log2 + statistics_ungappedLogK) / statistics_ungappedLambda;

	return floor(nominalScore);
}

// Convert an ungapped nominal score to a normalized score
float statistics_ungappedNominal2normalized(int4 nominalScore)
{
	float normalizedScore;

	normalizedScore = (((float)nominalScore * statistics_ungappedLambda) - statistics_ungappedLogK)
	                  / statistics_log2;

	return normalizedScore;
}

// Calculate the evalue for a given ungapped normalizedScore
double statistics_ungappedCalculateEvalue(float normalizedScore)
{
	return statistics_searchSpaceSize / pow(2, normalizedScore);
}

// Given an evalue (such as a cutoff) calculate the minimum ungapped nominal score needed to attain it
int4 statistics_ungappedEvalue2nominal(double evalue)
{
	double normalizedScore;

	normalizedScore = log(statistics_searchSpaceSize / evalue) / (double)statistics_log2;
	return ceil((statistics_log2 * normalizedScore + statistics_ungappedLogK)
	       / statistics_ungappedLambda);
}

// Convert a gapped nominal score to a normalized score
float statistics_gappedNominal2normalized(int4 nominalScore)
{
	float normalizedScore;

	normalizedScore = (((float)nominalScore * statistics_gappedParams.lambda)
                    - statistics_gappedLogK) / statistics_log2;

	return normalizedScore;
}

// Given a normalized score calculate the lowest gapped nominal score needed to achieve it
int4 statistics_gappedNormalized2nominal(double normalizedScore)
{
	return ceil(((statistics_log2 * normalizedScore + statistics_gappedLogK)
	       / statistics_gappedParams.lambda) - 0.5);
}

// Calculate the evalue for a given gapped normalizedScore
double statistics_gappedCalculateEvalue(float normalizedScore)
{
	return statistics_searchSpaceSize / pow(2, normalizedScore);
}

// Given an evalue (such as a cutoff) calculate the minimum gapped nominal score needed to attain it
int4 statistics_gappedEvalue2nominal(double evalue)
{
	double normalizedScore;

	normalizedScore = log(statistics_searchSpaceSize / evalue) / (double)statistics_log2;

//    printf("[%f,%f,%f,%f, %f]\n", normalizedScore, statistics_log2 * normalizedScore, statistics_log2 * normalizedScore + statistics_gappedLogK,(statistics_log2 * normalizedScore + statistics_gappedLogK)
//	       / (double)statistics_gappedParams.lambda, statistics_gappedParams.lambda);

	return ceil((statistics_log2 * normalizedScore + statistics_gappedLogK)
	       / statistics_gappedParams.lambda);
}

// Calculate minimum nominal score required to trigger gapping for nucleotide searches
int4 statistics_ungappedNucleotideTrigger(struct PSSMatrix PSSMatrix)
{
	double evalue, normalizedScore;
	int4 averageSubjectLength;

    evalue = 0.025;
    averageSubjectLength = statistics_databaseSize / statistics_numberOfSequences;

	if (averageSubjectLength > PSSMatrix.length)
	{
		normalizedScore = log(PSSMatrix.length * averageSubjectLength / evalue)
		                / (double)statistics_log2;
	}
	else
	{
		normalizedScore = log(averageSubjectLength * averageSubjectLength / evalue)
		                / (double)statistics_log2;
	}

	return ceil((statistics_log2 * normalizedScore + statistics_ungappedLogK)
	       / statistics_ungappedLambda);
}

/***************************************************************************************/
/* For each matrix, gapStart, gapExtend - precomputed gapped lambda, K, H, alpha, beta */
/***************************************************************************************/

#define NUM_PRECOMPUTED_PARAMS 60
static const struct statisticalParameters statistics_precomputedParams[NUM_PRECOMPUTED_PARAMS] =
{
    {"BLOSUM45", 13, 3, 0.207, 0.049, 0.14, 1.5, -22},
    {"BLOSUM45", 12, 3, 0.199, 0.039, 0.11, 1.8, -34},
    {"BLOSUM45", 11, 3, 0.190, 0.031, 0.095, 2.0, -38},
    {"BLOSUM45", 10, 3, 0.179, 0.023, 0.075, 2.4, -51},
    {"BLOSUM45", 16, 2, 0.210, 0.051, 0.14, 1.5, -24},
    {"BLOSUM45", 15, 2, 0.203, 0.041, 0.12, 1.7, -31},
    {"BLOSUM45", 14, 2, 0.195, 0.032, 0.10, 1.9, -36},  /* Recommended */
    {"BLOSUM45", 13, 2, 0.185, 0.024, 0.084, 2.2, -45},
    {"BLOSUM45", 12, 2, 0.171, 0.016, 0.061, 2.8, -65},
    {"BLOSUM45", 19, 1, 0.205, 0.040, 0.11, 1.9, -43},
    {"BLOSUM45", 18, 1, 0.198, 0.032, 0.10, 2.0, -43},
    {"BLOSUM45", 17, 1, 0.189, 0.024, 0.079, 2.4, -57},
    {"BLOSUM45", 16, 1, 0.176, 0.016, 0.063, 2.8, -67},
    {"BLOSUM50", 13, 3, 0.212, 0.063, 0.19, 1.1, -16},
    {"BLOSUM50", 12, 3, 0.206, 0.055, 0.17, 1.2, -18},
    {"BLOSUM50", 11, 3, 0.197, 0.042, 0.14, 1.4, -25},
    {"BLOSUM50", 10, 3, 0.186, 0.031, 0.11, 1.7, -34},
    {"BLOSUM50", 9, 3, 0.172, 0.022, 0.082, 2.1, -48},
    {"BLOSUM50", 16, 2, 0.215, 0.066, 0.20, 1.05, -15},
    {"BLOSUM50", 15, 2, 0.210, 0.058, 0.17, 1.2, -20},
    {"BLOSUM50", 14, 2, 0.202, 0.045, 0.14, 1.4, -27},
    {"BLOSUM50", 13, 2, 0.193, 0.035, 0.12, 1.6, -32}, /* Recommended */
    {"BLOSUM50", 12, 2, 0.181, 0.025, 0.095, 1.9, -41},
    {"BLOSUM50", 19, 1, 0.212, 0.057, 0.18, 1.2, -21},
    {"BLOSUM50", 18, 1, 0.207, 0.050, 0.15, 1.4, -28},
    {"BLOSUM50", 17, 1, 0.198, 0.037, 0.12, 1.6, -33},
    {"BLOSUM50", 16, 1, 0.186, 0.025, 0.10, 1.9, -42},
    {"BLOSUM50", 15, 1, 0.171, 0.015, 0.063, 2.7, -76},
    {"BLOSUM62", 11, 2, 0.297, 0.082, 0.27, 1.1, -10},
    {"BLOSUM62", 10, 2, 0.291, 0.075, 0.23, 1.3, -15},
    {"BLOSUM62", 9, 2, 0.279, 0.058, 0.19, 1.5, -19},
    {"BLOSUM62", 8, 2, 0.264, 0.045, 0.15, 1.8, -26},
    {"BLOSUM62", 7, 2, 0.239, 0.027, 0.10, 2.5, -46},
    {"BLOSUM62", 6, 2, 0.201, 0.012, 0.061, 3.3, -58},
    {"BLOSUM62", 13, 1, 0.292, 0.071, 0.23, 1.2, -11},
    {"BLOSUM62", 12, 1, 0.283, 0.059, 0.19, 1.5, -19},
    {"BLOSUM62", 11, 1, 0.267, 0.041, 0.14, 1.9, -30},  /* Recommended */
    {"BLOSUM62", 10, 1, 0.243, 0.024, 0.10, 2.5, -44},
    {"BLOSUM62", 9, 1, 0.206, 0.010, 0.052, 4.0, -87},
    {"BLOSUM80", 25, 2, 0.342, 0.17, 0.66, 0.52, -1.6},
    {"BLOSUM80", 13, 2, 0.336, 0.15, 0.57, 0.59, -3},
    {"BLOSUM80", 9, 2, 0.319, 0.11, 0.42, 0.76, -6},
    {"BLOSUM80", 8, 2, 0.308, 0.090, 0.35, 0.89, -9},
    {"BLOSUM80", 7, 2, 0.293, 0.070, 0.27, 1.1, -14},
    {"BLOSUM80", 6, 2, 0.268, 0.045, 0.19, 1.4, -19},
    {"BLOSUM80", 11, 1, 0.314, 0.095, 0.35, 0.90, -9},
    {"BLOSUM80", 10, 1, 0.299, 0.071, 0.27, 1.1, -14},  /* Recommended */
    {"BLOSUM80", 9, 1, 0.279, 0.048, 0.20, 1.4, -19},
    {"PAM30", 7, 2, 0.305, 0.15, 0.87, 0.35, -3},
    {"PAM30", 6, 2, 0.287, 0.11, 0.68, 0.42, -4},
    {"PAM30", 5, 2, 0.264, 0.079, 0.45, 0.59, -7},
    {"PAM30", 10, 1, 0.309, 0.15, 0.88, 0.35, -3},  /* Recommended */
    {"PAM30", 9, 1, 0.294, 0.11, 0.61, 0.48, -6},
    {"PAM30", 8, 1, 0.270, 0.072, 0.40, 0.68, -10},
    {"PAM70", 8, 2, 0.301, 0.12, 0.54, 0.56, -5},
    {"PAM70", 7, 2, 0.286, 0.093, 0.43, 0.67, -7},
    {"PAM70", 6, 2, 0.264, 0.064, 0.29, 0.90, -12},
    {"PAM70", 11, 1, 0.305, 0.12, 0.52, 0.59, -6},
    {"PAM70", 10, 1, 0.291, 0.091, 0.41, 0.71, -9},
    {"PAM70", 9, 1, 0.270, 0.060, 0.28, 0.97, -14},  /* Recommended */
};

// Given a scoring matrix and gap penalties, looks up and returns the gapped
// precomputed lambda, k, h, alpha, and beta values
struct statisticalParameters statistics_lookupPrecomputedParams(
       char* matrix, int2 startGap, int2 extendGap)
{
	int4 count = 0;
	struct statisticalParameters parameters;

    // Search through each entry in the table
	while (count < NUM_PRECOMPUTED_PARAMS)
    {
    	parameters = statistics_precomputedParams[count];

        // If entry matches, return the details
		if (strcmp(parameters.matrix, matrix) == 0 && parameters.startGap == startGap &&
            parameters.extendGap == extendGap)
        {
        	return statistics_precomputedParams[count];
        }
    	count++;
    }

    // No matching entries found, error
    fprintf(stderr, "Error: Matrix %s with gap penalties %d,%d not supported\n",
            matrix, startGap, extendGap);
    exit(-1);
}

/*************************************************/
/* Code for calculating ungapped lambda, K and H */
/*************************************************/

void statistics_calculateUngappedKarlinParameters(struct PSSMatrix PSSMatrix)
{
	int4 highest, lowest;
	double* scoreProbabilities;
	int4 queryPosition;
	unsigned char subjectCode;
	double probability, sum = 0;
	int4 score;
    int4 numRegularLettersInQuery = 0;
	static float* alphabetFrequencies;

    if (encoding_alphabetType == encoding_nucleotide)
    {
        alphabetFrequencies = Nucleotide_prob;
    }
    else
    {
    	// Use Robinson&Robinson frequencies for protein alphabet
    	alphabetFrequencies = Robinson_prob;
    }

    highest = PSSMatrix.highestValue;
	lowest = PSSMatrix.lowestValue;

	// Initialize array to hold probability values for range of possible scores
	scoreProbabilities = (double*)global_malloc(sizeof(double) * (highest - lowest + 1));
	scoreProbabilities -= lowest;
	score = lowest;
	while (score <= highest)
	{
		scoreProbabilities[score] = 0.0;
		score++;
	}

    // Determine the number of regular letters (ie. non-ambigious codes) in the query
	queryPosition = 0;
	while (queryPosition < statistics_querySize)
	{
//    	if (PSSMatrix.queryCodes[queryPosition] < encoding_numRegularLetters)
    	if (PSSMatrix.queryCodes[queryPosition] != encoding_numRegularLetters)
        {
			numRegularLettersInQuery++;
        }
        queryPosition++;
	}

	// For each position in the query PSSMatrix that does not represent an ambigious code
	queryPosition = 0;
	while (queryPosition < statistics_querySize)
	{
//    	if (PSSMatrix.queryCodes[queryPosition] < encoding_numRegularLetters)
    	if (PSSMatrix.queryCodes[queryPosition] != encoding_numRegularLetters)
        {
            // For each regular amino-acid subject code
            subjectCode = 0;
            while (subjectCode < encoding_numRegularLetters)
            {
                // Calculate probability that, if a residue was randomly chosen from the subject database
                // and a position was randomly chosen from query, they would be subjectCode and queryPosition
                // respectively
                probability = alphabetFrequencies[subjectCode] / (float)numRegularLettersInQuery;

                // Calculate score of aligning query position to this subject residue
                score = PSSMatrix.matrix[queryPosition][subjectCode];

                // Add to set of probabilities
                scoreProbabilities[score] += probability;

                subjectCode++;
            }
		}
        queryPosition++;
	}

	// Calculate residue frequency normalizing value
	subjectCode = 0;
	while (subjectCode < encoding_numRegularLetters)
	{
		sum += alphabetFrequencies[subjectCode];
		subjectCode++;
	}

	// Normalized probabilities by dividing them by the sum of the robinson frequency values
	score = lowest;
	while (score <= highest)
	{
		scoreProbabilities[score] /= sum;
		score++;
	}

	// Calculate the Lambda, H and K values for these score probabilities
	// See karlin.c for implementation
	BlastKarlinBlkCalc(scoreProbabilities, lowest, highest);

	statistics_ungappedH = BlastKarlin_H;
	statistics_ungappedK = BlastKarlin_K;
	statistics_ungappedLambda = BlastKarlin_lambda;

    if (statistics_ungappedK <= 0.0001)
    {
    	fprintf(stderr, "Error: unable to calculate statistical significance for given scoring scheme\n");
        fflush(stderr);
        exit(-1);
    }

	scoreProbabilities += lowest;
    free(scoreProbabilities);
}

// Calculate length adjustment (aka. "effective HSP length") by iteratively
// applying equation from "Local Alignment Statistics" Altschul,Gish
// using ungapped K & H values
float statistics_calcLengthAdjust(int4 numberOfSequences)
{
	int4 count;
	float lengthAdjust;
	float minimumQueryLength;

	// Query can not be shorter than 1/K
	minimumQueryLength = (float)1 / statistics_ungappedK;

	count = 0;
	lengthAdjust = 0;
	// Calculate length adjustment value
	while (count < 5)
	{
		lengthAdjust = (statistics_ungappedLogK +
			log((statistics_querySize - lengthAdjust) *
			    (statistics_databaseSize - numberOfSequences *
				 lengthAdjust))) / statistics_ungappedH;

		// Length adjust can not shorten query to less than minimum query length
		if (lengthAdjust > statistics_querySize - minimumQueryLength)
		{
			lengthAdjust = statistics_querySize - minimumQueryLength;
		}

		count++;
	}

	return lengthAdjust;
}
