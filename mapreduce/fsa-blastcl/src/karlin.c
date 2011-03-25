// karlin.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// This code has been copied from blastkar.c and modified to make it work standalone
// without using implementations of mathematical functions in ncbimath.c and other
// definitions in the NCBI BLAST toolkit.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "blast.h"
#include "karlin.h"

/**************** Statistical Significance Parameter Subroutine ****************

    Version 1.0     February 2, 1990
    Version 1.2     July 6,     1990

    Program by:     Stephen Altschul

    Address:        National Center for Biotechnology Information
                    National Library of Medicine
                    National Institutes of Health
                    Bethesda, MD  20894

    Internet:       altschul@ncbi.nlm.nih.gov

See:  Karlin, S. & Altschul, S.F. "Methods for Assessing the Statistical
    Significance of Molecular Sequence Features by Using General Scoring
    Schemes,"  Proc. Natl. Acad. Sci. USA 87 (1990), 2264-2268.

    Computes the parameters lambda and K for use in calculating the
    statistical significance of high-scoring segments or subalignments.

    The scoring scheme must be integer valued.  A positive score must be
    possible, but the expected (mean) score must be negative.

    A program that calls this routine must provide the value of the lowest
    possible score, the value of the greatest possible score, and a pointer
    to an array of probabilities for the occurence of all scores between
    these two extreme scores.  For example, if score -2 occurs with
    probability 0.7, score 0 occurs with probability 0.1, and score 3
    occurs with probability 0.2, then the subroutine must be called with
    low = -2, high = 3, and pr pointing to the array of values
    { 0.7, 0.0, 0.1, 0.0, 0.0, 0.2 }.  The calling program must also provide
    pointers to lambda and K; the subroutine will then calculate the values
    of these two parameters.  In this example, lambda=0.330 and K=0.154.

    The parameters lambda and K can be used as follows.  Suppose we are
    given a length N random sequence of independent letters.  Associated
    with each letter is a score, and the probabilities of the letters
    determine the probability for each score.  Let S be the aggregate score
    of the highest scoring contiguous segment of this sequence.  Then if N
    is sufficiently large (greater than 100), the following bound on the
    probability that S is greater than or equal to x applies:

            P( S >= x )   <=   1 - exp [ - KN exp ( - lambda * x ) ].

    In other words, the p-value for this segment can be written as
    1-exp[-KN*exp(-lambda*S)].

    This formula can be applied to pairwise sequence comparison by assigning
    scores to pairs of letters (e.g. amino acids), and by replacing N in the
    formula with N*M, where N and M are the lengths of the two sequences
    being compared.

    In addition, letting y = KN*exp(-lambda*S), the p-value for finding m
    distinct segments all with score >= S is given by:

                           2             m-1           -y
            1 - [ 1 + y + y /2! + ... + y   /(m-1)! ] e

    Notice that for m=1 this formula reduces to 1-exp(-y), which is the same
    as the previous formula.

*******************************************************************************/
float BlastKarlin_lambda;
float BlastKarlin_K;
float BlastKarlin_H;

double BlastKarlinLHtoK(double* sfp, int4 min, int4 max, double lambda, double H);
double BlastKarlinLambdaNR(double* sfp, int4 min, int4 max);
double BlastKarlinLtoH(double* sfp, int4 min, int4 max, double lambda);

void BlastKarlinBlkCalc(double* scoreProbabilities, int4 min, int4 max)
{
	/* Calculate the parameter Lambda */
	BlastKarlin_lambda = BlastKarlinLambdaNR(scoreProbabilities, min, max);

	/* Calculate H */
	BlastKarlin_H = BlastKarlinLtoH(scoreProbabilities, min, max, BlastKarlin_lambda);

	/* Calculate K */
	BlastKarlin_K = BlastKarlinLHtoK(scoreProbabilities, min, max, BlastKarlin_lambda,
	                                 BlastKarlin_H);
}

#define DIMOFP0	(iter*range + 1)
#define DIMOFP0_MAX (BLAST_KARLIN_K_ITER_MAX*BLAST_SCORE_RANGE_MAX+1)
/****************************************************************************
For more accuracy in the calculation of K, set K_SUMLIMIT to 0.00001.
For high speed in the calculation of K, use a K_SUMLIMIT of 0.001
Note:  statistical significance is often not greatly affected by the value
of K, so high accuracy is generally unwarranted.
*****************************************************************************/
/* K_SUMLIMIT_DEFAULT == sumlimit used in BlastKarlinLHtoK() */
#define BLAST_KARLIN_K_SUMLIMIT_DEFAULT 0.01
#define BLAST_KARLIN_K_ITER_MAX 100
/*
SCORE_MIN is (-2**31 + 1)/2 because it has been observed more than once that
a compiler did not properly calculate the quantity (-2**31)/2.  The list
of compilers which failed this operation have at least at some time included:
NeXT and a version of AIX/370's MetaWare High C R2.1r.
For this reason, SCORE_MIN is not simply defined to be LONG_MIN/2.
*/
#define BLAST_SCORE_1MIN (-10000)
#define BLAST_SCORE_1MAX ( 1000)
#define BLAST_SCORE_RANGE_MAX	(BLAST_SCORE_1MAX - BLAST_SCORE_1MIN)
/* Initial guess for the value of Lambda in BlastKarlinLambdaNR */
#define BLAST_KARLIN_LAMBDA0_DEFAULT    0.5
/* LAMBDA_ACCURACY_DEFAULT == accuracy to which Lambda should be calc'd */
#define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT    (1.e-5)
/* LAMBDA_ITER_DEFAULT == no. of iterations in LambdaBis = ln(accuracy)/ln(2)*/
#define BLAST_KARLIN_LAMBDA_ITER_DEFAULT        17
#define ABS(a)	((a)>=0?(a):-(a))

long gcd(long a, long b);
int4 BlastScoreChk(int4 lo, int4 hi);

double BlastKarlinLHtoK(double* sfp, int4 min, int4 max, double lambda, double H)
{
#ifndef BLAST_KARLIN_STACKP
	double *P0 = NULL;
#else
	double P0 [DIMOFP0_MAX];
#endif
	int4	low;	/* Lowest score (must be negative) */
	int4	high;	/* Highest score (must be positive) */
	double	K;			/* local copy of K */
	double	ratio;
	int4		i, j;
	int4	range, lo, hi, first, last;
	double	sum;
	double	Sum, av, oldsum, oldsum2;
	int4		iter;
	double	sumlimit;
	double	*p, *ptrP, *ptr1, *ptr2, *ptr1e;
	double	etolami, etolam;

	if (lambda <= 0. || H <= 0.) {
		return -1.;
	}

	low = min;
	high = max;
	range = high - low;

	av = H/lambda;
	etolam = exp((double)lambda);
	if (low == -1 || high == 1) {
		K = av;
		return K * (1.0 - 1./etolam);
	}

	sumlimit = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;

	iter = BLAST_KARLIN_K_ITER_MAX;

	if (DIMOFP0 > DIMOFP0_MAX) {
		return -1.;
	}
#ifndef BLAST_KARLIN_STACKP
	P0 = (double *)global_malloc(DIMOFP0 * sizeof(*P0));
	if (P0 == NULL)
		return -1.;
#else
	memset((CharPtr)P0, 0, DIMOFP0*sizeof(P0[0]));
#endif

	Sum = 0.;
	lo = hi = 0;
	p = &sfp[low];
	P0[0] = sum = oldsum = oldsum2 = 1.;
    for (j = 0; j < iter && sum > sumlimit; Sum += sum /= ++j) {
        first = last = range;
		lo += low;
		hi += high;
        for (ptrP = P0+(hi-lo); ptrP >= P0; *ptrP-- =sum) {
            ptr1 = ptrP - first;
            ptr1e = ptrP - last;
            ptr2 = p + first;
            for (sum = 0.; ptr1 >= ptr1e; )
                sum += *ptr1--  *  *ptr2++;
            if (first)
                --first;
            if (ptrP - P0 <= range)
                --last;
        }
		etolami = pow((double)etolam, lo - 1);
        for (sum = 0., i = lo; i != 0; ++i) {
			etolami *= etolam;
            sum += *++ptrP * etolami;
		}
        for (; i <= hi; ++i)
            sum += *++ptrP;
		oldsum2 = oldsum;
		oldsum = sum;
    }

	/* Terms of geometric progression added for correction */
	ratio = oldsum / oldsum2;
	if (ratio >= (1.0 - sumlimit*0.001)) {
		K = -1.;
		goto CleanUp;
	}
	sumlimit *= 0.01;
	while (sum > sumlimit) {
		oldsum *= ratio;
		Sum += sum = oldsum / ++j;
	}

/* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
Karlin&Altschul (1990) */
    	for (i = 1, j = -low; i <= range && j > 1; ++i)
        	if (p[i])
            		j = gcd(j, i);

	if (j*etolam > 0.05)
	{
		etolami = pow((double)etolam, -j);
    		K = j*exp((double)-2.0*Sum) / (av*(1.0 - etolami));
	}
	else
	    K = -j*exp((double)-2.0*Sum) / (av*(exp(-j*(double)lambda) - 1));

CleanUp:
#ifndef BLAST_KARLIN_K_STACKP
	if (P0 != NULL)
		free(P0);
#endif
	return K;
}

/*
	BlastKarlinLambdaBis

	Calculate Lambda using the bisection method (slow).
*/
double BlastKarlinLambdaBis(double* sfp, int4 min, int4 max)
{
	double* sprob;
	double lambda, up, newval;
	int4	i, low, high;
	int4		j;
	double sum, x0, x1;

	low = min;
	high = max;
	if (BlastScoreChk(low, high) != 0)
		return -1.;

	sprob = sfp;
	up = BLAST_KARLIN_LAMBDA0_DEFAULT;
	for (lambda=0.; ; ) {
		up *= 2;
		x0 = exp((double)up);
		x1 = pow((double)x0, low - 1);
		if (x1 > 0.) {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * (x1 *= x0);
		}
		else {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * exp(up * i);
		}
		if (sum >= 1.0)
			break;
		lambda = up;
	}

	for (j=0; j<BLAST_KARLIN_LAMBDA_ITER_DEFAULT; ++j) {
		newval = (lambda + up) / 2.;
		x0 = exp((double)newval);
		x1 = pow((double)x0, low - 1);
		if (x1 > 0.) {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * (x1 *= x0);
		}
		else {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * exp(newval * i);
		}
		if (sum > 1.0)
			up = newval;
		else
			lambda = newval;
	}
	return (lambda + up) / 2.;
}

/******************* Fast Lambda Calculation Subroutine ************************
	Version 1.0	May 16, 1991
	Program by:	Stephen Altschul

	Uses Newton-Raphson method (fast) to solve for Lambda, given an initial
	guess (lambda0) obtained perhaps by the bisection method.
*******************************************************************************/
double BlastKarlinLambdaNR(double* sfp, int4 min, int4 max)
{
	int4	low;			/* Lowest score (must be negative)  */
	int4	high;			/* Highest score (must be positive) */
	int4	j;
	int4	i;
	double *sprob;
	double lambda0, sum, slope, temp, x0, x1, amt;

	low = min;
	high = max;
	if (BlastScoreChk(low, high) != 0)
		return -1.;

	lambda0 = BLAST_KARLIN_LAMBDA0_DEFAULT;

	/* Calculate lambda */
	sprob = sfp;
	for (j=0; j<20; ++j) { /* limit of 20 should never be close-approached */
		sum = -1.0;
		slope = 0.0;
		if (lambda0 < 0.01)
			break;
		x0 = exp((double)lambda0);
		x1 = pow((double)x0, low - 1);
		if (x1 == 0.)
			break;
		for (i=low; i<=high; ++i) {
			sum += (temp = sprob[i] * (x1 *= x0));
			slope += temp * i;
		}
		lambda0 -= (amt = sum/slope);
		if (ABS(amt/lambda0) < BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT) {
			/*
			Does it appear that we may be on the verge of converging
			to the ever-present, zero-valued solution?
			*/
			if (lambda0 > BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT)
				return lambda0;
			break;
		}
	}
	return BlastKarlinLambdaBis(sfp, min, max);
}


/*
	BlastKarlinLtoH

	Calculate H, the relative entropy of the p's and q's
*/
double BlastKarlinLtoH(double* sfp, int4 min, int4 max, double lambda)
{
	int4 score;
	double av, etolam, etolami;

	if (lambda < 0.) {
		return -1.;
	}
	if (BlastScoreChk(min, max) != 0)
		return -1.;

	etolam = exp((double)lambda);
	etolami = pow((double)etolam, min - 1);
	if (etolami > 0.)
	{
	    av = 0.0;
	    for (score=min; score<=max; score++)
   			av += sfp[score] * score * (etolami *= etolam);
	}
	else
	{
	    av = 0.0;
	    for (score=min; score<=max; score++)
   			av += sfp[score] * score * exp(lambda * score);
	}

   	return lambda * av;
}

long gcd(long a, long b)
{
	long c;

	b = ABS(b);
	if (b > a)
		c=a, a=b, b=c;

	while (b != 0) {
		c = a%b;
		a = b;
		b = c;
	}
	return a;
}

int4 BlastScoreChk(int4 lo, int4 hi)
{
	if (lo >= 0 || hi <= 0 ||
			lo < BLAST_SCORE_1MIN || hi > BLAST_SCORE_1MAX)
		return 1;

	if (hi - lo > BLAST_SCORE_RANGE_MAX)
		return 1;

	return 0;
}

/**
 * Computes the adjustment to the lengths of the query and database sequences
 * that is used to compensate for edge effects when computing evalues.
 *
 * The length adjustment is an integer-valued approximation to the fixed
 * point of the function
 *
 *    f(ell) = beta +
 *               (alpha/lambda) * (log K + log((m - ell)*(n - N ell)))
 *
 * where m is the query length n is the length of the database and N is the
 * number of sequences in the database. The values beta, alpha, lambda and
 * K are statistical, Karlin-Altschul parameters.
 * 
 * The value of the length adjustment computed by this routine, A, 
 * will always be an integer smaller than the fixed point of
 * f(ell). Usually, it will be the largest such integer.  However, the
 * computed length adjustment, A, will also be so small that
 *
 *    K * (m - A) * (n - N * A) > min(m,n).
 *
 * Moreover, an iterative method is used to compute A, and under
 * unusual circumstances the iterative method may not converge. 
 *
 * @param K      the statistical parameter K
 * @param logK   the natural logarithm of K
 * @param alpha_d_lambda    the ratio of the statistical parameters
 *                          alpha and lambda (for ungapped alignments, the
 *                          value 1/H should be used)
 * @param beta              the statistical parameter beta (for ungapped
 *                          alignments, beta == 0)
 * @param query_length      the length of the query sequence
 * @param db_length         the length of the database
 * @param db_num_seq        the number of sequences in the database
 * @param length_adjustment the computed value of the length adjustment [out]
 *
 * @return   0 if length_adjustment is known to be the largest integer less
 *           than the fixed point of f(ell); 1 otherwise.
 */

#define maximum(a,b) ((a > b) ? a : b)

int4 BlastComputeLengthAdjustment(float K,
                             float logK,
                             float alpha_d_lambda,
                             float beta,
                             int4 query_length,
                             uint4 db_length,
                             int4 db_num_seqs,
                             int4 *length_adjustment)
{
    int4 i;                     /* iteration index */
    const int4 maxits = 20;     /* maximum allowed iterations */
    double m = query_length, n = db_length, N = db_num_seqs;

    double ell;            /* A float value of the length adjustment */
    double ss;             /* effective size of the search space */
    double ell_min = 0, ell_max;   /* At each iteration i,
                                         * ell_min <= ell <= ell_max. */
    char converged = 0;       /* True if the iteration converged */
    double ell_next = 0;   /* Value the variable ell takes at iteration
                                 * i + 1 */
    /* Choose ell_max to be the largest nonnegative value that satisfies
     *
     *    K * (m - ell) * (n - N * ell) > max(m,n)
     *
     * Use quadratic formula: 2 c /( - b + sqrt( b*b - 4 * a * c )) */
    { /* scope of a, mb, and c, the coefficients in the quadratic formula
       * (the variable mb is -b) */
        double a  = N;
        double mb = m * N + n;
        double c  = n * m - maximum(m, n) / K;

        if(c < 0) {
            *length_adjustment = 0;
            return 1;
        } else {
            ell_max = 2 * c / (mb + sqrt(mb * mb - 4 * a * c));
        }
    } /* end scope of a, mb and c */

    for(i = 1; i <= maxits; i++) {      /* for all iteration indices */
        double ell_bar;    /* proposed next value of ell */
        ell      = ell_next;
        ss       = (m - ell) * (n - N * ell);
        ell_bar  = alpha_d_lambda * (logK + log(ss)) + beta;
        if(ell_bar >= ell) { /* ell is no bigger than the true fixed point */
            ell_min = ell;
            if(ell_bar - ell_min <= 1.0) {
                converged = 1;
                break;
            }
            if(ell_min == ell_max) { /* There are no more points to check */
                break;
            }
        } else { /* else ell is greater than the true fixed point */
            ell_max = ell;
        }
        if(ell_min <= ell_bar && ell_bar <= ell_max) {
          /* ell_bar is in range. Accept it */
            ell_next = ell_bar;
        } else { /* else ell_bar is not in range. Reject it */
            ell_next = (i == 1) ? ell_max : (ell_min + ell_max) / 2;
        }
    } /* end for all iteration indices */
    if(converged) { /* the iteration converged */
        /* If ell_fixed is the (unknown) true fixed point, then we
         * wish to set (*length_adjustment) to floor(ell_fixed).  We
         * assume that floor(ell_min) = floor(ell_fixed) */
        *length_adjustment = (int4) ell_min;
        /* But verify that ceil(ell_min) != floor(ell_fixed) */
        ell = ceil(ell_min);
        if( ell <= ell_max ) {
          ss = (m - ell) * (n - N * ell);
          if(alpha_d_lambda * (logK + log(ss)) + beta >= ell) {
            /* ceil(ell_min) == floor(ell_fixed) */
            *length_adjustment = (int4) ell;
          }
        }
    } else { /* else the iteration did not converge. */
        /* Use the best value seen so far */
        *length_adjustment = (int4) ell_min;
    }

    return converged ? 0 : 1;
}
