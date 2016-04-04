/***************************************************************************
 * This subclass of object implements a general probability distribution   *
 * function.                                                               *
 *                                                                         *
 * File: /User/frank/Objc_Classes/Probability/Probability.h                *
 *                                                                         *
 * Revision History:                                                       *
 *  1. 05/04/92 - Started                                                  *
 *  2. 05/18/92 - Add dimensions                                           *
 *  3. 07/09/92 - Add gamma, log_gamma, log_normal, and correlation        *
 *                coefficient per Bill Kushner.                            *
 *  4. 11/19/92 - Add GAUSS_MARKOV, distributionFunction -> densityFunction*
 *  5. 02/10/94 - Add Rician, Rayleigh, noncentral chi squared              *
 *  6. 03/24/94 - Add normal01                                              *
 ***************************************************************************/

#import <objc/Object.h>

#define NORMAL              0
#define UNIFORM             1
#define EXPONENTIAL_DIST    2
#define ERLANG_DIST         3               /* special case of Gamma dist   */
#define LOG_NORMAL          4
#define GAUSS_MARKOV        5
#define RICIAN              6
#define RAYLEIGH            7
#define NON_CENTRAL_CHI     8               /* non central chi squared      */
#define LOG_RICIAN          9



@interface Probability:Object
{
  id     randomObject;                      /* random # generator               */
  double *mean;                             /* mean of the distribution         */
  double *variance;                         /* variance                         */
  double *stdDev;                           /* standard deviation               */
  double rho;                               /* correlation coefficient          */
  double probability;                       /* Probability of occurrence        */
  double noncentralParameter;               /* NA**2, for chi squared           */
  int    distribution;                      /* distribution type, see above     */
  int    dimensions;                        /* # of dimensions in distribution  */
  int   degrees;                            /* Degrees of freedom in noncent    */
  double *doubleOutputs;                    /* Used for array outputs           */
}


/*********************************
 * Initialization:               *
 *********************************/
- init;
- initMeanVariance;
- initializeSeeds;

/*******************************
 * These methods set parameters:*
 *******************************/
- setMean:(double*)newMean;
- setVariance: (double*)newVariance;
- setCorrelationCoefficient: (double)newRho;
- setProbDistribution: (int)newDistribution;
- setProbability: (double)newProbability;
- setDimensions: (int)newDimension;
- setNoncentralParameter:(double)parameter;
- setDegreesFreedom:(int)degree;

/*******************************
 * These methods get parameters*
 *******************************/
- (double *)getMean;
- (double *)getVariance;
- (int)   getProbDistribution;
- (double)getProbability;
- (int) getDimensions;
- (int) degrees;

/*******************************
 * These methods get outputs:  *
 *******************************/
- (double *)newSample;
- (double  )densityFunctionAtX: (double *)x;
- (double  )normal01;

@end