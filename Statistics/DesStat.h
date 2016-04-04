/************************************************************************
 * This subclass of object implements a statistical capture object for  *
 * calculating statistics.                                              *
 *                                                                      *
 * File:DesStat.h                                                       *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/26/92  - Started                                              *
 ************************************************************************/
#import <objc/Object.h>


@interface DesStat:Object
{
    int         dim;                /* dimension of data */
    int         arraySize;      /* symmetric storage size */
    long int    count;
    long int    countAtMean;    /* dates last cal'c of mean */
    long int    countAtStdDev;  /* dates last cal'c of std dev */
    long int    countAtCovar;
    long int    countAtNormCovar;
    long double *vecSum;
    long double *arraySum;
    float       *mean;      /* mean */
    float       *stdDev;    /* standard deviation */
    float       *covar;     /* covariance matrix */
    float       *normCovar; /* normalized covariance matrix */
}


/**********************
 * Initialization:    *
 **********************/
- init;
- initWithDimension:(int)dim;

/**********************
 * Setting parameters *
 **********************/
- resetWithDimension:(int)dim;
- reset;

/**********************
 * Getting parameters:*
 **********************/
- (int)getDim;
- (long)getCountAtMean;
- (long)getCountAtStdDev;
- (long)getCountAtCovar;
- (long)getCountAtNormCovar;
- (long)getCount;
- (long double *)getVecSum;
- (long double *)getArraySum;
- (float *)getMean;
- (float *)getStdDev;
- (float *)getCovar;
- (float *)getNormCovar;

/**********************
 * Accumulating       *
 * statistics:        *
 **********************/
- addData: (long double *)pData;


@end