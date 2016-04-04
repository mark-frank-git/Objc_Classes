/************************************************************************
 * This subclass of object implements a statistical capture object for  *
 * calculating statistics.                                              *
 *                                                                      *
 * File:DesStat.m                                                       *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/26/92  - Started                                              *
 *  2. 08/27/92  - Use symmetric storage for arraySum                   *
 ************************************************************************/
 
#import "DesStat.h"
#import <stdlib.h>
#import <stdio.h>
#import <math.h>

#define ABS(a)    ((a) >= 0 ? (a) : (-a))
#define SIGN(a,b) ( ((b)>=0 ) ? ABS((a)) : -ABS((a)) )


   
@implementation DesStat


- init
{
  [super init];
  vecSum    = NULL;
  arraySum  = NULL;
  mean      = NULL;
  stdDev    = NULL;
  covar     = NULL;
  normCovar = NULL;
   return self;
}

- initArrays
{
  if(vecSum != NULL)
  {
    free(vecSum);
    free(arraySum);
    free(mean);
    free(stdDev);
    free(covar);
    free(normCovar);
  }
//
// initialize instance variables:
//
  arraySize     = (dim*(dim+1))/2;
  count         = 0;
  countAtMean   = -1;
  countAtStdDev = -1;
  countAtCovar  = -1;
  countAtNormCovar  = -1;
  vecSum        = (long double *)calloc(dim,sizeof(long double));     
  arraySum      = (long double *)calloc(arraySize,sizeof(long double));     
  mean          = (float *)malloc(dim*sizeof(float));     
  stdDev        = (float *)malloc(dim*sizeof(float));     
  covar         = (float *)malloc(dim*dim*sizeof(float));     
  normCovar     = (float *)malloc(dim*dim*sizeof(float));
  
  return self; 
}

- initWithDimension:(int)dimension
{
  
  [self init];
  dim           = dimension;
  [self initArrays];
  
  return self;
}

/*#########################*
 * The following methods   *
 * reset parameters:       *
 *#########################*/
- resetWithDimension:(int)dimension
{
  dim           = dimension;
  [self initArrays];
  
  return self;
}

- reset
{ int       i;

  count         = 0;
  countAtMean   = -1;
  countAtStdDev = -1;
  countAtCovar  = -1;
  countAtNormCovar  = -1;
  for( i=0; i<dim; i++)
    vecSum[i] = 0.;
  for( i=0; i<arraySize; i++)
    arraySum[i] = 0.;
  
  return self;
}

/*#########################*
 * The following methods   *
 * return parameters:      *
 *#########################*/
- (int)getDim    { return dim; }
- (long)getCountAtMean      { return countAtMean; }
- (long)getCountAtStdDev    { return countAtStdDev; }
- (long)getCountAtCovar     { return countAtCovar; }
- (long)getCountAtNormCovar { return countAtNormCovar; }
- (long)getCount    { return count; }
- (long double *)getVecSum { return vecSum; }
- (long double *)getArraySum { return arraySum; }

- (float *)getMean
{
  int       i;
  long double dcount;

  if( count == 0 ) return NULL;

  dcount = (long double)count;

  if( countAtMean != count )
  {
    for( i=0; i<dim; i++ )
      mean[i] = (float) (vecSum[i]/dcount);
    countAtMean = count;
  }
  return mean;
}

- (float *)getStdDev
{   
  int  i, j;
  long double dcount;

  if( count == 0 ) return NULL;

  dcount = (long double)count;

  if( countAtStdDev != count )
  { if( countAtMean != count )
      [self getMean];
    j = 0;
    for( i=0; i<dim; i++ )
    {
      stdDev[i]=sqrt( arraySum[j]/dcount - mean[i]*mean[i]);
      j += dim-i;
    }
    countAtStdDev = count;
  }
  return stdDev;
}

- (float *)getCovar
{   
  int  i,j,k;
  long double dcount;

  if( count == 0 ) return NULL;

  dcount = (long double)count;
    
  if( countAtCovar != count )
  { if( countAtMean != count )
      [self getMean];
    k = 0;
    for( i=0; i<dim; i++ )
    {
      for( j=0; j<dim; j++ )
      {
        if( i <= j )
          covar[i*dim+j] = arraySum[k++]/dcount - mean[i]*mean[j];
        else
          covar[i*dim+j] = covar[j*dim+i];
      }
    }
    countAtCovar = count;
  }
  return covar;
}

- (float *)getNormCovar
{   
  int  i,j;
  long double dcount;

  if( count == 0 ) return NULL;

  dcount = (long double)count;
    
  if( countAtNormCovar != count )
  { if( countAtCovar != count )
      [self getCovar];
    if( countAtStdDev != count )
    { for( i=0; i<dim; i++ )
        stdDev[i]=sqrt( covar[i*dim+i] );
      countAtStdDev = count;
    }
    for( i=0; i<dim; i++ )
      for( j=0; j<dim; j++ )
        normCovar[i*dim+j] = covar[i*dim+j]/( stdDev[i]*stdDev[j] );
    countAtNormCovar = count;
  }
  return normCovar;
}

/*#########################*
 * Add in new data:        *
 *                         *
 * Note: Only upper tri-   *
 * angular part of matrix  *
 * is used.Also, k=i*dim+j *
 *#########################*/
- addData: (long double *)pData
{   
  int i,j,k;

  k = 0;
  for( i=0; i<dim; i++ )
  {
    vecSum[i] += pData[i];
    for( j=i; j<dim; j++ )
      arraySum[k++] += pData[i]*pData[j];
  }
  count++;
    
  return self;
}

@end