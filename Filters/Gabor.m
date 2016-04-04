/************************************************************************
 * This subclass of object implements a discrete Gabor filter as        *
 * described in:                                                        *
 *                                                                      *
 * J. Wexler and S. Raz, "Discrete Gabor expansions," Signal Processing,*
 *  vol. 21, no. 3, pp. 207-220,  Nov. 1990.                            *
 *                                                                      *
 * File: Gabor.m                                                        *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/14/93  - Started                                              *
 ************************************************************************/

/**************************
 * Include files:         *
 **************************/
#import "Gabor.h"
#import "ComplexMatrix.h"
#import <math.h>
#import <stdio.h>
#import <stdlib.h>

@implementation Gabor

- initWithType:(int)type length:(int)length width:(int)width shift:(int)shift
{
  int n;
  [super init];

  windowType    = type;
  windowLength  = length;
  effectiveWidth= width;
  shiftParameter= shift;
  maximumN      = shift/2;
  if(shift>0)
    maximumM    = length/shift/2;
  else
    maximumM    = 0;

  n = windowLength;
  [self allocateMemory];
  [self calculateWindow];
  [self calculateBiorthogonal];
  
  return self;
}

/*##############################*
 * Allocate memory needed:  *
 *##############################*/
- allocateMemory
{
  int n;

  n = windowLength;
  complexMatrix        = [[ComplexMatrix alloc] initWithRows:n columns:n];
  window               = (double *)malloc(n*sizeof(double));
  shiftedWindow        = (double *)malloc(n*sizeof(double));
  biorthogonalFunction = (double *)malloc(n*sizeof(double));
  reconstructedSignal  = (double *)malloc(n*sizeof(double));
  expansionCoefficients = (COMPLEX *)malloc(n*sizeof(COMPLEX));
  return self;
}

/*##############################*
 * Free up memory       *
 *##############################*/
- freeMemory
{
  [complexMatrix free];
  free(window);
  free(shiftedWindow);
  free(biorthogonalFunction);
  free(expansionCoefficients);
  free(reconstructedSignal);
  window = biorthogonalFunction = NULL;
  expansionCoefficients = NULL;
  return self;
}

/*##############################*
 * Free myself:         *
 *##############################*/
- free
{
  [self freeMemory];
  [super free];
  return self;
}

/*##############################*
 * Calculate the window function*
 *##############################*/
- calculateWindow
{
  int i, j, n, start;
  double scale, center, sum;

  scale = sqrt(1./effectiveWidth);
  switch(windowType)
  {
    case RECT_WINDOW:
    default:
      for(i=0; i<windowLength; i++)
        window[i] = 0.;
      start = (windowLength-effectiveWidth)/2;
      for(i=start; i<start+effectiveWidth; i++)
        window[i] = scale;
      break;
    case GAUSSIAN_WINDOW:
      center = (windowLength - 1.)/2.;
      scale *= sqrt(SQRT2);
      for(i=0; i<windowLength; i++)
        window[i] = scale*exp(-PI*(i-center)*(i-center)/
                    effectiveWidth/effectiveWidth);
      break;
    case RAISED_COS:
      for(i=0; i<windowLength; i++)
        window[i] = 0.;
      start = (windowLength - 1)/2 - effectiveWidth + 2;
      j = 0;
      n = 1.8*effectiveWidth;
      sum = 0.;
      for(i=start; i<start+n; i++)
      {
        window[i] = (1.-cos(TWOPI*j/(n-1)));
        window[i] = MAX(0., window[i]);
        sum += window[i]*window[i];
        j++;
      }
      scale = 1./sqrt(sum);
      for(i=start; i<start+n; i++)
        window[i] *= scale;
      break;
  }
  return self;
}

/*##############################*
 * Calculate the Biorthogonal   *
 * function.            *
 *##############################*/
- calculateBiorthogonal
{
  int i, j, k, m, n, index, m_max;
  COMPLEX *b, *input, *solution;
  COMPLEX w, cexp, ctemp;

  b = (COMPLEX *)malloc(windowLength*sizeof(COMPLEX));
  input = (COMPLEX *)malloc(windowLength*windowLength*sizeof(COMPLEX));

  if(shiftParameter>0 && windowLength>shiftParameter)
    m_max = windowLength/shiftParameter;
  else
    return self;
  i=j=0;
  for(m=0; m<m_max; m++)
  {
    for(n=0; n<shiftParameter; n++)
    {
      for(k=0; k<windowLength; k++)
      {
        index = (k+m*shiftParameter)%(windowLength);
        CMPLX(ctemp, window[index], 0.);
        CMPLX(cexp, 0., TWOPI*n*k/shiftParameter);
        CEXP(w, cexp);
        CMULT(input[i], ctemp, w);
        i++;
      }
      if((n==0) && (m==0))
      {
        CMPLX(b[j],1.,0.);
      }
      else
      {
        CMPLX(b[j],0.,0.);
      }
      j++;
    }
  }
  [complexMatrix setMatrix:input];
  [complexMatrix linearSolve:b];
  solution = [complexMatrix xSolution];
  for(i=0; i<windowLength; i++)
    biorthogonalFunction[i] = solution[i].x;

  free(b);
  free(input);
  return self;  
} 

/*##############################*
 * Calculate the expansion  *
 * coefficients for an input    *
 * sequence.            *
 *##############################*/
- calculateExpansionCoefficientsFor:(double *)inputX
{
  int i, k, m, n, index, m_max;
  double  temp;
  COMPLEX w, cexp, ctemp, ctemp1;

  if(biorthogonalFunction == NULL)
    return self;

  if(shiftParameter>0 && windowLength>shiftParameter)
    m_max = windowLength/shiftParameter;
  else
    return self;
  i = 0;
  for(m=0; m<m_max; m++)
  {
    for(n=0; n<shiftParameter; n++)
    {
      CMPLX(expansionCoefficients[i], 0., 0.);
      for(k=0; k<windowLength; k++)
      {
        index = (k-m*shiftParameter);
        while(index<0)
          index += windowLength;
        index %= windowLength;
        temp  = inputX[k]*biorthogonalFunction[index];
        CMPLX(ctemp, temp, 0.);                        /* x[k]*gamma(index) */
        CMPLX(cexp, 0., TWOPI*n*k/shiftParameter);     
        CEXP(w, cexp);                                 /* w(nk)             */
        CMULT(ctemp1, w, ctemp);                       /* x[k]gammamn(k)    */
        expansionCoefficients[i].x += ctemp1.x;
        expansionCoefficients[i].y -= ctemp1.y;
      }
      i++;
    }
  }
  return self;  
} 

/*##############################*
 * Calculate the expansion  *
 * coefficients for an input    *
 * sequence.            *
 *##############################*/
- calculateReducedCoefficientSetFor: (double *)inputX
{
  int i, k, m, n, m_max, m_length, n_max;
  int c_index, b_index, m_symmetric, n_symmetric;
  double  temp;
  COMPLEX w, cexp, ctemp, ctemp1;

  if(biorthogonalFunction == NULL)
    return self;

  if(shiftParameter>0 && windowLength>shiftParameter)
    m_length = windowLength/shiftParameter;
  else
    return self;
  m_max = MIN(m_length/2, maximumM);
  n_max = MIN(shiftParameter/2, maximumN);
  for(m=0; m<=m_max; m++)
  {
    c_index = m*shiftParameter;
    for(n=0; n<=n_max; n++)
    {
      i = n + c_index;
      CMPLX(expansionCoefficients[i], 0., 0.);
      for(k=0; k<windowLength; k++)
      {
        b_index = (k-m*shiftParameter);
        while(b_index<0)
          b_index += windowLength;
        temp  = inputX[k]*biorthogonalFunction[b_index];
        CMPLX(ctemp, temp, 0.);                        /* x[k]*gamma(index) */
        CMPLX(cexp, 0., TWOPI*n*k/shiftParameter);     
        CEXP(w, cexp);                                 /* w(nk)             */
        CMULT(ctemp1, w, ctemp);                       /* x[k]gammamn(k)    */
        expansionCoefficients[i].x += ctemp1.x;
        expansionCoefficients[i].y -= ctemp1.y;
      }
    }
    c_index = m*shiftParameter;
    for(n=1; (n<=n_max)&&(n<shiftParameter/2); n++)
    {
      n_symmetric = shiftParameter-n;
      i           = c_index + n_symmetric;
      CMPLX(expansionCoefficients[i], 0., 0.);
      for(k=0; k<windowLength; k++)
      {
        b_index = (k-m*shiftParameter);
        while(b_index<0)
          b_index += windowLength;
        temp  = inputX[k]*biorthogonalFunction[b_index];
        CMPLX(ctemp, temp, 0.);                        /* x[k]*gamma(index) */
        CMPLX(cexp, 0., TWOPI*n_symmetric*k/shiftParameter);     
        CEXP(w, cexp);                                 /* w(nk)             */
        CMULT(ctemp1, w, ctemp);                       /* x[k]gammamn(k)    */
        expansionCoefficients[i].x += ctemp1.x;
        expansionCoefficients[i].y -= ctemp1.y;
      }
    }
  }
  for(m=1; m<m_max; m++)
  {
    m_symmetric = m_length - m;
    c_index = m_symmetric*shiftParameter;
    for(n=0; n<=n_max; n++)
    {
      i = c_index + n;
      CMPLX(expansionCoefficients[i], 0., 0.);
      for(k=0; k<windowLength; k++)
      {
        b_index = (k-m_symmetric*shiftParameter);
        while(b_index<0)
          b_index += windowLength;
        temp  = inputX[k]*biorthogonalFunction[b_index];
        CMPLX(ctemp, temp, 0.);                        /* x[k]*gamma(index) */
        CMPLX(cexp, 0., TWOPI*n*k/shiftParameter);     
        CEXP(w, cexp);                                 /* w(nk)             */
        CMULT(ctemp1, w, ctemp);                       /* x[k]gammamn(k)    */
        expansionCoefficients[i].x += ctemp1.x;
        expansionCoefficients[i].y -= ctemp1.y;
      }
    }
    c_index = m_symmetric*shiftParameter;
    for(n=1; (n<=n_max)&&(n<shiftParameter/2); n++)
    {
      n_symmetric = shiftParameter-n;
      i = c_index + n_symmetric;
      CMPLX(expansionCoefficients[i], 0., 0.);
      for(k=0; k<windowLength; k++)
      {
        b_index = (k-m_symmetric*shiftParameter);
        while(b_index<0)
          b_index += windowLength;
        temp  = inputX[k]*biorthogonalFunction[b_index];
        CMPLX(ctemp, temp, 0.);                        /* x[k]*gamma(index) */
        CMPLX(cexp, 0., TWOPI*n_symmetric*k/shiftParameter);     
        CEXP(w, cexp);                                 /* w(nk)             */
        CMULT(ctemp1, w, ctemp);                       /* x[k]gammamn(k)    */
        expansionCoefficients[i].x += ctemp1.x;
        expansionCoefficients[i].y -= ctemp1.y;
      }
    }
  }
  return self;  
} 

/*##############################*
 * Reconstruct the signal from  *
 * the expansion coefficients.  *
 *##############################*/
- reconstructSignal
{
  int i, k, m, n, w_index, c_index, m_length, m_max, n_max;
  int m_symmetric, n_symmetric;
  double sum;
  COMPLEX w, cexp, ctemp, ctemp1, hkmn;

  if(biorthogonalFunction == NULL)
    return self;

  if(shiftParameter>0 && windowLength>shiftParameter)
    m_length = windowLength/shiftParameter;
  else
    return self;
  m_max = MIN(m_length/2, maximumM);
  n_max = MIN(shiftParameter/2, maximumN);
/* We use the symmetry of the expansion coefficients below */
  for(k=0; k<windowLength; k++)
  {
    sum = 0.;
    for(m=0; m<=m_max; m++)
    {
      c_index = m*shiftParameter;
      w_index = (k-m*shiftParameter);
      while(w_index<0)
        w_index += windowLength;
      CMPLX(hkmn, window[w_index], 0.);              /* h(k-mN)      */
      for(n=0; n<=n_max; n++)
      {
        CMPLX(cexp, 0., TWOPI*n*k/shiftParameter);
        CEXP(w, cexp);                               /* w(nk)        */
        CMULT(ctemp1, hkmn, w);                      /* hmn(k)       */
        i = c_index + n;
        CMULT(ctemp, ctemp1, expansionCoefficients[i]); /* amn*hmn(k) */
        sum    += ctemp.x;
        if(n!=0 && n!= shiftParameter/2)
        {
          n_symmetric = shiftParameter - n;
          CMPLX(cexp, 0., TWOPI*n_symmetric*k/shiftParameter);
          CEXP(w, cexp);                              /* w(nk)        */
          CMULT(ctemp1, hkmn, w);                     /* hmn(k)       */
          i = c_index + n_symmetric;
          CMULT(ctemp, ctemp1, expansionCoefficients[i]);/* amn*hmn(k) */
          sum    += ctemp.x;
        }
      }
    }
    for(m=1; (m<m_max && m<m_length/2); m++)
    {
      m_symmetric = m_length - m;
      c_index = m_symmetric*shiftParameter;
      w_index = (k-m_symmetric*shiftParameter);
      while(w_index<0)
        w_index += windowLength;
      CMPLX(hkmn, window[w_index], 0.);              /* h(k-mN)      */
      for(n=0; n<=n_max; n++)
      {
        CMPLX(cexp, 0., TWOPI*n*k/shiftParameter);
        CEXP(w, cexp);                               /* w(nk)        */
        CMULT(ctemp1, hkmn, w);                      /* hmn(k)       */
        i = c_index + n;
        CMULT(ctemp, ctemp1, expansionCoefficients[i]); /* amn*hmn(k) */
        sum    += ctemp.x;
        if(n!=0 && n!= shiftParameter/2)
        {
          n_symmetric = shiftParameter - n;
          CMPLX(cexp, 0., TWOPI*n_symmetric*k/shiftParameter);
          CEXP(w, cexp);                              /* w(nk)        */
          CMULT(ctemp1, hkmn, w);                     /* hmn(k)       */
          i = c_index + n_symmetric;
          CMULT(ctemp, ctemp1, expansionCoefficients[i]);/* amn*hmn(k) */
          sum    += ctemp.x;
        }
      }
    }
    reconstructedSignal[k] = sum;
  }
  return self;  
}

/*########################*
 * Zero out the expansion *
 * Coefficients       *
 *########################*/
- resetCoefficients;
{
  int i;
  for(i=0; i<windowLength; i++)
  {
    CMPLX(expansionCoefficients[i], 0., 0.);
  }
  return self;
}

/*########################*
 * Set instance variable  *
 *########################*/
- setWindowType: (int)type
{
  windowType = type;
  [self calculateWindow];
  [self calculateBiorthogonal];
  return self;
}
- setShift: (int)shift
{
  shiftParameter = shift;
  [self calculateWindow];
  [self calculateBiorthogonal];
  return self;
}
- setWindowLength: (int)length
{
  windowLength = length;
  [self freeMemory];
  [self allocateMemory];
  [self calculateWindow];
  [self calculateBiorthogonal];
  return self;
}
- setWindowWidth: (int)width
{
  effectiveWidth = width;
  [self calculateWindow];
  [self calculateBiorthogonal];
  return self;
}
- setMaximumM: (int)m
{
  maximumM = m;
  return self;
}
- setMaximumN: (int)n
{
  maximumN = n;
  return self;
}

/*#######################*
 * Get instance variables*
 *#######################*/
- (int)windowType               { return windowType;}
- (int)windowLength             { return windowLength;}
- (int)windowWidth              { return effectiveWidth;}
- (int)windowShift      { return shiftParameter;}
- (int)mValue { return ((shiftParameter>0)?(windowLength/shiftParameter):1);}
- (int)maximumM                 { return maximumM;}
- (int)maximumN         { return maximumN;}
- (double *)windowFunction      { return window;}
- (double *)biorthogonalFunction{ return biorthogonalFunction;}
- (double *)reconstructedSignal { return reconstructedSignal;}
- (COMPLEX *)expansionCoefficients { return expansionCoefficients;}
- (double *)shiftedWindow:(int)shift
{
  int i, index;

  for(i=0; i<windowLength; i++)
  {
    index = (i-shift);
    while(index<0)
      index += windowLength;
    index %= windowLength;
    shiftedWindow[i] = window[index];
  }
  return shiftedWindow;
}

@end