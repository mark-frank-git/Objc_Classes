/************************************************************************
 * This subclass of object implements a polynomial object.  Three       *
 * polynomials are contained in this object: a, b, c.  Most of the      *
 * operations take place on the a polynomial, with the b and c poly-    *
 * nomials used for multiplying, etc.                                   *
 *                                                                      *
 * File:Polynomial.m                                                    *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/18/92  - Started                                              *
 ************************************************************************/

#import "Polynomial.h"
#import <stdlib.h>

#define MAX(a,b)  ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)    ((a) >= 0 ? (a) : (-a))

void   jt(COMPLEX *coeff, int degree, COMPLEX *ans, int scale, int *fail);

@implementation Polynomial

- init
{
  [super init];
  
// init instance variables:
  aOrder = bOrder = cOrder  = 0;
  aLength = bLength = cLength   = 0;
  a = b = c             = NULL;
  aRoots            = NULL;
  
  return self;
}


- initWithAPoly: (double *)aPoly order:(int)n
{
  int i;
  
  [self init];
  
  aOrder  = n;
  aOrder  = MAX(aOrder, 0);
  aLength = 2*n+1;
  a   = (double *)malloc(aLength*sizeof(double));
  for(i=0; i<n+1; i++)
    a[i] = aPoly[i];

  return self;
}

/*##############################*
 * Set a new a polynomial:  *
 *##############################*/
- setAPoly: (double *)aPoly order:(int)n
{
  int i;
  
  aOrder = n;
  aOrder  = MAX(aOrder, 0);
  if((n+1)>aLength)
  {
    if(a != NULL)
      free(a);
    aLength = 2*n+1;
    a = (double *)malloc(aLength*sizeof(double));
  }
  for(i=0; i<(n+1); i++)
    a[i] = aPoly[i];
  
  return self;
}

/*##############################*
 * Set a new b polynomial:  *
 *##############################*/
- setBPoly: (double *)bPoly order:(int)n
{
  int i;
  
  bOrder = n;
  bOrder  = MAX(bOrder, 0);
  if((n+1)>bLength)
  {
    if(b != NULL)
      free(b);
    bLength = 2*n+1;
    b = (double *)malloc(bLength*sizeof(double));
  }
  for(i=0; i<(n+1); i++)
    b[i] = bPoly[i];
  
  return self;
}

/*##############################*
 * Get parameters:      *
 *##############################*/
- (double *)getAPoly {return a;}
- (double *)getBPoly {return b;}
- (double *)getCPoly {return c;}
- (int)getAOrder     {return aOrder;}
- (int)getBOrder     {return bOrder;}
- (int)getCOrder     {return cOrder;}


/*##############################*
 * Copy polynomial in c to a    *
 *##############################*/
- copyCToA
{
  int i, n;
  
  n      = cOrder;
  aOrder = n;
  if((n+1)>aLength)
  {
    if(a != NULL)
      free(a);
    aLength = 2*n+1;
    a = (double *)malloc(aLength*sizeof(double));
  }
  for(i=0; i<(n+1); i++)
    a[i] = c[i];
  
  return self;
}

/*##############################*
 * Multiply polynomial in a by  *
 * polynomial in b, store result*
 * in c.            *
 *##############################*/
- multABToC
{
  int i, j, k, n;
  

  cOrder = aOrder+bOrder;
  n      = cOrder;
  if((n+1)>cLength)
  {
    if(c != NULL)
      free(c);
    cLength = 2*n+1;
    c = (double *)malloc(cLength*sizeof(double));
  }

  for(i=0; i<=n; i++)
  {
    c[i] = 0.;
    for(j=0; j<=aOrder; j++)
    {
      for(k=0; k<=bOrder; k++)
      {
        if(j+k == i)
          c[i] += a[j]*b[k];
      }
    }
  }
  return self;
}

/*##############################*
 * Reverse polynomial in a. *
 *##############################*/
- reverseA
{
  int i, n;
  double temp;
  
  n = aOrder;
  for(i=0; i<(aOrder+1)/2; i++)
  {
    temp = a[i];
    a[i] = a[n];
    a[n--] = temp;
  }
  return self;
}

/*##############################*
 * Shift the polynomial in a    *
 * right, and fill in with zeros*
 *##############################*/
- shiftARightBy: (int)shift
{
  int i, n, length;
  
  n = shift;
  length = aOrder+1;
  for(i=0; i<shift; i++)
  {
    if(n<length)
      a[n++] = a[i];
    a[i] = 0.;
  }
  return self;
}

/*##############################*
 * Evaluate the polynomial in a *
 * at a complex point.      *
 *##############################*/
- (COMPLEX)evaluateAAtComplexPoint:(COMPLEX)point
{
  int     i;
  COMPLEX result, cprod, ctemp, csum; 

/****************************
 * Evaluate polynomial using*
 * Horner's method:         *
 ****************************/
  if(aOrder>0)
  {
   CMPLX(ctemp,a[aOrder],0.);
   CMULT(cprod,ctemp,point);  
   for(i=(aOrder-1); i>0; i--)
   {                  
     CMPLX(ctemp,a[i],0.);
     CADD(csum,ctemp,cprod);
     CMULT(cprod,csum,point); 
   } 
   CMPLX(ctemp,a[0],0.);    
   CADD(result,cprod,ctemp);
  }
  else
    CMPLX(result, a[0], 0.);  
  return result;
}


#define REAL_SCALE  1.e5
/*###############################*
 * Find the roots of the poly-  *
 * nomial in a.         *
 * return roots in the form     *
 * c[0], c[1] .. r[0], r[1] ... *
 * c*[0], c*[1] ...             *
 * where the c's are the complex*
 * roots, and the r's are the   *
 * real roots.                  *
 *##############################*/
- (COMPLEX *)rootsOfA
{
  int i, j, fail, scale;
  int number_complex, number_real;
  BOOL    new_complex_root;
  COMPLEX *complex_a, *complex_roots, *real_roots;
  
  complex_a = (COMPLEX *)malloc((aOrder+1)*sizeof(COMPLEX));
  j = aOrder;
  for(i=0; i<=aOrder; i++)
    CMPLX(complex_a[i], a[j--], 0.);
    
  if(aRoots != NULL)
    free(aRoots);
  aRoots = (COMPLEX *)malloc(aOrder*sizeof(COMPLEX));

/*******************************
 * Use Jenkins-Traub, scale    *
 * sets the accuracy.          *
 *******************************/
  scale = 2;
  fail = 1;
  while(fail || scale<10000)
  {
    scale *= 2;
    jt(complex_a, aOrder, aRoots, scale, &fail);
  }
  free(complex_a);
  
/*******************************
 * Now, re-order roots:        *
 *******************************/
  real_roots    = (COMPLEX *)malloc(aOrder*sizeof(COMPLEX));
  complex_roots = (COMPLEX *)malloc(aOrder*sizeof(COMPLEX));
  number_complex = number_real = 0;
  for(i=0; i<aOrder; i++)
  {
    if( ABS(aRoots[i].x) > REAL_SCALE*ABS(aRoots[i].y) )
    {
      real_roots[number_real] = aRoots[i];
      number_real++;
    }
    else
    {
      new_complex_root = YES;
      for(j=0; j<number_complex; j++)
      {
        if( ABS(aRoots[i].x - complex_roots[j].x) <=
                            ABS(aRoots[i].x/REAL_SCALE) )
    {
      new_complex_root = NO;
      break;
    }
      }
      if( new_complex_root )
      {
        complex_roots[number_complex] = aRoots[i];
        number_complex++;
      }
    }
  }
  number_complex = aOrder - number_real;    /* error check */
  number_complex /= 2;
  for(i=0; i<number_complex; i++)
    aRoots[i] = complex_roots[i];
  for(j=0; j<number_real; j++)
  {
    aRoots[i] = real_roots[j];
    i++;
  }
  for(j=0; j<number_complex; j++)   /* conjugate roots */
  {
    aRoots[i].x = complex_roots[j].x;
    aRoots[i++].y = -complex_roots[j].y;
  }
  free(real_roots);
  free(complex_roots);
      
  return aRoots;
}


@end