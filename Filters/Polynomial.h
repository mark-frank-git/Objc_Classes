/************************************************************************
 * This subclass of object implements a polynomial object.  Three       *
 * polynomials are contained in this object: a, b, c.  Most of the      *
 * operations take place on the a polynomial, with the b and c poly-    *
 * nomials used for multiplying, etc.                                   *
 *                                                                      *
 * File:Polynomial.h                                                    *
 *                                                                      *
 * Note: for p(x) = a[n]*x^n + a[n-11]*x^(n-1) + . . . + a[0]           *
 *       p(x) has order n, but the length of the array is n+1.          *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 08/18/92  - Started                                              *
 ************************************************************************/
#import <objc/Object.h>
#import "c_headers.h"

@interface Polynomial:Object
{
  int   aOrder, bOrder, cOrder;                     /* orders of the polynomials    */
  int   aLength, bLength, cLength;                  /* actual lengths of the arrays */
  double *a, *b, *c;                                /* polynomials          */

  COMPLEX *aRoots;                                  /* Complex roots of a       */
}

/**********************
 * Initialization:    *
 **********************/
- initWithAPoly: (double *)aPoly order:(int)n;

/**********************
 * Set parameters:    *
 **********************/
- setAPoly: (double *)aPoly order:(int)n;
- setBPoly: (double *)bPoly order:(int)n;

/**********************
 * Get parameters:    *
 **********************/
- (double *)getAPoly;
- (double *)getBPoly;
- (double *)getCPoly;
- (int)getAOrder;
- (int)getBOrder;
- (int)getCOrder;

/**********************
 * Calculating:       *
 **********************/
- copyCToA;
- multABToC;
- reverseA;
- shiftARightBy:(int)shift;
- (COMPLEX)evaluateAAtComplexPoint:(COMPLEX)point;
- (COMPLEX *)rootsOfA;


@end