/************************************************************************
 * This subclass of object implements a discrete Gabor filter as        *
 * described in:                                                        *
 *                                                                      *
 * J. Wexler and S. Raz, "Discrete Gabor expansions," Signal Processing,*
 *  vol. 21, no. 3, pp. 207-220,  Nov. 1990.                            *
 *                                                                      *
 * File: Gabor.h                                                        *
 *                                                                      *
 * Revision history:                                                    *
 *  1. 04/14/93  - Started                                              *
 ************************************************************************/

#import <objc/Object.h>
#import "c_headers.h"

#define RECT_WINDOW 0       /* windowType       */
#define GAUSSIAN_WINDOW 1
#define RAISED_COS  2       /* Hanning      */


@interface Gabor:Object
{
  id    complexMatrix;          /* For calculating biorthog.    */
  int   windowType;             /* Window function      */
  int   windowLength;           /* Total length         */
  int   effectiveWidth;         /* Window effective width   */
  int   shiftParameter;         /* Window shift         */
  int   maximumM, maximumN;     /* # of expansion coeffs to use */
  double *window, *biorthogonalFunction;/* Window and its biorth. fn    */
  double *reconstructedSignal;
  double *shiftedWindow;
  COMPLEX *expansionCoefficients;   /* Gabor expansion coefficients */
}

/**********************
 * Initialization:    *
 **********************/
- initWithType:(int)type length:(int)length width:(int)width shift:(int)shift;
- allocateMemory;
- freeMemory;
- free;

/**********************
 * Calculating:       *
 **********************/
- calculateWindow;
- calculateBiorthogonal;
- calculateExpansionCoefficientsFor:(double *)xInput;
- calculateReducedCoefficientSetFor: (double *)xInput;
- reconstructSignal;
- resetCoefficients;

/**********************
 * Set parameters:    *
 **********************/
- setWindowType:   (int)type;
- setWindowLength: (int)length;
- setWindowWidth:  (int)width;
- setShift:        (int)shift;
- setMaximumM:     (int)m;
- setMaximumN:     (int)n;

/**********************
 * Get parameters:    *
 **********************/
- (int)windowType;
- (int)windowLength;
- (int)windowWidth;
- (int)windowShift;
- (int)mValue;
- (int)maximumM;
- (int)maximumN;
- (double *)windowFunction;
- (double *)biorthogonalFunction;
- (double *)reconstructedSignal;
- (double *)shiftedWindow:(int)shift;
- (COMPLEX *)expansionCoefficients;

@end